/************************************************************************
Contains CRS_GENERIC implementations of procedures:

!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALGORITHMS
!!!!!!!!!!!! FOR CRS and CRS_GENERIC
!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING

  lar_allocate_SM_and_LV - to allocate space for stiffness matrix and load vector
  lar_initialize_SM_and_LV - to initialize stiffness matrix and/or load vector
  lar_get_storage - to compute storage of SM, LV and preconditioner
  lar_fill_assembly_table_int_ent - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
  lar_assemble_SM_and_LV_with_table - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
  lar_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  lar_allocate_preconditioner - to allocate space for preconditioner
  lar_fill_preconditioner - to fill preconditioner
  lar_free_preconditioner - to free space of preconditioner structure
  lar_free_SM_and_LV - to free space of matrix structure

  lar_create_solver_structures_accel - utility to create data structures on GPU
  lar_get_crs_data - utility to get CRS parameters and pointers

  lar_compute_residual - to compute the residual of the not preconditioned
	system of equations, v = ( b - Ax )
  lar_compute_preconditioned_residual_crs_generic - to compute the residual of the
	preconditioned system of equations, v = M^-1 * ( b - Ax )
		where M^-1 corresponds directly to the stored preconditioner matrix
  lar_perform_BJ_or_GS_iterations - to perform one iteration of block Gauss-
	Seidel or block Jacobi algorithm, v_out = v_in + M^-1 * ( b - A * v_in )
  lar_perform_rhsub - to perform forward reduction and back-substitution for ILU
		   preconditioning

------------------------------
History:
  10.2015 - Kazimierz Chlon
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<assert.h>
#include<omp.h>

// API for linear algebra routines supporting MKB solver
#include <modfem/ls_mkb/lah_intf.h>

// for debugging mpi
//#include "../../include/pch_intf.h"

/* internal information for the solver module */
#include "lah_crs_generic.h"
#include <modfem/lin_alg_intf.h>

/* GLOBAL VARIABLES */
 int   itv_nr_crs_generic_matrices=0;   /* the number of matrices managed by module */
 itt_crs_generic_matrices itv_crs_generic_matrices[LAC_MAX_MATRICES];
											 /* array of CRS matrices */


/*---------------------------------------------------------
  lar_allocate_SM_and_LV_crs_generic - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
int lar_allocate_SM_and_LV_crs_generic( // returns: matrix index in itv_crs_generic_matrices array
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]==0 - iblock is a ghost block with no rows associated with
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
	)
{
  //printf("\n\n CRS-generic Solver \n\n\n");

/* auxiliary variable */
  int i, j, iblock, bloff, nrblocks, nrdofbl, nroffbl;
  int nnz=0,row_nnz, rowi=0,rowj=0,neig=0,neigi,n, key;
/* v2 */
  int val_index = 0, col_block_index=0,block_ptr_index=0;

/*++++++++++++++++ executable statements ++++++++++++++++*/


// create structures for a new matrix
  itt_crs_generic_matrices *it_matrix;
  itv_nr_crs_generic_matrices++;

  if(itv_nr_crs_generic_matrices>=LAC_MAX_MATRICES){
	printf("Too much (%d) matrices requested in las_block_intf! (correct lah_block.h)\n",
	   itv_nr_crs_generic_matrices);
  }

  it_matrix = &itv_crs_generic_matrices[itv_nr_crs_generic_matrices-1];

  it_matrix->Max_SM_size = Max_SM_size;
  it_matrix->Nrblocks = Nrblocks;
  it_matrix->Nrdofgl = Nrdof_glob;
  int inter_Nrdof=0;
  int inter_nrbl=0;



  int jblock;
  for(iblock=0;iblock<Nrblocks;iblock++){
	if(Nroffbl[iblock]>0){ // ghost blocks have Nroffbl[index]<=0
	 inter_Nrdof+=Nrdofbl[iblock];
	 inter_nrbl++;
	 nnz += Nrdofbl[iblock]*Nrdofbl[iblock]; // diagonal block
	 for(jblock=0; jblock< Nroffbl[iblock]; jblock++){ // off diagonal blocks
	  int nrdofbl = Nrdofbl[L_offbl[iblock][jblock]];
	  nnz+=Nrdofbl[iblock]*nrdofbl;
	  //printf("iblock=%d, jblock %d, nrdofbl %d  \n", iblock,jblock,nrdofbl);
	 }
	}
  }
  it_matrix->Nnz = nnz;
  it_matrix->Nrdof_internal = inter_Nrdof;


  //allocate crs structure
  it_matrix->crs_row_ptr = (int*)malloc( (inter_Nrdof+1)*sizeof(int));
  if( it_matrix->crs_row_ptr==NULL ) {
	printf("Not enough memory for allocating crs structure for row_ptr vector\n");
	exit(-1);
  }

/*
#ifdef PARALLEL
  printf("\n\n ### pcv_my_proc_id=%d nnz=%d Nrdof_glob=%d inter_Nrdof=%d Nrblocks=%d\n",pcv_my_proc_id,nnz,Nrdof_glob,inter_Nrdof_glob,Nrblocks);
#endif
  printf("\n\n ### nnz=%d Nrdof_glob=%d inter_Nrdof=%d inter_nrbl=%d\n",nnz,Nrdof_glob,inter_Nrdof,inter_nrbl);
/**/

  it_matrix->crs_col_ind = (int*)malloc(nnz*sizeof(int));
  if( it_matrix->crs_col_ind==NULL ) {
	printf("Not enough memory for allocating crs structure for col_ind vector\n");
	exit(-1);
  }

  it_matrix->crs_val= (double*)malloc(nnz*sizeof(double));
  if( it_matrix->crs_val==NULL ) {
	printf("Not enough memory for allocating crs structure for val vector\n");
	exit(-1);
  }


  it_matrix->rhs = (double*)malloc( (Nrdof_glob)*sizeof(double));
  if( it_matrix->rhs==NULL ) {
	printf("Not enough memory for allocating crs structure for rhs vector\n");
	exit(-1);
  }


  it_matrix->posg = (int*)malloc( (Nrblocks+1)*sizeof(int));
  if( it_matrix->posg==NULL ) {
	printf("Not enough memory for allocating crs structure for posg vector\n");
	exit(-1);
  }


  for(iblock=0;iblock<Nrblocks;iblock++){
	it_matrix->posg[iblock] = Posglob[iblock];
  }
  it_matrix->posg[Nrblocks] = it_matrix->Nrdofgl;

  //filling crs_row and crs_col structure

#pragma omp parallel for default(none) firstprivate(it_matrix, Nrdof_glob)
  for(iblock=0;iblock<Nrdof_glob;iblock++)it_matrix->rhs[iblock]=0.0;
#pragma omp parallel for default(none) firstprivate(it_matrix, nnz)
  for(iblock=0;iblock<nnz;iblock++)it_matrix->crs_val[iblock]=0.0;


  //it_matrix->block_ptr=it_matrix->crs_row_ptr;



  // for(iblock=0;iblock<Nrblocks;iblock++){
  // for(neigi=0;neigi<it_matrix->nr_neig[iblock];neigi++)
  // printf("iblock=%d, %d\n",iblock,L_offbl[iblock][neigi]);}


  //v3
  col_block_index = 0;
  it_matrix->crs_row_ptr[block_ptr_index]=col_block_index;
  block_ptr_index++;
  for(iblock=0;iblock<inter_nrbl;iblock++){
	//printf("(%d %d) ",col_block_index,val_index);

	//col_block_index+=it_matrix->nrdofbl[iblock];
	for(rowj=0;rowj<Nrdofbl[iblock];++rowj){
	  //printf("ww\n");
	  for(rowi=0;rowi<Nrdofbl[iblock];++rowi,++col_block_index){
	//it_matrix->crs_col_ind[col_block_index]=iblock*Nrdofbl[iblock]+rowi; // ERROR - corrected 2018
	it_matrix->crs_col_ind[col_block_index]=it_matrix->posg[iblock]+rowi;
	  }
	  //printf("dd\n");
	  for(neigi=0;neigi<Nroffbl[iblock];neigi++){
	for(rowi=0;rowi<Nrdofbl[L_offbl[iblock][neigi]];++rowi,++col_block_index){
	  it_matrix->crs_col_ind[ col_block_index]= it_matrix->posg[L_offbl[iblock][neigi]]+rowi;
	}

	  }
	  it_matrix->crs_row_ptr[block_ptr_index++]=col_block_index;
	}


  }
  it_matrix->crs_row_ptr[inter_Nrdof]=nnz;


  //printf("\n\nblock_ptr_index= %d\n\n",block_ptr_index);

  //it_matrix->crs_row_ptr[block_ptr_index]=col_block_index+1;

  //it_matrix->crs_col_ind[col_block_index]=0;


  /*kbw
//if(pcv_my_proc_id==1){
	for(iblock=0,i=0;iblock<nnz;++i){
	for(j=it_matrix->crs_row_ptr[i];j<it_matrix->crs_row_ptr[i+1];j++,iblock++)
	printf("(%d %d)",j,it_matrix->crs_col_ind[j]);
	//printf("(%d %d)",i,iblock);
	printf("\n");
	}
//}
	/**/

  /******** Sorting crs_col_ind - sir_sort_short line 1072 sis_mkb_intf.c ??? *********/

  for(j=0; j<inter_Nrdof;j++){
	n=it_matrix->crs_row_ptr[j+1];
	for(i=it_matrix->crs_row_ptr[j]+1;i<n;i++){
	  key=it_matrix->crs_col_ind[i];
	  rowi=i-1;
	  while(rowi>=it_matrix->crs_row_ptr[j] && it_matrix->crs_col_ind[rowi]>key){
		it_matrix->crs_col_ind[rowi+1] = it_matrix->crs_col_ind[rowi]; //przesuwanie elementów
		--rowi;
	  }
	  it_matrix->crs_col_ind[rowi+1] = key;
	}
  }

  /*kbw
//if(pcv_my_proc_id==1){
  printf("\nAfter sorting\n");
  for(iblock=0,i=0;iblock<nnz;++i){
	for(j=it_matrix->crs_row_ptr[i];j<it_matrix->crs_row_ptr[i+1];j++,iblock++)
	  printf("(%d %d)",j,it_matrix->crs_col_ind[j]);
	//printf("(%d %d)",i,iblock);
	printf("\n");
  }
  getchar();getchar();
//}
	/**/

/*ok_kbw*/
  printf("\nAllocated CRS-generic matrix %d: n = %d, nnz = %d\n",
	 itv_nr_crs_generic_matrices-1, it_matrix->Nrdofgl, it_matrix->Nnz);
/*kew*/

  return(itv_nr_crs_generic_matrices-1);
}

/**---------------------------------------------------------
  lar_create_solver_structures_accel_crs_generic - utility to create data structures on GPU
---------------------------------------------------------*/
int lar_create_solver_structures_accel_crs_generic(
  int Matrix_id   /* in: matrix ID */
)
{
  int info=1;


  // call OpenCL routine
#ifdef OPENCL_GPU


  //  pass created arrays to tmr routine ???
  /* tmr_ocl_create_solver_structures_crs(  ??? */
  /* tmr_create_solver_structures_accel_crs(  */
  /* 					 int Nrdof_glob, */
  /* 					 int nnz, */
  /* 					 int* crs_col_ind, */
  /* 					 int* crs_row_ptr, */
  /* 					 double* crs_val, */
  /* 					 double* rhs */
  /* 					 itd, itp ); */




#endif

  return(info);
}

/**---------------------------------------------------------
lar_get_crs_data_crs_generic - utility to get CRS parameters and pointers
---------------------------------------------------------*/
int lar_get_crs_data_crs_generic( /* returns: >0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int* Nrdof_glob_p, /* output pointers: nrdof, nnz etc. - self explanatory */
  int* Nnz_p,
  int** Crs_col_ind_p,
  int** Crs_row_ptr_p,
  double** Crs_val_p,
  double** Rhs_p
			  )
{

  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];

  *Nrdof_glob_p = it_matrix->Nrdofgl;
  *Nnz_p = it_matrix->Nnz;
  *Crs_col_ind_p = it_matrix->crs_col_ind;
  *Crs_row_ptr_p = it_matrix->crs_row_ptr;
  *Crs_val_p = it_matrix->crs_val;
  *Rhs_p = it_matrix->rhs;

  return(0);
}

/**--------------------------------------------------------
  lar_initialize_SM_and_LV_crs_generic - to initialize stiffness matrix and/or load vector
---------------------------------------------------------*/
int lar_initialize_SM_and_LV_crs_generic(
  int Matrix_id,   /* in: matrix ID */
  int Scope    /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  ){
	printf("\n\n CRS generic\n\n");

	itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];

	int nrdof_glob = it_matrix->Nrdofgl;
	int nrdof_internal =it_matrix->Nrdof_internal;
	int iblock;

#pragma omp parallel for default(none) firstprivate(it_matrix, nrdof_glob)
	for(iblock=0;iblock<nrdof_glob;iblock++) {
	  it_matrix->rhs[iblock]=0.0;
	}

	if(Scope==LAC_SCOPE_SM_AND_LV){
	  // parallel loop is the same as in matrix-vector product to properly assign
	  // parts of it_matrix->crs_val to threads (for NUMA architectures)
#pragma omp parallel for default(none) firstprivate(it_matrix, nrdof_glob,nrdof_internal)
	  for(iblock=0;iblock<nrdof_internal;iblock++) {
	int icrs=0;
	for(icrs=it_matrix->crs_row_ptr[iblock];
		icrs<it_matrix->crs_row_ptr[iblock+1];
		icrs++) it_matrix->crs_val[icrs]=0.0;

	  // quick and dirty
/* #pragma omp parallel for default(none) firstprivate(it_matrix, nrdof_glob) */
/*       for(iblock=0;iblock<it_matrix->Nnz;iblock++) it_matrix->crs_val[iblock]=0.0; */

	  }
	}

/*ok_kbw*/
  printf("\n\nInitialized CRS-generic matrix %d: n = %d, nnz = %d, half bandwidth = %d\n",
	 Matrix_id, it_matrix->Nrdofgl, it_matrix->Nnz, it_matrix->Half_bandwidth);
  printf("The average number of non-zeros in a single row: %d\n",
	 it_matrix->Nnz/it_matrix->Nrdofgl);
/*kew*/
	return 0;
  }

/**--------------------------------------------------------
  lar_get_storage_crs_generic - to compute storage of SM, LV and preconditioner
---------------------------------------------------------*/
double lar_get_storage_crs_generic( /* returns: storage in MB */
  int Matrix_id   /* in: matrix ID */
				   )
{
	return -1.0;
}


/**-----------------------------------------------------------
  lar_fill_assembly_table_int_ent_crs_generic - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
int lar_fill_assembly_table_int_ent_crs_generic(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int *Assembly_table_int_ent /* part of the global assembly table */
){

/* get pointer to the level structure */
  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];

/* compute local stiffness matrix nrdof */
  int nrdof=0; int iblock;
  int max_smallbl_size=0; int smallbl_size;
  for(iblock=0;iblock<Nr_dof_bl;iblock++){
	//nrdof += L_bl_nrdof[iblock];
	smallbl_size = it_matrix->posg[L_bl_id[iblock]+1] - it_matrix->posg[L_bl_id[iblock]];
	nrdof += smallbl_size;
	if(max_smallbl_size<smallbl_size)max_smallbl_size=smallbl_size;
  }

/*kbw
  printf("In lar_crs_fill_assembly_table_int_ent: matrix_id %d\n",
	 Matrix_id);
  printf("nr_dof_ent_loc %d, position %u, nrdof %d\n",
	 Nr_dof_bl, Assembly_table_int_ent, nrdof);
  int ibl;
  for(ibl=0;ibl<Nr_dof_bl; ibl++) printf("bl_id %d (global %d)\n",
					 L_bl_id[ibl], it_matrix->posg[L_bl_id[ibl]]);
  printf("\n");
  getchar();getchar();
/*kew*/

  int idofbl;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {
	int bl_id_y = it_matrix->posg[L_bl_id[idofbl]];
	//int nrdof_i = it_matrix->posg[L_bl_id[idofbl]+1] - bl_id_y;

	//int posloc_j=0;
	int jdofbl;
	for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) {
	  int bl_id_x = it_matrix->posg[L_bl_id[jdofbl]];
	  //int nrdof_j = it_matrix->posg[L_bl_id[jdofbl]+1] - bl_id_x;

	  //int idof;
	  //for(idof=0;idof<nrdof_i;idof++){

		//int jaux = posloc_j+(posloc_i+idof)*nrdof;
	//int jdof;
		//for(jdof=0;jdof<nrdof_j;jdof++) {

	  int jaux = jdofbl+idofbl*Nr_dof_bl;

/*kbw
	  printf("entry %d (idofbl %d, bl_id_y %d, jdofbl %d bl_id_x %d)\n",
		 jaux, idofbl, bl_id_y, jdofbl, bl_id_x);
/*kew*/

#ifdef PARALLEL
  if(bl_id_x<it_matrix->Nrdof_internal){
#endif

		  //int jcrs=bl_id_x+jdof;
	  //int icrs;
		  //for(icrs=it_matrix->crs_row_ptr[jcrs];
	   //   icrs<it_matrix->crs_row_ptr[jcrs+1];
	   //   icrs++){

		  //  if((it_matrix->crs_col_ind[icrs])==(bl_id_y+idof)){

	  int icrs,index,offset=0;
	  for(icrs=it_matrix->crs_row_ptr[bl_id_x];
	  icrs<it_matrix->crs_row_ptr[bl_id_x+1];
	  ++icrs,++offset){

	if(it_matrix->crs_col_ind[icrs]==bl_id_y){

/*kbw
		  printf("found for entry %d position in CRS %d\n", jaux, icrs);
/*kew*/

		 // Assembly_table_int_ent[jaux+jdof]=icrs;

		//  break;
	   // }
	 // }

if(max_smallbl_size>1){
	// for max_smallbl_size>1 we need two information:
	//		offset - index of nonzero value in crs_row (on first 15bits)
	//		bl_id_x - row index (on last 16bits)
	//printf("(soffset=%d index=%d) ",offset,bl_id_y);
	if(bl_id_y > 65535 || offset > 32767) {
	  printf("Fatal error num_row or number of elements in row is bigger than 2^16 (65536)\n");
	  printf("Please turn off assemble table\n");
	  exit(-1);
	}
	  offset=(offset<<16);
	  index=offset+bl_id_x;


	  Assembly_table_int_ent[jaux]=index;
}else
	  Assembly_table_int_ent[jaux]=icrs;

	  break;
		}
	  }
	 //posloc_j += nrdof_j;

#ifdef PARALLEL
  }else
	 Assembly_table_int_ent[jaux]=-1;
#endif
	}
  }

/*<<<<<<<<<< RHS part of Assembly_table exchanged for local_to_global array >>>>>>>>>>>>*/
/*   int jdofbl; */
/*   for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) { */
/*     int bl_id_x=it_matrix->posg[L_bl_id[jdofbl]]; */

/*     int jaux = Nr_dof_bl*Nr_dof_bl+jdofbl; */

/* /\*kbw */
/* 	  printf("filling entry %d for RHS block %d (block %d)\n",  */
/* 		 jaux, jdofbl, bl_id_x); */
/* /\*kew*\/ */

/*     Assembly_table_int_ent[jaux]=bl_id_x; */

/*   } */

/*kbw
	  getchar();
/*kew*/

  return(1);
}

/*------------------------------------------------------------
  lar_assemble_SM_and_LV_with_table_crs_generic - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_with_table_crs_generic(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Scope,    /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  int* Local_to_global_int_ent, /* part of the global local_to_global table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
  )
{

  //printf("\n\n CRS-assembling with table \n\n\n");

  /* get pointer to the matrix structure */
  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];
  int idof,jdof;
  int smallbl_size, max_smallbl_size=0;

/* compute local stiffness matrix nrdof */
  int nrdof=0; int iblock;
  for(iblock=0;iblock<Nr_dof_bl;iblock++){
	//nrdof += L_bl_nrdof[iblock];
	smallbl_size = it_matrix->posg[Local_to_global_int_ent[iblock]+1] -
				   it_matrix->posg[Local_to_global_int_ent[iblock]];

	/*jbw
	printf("iblock = %d / Nr_dof_bl = %d / smallbl_size =%d\n",
	   iblock,Nr_dof_bl,smallbl_size);
	/*jbw*/
	nrdof +=smallbl_size;
	if(max_smallbl_size<smallbl_size)max_smallbl_size=smallbl_size;
  }

  /*kbw
  printf("In lar_crs_assembly: matrix_id %d, nr_dof_ent_loc %d, nrdofs %d\n",
	 Matrix_id, Nr_dof_bl, nrdof);
  int ibl;
  for(ibl=0;ibl<Nr_dof_bl; ibl++){
	printf("bl_id %d (global %d) nrdof %d\n",
	   Local_to_global_int_ent[ibl],
	   it_matrix->posg[Local_to_global_int_ent[ibl]],
	   it_matrix->posg[Local_to_global_int_ent[ibl]+1] - it_matrix->posg[Local_to_global_int_ent[ibl]]);
  }
  printf("\n");
  getchar();
/*kew*/

  int posloc_i=0;
  int idofbl;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {
	int bl_id_y = it_matrix->posg[Local_to_global_int_ent[idofbl]];
	int nrdof_i = it_matrix->posg[Local_to_global_int_ent[idofbl]+1] - bl_id_y;

	int posloc_j=0;
	int jdofbl;
	for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) {
	  int bl_id_x = it_matrix->posg[Local_to_global_int_ent[jdofbl]];
	  int nrdof_j = it_matrix->posg[Local_to_global_int_ent[jdofbl]+1] - bl_id_x;

	  int jaux = jdofbl+idofbl*Nr_dof_bl;

	  int icrs = Assembly_table_int_ent[jaux];

#ifdef PARALLEL
  if(icrs!=-1){
#endif

if(max_smallbl_size>1){  //vector problem


	  int offset=(icrs>>16);
	  int index=(icrs-(offset<<16));
	  //printf("(roffset=%d index=%d)",offset,index);

for(idof=0;idof<nrdof_i;++idof){
	jaux = posloc_j+(posloc_i+idof)*nrdof;
	for(jdof=0;jdof<nrdof_j;++jdof){

		icrs=it_matrix->crs_row_ptr[index+jdof];
/*kbw
	  if(fabs(Stiff_mat[jaux+jdof])>1.e-6){
		printf("found for entry %d position in CRS %d\n", icrs+offset+jdof, icrs);
	  }
/*kew*/
/*kbw
	  if(fabs(Stiff_mat[jaux+jdof])>1.e-6){
	  printf("filling with %lf: before %lf\n",
		 Stiff_mat[jaux+idof], it_matrix->crs_val[icrs+offset+jdof]);
}
/*kew*/
	   it_matrix->crs_val[icrs+offset+idof] += Stiff_mat[jaux+jdof];
/*kbw
	  if(fabs(Stiff_mat[jaux+jdof])>1.e-6){
	  printf("filling with %lf: after %lf\n",
		 Stiff_mat[jaux+idof], it_matrix->crs_val[icrs+offset+jdof]);
}
/*kew*/
	}
}
  }else //max_smallbl_size==1 , scalar problem
	it_matrix->crs_val[icrs] += Stiff_mat[jaux];

#ifdef PARALLEL
  }
#endif

	posloc_j += nrdof_j;
	}


	for (idof=0;idof<nrdof_i;idof++) {

/*kbw
	  if(fabs(Rhs_vect[posloc_i+idof])>1.e-6){
	printf("filling entry %d for RHS block %d\n",
		   bl_id_y+idof, posloc_i+idof);
	  }
/*kew*/
/*kbw
	  if(fabs(Rhs_vect[posloc_i+idof])>1.e-6){
	  printf("filling with %lf: before %lf\n",
		 Rhs_vect[posloc_i+idof], it_matrix->rhs[bl_id_y+idof]);
}
/*kew*/

	  it_matrix->rhs[bl_id_y+idof] += Rhs_vect[posloc_i+idof];

/*kbw
	  if(fabs(Rhs_vect[posloc_i+idof])>1.e-6){
	  printf("filling with %lf: after %lf\n",
		 Rhs_vect[posloc_i+idof], it_matrix->rhs[bl_id_y+idof]);
}
/*kew*/

	}
	posloc_i += nrdof_i;
  }

/*kbw
  getchar();
/*kew*/

  return(1);
}

/*------------------------------------------------------------
  lar_assemble_SM_and_LV_crs_generic - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_crs_generic(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Scope,    /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
  )
{

  /* pointers to solver structure and stiffness matrix blocks' info */
  //itt_blocks *block_i, *block_j;	/* to simplify typing */

  /* auxiliary variables */
  int iblock, jblock, nrdof_i, nrdof_j, posbl_i, posbl_j;
  int posloc_i, posloc_j, nrdof, nrdof_glob, ibl;
  int i,j,k, iaux, jaux;
  int jcrs,icrs,bl_id_y,bl_id_x,idofbl,jdofbl,val_index,idof,jdof,ineig;

  //printf("\n\n CRS-assembling \n\n\n");

/* get pointer to the level structure */
  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];
  int nrdof_internal =it_matrix->Nrdof_internal;

   //printf("\n\n\n\n\nasembling %d !!!!!!!!!!! \n\n",it_matrix->crs_row_ptr[1]);

  /*kbw
  printf("In lar_crs_assembly: matrix_id %d, nr_dof_ent_loc %d\n",
	 Matrix_id, Nr_dof_bl);
  for(ibl=0;ibl<Nr_dof_bl; ibl++){
	printf("bl_id %d (global %d) nrdof %d\n",
	   L_bl_id[ibl], it_matrix->posg[L_bl_id[ibl]], L_bl_nrdof[ibl]);
  }
  printf("\n");
  getchar();
/*kew*/

   // //******************CRS************************************//
// /* compute local stiffness matrix nrdof */
  nrdof=0;
  for(iblock=0;iblock<Nr_dof_bl;iblock++){
	nrdof += L_bl_nrdof[iblock];
  }


  posloc_i=0;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {
	bl_id_y=it_matrix->posg[L_bl_id[idofbl]];
	nrdof_i = L_bl_nrdof[idofbl];

	posloc_j=0;
	for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) {
	  bl_id_x=it_matrix->posg[L_bl_id[jdofbl]];
	  nrdof_j = L_bl_nrdof[jdofbl];
 #ifdef PARALLEL
  if(bl_id_x<nrdof_internal)
#endif
	  for(idof=0;idof<nrdof_i;idof++){

		jaux = posloc_j+(posloc_i+idof)*nrdof;
		for(jdof=0;jdof<nrdof_j;jdof++) {

/*kbw
	  printf("entry %d (idofbl %d, bl_id_y %d, jdofbl %d bl_id_x %d)\n",
		 jaux, idofbl, bl_id_y, jdofbl, bl_id_x);
/*kew*/

		  jcrs=bl_id_x+jdof;
		  for(icrs=it_matrix->crs_row_ptr[jcrs];
		  icrs<it_matrix->crs_row_ptr[jcrs+1];
		  icrs++){
			//if(icrs>=it_matrix->Nnz){
		//  printf("\n\n icrs=%d > %d=Nnz",icrs,it_matrix->Nnz);
		//  exit(0);
		//}
			if((it_matrix->crs_col_ind[icrs])==(bl_id_y+idof)){

/*kbw
		  if(fabs(Stiff_mat[jaux+jdof])>1.e-6){
		printf("found for entry %d position in CRS %d\n", jaux+jdof, icrs);
		  }
/*kew*/
/*kbw
		  if(fabs(Stiff_mat[jaux+jdof])>1.e-6){
		printf("filling with %lf: before %lf\n",
			   Stiff_mat[jaux+jdof], it_matrix->crs_val[icrs]);
		  }
/*kew*/

		  it_matrix->crs_val[icrs] += Stiff_mat[jaux+jdof];

/*kbw
		  if(fabs(Stiff_mat[jaux+jdof])>1.e-6){
		printf("filling with %lf: after %lf\n",
			   Stiff_mat[jaux+jdof], it_matrix->crs_val[icrs]);
		  }
/*kew*/

		  break;
		}

		  }
		}
	  }
	  posloc_j += nrdof_j;

/*kbw
  getchar();
/*kew*/

	}

	for (idof=0;idof<nrdof_i;idof++) {

/*kbw
	  if(fabs(Rhs_vect[posloc_i+idof])>1.e-6){
	printf("filling entry %d for RHS block %d (block %d)\n",
		   bl_id_y+idof, posloc_i+idof);
	  }
/*kew*/
/*kbw
	  if(fabs(Rhs_vect[posloc_i+idof])>1.e-6){
	printf("filling with %lf: before %lf\n",
		   Rhs_vect[posloc_i+idof], it_matrix->rhs[bl_id_y+idof]);
	  }
/*kew*/

	  /* assemble right hand side block's entries */
	  it_matrix->rhs[bl_id_y+idof] += Rhs_vect[posloc_i+idof];

/*kbw
	  if(fabs(Rhs_vect[posloc_i+idof])>1.e-6){
	printf("filling with %lf: after %lf\n",
		   Rhs_vect[posloc_i+idof], it_matrix->rhs[bl_id_y+idof]);
	  }
/*kew*/

	}
	posloc_i += nrdof_i;
  }

/*kbw
  getchar();
/*kew*/

  // printf("\n\n");
  // for(iblock=0;iblock<56;iblock++)printf("%d-%lf\n",it_matrix->crs_col_ind[iblock],it_matrix->crs_val[iblock]);
  // printf("\nFFF\n");
  //for(iblock=0;iblock<it_matrix->Nrdofgl;iblock++)printf("%lf  ",it_matrix->rhs[iblock]);
  return(1);
}

/**--------------------------------------------------------
lar_allocate_preconditioner_crs_generic - to allocate space for preconditioner
---------------------------------------------------------*/
int lar_allocate_preconditioner_crs_generic( /* returns:   >0 number of diagonal blocks */
						  /*	       <=0 - error */
  int Matrix_id,   /* in: matrix ID */
  int Precon,      /* in: type of preconditioner (lah_block.h line circa 45) */
  int ILU_k // in: for ILU(k) - k;
  )
{

  //
  //!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALL ALGORITHMS
  //!!!!!!!!!!!! FOR CRS and CRS_GENERIC
  //!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING
  //


  int i,Nrdofgl,kk;
  int block_row_gl, sm_block, block_col_gl;

  itt_crs_generic_matrices *it_matrix;
  it_matrix = &itv_crs_generic_matrices[Matrix_id];
  Nrdofgl = it_matrix->Nrdofgl;
  it_matrix->Precon=Precon;
  int nrdof_internal =it_matrix->Nrdof_internal;

  //printf("\nCRS_lar_allocate_preconditioner\n");

  it_matrix->diag_ptr = (int*)malloc( nrdof_internal*sizeof(int));
  if( it_matrix->diag_ptr==NULL ) {
	printf("Not enough memory for allocating crs structure in preconditioner for diag_ptr vector\n");
	exit(-1);
  }

  it_matrix->diag_precon = (double*)malloc( nrdof_internal*sizeof(double));
  if( it_matrix->diag_precon==NULL ) {
	printf("Not enough memory for allocating crs structure in preconditioner for diag_precon vector\n");
	exit(-1);
  }



  for (block_row_gl = 0;  block_row_gl < nrdof_internal; block_row_gl++){

	for(sm_block = it_matrix->crs_row_ptr[ block_row_gl ];
		sm_block < it_matrix->crs_row_ptr[ block_row_gl+1 ];
		sm_block++){

	  block_col_gl = it_matrix->crs_col_ind[ sm_block ];

	  if(block_col_gl == block_row_gl){
		it_matrix->diag_ptr[block_row_gl] = sm_block;
	  }
	}
  }

}


/**--------------------------------------------------------
  lar_fill_preconditioner_crs_generic - to fill preconditioner
---------------------------------------------------------*/
int lar_fill_preconditioner_crs_generic(
  int Matrix_id   /* in: matrix ID */
	)
{

  //
  //!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALL ALGORITHMS
  //!!!!!!!!!!!! FOR CRS and CRS_GENERIC
  //!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING
  //

  int i,Nrdofgl,diag_index;
  double *pivots;


  itt_crs_generic_matrices *it_matrix;
  it_matrix = &itv_crs_generic_matrices[Matrix_id];
  Nrdofgl = it_matrix->Nrdofgl;
  int nrdof_internal =it_matrix->Nrdof_internal;

  //printf("\nCRS_lar_fill_preconditioner\n");

  for(i=0;i<nrdof_internal;++i){
	diag_index=it_matrix->diag_ptr[i];
	it_matrix->diag_precon[i]=1.0/it_matrix->crs_val[diag_index];

  }
  return 0;
}

/**--------------------------------------------------------
  lar_free_preconditioner_crs_generic - to free space for a block structure
---------------------------------------------------------*/
int lar_free_preconditioner_crs_generic(
  int Matrix_id   /* in: matrix ID */
  )
{
  //
  //!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALL ALGORITHMS
  //!!!!!!!!!!!! FOR CRS and CRS_GENERIC
  //!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING
  //

  itt_crs_generic_matrices *it_matrix;
  it_matrix = &itv_crs_generic_matrices[Matrix_id];

  free(it_matrix->diag_ptr);
  free(it_matrix->diag_precon);
  return 0;
}


/*---------------------------------------------------------
lar_compute_residual_crs_generic - to compute the residual of the system of equations,
	v = ( b - Ax ) (used also to compute the product v = -Ax)
---------------------------------------------------------*/
void lar_compute_residual_crs_generic(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: indicator whether to use RHS */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Ndof, 	/* in: number of unknowns (components of x) */
  double* X, 	/* in: input vector (may be NULL if Ini_zero==0) */
  double* B,	/* in:  the rhs vector, if NULL take rhs */
				/*      from block data structure ( B is not taken */
				/*      into account if Use_rhs!=1) */
  double* V 	/* out: v = b-Ax */
	)
{

  //
  //!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALL ALGORITHMS
  //!!!!!!!!!!!! FOR CRS and CRS_GENERIC
  //!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING
  //

/* constants */
  int ione = 1;
  double done = 1.;

  int i,j,k;
  int Nrdofgl, *z;
  double sum;
  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];
  Nrdofgl = it_matrix->Nrdofgl;
  int nrdof_internal =it_matrix->Nrdof_internal;
  assert(Ndof==it_matrix->Nrdofgl);

  //printf("\nCRS_lar_compute_residual\n");


  // for(i=0;i<56;i++)printf("%d-%lf  ",it_matrix->crs_col_ind[i],it_matrix->crs_val[i]);

  // printf("\n\n ************* in lar_compute_residual_crs  %d %d **************\n\n",Ndof, Nrdofgl);


#pragma omp parallel for default(none) firstprivate(Ndof,it_matrix,V)
  for(i=0;i<Ndof;i++) V[i]=0.0;

  /***** Standard */

  if(!Ini_zero){
#pragma omp parallel for default(none) firstprivate(Ndof,it_matrix,V,X,nrdof_internal)
	for(i=0;i<nrdof_internal;++i){
	  int k;
	  for(k=it_matrix->crs_row_ptr[i];k<it_matrix->crs_row_ptr[i+1];++k){
	V[i] -= it_matrix->crs_val[k]*X[it_matrix->crs_col_ind[k]];
	  }
	}
  }
  if(Use_rhs==1){
	if(B==NULL){
#pragma omp parallel for default(none) firstprivate(Ndof, V, it_matrix)
	  for(i=0;i<Ndof;i++) V[i] += it_matrix->rhs[i];
	  //      daxpy_(&Ndof, &done, it_matrix->rhs, &ione, V, &ione);
	}
	else{
#pragma omp parallel for default(none) firstprivate(Ndof, V, B)
	  for(i=0;i<Ndof;i++) V[i] += B[i];
	  //      daxpy_(&Ndof, &done, B, &ione, V, &ione);
	}
  }
  /*****/


  // printf("\nRozwiązanie:\n\n");
  // for(i=0;i<Ndof;i++)printf("%lf ",V[i]);

}

/**--------------------------------------------------------
lar_compute_preconditioned_residual_crs_generic - to compute the residual of the
	preconditioned system of equations, v = M^-1 * ( b - Ax )
		where M^-1 corresponds directly to the stored preconditioner matrix
---------------------------------------------------------*/
void lar_compute_preconditioned_residual_crs_generic (
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: indicator whether to use RHS */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Ndof, 	/* in: number of unknowns (components of x) */
  double* X, 	/* in: initial guess vector */
  double* B,	/* in:  the rhs vector, if NULL take rhs */
				/*      from block data structure */
  double* V 	/* out: preconditioned residual, v = M^-1*(b-Ax) */
								  )
  {

  //
  //!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALL ALGORITHMS
  //!!!!!!!!!!!! FOR CRS and CRS_GENERIC
  //!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING
  //

}

/**--------------------------------------------------------
lar_perform_BJ_or_GS_iterations_crs_generic - to perform one iteration of block Gauss-Seidel
	or block Jacobi algorithm:  v_out = v_in + M^-1 * ( b - A * v_in )
		where M^-1 results from stored preconditioner matrix and the algorithm
---------------------------------------------------------*/
void lar_perform_BJ_or_GS_iterations_crs_generic(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: 0 - no rhs, 1 - with rhs */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Nr_prec,  /* in: number of preconditioner iterations */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* in,out: vector of unknowns updated */
				/* during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
				/*      from block data structure */
	)
{
  //
  //!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALL ALGORITHMS
  //!!!!!!!!!!!! FOR CRS and CRS_GENERIC
  //!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING
  //
  itt_crs_generic_matrices *it_matrix;
  it_matrix = &itv_crs_generic_matrices[Matrix_id];

  int nrdof_internal =it_matrix->Nrdof_internal;

  int i_loop;		/* counter for loops for GS */
  int ione = 1;
  double done = 1.;
  int i,ret;

  int block_row_gl, sm_block, block_col_gl,block_val_ind;

  //printf("\nCRS_lar_perform_BJ_or_GS_iterations\n");
  //fflush(stdout);


  /*kbw

//if(pcv_my_proc_id==1){
#define SAVE_LOG_TO_FILE
#ifdef SAVE_LOG_TO_FILE

  char filename[256];
  const char *program_name = getenv("_");

  sprintf(filename,"CRS_GENERIC_Solution.txt",program_name);

  static int file_generic_opened_by_program = 0;
  FILE *file = NULL;
  if(file_generic_opened_by_program == 0) {
	file = fopen(filename,"w+");
	file_generic_opened_by_program = 1;
  } else {
	file = fopen(filename,"a+");
  }
#else
  FILE *file = stdout;
#endif

  int row_ind;
  fprintf(file,"\n\n CRS on entry: row, column, value\n");
  for(i=0;i<nrdof_internal;++i){
	for(row_ind=it_matrix->crs_row_ptr[i];row_ind<it_matrix->crs_row_ptr[i+1];row_ind++){
	  fprintf(file,"%6d%6d%20.10lf\n",i,it_matrix->crs_col_ind[row_ind],it_matrix->crs_val[row_ind]);
	}
  }

  fprintf(file,"\n RHS on entry: index, value\n");
  for(i=0;i<Ndof;++i){
	fprintf(file,"%6d%20.10lf\n",i,it_matrix->rhs[i]);
  }

  fprintf(file,"\n SOLUTION [before solve] on entry: index, value\n");
  for(i=0;i<Ndof;++i) {
	fprintf(file,"%6d%20.10lf\n",i,V[i]);
  }

#ifdef SAVE_LOG_TO_FILE
  fclose(file);
#endif
//}//if(pcv_my_proc_id==1)

  // Exit - uncomment if needed for testing
  //exit(0);

/*kew*/



  /* loop over loops */
  for(i_loop=0;i_loop<Nr_prec;i_loop++){

	/* loop over all blocks */
#pragma omp parallel default(none) firstprivate(V,B,it_matrix,i_loop,Use_rhs,Ini_zero,Ndof,ione,done,ret, nrdof_internal) \
  private(block_row_gl, sm_block, block_col_gl, block_val_ind, i)
	{
	  // each thread gets its copy of input vector
	  double *vtemp_GS_BJ;		/* temporary vector*/
	  vtemp_GS_BJ = (double*)malloc(Ndof*sizeof(double));

#ifdef _OPENMP
	  int my_id = omp_get_thread_num();
	  int portion = nrdof_internal/omp_get_num_threads() + 1;
#else
	  int my_id = 0;
	  int portion = nrdof_internal;
#endif


	  int my_first_x = my_id*portion - it_matrix->Half_bandwidth - 10;
	  if(my_first_x<0) my_first_x = 0;
	  int my_last_x = (my_id+1)*portion + it_matrix->Half_bandwidth + 10;
	  if(my_last_x>nrdof_internal) my_last_x=nrdof_internal;
	  // initialize first it_matrix->Half_bandwidth in lar_allocate_SM_and_LV
	  //for(i=my_first_x; i<my_last_x; i++) vtemp_GS_BJ[i]=V[i];
	  for(i=0;i<Ndof;i++) vtemp_GS_BJ[i]=V[i];



	  //#pragma omp barrier

	  int my_first = my_id*portion;
	  if(my_first<0) my_first = 0;
	  int my_last = (my_id+1)*portion;
	  if(my_last>nrdof_internal) my_last=nrdof_internal;

	  //#pragma omp for
	  //for (block_row_gl = 0; block_row_gl < nrdof_internal; block_row_gl++){
	// initialize first it_matrix->Half_bandwidth in lar_allocate_SM_and_LV
	  for (block_row_gl = my_first; block_row_gl < my_last; block_row_gl++){

		double vloc=0.0;

	//Lx
	if(it_matrix->Precon==BLOCK_JACOBI){
	  for(sm_block = it_matrix->crs_row_ptr[block_row_gl];
		  sm_block < it_matrix->diag_ptr[block_row_gl]; sm_block++){

		block_col_gl = it_matrix->crs_col_ind[sm_block];

		assert(block_row_gl != block_col_gl);

			vloc -= it_matrix->crs_val[ sm_block] * V[ block_col_gl];
	  }
	}
	else{
	  for(sm_block = it_matrix->crs_row_ptr[block_row_gl];
		  sm_block < it_matrix->diag_ptr[block_row_gl]; sm_block++){

		block_col_gl = it_matrix->crs_col_ind[sm_block];

		assert(block_row_gl != block_col_gl);

		vloc -= it_matrix->crs_val[ sm_block] * vtemp_GS_BJ[ block_col_gl];
	  }
		}
	  //Ux
		for(sm_block = it_matrix->diag_ptr[block_row_gl]+1;
			sm_block < it_matrix->crs_row_ptr[block_row_gl+1];sm_block++){

		  block_col_gl = it_matrix->crs_col_ind[sm_block];

		  assert(block_row_gl != block_col_gl);

		  vloc -= it_matrix->crs_val[ sm_block ] * V[ block_col_gl];
		}

	// now RHS: b - Ax

		if(Use_rhs==1){
		  if(B==NULL){
			vloc+=it_matrix->rhs[block_row_gl];
		  }
		  else{
			vloc+=B[block_row_gl];
		  }
		}

		vloc*=it_matrix->diag_precon[block_row_gl];
		vtemp_GS_BJ[block_row_gl]=vloc;

	  } // end of parallel loop over rows

#pragma omp barrier

	  //#pragma omp for
	  //for (block_row_gl = 0;  block_row_gl < nrdof_internal; block_row_gl++){
	  for (block_row_gl = my_first;  block_row_gl < my_last; block_row_gl++){
		V[block_row_gl]=vtemp_GS_BJ[block_row_gl];
	  }
	  free(vtemp_GS_BJ);

	}// end parallel region
  } //loop over preconditioner iterations


  /*kbw
{

	printf("V on leaving\n");
	for(i=0;i<Ndof;i++) printf("%20.15lf",V[i]);
	//for(i=0;i<Ndof;i++) printf("%10.6lf",V[i]);
	printf("\n");
   getchar();
}
/*kew*/

}

/**--------------------------------------------------------
lar_perform_rhsub_crs_generic - to perform forward reduction and back-substitution for ILU
		   preconditioning
---------------------------------------------------------*/
void lar_perform_rhsub_crs_generic(
  int Matrix_id,   /* in: matrix ID */
  int Ndof,	   /* in: number of unknowns (components of v*) */
  double* V,	   /* in,out: vector of unknowns updated */
				   /* during the loop over subdomains */
  double* B	   /* in:  the rhs vector, if NULL take rhs */
				   /*      from block data structure */
	)
{
  int i,j,Nrdofgl;
  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];
  Nrdofgl = it_matrix->Nrdofgl;

  int  *z;
  double sum;

  //
  //!!!!!!!!!!!! THERE SHOULD BE ONE IMPLEMENTATION OF ALL ALGORITHMS
  //!!!!!!!!!!!! FOR CRS and CRS_GENERIC
  //!!!!!!!!!!!! THEY SHOULD DIFFER ONLY IN CREATION/FILLING/DESTROYING
  //


  //printf("\nCRS_lar_perform_rhsub\n");


// z=(int*)malloc(Nrdofgl*sizeof(int));
// for(i=0;i<Nrdofgl;i++)z[i]=0.0;
// for(i=0;i<Nrdofgl;i++){
  // sum=0.0;
  // for(j=it_matrix->crs_row_ptr[i];j<it_matrix->diag_ptr[i];j++)
	// sum+=it_matrix->crs_val[j]*z[it_matrix->crs_col_ind[j]];
  // //printf("%d %lf\n",i,B[i]);
  // if(B==NULL)z[i]=it_matrix->pivots[i]*(it_matrix->rhs[i]-sum);
  // else z[i]=it_matrix->pivots[i]*(B[i]-sum);
// }
// for(i=Nrdofgl-1;i>=0;--i){
  // sum=0.0;
  // for(j=it_matrix->diag_ptr[i]+1;j<it_matrix->crs_row_ptr[i+1];j++){
	// sum+=it_matrix->crs_val[j]*V[it_matrix->crs_col_ind[j]];
	// V[i]=z[i]-it_matrix->pivots[i]*sum;
  // }
// }


}



/*---------------------------------------------------------
  lar_free_SM_and_LV_crs_generic - to free space for a block structure
---------------------------------------------------------*/
int lar_free_SM_and_LV_crs_generic(
  int Matrix_id   /* in: matrix ID */
  )
{
  // to keep matrix management as simple as possible it is required that
  // the matrices are freed in the reverse order wrt creating (as in heap)
  if(Matrix_id!=itv_nr_crs_generic_matrices-1){
	printf("Matrix destroyed %d (from %d) is not the last in lar_free_SM_and_LV_bcrs!!!\n",
	   Matrix_id, itv_nr_crs_generic_matrices);
	exit(-1);
  }


//printf("\n\n CRS-free \n\n");
  /* pointer to solver structure */
	itt_crs_generic_matrices *it_matrix;
	it_matrix = &itv_crs_generic_matrices[Matrix_id];
	  free(it_matrix->crs_val);
	  free(it_matrix->crs_col_ind);
	  free(it_matrix->crs_row_ptr);
	  free(it_matrix->rhs);
	  free(it_matrix->posg);
	itv_nr_crs_generic_matrices--;
  //printf("\n\n CRS-free-end \n\n");
  return(0);
}

/*---------------------------------------------------------
  lar_get_SM_and_LV_crs_from_crs_generic - convert crs_generic matrix format to crs
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_crs_from_crs_generic( // returns: flag - 0 if not needed delete crs matrices
  int Matrix_id, /* in: matrix ID */
  int offset, /* in: offset in crs_row and crs_col matrix */
  int** crs_row,	   /* out: matrix of rows in crs */
  int** crs_col,	   /* out: matrix of column in crs */
  double** crs_val,	   /* out: matrix of value in crs */
  double** rhs	   /* out: rhs vector */
			   )
{
  int i;
  int info;

  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];
  *rhs=it_matrix->rhs;
  if(offset==0){
	*crs_row=it_matrix->crs_row_ptr;
	*crs_col=it_matrix->crs_col_ind;
	*crs_val=it_matrix->crs_val;
	info=0;
  }else{
	// allocate new crs structure
	if(*crs_row==NULL){
	  *crs_row = (int*)malloc( (it_matrix->Nrdofgl+1)*sizeof(int));
	  if( *crs_row==NULL ) {
		printf("Not enough memory for allocating crs structure for row_ptr vector\n");
		exit(-1);
	  }
	}
	if(*crs_col==NULL){
	  *crs_col = (int*)malloc(it_matrix->Nnz*sizeof(int));
	  if( *crs_col==NULL ) {
		printf("Not enough memory for allocating crs structure for col_ind vector\n");
		exit(-1);
	  }
	}
	if(*crs_val==NULL){
	  *crs_val= (double*)malloc(it_matrix->Nnz*sizeof(double));
	  if( *crs_val==NULL ) {
		printf("Not enough memory for allocating crs structure for val vector\n");
		exit(-1);
	  }
	}


	//insert value to new crs structure

	#pragma omp parallel for default(none) private(i) shared(crs_row,it_matrix,offset)
	for(i=0;i<=it_matrix->Nrdofgl;i++)(*crs_row)[i]=it_matrix->crs_row_ptr[i]+offset;
	#pragma omp parallel for default(none) private(i) shared(crs_col,it_matrix,offset)
	for(i=0;i<it_matrix->Nnz;i++)(*crs_col)[i]=it_matrix->crs_col_ind[i]+offset;
	#pragma omp parallel for default(none) private(i) shared(crs_val,it_matrix,offset)
	for(i=0;i<it_matrix->Nnz;i++)(*crs_val)[i]=it_matrix->crs_val[i];

	info=1;
  }
  return(info);
}


/*---------------------------------------------------------
  lar_get_SM_and_LV_coo_from_crs_generic - convert crs matrix format to coo for external solver
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_coo_from_crs_generic( // returns: flag - 0 if not needed delete coo matrices
  int Matrix_id, /* in: matrix ID */
  int* nnz,		/* out: Number of nonzero elements in matrix */
  int offset, /* in: offset in coo_row and coo_col matrix */
  int** coo_row,	   /* out: matrix of rows in coo */
  int** coo_col,	   /* out: matrix of column in coo */
  double** coo_val,	   /* out: matrix of value in coo */
  double** rhs	   /* out: rhs vector */
			   )
{

  itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];
  *rhs=it_matrix->rhs;

  int info;
  int i,j;

  if(offset==0){
	*coo_col=it_matrix->crs_col_ind;
	info=0;
  }else{
	if(*coo_col==NULL){
	  *coo_col = (int*)malloc(it_matrix->Nnz*sizeof(int));
	  if( *coo_col==NULL ) {
		printf("Not enough memory for allocating crs structure for col_ind vector\n");
		exit(-1);
	  }
	}
	#pragma omp parallel for default(none) private(i) shared(coo_col,it_matrix,offset)
	for(i=0;i<it_matrix->Nnz;i++)(*coo_col)[i]=it_matrix->crs_col_ind[i]+offset;

	info=1;
  }

  if(*coo_row==NULL){
	*coo_row = (int*)malloc(it_matrix->Nnz*sizeof(int));
	if( *coo_row==NULL ) {
	  printf("Not enough memory for allocating crs structure for row_ptr vector\n");
	  exit(-1);
	}
  }
  if(*coo_val==NULL){
	*coo_val= (double*)malloc(it_matrix->Nnz*sizeof(double));
	if( *coo_val==NULL ) {
	  printf("Not enough memory for allocating crs structure for val vector\n");
	  exit(-1);
	}
  }



  #pragma omp parallel for default(none) private(i,j) shared(coo_row,it_matrix,offset)
  for(i=0;i<it_matrix->Nrdofgl;i++){
	  for(j=it_matrix->crs_row_ptr[i];j<it_matrix->crs_row_ptr[i+1];j++)
		(*coo_row)[j]=i+offset;
  }
  #pragma omp parallel for default(none) private(i) shared(coo_val,it_matrix,offset)
  for(i=0;i<it_matrix->Nnz;i++)(*coo_val)[i]=it_matrix->crs_val[i];

  *nnz=it_matrix->Nnz;

  //printf("\n\n    !!!!     Test %d\n",*nnz);
  return info;
}
