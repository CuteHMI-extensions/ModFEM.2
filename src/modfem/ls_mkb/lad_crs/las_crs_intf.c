
/************************************************************************
Contains CRS implementations of procedures:

  lar_allocate_SM_and_LV_crs - to allocate space for stiffness matrix and load vector
  lar_initialize_SM_and_LV_crs - to initialize stiffness matrix and/or load vector
  lar_get_storage_crs - to compute storage of SM, LV and preconditioner
  lar_fill_assembly_table_int_ent_crs - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
  lar_assemble_SM_and_LV_with_table_crs - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
  lar_assemble_SM_and_LV_crs - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  lar_allocate_preconditioner_crs - to allocate space for preconditioner
  lar_fill_preconditioner_crs - to fill preconditioner
  lar_free_preconditioner_crs - to free space of preconditioner structure
  lar_free_SM_and_LV_crs - to free space of matrix structure

  lar_create_solver_structures_accel_crs - utility to create data structures on GPU
  lar_get_crs_data_crs - utility to get CRS parameters and pointers

  lar_compute_residual_crs - to compute the residual of the not preconditioned
	system of equations, v = ( b - Ax )
  lar_perform_BJ_or_GS_iterations_crs - to perform one iteration of block Gauss-
	Seidel or block Jacobi algorithm, v_out = v_in + M^-1 * ( b - A * v_in )
  lar_perform_rhsub_crs - to perform forward reduction and back-substitution for ILU
		   preconditioning

TODO
  lar_compute_preconditioned_residual_crs - to compute the residual of the
	preconditioned system of equations, v = M^-1 * ( b - Ax )
		where M^-1 corresponds directly to the stored preconditioner matrix
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
#include "./lah_crs.h"

#include <modfem/lin_alg_intf.h>

#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)

  #include "tmh_intf.h"

  // interface of opencl implementation
  #include "../../tmd_opencl/tmh_ocl.h"

#endif

/* GLOBAL VARIABLES */
int   itv_nr_crs_matrices=0;  /* the number of matrices managed by module */
itt_crs_matrices itv_crs_matrices[LAC_MAX_MATRICES];  /* array of CRS matrices */


/*---------------------------------------------------------
  lar_allocate_SM_and_LV_crs - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
int lar_allocate_SM_and_LV_crs( // returns: matrix index in itv_crs_matrices array
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]==0 - iblock is a ghost block with no rows associated with
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
	)
{
  //printf("\n\n CRS-allocate \n\n\n");

/* auxiliary variable */
  int i, j, iblock;

/*++++++++++++++++ executable statements ++++++++++++++++*/


// create structures for a new matrix
  itt_crs_matrices *it_matrix;
  itv_nr_crs_matrices++;
  int inter_Nrdof_glob=0;

  if(itv_nr_crs_matrices>=LAC_MAX_MATRICES){
	printf("Too much (%d) matrices requested in las_block_intf! (correct lah_block.h)\n",
	   itv_nr_crs_matrices);
  }

  it_matrix = &itv_crs_matrices[itv_nr_crs_matrices-1];

  it_matrix->Max_SM_size = Max_SM_size;

/*kbw
  printf("\nEntering lar_allocate_SM_and_LV_crs: max_sm_size %d, nrdof_glob %d\n",
	 Max_SM_size, Nrdof_glob);
  printf("\n");
// !!! offset 1 numbering of blocks moved to lad_block !!!
  printf("block_id,\tnroffbl\t\tneighbors\n");
  for(i=0; i<Nrdof_glob; i++){
	printf("%d\t\t%d\t\t%d\t",i, Nroffbl[i]);
	for(j=0;j<Nroffbl[i]; j++){
	  printf("%10d",L_offbl[i][j]);
	}
	printf("\n");
  }
/*kew*/

  //allocate crs structure

	int nnz=0;
  for(iblock=0;iblock<Nrdof_glob;iblock++){
	if(Nroffbl[iblock]>0){ // ghost blocks have Nroffbl[index]<=0
	  nnz+=(Nroffbl[iblock]+1);
	  ++inter_Nrdof_glob;
	  //printf("iblock=%d, Nroffbl %d\n",iblock,Nroffbl[iblock]);
	}
  }
  it_matrix->Nnz = nnz;
  it_matrix->Nrdof_internal = inter_Nrdof_glob;
  it_matrix->Nrdofgl=Nrdof_glob;
  //printf("\n\n nnz=%d inter_Nrdof_glob=%d Nrdof_glob=%d \n",nnz,inter_Nrdof_glob,Nrdof_glob);

#ifdef PARALLEL
//  printf("\n\n ### pcv_my_proc_id=%d nnz=%d Nrdof_glob=%d inter_Nrdof_glob=%d\n",pcv_my_proc_id,nnz,Nrdof_glob,inter_Nrdof_glob);
#endif

  it_matrix->crs_row_ptr = (int*)malloc( (Nrdof_glob+1)*sizeof(int));
  if( it_matrix->crs_row_ptr==NULL ) {
	printf("Not enough memory for allocating crs structure for row_ptr vector\n");
	exit(-1);
  }

  it_matrix->crs_col_ind = (int*)malloc((nnz)*sizeof(int));
  if( it_matrix->crs_col_ind==NULL ) {
	printf("Not enough memory for allocating crs structure for col_ind vector\n");
	exit(-1);
  }

  it_matrix->crs_val= (double*)malloc((nnz)*sizeof(double));
  if( it_matrix->crs_val==NULL ) {
	printf("Not enough memory for allocating crs structure for val vector\n");
	exit(-1);
  }


  it_matrix->rhs = (double*)malloc( (Nrdof_glob)*sizeof(double));
  if( it_matrix->rhs==NULL ) {
	printf("Not enough memory for allocating crs structure for rhs vector\n");
	exit(-1);
  }


  //filling crs_row and crs_col structure
  // loops are parallelized to spread matrices in memory for NUMA architectures
#pragma omp parallel for default(none) firstprivate(it_matrix, Nrdof_glob)
  for(iblock=0;iblock<Nrdof_glob;iblock++)it_matrix->rhs[iblock]=0.0;
#pragma omp parallel for default(none) firstprivate(it_matrix, nnz)
  for(iblock=0;iblock<nnz;iblock++)it_matrix->crs_val[iblock]=0.0;


  int bandwidth = 0;

  int neigi;

/*kbw
  for(iblock=0;iblock<my_inter_Nrdof_glob;iblock++){
	for(neigi=0;neigi<Nroffbl[iblock];neigi++){
	  printf("iblock=%d, %d\n",iblock,L_offbl[iblock][neigi]);
	}
  }
//*kew*/


  //v3
  int val_index = 0, col_block_index=0, block_ptr_index=0;

  it_matrix->crs_row_ptr[block_ptr_index]=col_block_index;
  block_ptr_index++;

  for(iblock=0;iblock<inter_Nrdof_glob;iblock++){

	// version for diagonal entries first in a row
	//it_matrix->crs_col_ind[col_block_index]=iblock;
	//++col_block_index;

	int inserted=0;
	for(neigi=0;neigi<Nroffbl[iblock];neigi++){

	  // checking assumptions
	  // list of neighbors is already sorted (sis_mkb_intf.c, line circa 1111)
	  if(neigi>0){
	//assert(L_offbl[iblock][neigi-1]<L_offbl[iblock][neigi]);
	if(L_offbl[iblock][neigi-1]>=L_offbl[iblock][neigi]){
	  printf("not sorted list of neighbors in lar_allocate_SM_and_LV_crs\n");
	  exit(-1);
	}
	  }

	  if(abs(L_offbl[iblock][neigi] - iblock)>bandwidth) {
		bandwidth=abs(L_offbl[iblock][neigi] - iblock);
	  }

	  // version for diagonal entries first in a row
	  //it_matrix->crs_col_ind[ col_block_index] = L_offbl[iblock][neigi];
	  //++col_block_index;

	  // version for diagonal entries inserted in the middle
	  if(inserted==0){

		if(L_offbl[iblock][neigi] < iblock){
		  it_matrix->crs_col_ind[ col_block_index] = L_offbl[iblock][neigi];
		  col_block_index++;
		}
		else{

		  it_matrix->crs_col_ind[col_block_index]=iblock;
		  col_block_index++;

		  it_matrix->crs_col_ind[ col_block_index] = L_offbl[iblock][neigi];
		  col_block_index++;

		  inserted = 1;

		}

	  }
	  else{

		it_matrix->crs_col_ind[ col_block_index] = L_offbl[iblock][neigi];
		col_block_index++;

	  }

	} // end loop over neighbor blocks

	// if diagonal was the last
	if(inserted==0){
	  it_matrix->crs_col_ind[col_block_index]=iblock;
	  col_block_index++;
	}

	it_matrix->crs_row_ptr[block_ptr_index]=col_block_index;
	block_ptr_index++;


  }

  //last diagonal entries
  //it_matrix->crs_col_ind[col_block_index]=iblock-1;
  //col_block_index++;
  //printf("\n nnz=%d col_block_index=%d it_matrix->Nrdofgl=%d it_matrix->Nnz=%d inter_Nrdof_glob=%d Nrdof_glob=%d\n\n",nnz,col_block_index,it_matrix->Nrdofgl,it_matrix->Nnz,inter_Nrdof_glob,Nrdof_glob);
  assert(nnz==col_block_index);
  it_matrix->crs_row_ptr[inter_Nrdof_glob]=nnz;

  it_matrix->Half_bandwidth = bandwidth;

/*ok_kbw*/
  printf("\nAllocated CRS matrix %d: nglob = %d, nloc=%d nnz = %d, half bandwidth = %d\n",
	 itv_nr_crs_matrices-1, it_matrix->Nrdofgl, it_matrix->Nrdof_internal, it_matrix->Nnz, it_matrix->Half_bandwidth);
/*kew*/

  //printf("\n\nblock_ptr_index= %d\n\n",block_ptr_index);

  //it_matrix->crs_row_ptr[block_ptr_index]=col_block_index+1;

  //it_matrix->crs_col_ind[col_block_index]=0;


  /*kbw
	for(iblock=0,i=0;iblock<nnz;++i){
	for(j=it_matrix->crs_row_ptr[i];j<it_matrix->crs_row_ptr[i+1];j++,iblock++)
	printf("(%d %d)",j,it_matrix->crs_col_ind[j]);
	//printf("(%d %d %lf)",j,it_matrix->crs_col_ind[j],it_matrix->crs_val[j]);
	//printf("(%d %d)",i,iblock);
	printf("\n");
	}
	/**/


  return(itv_nr_crs_matrices-1);
}

/**---------------------------------------------------------
  lar_create_solver_structures_accel_crs - utility to create data structures on GPU
---------------------------------------------------------*/
int lar_create_solver_structures_accel_crs(
  int Matrix_id   /* in: matrix ID */
)
{
  int info=1;
  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];

  //  pass created arrays to tmr routine ???
  /* tmr_ocl_create_solver_structures_crs(  ??? */
  /* tmr_ocl_create_solver_structures_accel_crs(  */
  /* 					 int Nrdof_glob, */
  /* 					 int nnz, */
  /* 					 int* crs_col_ind, */
  /* 					 int* crs_row_ptr, */
  /* 					 double* crs_val, */
  /* 					 double* rhs */
  /* 					 itd, itp ); */



  // call OpenCL routine
#ifdef OPENCL_GPU
  info = tmr_ocl_create_solver_structures_crs(
					 it_matrix->Nrdofgl,
					 it_matrix->Nnz,
					 it_matrix->crs_col_ind,
					 it_matrix->crs_row_ptr,
					 it_matrix->crs_val,
					 it_matrix->rhs
					 );

#endif

  return(info);
}

/**---------------------------------------------------------
lar_get_crs_data_crs - utility to get CRS parameters and pointers
---------------------------------------------------------*/
int lar_get_crs_data_crs( /* returns: >0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int* Nrdof_glob_p, /* output pointers: nrdof, nnz etc. - self explanatory */
  int* Nnz_p,
  int** Crs_col_ind_p,
  int** Crs_row_ptr_p,
  double** Crs_val_p,
  double** Rhs_p
			  )
{

  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];

  *Nrdof_glob_p = it_matrix->Nrdofgl;
  *Nnz_p = it_matrix->Nnz;
  *Crs_col_ind_p = it_matrix->crs_col_ind;
  *Crs_row_ptr_p = it_matrix->crs_row_ptr;
  *Crs_val_p = it_matrix->crs_val;
  *Rhs_p = it_matrix->rhs;

  return(0);
}

/*---------------------------------------------------------
  lar_initialize_SM_and_LV_crs - to initialize stiffness matrix and/or load vector
---------------------------------------------------------*/
int lar_initialize_SM_and_LV_crs(
  int Matrix_id,   /* in: matrix ID */
  int Scope    /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  ){

	itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];

	int nrdof_glob = it_matrix->Nrdofgl;
	int nrdof_internal =it_matrix->Nrdof_internal;

	// zero SM and LV!!!!!!!!!!!!!!!!!!!!
	int iblock;

	// parallel loop is the same as in matrix-vector product to properly assign
	// parts of it_matrix->crs_val to threads (for NUMA architectures)
#pragma omp parallel for default(none) firstprivate(it_matrix, nrdof_glob)
	for(iblock=0;iblock<nrdof_glob;iblock++) {
	  it_matrix->rhs[iblock]=0.0;
	}

	if(Scope==LAC_SCOPE_SM_AND_LV){
#pragma omp parallel for default(none) firstprivate(it_matrix, nrdof_glob,nrdof_internal)
	  for(iblock=0;iblock<nrdof_internal;iblock++) {
	int icrs;
	for(icrs=it_matrix->crs_row_ptr[iblock];
		icrs<it_matrix->crs_row_ptr[iblock+1];
		icrs++) {
	  it_matrix->crs_val[icrs]=0.0;
	}
	  }
	}

/*ok_kbw*/
  printf("\n\nInitialized CRS matrix %d: n = %d, nnz = %d, half bandwidth = %d\n",
	 Matrix_id, it_matrix->Nrdofgl, it_matrix->Nnz, it_matrix->Half_bandwidth);
  printf("The average number of non-zeros in a single row: %d\n",
	 it_matrix->Nnz/it_matrix->Nrdofgl);
/*kew*/

	return 0;
  }

/*---------------------------------------------------------
  lar_get_storage_crs - to compute storage of SM, LV and preconditioner
---------------------------------------------------------*/
double lar_get_storage_crs( /* returns: storage in MB */
  int Matrix_id   /* in: matrix ID */
			   )
{
	return -1;
}


/**-----------------------------------------------------------
  lar_fill_assembly_table_int_ent_crs - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
int lar_fill_assembly_table_int_ent_crs(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int* Assembly_table_int_ent /* part of the global assembly table */
  )
{

/* get pointer to the level structure */
  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];

/*kbw
  printf("In lar_crs_fill_assembly_table_int_ent: matrix_id %d\n",
	 Matrix_id);
  printf("nr_dof_ent_loc %d, position %u\n",
	 Nr_dof_bl, Assembly_table_int_ent);
  int ibl;
  for(ibl=0;ibl<Nr_dof_bl; ibl++) printf("bl_id %d ", L_bl_id[ibl]);
  printf("\n");
/*kew*/

  int idofbl;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {
	int bl_id_y=L_bl_id[idofbl];

	int jdofbl;
	for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) {
	  int bl_id_x=L_bl_id[jdofbl];

	  int jaux = jdofbl+idofbl*Nr_dof_bl;

/*kbw
	  printf("entry %d (idofbl %d, bl_id_y %d, jdofbl %d bl_id_x %d)\n",
		 jaux, idofbl, bl_id_y, jdofbl, bl_id_x);
/*kew*/



#ifdef DEBUG_LSM
	  int success=0;
#endif


#ifdef PARALLEL
  if(bl_id_x<it_matrix->Nrdof_internal){
#endif


	  int icrs;
	  for(icrs=it_matrix->crs_row_ptr[bl_id_x];
	  icrs<it_matrix->crs_row_ptr[bl_id_x+1];
	  icrs++){

	if(it_matrix->crs_col_ind[icrs]==bl_id_y){

/*kbw
	  printf("found for entry %d position in CRS %d\n", jaux, icrs);
/*kew*/
#ifdef DEBUG_LSM
	  success=1;
#endif

	  Assembly_table_int_ent[jaux]=icrs;

	  break;
		}
	  }

#ifdef PARALLEL
  }else
	Assembly_table_int_ent[jaux]=-1;
#endif


#ifdef DEBUG_LSM
	  if(success==0){
	printf("NOT found for entry %d (local %d, %d; global %d, %d) position in CRS\n",
		   jaux, idofbl, jdofbl, bl_id_y, bl_id_x);
	getchar(); getchar();
	  }
#endif

	}
  }

/*<<<<<<<<<< RHS part of Assembly_table exchanged for local_to_global array >>>>>>>>>>>>*/
/*   int jdofbl; */
/*   for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) { */
/*     int bl_id_x=L_bl_id[jdofbl]; */

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
  lar_assemble_SM_and_LV_with_table_crs - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_with_table_crs(
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
  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];

  //printf("\n\n\n\n\nAssembly_table_int_ent=%d ->",Assembly_table_int_ent[0]);

  int idofbl;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {

	int jdofbl;
	for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) {

	  int jaux = jdofbl+idofbl*Nr_dof_bl;

	  int icrs = Assembly_table_int_ent[jaux];
#ifdef PARALLEL
  if(icrs!=-1){
#endif


/*kbw
		  if(icrs==0){
		  printf("found for entry %d position in CRS %d Nr_dof_bl=%d\n", jaux, icrs,Nr_dof_bl);
		  printf("filling with %lf: before %lf\n",
			  Stiff_mat[jaux], it_matrix->crs_val[icrs]);getchar();}

/*kew*/
	  it_matrix->crs_val[icrs] += Stiff_mat[jaux];
/*kbw
		  printf("filling with %lf: after %lf\n",
			 Stiff_mat[jaux], it_matrix->crs_val[icrs]);
/*kew*/


#ifdef PARALLEL
  }
#endif


	}


  }

/*<<<<<<<<<< RHS part of Assembly_table exchanged for local_to_global array >>>>>>>>>>>>*/
//  int jaux = Nr_dof_bl*Nr_dof_bl;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {

	//int irhs = Assembly_table_int_ent[jaux+idofbl];
	int irhs = Local_to_global_int_ent[idofbl];
/*kbw
	  printf("filling entry %d for RHS block %d\n",
		 irhs, idofbl);
	  printf("filling with %lf: before %lf\n",
		 Rhs_vect[idofbl], it_matrix->rhs[irhs]);
/*kew*/
	it_matrix->rhs[irhs] += Rhs_vect[idofbl];

/*kbw
	  printf("filling with %lf: before %lf\n",
		 Rhs_vect[idofbl], it_matrix->rhs[irhs]);
/*kew*/
  }

/*kbw
  getchar();
/*kew*/

  return(1);
}

/*------------------------------------------------------------
  lar_assemble_SM_and_LV_crs - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_crs(
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
  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];

   //printf("\n\n\n\n\nasembling %d !!!!!!!!!!! \n\n",it_matrix->crs_row_ptr[1]);

  /*kbw
  printf("In lar_crs_assembly: matrix_id %d, nr_dof_ent_loc %d\n",
	 Matrix_id, Nr_dof_bl);
  for(ibl=0;ibl<Nr_dof_bl; ibl++) printf("bl_id %d \n",  L_bl_id[ibl]);
  printf("\n");
/*kew*/

   // //******************CRS************************************//
// /* compute local stiffness matrix nrdof */
  nrdof=0;
  for(iblock=0;iblock<Nr_dof_bl;iblock++){
	nrdof += L_bl_nrdof[iblock];
  }


  posloc_i=0;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {
	bl_id_y=L_bl_id[idofbl];
	nrdof_i = L_bl_nrdof[idofbl];

	posloc_j=0;
	for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) {
	  bl_id_x=L_bl_id[jdofbl];
	  nrdof_j = L_bl_nrdof[jdofbl];

	  for(idof=0;idof<nrdof_i;idof++){

		jaux = posloc_j+(posloc_i+idof)*nrdof;
		for(jdof=0;jdof<nrdof_j;jdof++) {

/*kbw
	  printf("entry %d (idofbl %d, bl_id_y %d, jdofbl %d bl_id_x %d)\n",
		 jaux, idofbl, bl_id_y, jdofbl, bl_id_x);
/*kew*/

		  jcrs=bl_id_x+jdof;

#ifdef PARALLEL
  if(jcrs<it_matrix->Nrdof_internal)
#endif
		  for(icrs=it_matrix->crs_row_ptr[jcrs];icrs<it_matrix->crs_row_ptr[jcrs+1];icrs++){
			//if(icrs>=it_matrix->Nnz){
		//  printf("\n\n icrs=%d > %d=Nnz",icrs,it_matrix->Nnz);
		//  exit(0);
		//}
			if((it_matrix->crs_col_ind[icrs])==(bl_id_y+idof)){

/*kbw
if(icrs==0){
		  printf("found for entry %d position in CRS %d\n", jaux, icrs);
		  printf("filling with %lf: before %lf\n",
			 Stiff_mat[jaux+jdof], it_matrix->crs_val[icrs]);}
/*kew*/

		  it_matrix->crs_val[icrs] += Stiff_mat[jaux+jdof];

/*kbw
if(icrs==0){
		  printf("filling with %lf: after %lf\n",
			 Stiff_mat[jaux+jdof], it_matrix->crs_val[icrs]);}
/*kew*/

		  break;
		}

		  }
		}
	  }
		posloc_j += nrdof_j;
	}

	for (idof=0;idof<nrdof_i;idof++) {

/*kbw
	  printf("filling entry %d for RHS block %d (block %d)\n",
		 bl_id_y+idof, posloc_i+idof);
	  printf("filling with %lf: before %lf\n",
		 Rhs_vect[posloc_i+idof], it_matrix->rhs[bl_id_y+idof]);
/*kew*/

	  /* assemble right hand side block's entries */
	  it_matrix->rhs[bl_id_y+idof] += Rhs_vect[posloc_i+idof];

/*kbw
	  printf("filling with %lf: after %lf\n",
		 Rhs_vect[posloc_i+idof], it_matrix->rhs[bl_id_y+idof]);
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



  /*---------------------------------------------------------
lar_allocate_preconditioner_crs - to create preconditioner blocks corresponding
					   to small subdomains of neighboring elements
---------------------------------------------------------*/
int lar_allocate_preconditioner_crs( /* returns:   >0 number of diagonal blocks */
						  /*	       <=0 - error */
  int Matrix_id,   /* in: matrix ID */
  int Precon,      /* in: type of preconditioner (lah_block.h line circa 45) */
  int ILU_k // in: for ILU(k) - k;
  )
{
  int i,Nrdofgl,kk;
  int block_row_gl, sm_block, block_col_gl;

  itt_crs_matrices *it_matrix;
  it_matrix = &itv_crs_matrices[Matrix_id];
  it_matrix->Precon = Precon;

  Nrdofgl = it_matrix->Nrdofgl;
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

/*---------------------------------------------------------
  lar_fill_preconditioner_crs - to fill preconditioner (here by factorizing
							diagonal blocks or by block ILU(0) factorization)
---------------------------------------------------------*/
int lar_fill_preconditioner_crs(
  int Matrix_id   /* in: matrix ID */
	)
{
  int i,Nrdofgl,diag_index;
  double *pivots;


  itt_crs_matrices *it_matrix;
  it_matrix = &itv_crs_matrices[Matrix_id];
  Nrdofgl = it_matrix->Nrdofgl;
  int nrdof_internal =it_matrix->Nrdof_internal;

  //printf("\nCRS_lar_fill_preconditioner\n");

  for(i=0;i<nrdof_internal;++i){
	diag_index=it_matrix->diag_ptr[i];
	it_matrix->diag_precon[i]=1.0/it_matrix->crs_val[diag_index];

  }


  return 0;
}

/*---------------------------------------------------------
  lar_free_preconditioner_crs - to free space for a block structure
---------------------------------------------------------*/
int lar_free_preconditioner_crs(
  int Matrix_id   /* in: matrix ID */
  )
{
  itt_crs_matrices *it_matrix;
  it_matrix = &itv_crs_matrices[Matrix_id];

  free(it_matrix->diag_ptr);
  free(it_matrix->diag_precon);

  return 0;
}


  /*---------------------------------------------------------
  lar_free_SM_and_LV_crs - to free space for a block structure
---------------------------------------------------------*/
int lar_free_SM_and_LV_crs(
  int Matrix_id   /* in: matrix ID */
  )
{
  // to keep matrix management as simple as possible it is required that
  // the matrices are freed in the reverse order wrt creating (as in heap)
  if(Matrix_id!=itv_nr_crs_matrices-1){
	printf("Matrix destroyed %d (from %d) is not the last in lar_free_SM_and_LV_crs!!!\n",
	   Matrix_id, itv_nr_crs_matrices);
	exit(-1);
  }

//printf("\n\n CRS-free \n\n");
  /* pointer to solver structure */
	itt_crs_matrices *it_matrix;
	it_matrix = &itv_crs_matrices[Matrix_id];
	  free(it_matrix->crs_val);
	  free(it_matrix->crs_col_ind);
	  free(it_matrix->crs_row_ptr);
	free(it_matrix->rhs);
	itv_nr_crs_matrices--;
  //printf("\n\n CRS-free-end \n\n");
  return(0);
}

/*---------------------------------------------------------
lar_compute_residual_crs - to compute the residual of the system of equations,
	v = ( b - Ax ) (used also to compute the product v = -Ax)
---------------------------------------------------------*/
void lar_compute_residual_crs(
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

/* constants */
  int ione = 1;
  double done = 1.;

  int i,j,k;

  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];
  int nrdof_internal =it_matrix->Nrdof_internal;
  assert(Ndof==it_matrix->Nrdofgl);

  //printf("\nCRS_lar_compute_residual\n");

	/*kbw
  //int i,j;
  int iblock;
   for(iblock=0,i=0;i<Ndof;++i)printf("%lf ",it_matrix->rhs[i]);
   printf("\n\n");


	for(iblock=0,i=0;iblock<it_matrix->Nnz;++i){
	for(j=it_matrix->crs_row_ptr[i];j<it_matrix->crs_row_ptr[i+1];j++,iblock++)
	//printf("(%d %d)",j,it_matrix->crs_col_ind[j]);
	printf("(%d %d %lf)",j,it_matrix->crs_col_ind[j],it_matrix->crs_val[j]);
	//printf("(%d %d)",i,iblock);
	printf("\n");
	}
	/**/


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
  //if(k>nrdof_internal)printf("\n\nOKI %d >%d SM=%lf LV=%lf\n\n",k,nrdof_internal,it_matrix->crs_val[k],X[it_matrix->crs_col_ind[k]]);
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


/*---------------------------------------------------------
lar_perform_BJ_or_GS_iterations_crs - to perform one iteration of block Gauss-Seidel algorithm
	(block Jacobi switched on by it_matrix->Precon==BLOCK_JACOBI)
		v_out = v_in + M^-1 * ( b - A * v_in )
---------------------------------------------------------*/
void lar_perform_BJ_or_GS_iterations_crs(
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
  itt_crs_matrices *it_matrix;
  it_matrix = &itv_crs_matrices[Matrix_id];

  int nrdof_internal =it_matrix->Nrdof_internal;

  int i_loop;		/* counter for loops for GS */
  int ione = 1;
  double done = 1.;
  int i,ret;

  int block_row_gl, sm_block, block_col_gl,block_val_ind;

  //printf("\nCRS_lar_perform_BJ_or_GS_iterations\n");


  /*kcw
  FILE *test_out;
  test_out = fopen("cpu_elems_crs.txt","w+");
  fprintf(test_out,"it_matrix.crs_val: \n");
  for(i=0;i<it_matrix->Nnz;++i){
	fprintf(test_out,"(%d: %lf) \n; ",i,it_matrix->crs_val[i]);
  }
  fclose(test_out);

  test_out = fopen("cpu_rhs.txt","w+");
  fprintf(test_out,"it_matrix.rhs: \n");
  for(i=0;i<Ndof;++i){
	fprintf(test_out,"(%d: %lf); \n ",i,it_matrix->rhs[i]);
  }
  fclose(test_out);

  //exit(0);
  /*kce*/

  /* loop over loops */
  for(i_loop=0;i_loop<Nr_prec;i_loop++){

	/* loop over all blocks */
#pragma omp parallel default(none) firstprivate(V,B,it_matrix,i_loop,Use_rhs,Ini_zero,Ndof,ione,done,ret,nrdof_internal) \
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
	  //for(i=my_first_x; i<my_last_x; i++) vtemp_GS_BJ[i]=V[i];
	  for(i=0;i<Ndof;i++) vtemp_GS_BJ[i]=V[i];

	  //#pragma omp barrier

	  int my_first = my_id*portion;
	  if(my_first<0) my_first = 0;
	  int my_last = (my_id+1)*portion;
	  if(my_last>nrdof_internal) my_last=nrdof_internal;

	  //#pragma omp for
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
	  //for (block_row_gl = 0;  block_row_gl < it_matrix->Nrdofgl; block_row_gl++){
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

/*---------------------------------------------------------
lar_perform_rhsub_crs - to perform forward reduction and back-substitution for ILU
		   preconditioning: v_out = M^-1 * b
---------------------------------------------------------*/
void lar_perform_rhsub_crs(
  int Matrix_id,   /* in: matrix ID */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* out: vector of unknowns updated */
		/*      during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	)
{
  int i,j,Nrdofgl;
  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];
  Nrdofgl = it_matrix->Nrdofgl;
  int nrdof_internal =it_matrix->Nrdof_internal;

  int  *z;
  double sum;

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
  lar_get_SM_and_LV_crs_from_crs - convert crs matrix format to crs
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_crs_from_crs( // returns: flag - 0 if not needed delete crs matrices
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

  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];
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
  lar_get_SM_and_LV_coo_from_crs - convert crs matrix format to coo for external solver
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_coo_from_crs( // returns: flag - 0 if not needed delete coo matrices
  int Matrix_id, /* in: matrix ID */
  int* nnz,		/* out: Number of nonzero elements in matrix */
  int offset, /* in: offset in coo_row and coo_col matrix */
  int** coo_row,	   /* out: matrix of rows in coo */
  int** coo_col,	   /* out: matrix of column in coo */
  double** coo_val,	   /* out: matrix of value in coo */
  double** rhs	   /* out: rhs vector */
			   )
{

  itt_crs_matrices *it_matrix = &itv_crs_matrices[Matrix_id];
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

//for(i=0;i<50;i++)printf("(%lf,%d,%d)\n",(*coo_val)[i],(*coo_row)[i],(*coo_col)[i]);

/*
  FILE *fp=fopen("mumps.txt", "w");
	int k=1;
	for(i=0;i<*nnz;i++){
		if(k!=(*coo_row)[i])fprintf(fp,"\n");
		fprintf(fp,"(%d,%d,%lf) ",(*coo_row)[i]-1,(*coo_col)[i]-1,(*coo_val)[i]);
		k= (*coo_row)[i];

	}exit(0);
*/

  return info;
}


