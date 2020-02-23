/*************************************************************
File contains BCRS implementations of  procedures:
  lar_allocate_SM_and_LV_bcrs - to allocate space for stiffness matrix and load vector
  lar_initialize_SM_and_LV_bcrs - to initialize stiffness matrix and/or load vector
  lar_get_storage_bcrs - to compute storage of SM, LV and preconditioner
  lar_fill_assembly_table_int_ent_bcrs - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
  lar_assemble_SM_and_LV_with_table_bcrs - to assemble entries to the global stiffness
						   matrix and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
  lar_assemble_SM_and_LV_bcrs - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  lar_allocate_preconditioner_bcrs - to allocate space for preconditioner
  lar_fill_preconditioner_bcrs - to fill preconditioner
  lar_free_preconditioner_bcrs - to free space for a bcrs structure
  lar_free_SM_and_LV_bcrs - to free space for a bcrs structure

lar_compute_residual_bcrs - to compute the residual of the not preconditioned
	system of equations, v = ( b - Ax )
lar_perform_BJ_or_GS_iterations_bcrs - to perform one iteration of block Gauss-Seidel
	or block Jacobi algorithm
	 v_out = v_in + M^-1 * ( b - A * v_in )
lar_perform_rhsub_bcrs - to perform forward reduction and back-substitution for ILU
		   preconditioning

------------------------------
History:
	02.2016 - Krzysztof Banas, initial version
	04.2016 - Kazimierz Chlon
	06.2016 - Krzysztof Banas, ILU(k) initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<assert.h>

#include<omp.h>

/* internal information for the BCRS module */
#include "./lah_bcrs.h"

// interface for linear algebra packages supporting MKB solver
#include <modfem/ls_mkb/lah_intf.h>

// for debugging mpi
#include <modfem/pch_intf.h>

/* LAPACK and BLAS procedures */
#include <modfem/lin_alg_intf.h>

#define TIME_TEST_ILU
#include <modfem/uth_system.h>

/* minimal number of dofs to make use of DGEMV reasonable */
#define MIN_DOF_DGEMV 10000
#define SMALL 1e-15 /* approximation of round-off error */
#define TOL 1e-9    /* default tolerance for different iterative processes */

#define MAX_BL_BL 20


/* GLOBAL VARIABLES */
int   itv_nr_bcrs_matrices=0;       /* the number of bcrs matrices */
itt_bcrs_matrices itv_bcrs_matrices[LAC_MAX_MATRICES];  /* array of matrices */


/************************ local utilities *****************/
#define LAC_END_OF_LIST -7
#define LAC_PUT_LIST_NOT_FOUND -6
#define LAC_PUT_LIST_FOUND -5
#define LAC_PUT_LIST_FULL -4

/*---------------------------------------------------------
lar_util_chk_list_bcrs - to check whether Num is on the list List
	with length Ll (filled with numbers and LAC_END_OF_LIST at the end)
---------------------------------------------------------*/
int lar_util_chk_list_bcrs(	/* returns: */
			/* >=0 - position on the list */
					/* LAC_PUT_LIST_NOT_FOUND - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
				);


/*---------------------------------------------------------
lar_util_put_list_bcrs - to put Num on the list List (with length Ll and filled with
					 meaningfull numbers and one or more LAC_END_OF_LIST at the end)
---------------------------------------------------------*/
int lar_util_put_list_bcrs( /* returns*/
		/*  >=0 - success, position on the list */
				/*  LAC_PUT_LIST_FOUND - not put - already on the list */
				/*  LAC_PUT_LIST_FULL - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
				);

int lar_sort_short(
		   int* A, // array
		   int p,  //first index
		   int k   // last index
		   );

//////////// DEFINITIONS OF PROCEDURES ///////////////

/*---------------------------------------------------------
  lar_allocate_SM_and_LV_bcrs - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
int lar_allocate_SM_and_LV_bcrs( // returns: matrix index in itv_matrices array
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Block_size,  /* in: size of SM blocks (-1 - non-constant size) */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]<=0 - iblock is a ghost block with no rows associated
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
	)
{

  printf("\nBCRS storage - lar_allocate_SM_and_LV_bcrs\n");
  assert(Block_size==BLS);

  int Nrblocks = Nrdof_glob / BLS;
  // Nrblocks - number of DOF blocks !!! (i.e. LINEAR blocks in the vector of unknowns)
  // DOF blocks correspond to square blocks in stiffness matrices (local and the global)
  // i.e. a single square block corresponds to one row DOF block (indexed by block_row_gl)
  // and one column DOF block (indexed by block_col_gl)
  // Indices block_row_gl, block_col_gl correspond to standard row_gl, col_gl CRS indices

  // assert: each Nrdofbl[iblock] = BLS - hence Nrdofbl[iblock] is not used
  // assert: each Posglob[iblock] = iblock*BLS - hence Posglob[iblock] is not used i.e.
  // DOF blocks are placed in the vector of unknowns in the same order as on the list L_offbl


  int nr_sm_blocks = 0; // number of square blocks in global SM
  int index, block_row_gl, block_ptr_index=0;
  int i,j,neigi,flag=0;

// create structures for a new matrix
  itt_bcrs_matrices *it_matrix;
  itv_nr_bcrs_matrices++;
  int inter_Nrbl=0;

  if(itv_nr_bcrs_matrices>=LAC_MAX_MATRICES){
	printf("Too much (%d) matrices requested in las_bcrs_intf!\n",
	   itv_nr_bcrs_matrices);
  }

  it_matrix = &itv_bcrs_matrices[itv_nr_bcrs_matrices-1];

  it_matrix->Max_SM_size = Max_SM_size;
  it_matrix->Nrblocks = Nrblocks;
  it_matrix->Nrdofgl = Nrdof_glob;


   //allocate bcrs structure
   //printf("   %d %d",Nrdof_glob,Nrblocks);getchar();getchar();


	int max_nrneig = 0;
  for (block_row_gl = 0;  block_row_gl < Nrblocks; block_row_gl++){

	if(Nroffbl[block_row_gl]>0){ // ghost blocks have Nroffbl[index]<=0

	  nr_sm_blocks += (Nroffbl[block_row_gl]+1); // off diagonal blocks + diagonal block
	  ++inter_Nrbl;
	}

	if(Nroffbl[block_row_gl]>max_nrneig) max_nrneig = Nroffbl[block_row_gl];

  }

  it_matrix->Nrbl_internal=inter_Nrbl;

//printf("\n\n\n\n   !!!  Nrblocks=%d  Nrdof_glob=%d  inter_Nrbl=%d  nr_sm_blocks=%d   \n\n\n\n",Nrblocks,Nrdof_glob,inter_Nrbl,nr_sm_blocks);

  it_matrix->block_row_ptr = (int*)malloc( (inter_Nrbl+1)*sizeof(int));
  if( it_matrix->block_row_ptr==NULL ) {
	printf("Not enough memory for allocating bcrs structure for row_ptr vector\n");
	exit(-1);
  }


  it_matrix->block_col_ind = (int*)malloc( (nr_sm_blocks+1)*sizeof(int));
  if( it_matrix->block_col_ind==NULL ) {
	printf("Not enough memory for allocating bcrs structure for col_ind vector\n");
	exit(-1);
  }

  it_matrix->block_val = (double*)malloc( (nr_sm_blocks*BLS*BLS)*sizeof(double));
  if( it_matrix->block_val==NULL ) {
	printf("Not enough memory for allocating bcrs structure for val vector\n");
	exit(-1);
  }

/*kbw
  printf("\n \n nr_sm_blocks %d, SM size %d \n \n",nr_sm_blocks,nr_sm_blocks*BLS*BLS);
//*kew*/

// IS THIS OPTIMAL FOR NUMA ARCHITECTURES?
  memset( it_matrix->block_val, 0, (nr_sm_blocks*BLS*BLS)*sizeof(double) );


  it_matrix->rhs = (double*)malloc( (Nrdof_glob)*sizeof(double));
  if( it_matrix->rhs==NULL ) {
	printf("Not enough memory for allocating crs structure for rhs vector\n");
	exit(-1);
  }
  memset( it_matrix->rhs, 0.0, (Nrdof_glob)*sizeof(double) );

  int bandwidth = 0;
  index = 0;
  for (block_row_gl = 0;  block_row_gl < inter_Nrbl; block_row_gl++){

	it_matrix->block_row_ptr[ block_row_gl ] = index;
	//printf("%d %d \n",block_row_gl,index);getchar();

	// all ghost blocks (identified by Nroffbl[block_row_gl] < 0) are skipped
   for( i=0,flag=0; i< Nroffbl[block_row_gl]; i++){


#ifdef DEBUG_LSM
	  // checking assumptions
	  // list of neighbors is already sorted (sis_mkb_intf.c, line circa 1111)
	  if(i>0){
	//assert(L_offbl[iblock][neigi-1]<L_offbl[iblock][neigi]);
	if(L_offbl[block_row_gl][i-1]>=L_offbl[block_row_gl][i]){
	  printf("not sorted list of neighbors in lar_allocate_SM_and_LV_crs\n");
	  exit(-1);
	}
	  }
#endif

	  if(abs(L_offbl[block_row_gl][i] - block_row_gl)>bandwidth) {
		bandwidth=abs(L_offbl[block_row_gl][i] - block_row_gl);
	  }


	  if(flag==0 && L_offbl[ block_row_gl ][i]>block_row_gl){
			it_matrix->block_col_ind[index] = block_row_gl;
			index ++;
			flag=1;
	  }
	  it_matrix->block_col_ind[index] = L_offbl[ block_row_gl ][i];
	  index ++;

	} //the end of loop over neighbors

	// if diagonal is the last entry (all neighbors have lower indices)
	if(flag==0 && Nroffbl[block_row_gl]>0 ){
	  it_matrix->block_col_ind[index] = block_row_gl;
	  index ++;
	}

  }
  it_matrix->block_row_ptr[inter_Nrbl]=index;

  assert(index==nr_sm_blocks);

  it_matrix->Nr_sm_blocks=nr_sm_blocks;
  it_matrix->Half_bandwidth = bandwidth;
  it_matrix->Max_nrneig = max_nrneig;


/*ok_kbw*/
  printf("Allocated BCRS matrix %d: nbl = %d/%d, n = %d, nr_sm_bl %d, nnz = %d, half bandwidth = %d\n",
	 itv_nr_bcrs_matrices-1,it_matrix->Nrbl_internal ,it_matrix->Nrblocks, it_matrix->Nrdofgl,
	 it_matrix->Nr_sm_blocks, it_matrix->Nr_sm_blocks*BLS*BLS,
	 it_matrix->Half_bandwidth);
/*kew*/

/*kbw
  printf("\n nr_sm_blocks %d %d %d ->%d %d %d\n\n",nr_sm_blocks,index,Nrblocks,Nrdofbl[0],Nrdofbl[1],Nrdofbl[2]);
   for(block_row_gl = 0;  block_row_gl < nr_sm_blocks; block_row_gl++)
	 printf("z=%d ",it_matrix->block_col_ind[block_row_gl]);
  printf("\n\n\n");
	 for(block_row_gl = 0;  block_row_gl <= Nrblocks; block_row_gl++)
	 printf("w=%d ",it_matrix->block_row_ptr[block_row_gl]);

  //getchar();
  //getchar();
//*kew*/

  return(itv_nr_bcrs_matrices-1);
}

/*---------------------------------------------------------
  lar_initialize_SM_and_LV_bcrs - to initialize stiffness matrix and/or load vector
---------------------------------------------------------*/
int lar_initialize_SM_and_LV_bcrs(
  int Matrix_id,   /* in: matrix ID */
  int Scope    /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  )

{

  int i;

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];

  // sequential version - it is possible that for NUMA architectures
  // all entries will be allocated in a single processor/core memory
  if(Scope==LAC_SCOPE_SM_AND_LV){
	memset( it_matrix->block_val, 0, (it_matrix->Nr_sm_blocks*BLS*BLS)*sizeof(double) );
  }
  memset( it_matrix->rhs, 0, (it_matrix->Nrdofgl)*sizeof(double) );

  //printf("\n\n lar_initialize_SM_and_LV %d\n\n", it_matrix->Nrdofgl);
  //getchar(); getchar();

return(1);
}

/*---------------------------------------------------------
  lar_get_storage_bcrs - to compute storage of SM, LV and preconditioner
---------------------------------------------------------*/
double lar_get_storage_bcrs( /* returns: storage in MB */
  int Matrix_id   /* in: matrix ID */
			   )
{
  // seem to be quite simple :)
	return -1;
}

/**-----------------------------------------------------------
  lar_fill_assembly_table_int_ent_bcrs - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
extern int lar_fill_assembly_table_int_ent_bcrs(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int *Assembly_table_int_ent /* part of the global assembly table */
  )
{
/* get pointer to the level structure */
  itt_bcrs_matrices *it_matrix = &itv_bcrs_matrices[Matrix_id];

/*kbw
  printf("In lar_fill_assembly_table_int_ent_bcrs: matrix_id %d\n",
	 Matrix_id);
  printf("nr_dof_ent_loc %d, position %u\n",
	 Nr_dof_bl, Assembly_table_int_ent);
  int ibl;
  for(ibl=0;ibl<Nr_dof_bl; ibl++) printf("bl_id %d ", L_bl_id[ibl]);
  printf("\n");
/*kew*/


#ifdef PARALLEL
  int nrbl_inter=it_matrix->Nrbl_internal;
#endif


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

#ifdef PARALLEL
  if(bl_id_y<nrbl_inter){
#endif


	  int icrs;
	  for(icrs=it_matrix->block_row_ptr[bl_id_y];
	  icrs<it_matrix->block_row_ptr[bl_id_y+1];
	  icrs++){

	if(it_matrix->block_col_ind[icrs]==bl_id_x){

/*kbw
	  printf("found for entry %d position in BCRS %d\n", jaux, icrs*BLS*BLS);
/*kew*/

	  Assembly_table_int_ent[jaux]=icrs*BLS*BLS;
	  break;
		}
	  }

#ifdef PARALLEL
  }else
	Assembly_table_int_ent[jaux]=-1;
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

/*     Assembly_table_int_ent[jaux]=bl_id_x*BLS; */

/*   } */

/*kbw
	  getchar();
/*kew*/

  return(1);
}


/*------------------------------------------------------------
  lar_assemble_SM_and_LV_with_table_bcrs - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
------------------------------------------------------------*/
extern int lar_assemble_SM_and_LV_with_table_bcrs(
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
  //printf("\n\n BCRS-assembling with table \n\n\n");

  /* get pointer to the matrix structure */
  itt_bcrs_matrices *it_matrix = &itv_bcrs_matrices[Matrix_id];

  int idofbl,k,l;
  for(idofbl=0;idofbl<Nr_dof_bl;idofbl++) {

	int jdofbl;
	for(jdofbl=0;jdofbl<Nr_dof_bl;jdofbl++) {

	  int jaux = jdofbl+idofbl*Nr_dof_bl;

	  int icrs = Assembly_table_int_ent[jaux];

#ifdef PARALLEL
  if(icrs!=-1){
#endif

/*kbw
		  printf("found for entry %d position in CRS %d\n", jaux, icrs);
		  printf("filling with %lf: before %lf\n",
			 Stiff_mat[jaux], it_matrix->crs_val[icrs]);
/*kew*/

	  int offset = (jdofbl*Nr_dof_bl * BLS + idofbl) * BLS;

	  for(k=0; k<BLS; k++){
		for(l=0; l<BLS; l++){


	  //(idofbl*Nr_dof_bl * BLS + jdofbl) * BLS + k*Nr_dof_bl * BLS + l =?= jaux*Nr_dof_bl * BLS* BLS + k*Nr_dof_bl * BLS + l
		  //it_matrix->block_val[icrs + k*BLS + l] += Stiff_mat[jaux*Nr_dof_bl * BLS* BLS + k*Nr_dof_bl * BLS + l];

	  // small blocks are stored by columns (same as input Stiff_mat!!!)
		  it_matrix->block_val[icrs + k*BLS + l] +=
		Stiff_mat[offset + k*Nr_dof_bl * BLS + l];

		  //if(icrs + k*BLS + l<26){printf("icrs=%d, valind=%d jaux=%d index=%d %12lf",icrs, icrs+ k*BLS + l,jaux,(jdofbl*Nr_dof_bl * BLS + idofbl) * BLS + k*Nr_dof_bl * BLS + l,Stiff_mat[(jdofbl*Nr_dof_bl * BLS + idofbl) * BLS + k*Nr_dof_bl * BLS + l]);getchar();}

		}
	  }

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
	int irhs = Local_to_global_int_ent[idofbl]*BLS;
/*kbw
	//if(!(fabs(it_matrix->rhs[38])<1.e15)){
	//if(irhs+Nr_dof_bl>it_matrix->Nrdofgl){
	if(irhs>35 && irhs<40){
	  printf("filling entry %d for RHS block %d\n",
		 irhs, idofbl);
	  printf("filling with %lf: before %lf\n",
		 Rhs_vect[idofbl*BLS], it_matrix->rhs[irhs]);
  getchar();
  getchar();
  }
/*kew*/

	for (l=0;l<BLS;l++) {
	  it_matrix->rhs[irhs+l] += Rhs_vect[idofbl*BLS+l];
	  //printf("%d rhs_as_ind=%d, rhs_ind=%d ",irhs,irhs+l,idofbl*BLS+l );getchar();
	}
/*kbw
  if(!(fabs(it_matrix->rhs[38])<1.e15)){
	  printf("filling with %lf: after %lf\n",
		 Rhs_vect[idofbl], it_matrix->rhs[irhs]);
  }
/*kew*/
  }

/*kbw
  getchar();
/*kew*/
/*kbw
  //if(!(fabs(it_matrix->rhs[38])<1.e15)){
  printf("RHS after assembly %d\n", irhs); int i;
	  for(i=0;i<40;i++) printf("%20.12lf",it_matrix->rhs[i]);
	  printf("\n");
	//getchar();
	//getchar();
	//}
/*kew*/

  return(1);
}


/*------------------------------------------------------------
  lar_assemble_SM_and_LV_bcrs - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
!!!!!!!!!******************************************!!!!!!!!!!!!!!
lar_assemble_SM_and_LV_bcrs allows for sending GHOST blocks on lists of blocks.
GHOST blocks  (indicated by block->Lngb==NULL) are blocks with associated block
structures, but without associated rows of GLOBAL stiffness matrix.
Such blocks are used for matrix-vector product with rows belonging to ordinary
(not ghost) blocks. In lar_assemble_SM_and_LV_bcrs we nevertheless assume that provided
Stiff_mat is square, i.e. there are entries for all blocks. Procedures
calling lar_assemble_SM_and_LV_bcrs must conform to this convention
!!!!!!!!!******************************************!!!!!!!!!!!!!!
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_bcrs(
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

  int nrdof = Nr_dof_bl * BLS; // !!!

  int block_row; // block_row is an index for subsequent LOCAL row blocks
  int block_val_ind=-1;
  int k,l;
  int block_col_gl,block_row_gl;

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];
  int nrbl_inter=it_matrix->Nrbl_internal;

  //printf("\n\n lar_assemble_SM_and_LV_bcrs start %d %d %d\n\n",nrbl_inter, Nr_dof_bl,nrdof);

  for( block_row = 0; block_row < Nr_dof_bl; block_row++){

	 block_row_gl = L_bl_id[ block_row ]; // global ID for local DOF=row block

#ifdef PARALLEL
  if(block_row_gl<nrbl_inter){
#endif

	int initial_index = it_matrix->block_row_ptr[ block_row_gl ]; // initial index in block_col_ind
	// (and indirectly val) for a given row of square blocks

	int final_index = it_matrix->block_row_ptr[ block_row_gl+1]; // final (+1) index in block_col_ind

	// globid=block_row_gl*it_matrix->Nrdofgl*BLS; - global index - not used

	// number of non-zero blocks in block_row (i.e. the row of blocks)
	int nr_neighbors = final_index - initial_index;

	int block_col; // block_col is an index for subsequent LOCAL column blocks

	for( block_col = 0; block_col < Nr_dof_bl; block_col++){

	  block_col_gl = L_bl_id[ block_col ];  // global ID for local DOF=column block


	  // globid=block_row_gl*BLS+block_col_gl*it_matrix->Nrdofgl*BLS; - global index - not used

	  int icrs;
	  for(icrs=initial_index; icrs<final_index; icrs++){

	if(it_matrix->block_col_ind[icrs]==block_col_gl){

/*kbw
	  printf("found for entry %d position in CRS %d\n", jaux, icrs*BLS*BLS);
/*kew*/

	  // the place of column block in the list of column blocks for a given row block
	  //icrs = lar_util_chk_list_bcrs( block_col_gl,
	  //				   &it_matrix->block_col_ind[initial_index],
	  //			   nr_neighbors);
	  //if(icrs==-1){printf("Element not found in lar_assemble_SM_and_LV_bcrs\n");exit(0);}
	  //  block_val_ind = (initial_index+icrs)*BLS*BLS;

	  block_val_ind = icrs*BLS*BLS;
	  break;

	}
	  }

	  for(k=0; k<BLS; k++){
	for(l=0; l<BLS; l++){


	  // small blocks are stored by columns (same as input Stiff_mat!!!)
	  it_matrix->block_val[ block_val_ind + k*BLS + l] +=
		Stiff_mat[ (block_col*nrdof + block_row) * BLS + k*nrdof + l];

	}
	  }
	}
#ifdef PARALLEL
	}
#endif

	for (l=0;l<BLS;l++) {
	  it_matrix->rhs[block_row_gl*BLS+l] += Rhs_vect[block_row*BLS+l];
	}

  }

  return 0;
}



/*---------------------------------------------------------
lar_allocate_preconditioner_bcrs - to create preconditioner blocks
---------------------------------------------------------*/
int lar_allocate_preconditioner_bcrs( /* returns:   >0 number of diagonal blocks */
						  /*	       <=0 - error */
  int Matrix_id,   /* in: matrix ID */
  int Precon,      /* in: type of preconditioner (lah_block.h line circa 45) */
  int ILU_k // in: for ILU(k) - k;
  )
{

/*kbw
  printf("Entering lar_allocate_preconditioner_bcrs: Matrix_id %d, Precon %d, ILU_k %d\n",
	 Matrix_id, Precon, ILU_k);
  //getchar();
  //getchar();
/*kew*/


  /* pointers to solver structure and stiffness matrix blocks' info */
  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];
  int nrbl_inter=it_matrix->Nrbl_internal;

  it_matrix->Precon = Precon;
  it_matrix->ILU_k = ILU_k;

/*kbw
  printf("lar_allocate_preconditioner_bcrs: Precon %d, ILU_k %d\n",
	 it_matrix->Precon, it_matrix->ILU_k);
  //getchar();
  //getchar();
/*kew*/

  if(Precon==BLOCK_JACOBI || Precon==BLOCK_GS){

	// Here: the simplest block-diagonal preconditioner
	// (diagonal preconditioner for CRS is obtained trivially)

	// auxiliary data structure for pointers to diagonal blocks in original matrix
	it_matrix->block_diag_ptr = (int*)malloc( (it_matrix->Nrblocks)*sizeof(int));
	if( it_matrix->block_diag_ptr==NULL ) {
	  printf("Not enough memory for allocating bcrs structure in preconditioner for block_diag_ptr vector\n");
	  exit(-1);
	}

	int block_row_gl;
	for (block_row_gl = 0;  block_row_gl < nrbl_inter; block_row_gl++){

	  int sm_block;
	  for(sm_block = it_matrix->block_row_ptr[ block_row_gl ];
	  sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];
	  sm_block++){

	int block_col_gl = it_matrix->block_col_ind[ sm_block ];

	if(block_col_gl == block_row_gl){

	  it_matrix->block_diag_ptr[block_row_gl] = sm_block;
	  //printf("(%d,%d) ",block_row_gl,sm_block);
	}

	  }

	}


	it_matrix->block_diag_ips = (int*)malloc( BLS*BLS*(it_matrix->Nrblocks)*sizeof(int));
	if( it_matrix->block_diag_ips==NULL ) {
	  printf("Not enough memory for allocating bcrs structure in preconditioner for block_diag_ips vector\n");
	  exit(-1);
	}

	//for (block_row_gl = 0;  block_row_gl < it_matrix->Nrblocks; block_row_gl++) printf("%d ",it_matrix->block_diag_ptr[block_row_gl]);
	//printf("\n\n");

	// 1 line for allocating preconditioner !!!!!!!!!
	it_matrix->block_diag_precon = (double*)malloc( (it_matrix->Nrblocks*BLS*BLS)*sizeof(double));
	if( it_matrix->block_diag_precon==NULL ) {
	  printf("Not enough memory for allocating bcrs structure in preconditioner for block_diag_precon vector\n");
	  exit(-1);
	}

  } // end if BLOCK_JACOBI or BLOCK_GS preconditioner
  else if( it_matrix->Precon==MULTI_ILU || it_matrix->Precon==BLOCK_ILU ){

#ifdef PARALLEL
  printf("\nWARNING: ILU is not adapted to parallel (distributed memory) execution. Exiting!!!\n");exit(-5);
#endif

	// first the simplest ILU(0)
	if(it_matrix->ILU_k==0){

	  int nrblocks = it_matrix->Nrblocks;
	  int nr_sm_blocks = it_matrix->Nr_sm_blocks;


	  // duplicate the original stiffness matrix structure
	  it_matrix->ilu_row_ptr = (int*)malloc( (nrblocks+1)*sizeof(int));
	  if( it_matrix->ilu_row_ptr==NULL ) {
	printf("Not enough memory for allocating bcrs structure for ilu_row_ptr vector\n");
	exit(-1);
	  }
	  memcpy(it_matrix->ilu_row_ptr, it_matrix->block_row_ptr, (nrblocks+1)*sizeof(int));

	  it_matrix->ilu_col_ind = (int*)malloc( (nr_sm_blocks+1)*sizeof(int));
	  if( it_matrix->ilu_col_ind==NULL ) {
	printf("Not enough memory for allocating bcrs structure for ilu_col_ind vector\n");
	exit(-1);
	  }
	  memcpy(it_matrix->ilu_col_ind, it_matrix->block_col_ind, (nr_sm_blocks+1)*sizeof(int));

	  it_matrix->ilu_val = (double*)malloc( (nr_sm_blocks*BLS*BLS)*sizeof(double));
	  if( it_matrix->ilu_val==NULL ) {
	printf("Not enough memory for allocating bcrs structure for ilu_val vector\n");
	exit(-1);
	  }

	  // auxiliary data structure for pointers to diagonal blocks in ILU matrix
	  it_matrix->ilu_diag_ptr = (int*)malloc( (it_matrix->Nrblocks)*sizeof(int));
	  if( it_matrix->ilu_diag_ptr==NULL ) {
	printf("Not enough memory for allocating bcrs structure in preconditioner for ilu_diag_ptr vector\n");
	exit(-1);
	  }

	  it_matrix->ILU_half_bandwidth = it_matrix->Half_bandwidth;
	  it_matrix->ILU_nr_sm_blocks = it_matrix->Nr_sm_blocks;

	}
	else{ // for ilu(k) k>0

	  int i,j,k;
	  int ILU_k = it_matrix->ILU_k;
	  int nrblocks = it_matrix->Nrblocks;

	  // create the graph for neighbors of blocks within subsequent rings
	  int* nrneig_ilu;
	  int** l_neig_ilu;

	  nrneig_ilu =  (int*)malloc( (it_matrix->Nrblocks)*sizeof(int));
	  l_neig_ilu =  (int**)malloc( (it_matrix->Nrblocks)*sizeof(int*));
	  int max_nrneig_ilu_ass = 0.5*(ILU_k+2)*(ILU_k+2)*(ILU_k+2)*it_matrix->Max_nrneig;

/*kbw
	  printf("\nILU(%d) for BCRS: max_nrneig %d, max_nrneig_ilu (assumed) %d\n",
		 it_matrix->ILU_k, it_matrix->Max_nrneig, max_nrneig_ilu_ass);
/*kew*/


	  for(i=0; i<it_matrix->Nrblocks; i++){
	l_neig_ilu[i] =  (int*)malloc( max_nrneig_ilu_ass*sizeof(int*));
	  }

	  for(i=0; i<it_matrix->Nrblocks; i++){
	nrneig_ilu[i] = 0;
	for(j=0;j<max_nrneig_ilu_ass; j++){
	  l_neig_ilu[i][j] = LAC_END_OF_LIST;
	}
	  }

	  int block_row_gl;
	  for (block_row_gl = 0;  block_row_gl < it_matrix->Nrblocks; block_row_gl++){

	int nr_neig = 0;

	// put all original neighbors of this block on the list
	int sm_block;
	for(sm_block = it_matrix->block_row_ptr[ block_row_gl ];
		sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];
		sm_block++){

	  l_neig_ilu[block_row_gl][nr_neig] = it_matrix->block_col_ind[ sm_block ];
	  nr_neig++;

	}

	assert(nr_neig== (it_matrix->block_row_ptr[ block_row_gl+1 ] -
			  it_matrix->block_row_ptr[ block_row_gl ] ));
	nrneig_ilu[block_row_gl] = nr_neig;


/*kbw
	//if(block_row_gl<100 && block_row_gl%10 == 0){
	if(block_row_gl<10){

	  printf("after ring 0, row block %d, list of %d neighbors:\n",
		 block_row_gl, nrneig_ilu[block_row_gl]);

	  for(i=0;i<nrneig_ilu[block_row_gl];i++){

		printf("%5d", l_neig_ilu[block_row_gl][i]);

	  }
	  printf("\n");

	}
/*kew*/

	  }

	  for (block_row_gl = 0;  block_row_gl < it_matrix->Nrblocks; block_row_gl++){

	// initial and final number of neighbors added for the last ring
	int initial_neig = 0;
	int final_neig = nrneig_ilu[block_row_gl];

	int i_iluk;
	for(i_iluk = 0; i_iluk < ILU_k; i_iluk++){

	  // for all neighbors added in the last ring (for i_iluk==0 for original neighbors)
	  int ineig;
	  for(ineig = initial_neig; ineig<final_neig; ineig++){

		// find block column number of the neighbor
		int block_row_neig = l_neig_ilu[block_row_gl][ineig];

		// for all original neighbors of this neighbor
		int sm_block;
		for(sm_block = it_matrix->block_row_ptr[ block_row_neig ];
		sm_block < it_matrix->block_row_ptr[ block_row_neig+1 ];
		sm_block++){

		  int block_col_gl = it_matrix->block_col_ind[ sm_block ];

		  int info = lar_util_put_list_bcrs(block_col_gl,
						l_neig_ilu[block_row_gl],
						max_nrneig_ilu_ass);
		  if(info>=0) nrneig_ilu[block_row_gl]++;
		  else if(info==LAC_PUT_LIST_FULL){
		printf("Too short list of neighbors in ILU(k) for bcrs.\n");
		printf("Correct max_nrneig_ilu_ass in lar_allocate_preconditioner_bcrs! Exiting.\n");
		exit(-1);
		  }

		}

	  } // end loop over neighbors added in the last ring

	  initial_neig = final_neig;
	  final_neig = nrneig_ilu[block_row_gl];

/*kbw
	  //if(block_row_gl<100 && block_row_gl%10 == 0){
	  if(block_row_gl<10){

		printf("after ring %d, row block %d, list of %d neighbors:\n",
		   i_iluk, block_row_gl, nrneig_ilu[block_row_gl]);

		for(i=0;i<nrneig_ilu[block_row_gl];i++){

		  printf("%5d", l_neig_ilu[block_row_gl][i]);

		}
		printf("\n");

	  }
/*kew*/

	}

	  }

	  // sort the lists of neighbors
	  for (block_row_gl = 0;  block_row_gl < it_matrix->Nrblocks; block_row_gl++){

	lar_sort_short(l_neig_ilu[block_row_gl], 0,
			   nrneig_ilu[block_row_gl]-1);

/*kbw
	//if(block_row_gl<100 && block_row_gl%10 == 0){
	//if(block_row_gl<10)
	{

	  printf("after sorting: row block %d, list of %d neighbors:\n",
		 block_row_gl, nrneig_ilu[block_row_gl]);

	  for(i=0;i<nrneig_ilu[block_row_gl];i++){

		printf("%5d", l_neig_ilu[block_row_gl][i]);

	  }
	  printf("\n");

	}
/*kew*/

	  }

	  // based on the list of neighbors create bcrs data structure


	  it_matrix->ilu_row_ptr = (int*)malloc( (nrblocks+1)*sizeof(int));
	  if( it_matrix->ilu_row_ptr==NULL ) {
	printf("Not enough memory for allocating bcrs structure for ilu_row_ptr vector\n");
	exit(-1);
	  }
	  memset(it_matrix->ilu_row_ptr, 0, (nrblocks+1)*sizeof(int));

	  int max_nrneig = 0;
	  int nr_sm_blocks = 0;
	  for (block_row_gl = 0;  block_row_gl < nrblocks; block_row_gl++){

	if(nrneig_ilu[block_row_gl]>0){ // ghost blocks are excluded

	  nr_sm_blocks += nrneig_ilu[block_row_gl];

	}

	if(nrneig_ilu[block_row_gl]>max_nrneig) max_nrneig = nrneig_ilu[block_row_gl];

	  }
	  it_matrix->ILU_nr_sm_blocks = nr_sm_blocks;

	  it_matrix->ilu_col_ind = (int*)malloc( (nr_sm_blocks+1)*sizeof(int));
	  if( it_matrix->ilu_col_ind==NULL ) {
	printf("Not enough memory for allocating bcrs structure for ilu_col_ind vector\n");
	exit(-1);
	  }
	  memset(it_matrix->ilu_col_ind, 0, (nr_sm_blocks+1)*sizeof(int));


	  it_matrix->ilu_val = (double*)malloc( (nr_sm_blocks*BLS*BLS)*sizeof(double));
	  if( it_matrix->ilu_val==NULL ) {
	printf("Not enough memory for allocating bcrs structure for ilu_val vector\n");
	exit(-1);
	  }

	  // auxiliary data structure for pointers to diagonal blocks in ILU matrix
	  it_matrix->ilu_diag_ptr = (int*)malloc( (it_matrix->Nrblocks)*sizeof(int));
	  if( it_matrix->ilu_diag_ptr==NULL ) {
	printf("Not enough memory for allocating bcrs structure in preconditioner for ilu_diag_ptr vector\n");
	exit(-1);
	  }

	  int bandwidth = 0;
	  int index = 0;
	  int ilu_row_gl;
	  for (ilu_row_gl = 0;  ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

	it_matrix->ilu_row_ptr[ ilu_row_gl ] = index;
	//printf("%d %d \n",ilu_row_gl,index);getchar();

	for( i=0; i< nrneig_ilu[ilu_row_gl]; i++){

#ifdef DEBUG_LSM
	  // checking assumptions
	  if(i>0){
		//assert(L_offbl[iblock][neigi-1]<L_offbl[iblock][neigi]);
		if(l_neig_ilu[ilu_row_gl][i-1]>=l_neig_ilu[ilu_row_gl][i]){
		  printf("not sorted list of neighbors in lar_allocate_SM_and_LV_crs\n");
		  exit(-1);
		}
	  }
#endif

	  if(abs(l_neig_ilu[ilu_row_gl][i] - ilu_row_gl)>bandwidth) {
		bandwidth=abs(l_neig_ilu[ilu_row_gl][i] - ilu_row_gl);
	  }

	  it_matrix->ilu_col_ind[index] = l_neig_ilu[ ilu_row_gl ][i];
	  index ++;

	} //the end of loop over neighbors

	  }
	  it_matrix->ilu_row_ptr[it_matrix->Nrblocks]=index;
	  it_matrix->ILU_half_bandwidth = bandwidth;

	  assert(index==nr_sm_blocks);

#ifdef DEBUG_LSM
	  // test
	  for (ilu_row_gl = 0;  ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

	int sm_block;
	for(sm_block = it_matrix->ilu_row_ptr[ ilu_row_gl ];
		sm_block < it_matrix->ilu_row_ptr[ ilu_row_gl+1 ];
		sm_block++){

	  int ineig = sm_block - it_matrix->ilu_row_ptr[ ilu_row_gl ];

	  if(ineig>=nrneig_ilu[ilu_row_gl]){
		printf("Error 125634 in setting BCRS for ILU(%d)\n", it_matrix->ILU_k);
	  }

	  if(it_matrix->ilu_col_ind[ sm_block ] != l_neig_ilu[ ilu_row_gl ][ineig]){
		printf("Error 125834 in setting BCRS for ILU(%d)\n", it_matrix->ILU_k);
	  }

	}
	  }

#endif

	}

	// ilu_diag_ptr - for ILU(0) and ILU(k>0)
	int ilu_row_gl;
#pragma omp parallel for default(none) firstprivate(it_matrix)
	for (ilu_row_gl = 0;  ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

	  int sm_block;
	  for(sm_block = it_matrix->ilu_row_ptr[ ilu_row_gl ];
	  sm_block < it_matrix->ilu_row_ptr[ ilu_row_gl+1 ];
	  sm_block++){

	int ilu_col_gl = it_matrix->ilu_col_ind[ sm_block ];

	if(ilu_col_gl == ilu_row_gl){

	  it_matrix->ilu_diag_ptr[ilu_row_gl] = sm_block;
	  //printf("(%d,%d) ",ilu_row_gl,sm_block);
	}

	  }

	}

#ifdef DEBUG_LSM
	for (ilu_row_gl = 0;  ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

	  int sm_block;
	  for(sm_block = it_matrix->ilu_row_ptr[ ilu_row_gl ];
	  sm_block < it_matrix->ilu_row_ptr[ ilu_row_gl+1 ];
	  sm_block++){

	if(it_matrix->ilu_col_ind[ it_matrix->ilu_diag_ptr[ilu_row_gl]] != ilu_row_gl){
	  printf("error 7111123 testing ilu bcrs matrices\n"); exit(-1);
	}

	  }
	}
#endif

  }
  else{

	printf("not implemented preconditioner %d in BCRS module\n",
	   it_matrix->Precon);

  }

  return(0);

}

/*---------------------------------------------------------
  lar_fill_preconditioner_bcrs - to fill preconditioner (here by factorizing
							diagonal blocks or by block ILU(0) factorization)
---------------------------------------------------------*/
int lar_fill_preconditioner_bcrs(
  int Matrix_id   /* in: matrix ID */
	)
{

#ifdef _OPENMP
  //#define TIME_TEST_ILU
#endif

#ifdef TIME_TEST_ILU
  double nr_oper = 0;
  double nr_access = 0;
  double exec_time = time_clock();

#endif

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];
  int nrbl_inter=it_matrix->Nrbl_internal;

  if(it_matrix->Precon==BLOCK_JACOBI || it_matrix->Precon==BLOCK_GS){

	int block_row_gl,k,l, blsl=BLS;
	int iaux=BLS;

	//printf("BCRS_lar_fill_preconditioner_bcrs %d\n\n",it_matrix->Nrblocks);
	// getchar();
	// getchar();

	int sm_block;

/*
	for (block_row_gl = 0;  block_row_gl < it_matrix->Nrblocks; ++block_row_gl){
	  for(sm_block = it_matrix->block_row_ptr[ block_row_gl ];
	  sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];
	  sm_block++){
	for(k=0; k<BLS; k++){
	  for(l=0; l<BLS; l++)
		printf("(%d,%20.15lf) %20.15lf\n",it_matrix->block_col_ind[sm_block],it_matrix->block_val[sm_block*BLS*BLS+k*BLS+l],it_matrix->rhs[block_row_gl*BLS+k]);
	}getchar();
	  }printf("\n");
	  printf("\n\n");
	}
	printf("\n\n\n");

 //*/


#pragma omp parallel for default(none) firstprivate(Matrix_id, it_matrix, blsl,nrbl_inter) private(block_row_gl,k,l,iaux)
	for (block_row_gl = 0;  block_row_gl < nrbl_inter; ++block_row_gl){

	  // // CRS
	  // it_matrix->block_diag_precon[block_row_gl*BLS*BLS] = 1.0/it_matrix->block_val[it_matrix->block_diag_ptr[block_row_gl]];
	  int info;
	  // BCRS
	  // rewrite diagonal block from global SM
	  for(k=0; k<BLS; k++){
	for(l=0; l<BLS; l++){

	  // it is possible to rewrite by rows or by columns !!!
	  it_matrix->block_diag_precon[block_row_gl*BLS*BLS + k*BLS + l] =
		it_matrix->block_val[ it_matrix->block_diag_ptr[block_row_gl]*BLS*BLS + k*BLS + l];

	}
	  }


	  /* LU decompose diagonal block - using LAPACK or better with manually optimized routine */

	  dgetrf_(&blsl,&blsl,&it_matrix->block_diag_precon[block_row_gl*BLS*BLS],
		  &blsl,&it_matrix->block_diag_ips[block_row_gl*BLS*BLS],&iaux);

	  //    printf(" %d",iaux);

	  // the best way is to invert the block and then use matrix-vector multiply
	  // in preconditioning

	} // end loop over diagonal blocks

#ifdef TIME_TEST_ILU
	nr_oper += it_matrix->Nrblocks * 2.0/3.0*BLS*BLS*BLS; // dgetrf
	nr_access += it_matrix->Nrblocks * (2*BLS*BLS+3*BLS*BLS); // rewriting and dgetrf
#endif


	/* printf("\n\n WWW\n"); */
	/* for (block_row_gl = 0;  block_row_gl < it_matrix->Nrblocks; ++block_row_gl){ */
	/*   for(k=0; k<BLS; k++){ */
	/* 	for(l=0; l<BLS; l++){ */
	/* 	  printf("(%d,%20.15lf) ",(block_row_gl*BLS*BLS + k*BLS + l),it_matrix->block_diag_precon[block_row_gl*BLS*BLS + k*BLS + l]);} */
	/* 	//printf("\n\n block_diag_ips=%d\n\n",it_matrix->block_diag_ips[block_row_gl*BLS*BLS+k]); */
	/*   } */
	/*   printf("\n"); */

	/* } */
	/* for (block_row_gl = 0;  block_row_gl < it_matrix->Nrblocks; ++block_row_gl){ */
	/*   for(sm_block = it_matrix->block_row_ptr[ block_row_gl ]; 	sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];	sm_block++){ */
	/* 	printf("(%d,%lf) ",it_matrix->block_col_ind[sm_block],it_matrix->block_val[sm_block]); */
	/*   } */
	/*   printf("\n"); */
	/* } */
	/* printf("\n\n\n"); */

  } // end if BLOCK_JACOBI or BLOCK_GS preconditioner
  else if( it_matrix->Precon==MULTI_ILU || it_matrix->Precon==BLOCK_ILU ){

	// first the simplest ILU(0)
	if(it_matrix->ILU_k==0){

	  // rewrite original matrix to ilu - in parallel (NUMA - OK)
	  //it_matrix->ilu_row_ptr[0]=it_matrix->block_row_ptr[0]; - done by memcpy

	  int block_row_gl;
#pragma omp parallel for default(none) firstprivate(it_matrix)
	  for(block_row_gl=0; block_row_gl < it_matrix->Nrblocks; block_row_gl++){

	//it_matrix->ilu_row_ptr[block_row_gl+1] = it_matrix->block_row_ptr[block_row_gl+1]; - done by memcpy
	int sm_block;
	for(sm_block = it_matrix->ilu_row_ptr[ block_row_gl ];
		sm_block < it_matrix->ilu_row_ptr[ block_row_gl+1 ];
		sm_block++){

	  //it_matrix->ilu_col_ind[ sm_block ] = it_matrix->block_col_ind[ sm_block ]; - done by memcpy
	  int block_col_gl = it_matrix->ilu_col_ind[ sm_block ];

	  int block_val_ind = sm_block*BLS*BLS;

	  int k;
	  for(k=0; k<BLS*BLS; k++){
		it_matrix->ilu_val[ block_val_ind + k] = it_matrix->block_val[ block_val_ind + k];
	  }

	} // end loop over horizontal blocks
	  } // end loop over row blocks


#ifdef DEBUG_LSM
	  //test
	  for(block_row_gl=0; block_row_gl < it_matrix->Nrblocks; block_row_gl++){

	if(it_matrix->block_row_ptr[block_row_gl] != it_matrix->ilu_row_ptr[block_row_gl] ||
	   it_matrix->block_row_ptr[block_row_gl+1] != it_matrix->ilu_row_ptr[block_row_gl+1]){
	  printf("error 1234 testing ilu bcrs matrices\n"); exit(-1);
	}

	if(it_matrix->block_col_ind[ it_matrix->ilu_diag_ptr[block_row_gl]] != block_row_gl){
	  printf("error 7123 testing ilu bcrs matrices\n"); exit(-1);
	}


	int sm_block;
	for(sm_block = it_matrix->block_row_ptr[ block_row_gl ];
		sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];
		sm_block++){

	  if(it_matrix->block_col_ind[ sm_block ] != it_matrix->ilu_col_ind[ sm_block ]){
		printf("error 123 testing ilu bcrs matrices\n"); exit(-1);
	  }

	  int block_val_ind = sm_block*BLS*BLS;
	  int k,l;

	  for(l=0; l<BLS; l++){
		for(k=0; k<BLS; k++){
		  if(it_matrix->block_val[ block_val_ind + k + l*BLS] !=
		 it_matrix->ilu_val[ block_val_ind + k + l*BLS]){
		printf("error 4123 testing ilu bcrs matrices\n"); exit(-1);
		  }
		}
	  }

	}
	  } // end test matrices
#endif

	} // end ILU(0)
	else{ // ILU(k>0)

	  // initialize preconditioner
	  memset(it_matrix->ilu_val, 0, (it_matrix->ILU_nr_sm_blocks*BLS*BLS)*sizeof(double));

	  int ilu_row_gl;
#pragma omp parallel for default(none) firstprivate(it_matrix)
	  for(ilu_row_gl=0; ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

	int sm_block_ilu;
	for(sm_block_ilu = it_matrix->ilu_row_ptr[ ilu_row_gl ];
		sm_block_ilu < it_matrix->ilu_row_ptr[ ilu_row_gl+1 ];
		sm_block_ilu++){

	  int ifound=0;
	  int sm_block;
	  for(sm_block = it_matrix->block_row_ptr[ ilu_row_gl ];
		  sm_block < it_matrix->block_row_ptr[ ilu_row_gl+1 ];
		  sm_block++){

		if(it_matrix->block_col_ind[sm_block]==it_matrix->ilu_col_ind[sm_block_ilu]){


		  int block_val_ind = sm_block*BLS*BLS;
		  int block_val_ind_ilu = sm_block_ilu*BLS*BLS;
		  int k,l;

		  for(l=0; l<BLS; l++){
		for(k=0; k<BLS; k++){

		  it_matrix->ilu_val[ block_val_ind_ilu + k + l*BLS] =
			it_matrix->block_val[ block_val_ind + k + l*BLS];

		}
		  }

		  ifound=1;
		  break;
		} // end if two blocks match

	  } // and loop over original neighbours

	} // and loop over ilu neighbours

	  } // and loop over row blocks

#ifdef DEBUG_LSM
	  //test
	  for(ilu_row_gl=0; ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

	int sm_block_ilu;
	for(sm_block_ilu = it_matrix->ilu_row_ptr[ ilu_row_gl ];
		sm_block_ilu < it_matrix->ilu_row_ptr[ ilu_row_gl+1 ];
		sm_block_ilu++){

	  // test 1 - original matrix
	  int ifound=0;
	  int sm_block;
	  for(sm_block = it_matrix->block_row_ptr[ ilu_row_gl ];
		  sm_block < it_matrix->block_row_ptr[ ilu_row_gl+1 ];
		  sm_block++){

		if(it_matrix->block_col_ind[sm_block]==it_matrix->ilu_col_ind[sm_block_ilu]){


		  int block_val_ind = sm_block*BLS*BLS;
		  int block_val_ind_ilu = sm_block_ilu*BLS*BLS;
		  int k,l;

		  for(l=0; l<BLS; l++){
		for(k=0; k<BLS; k++){
		  if(it_matrix->block_val[ block_val_ind + k + l*BLS] !=
			 it_matrix->ilu_val[ block_val_ind_ilu + k + l*BLS]){
			printf("error 444123 testing ilu bcrs matrices\n"); exit(-1);
		  }
		}
		  }

		  ifound=1;
		  break;
		} // end if two blocks match

	  } // and loop over original neighbours

	  if(ifound==0){
		// new entry (not present in original matrix)

		// 1. must be zero
		int block_val_ind_ilu = sm_block_ilu*BLS*BLS;
		int k,l;

		for(l=0; l<BLS; l++){
		  for(k=0; k<BLS; k++){
		if(fabs(it_matrix->ilu_val[ block_val_ind_ilu + k + l*BLS])>1.e-12){
		  printf("error 44412223 testing ilu bcrs matrices\n"); exit(-1);
		}
		  }
		}


		// 2. must have symmetric counterpart (structural symmetry)
		ifound = 0;
		// go into the row corresponding to the column
		int k_block = it_matrix->ilu_col_ind[sm_block_ilu];

		for(sm_block = it_matrix->ilu_row_ptr[ k_block ];
		sm_block < it_matrix->ilu_row_ptr[ k_block+1 ];
		sm_block++){

		  if(it_matrix->ilu_col_ind[sm_block]==ilu_row_gl){

		ifound = 1;

		// this new entry must also be zero
		int block_val_ind_ilu = sm_block*BLS*BLS;
		int k,l;
		for(l=0; l<BLS; l++){
		  for(k=0; k<BLS; k++){
			if(fabs(it_matrix->ilu_val[ block_val_ind_ilu + k + l*BLS])>1.e-12){

			  printf("original row %d, original column %d, block %d\n",
				 ilu_row_gl, k_block, sm_block_ilu);
			  printf("symmetric row %d, symmetric column %d, block %d\n",
				 k_block, it_matrix->ilu_col_ind[sm_block], sm_block);
			  printf("k %d, l %d, value %lf\n",
				 k, l, it_matrix->ilu_val[ block_val_ind_ilu + k + l*BLS]);

			  printf("error 444111222333 testing ilu bcrs matrices\n"); exit(-1);
			}
		  }
		}

		break;
		  } // end if found symmetric counterpart

		} // end loop over small block in the row corresponding to the column

		if(ifound == 0){

		  // if symmetric counterpart not found
		  printf("error 544111222333 testing ilu bcrs matrices\n"); exit(-1);


		}


	  }

	} // and loop over ilu neighbours


	  } // end test matrices
#endif


	}

#ifdef TIME_TEST_ILU
	  nr_access += it_matrix->ILU_nr_sm_blocks * (2*BLS*BLS); // rewriting
#endif

/*kbw
	  printf("BEFORE ILU DECOMPOSITION:\n");
	  {
	  int ilu_row_gl;
	  for(ilu_row_gl=0; ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

		int sm_block;
		for(sm_block = it_matrix->ilu_row_ptr[ ilu_row_gl ];
			sm_block < it_matrix->ilu_row_ptr[ ilu_row_gl+1 ];
			sm_block++){

		  int ilu_col_gl = it_matrix->ilu_col_ind[ sm_block ];

		  int ilu_val_ind = sm_block*BLS*BLS;

		  printf("row %d, column %d, small_block %d:\n",
			 ilu_row_gl, ilu_col_gl, sm_block);

		  int i;
		  for(i=0;i<BLS*BLS;i++){
			printf("%25.15lf", it_matrix->ilu_val[ilu_val_ind+i]);
		  }
		  printf("\n");
	  }
	  getchar();
	}
	  }
/*kew*/



/* PERFORM INCOMPLETE LU FACTORIZATION */
	int version=2;

#ifdef PARALLEL
version=0;
#endif

	printf("Starting factorization - choose the version (current %d)\n", version);
	//    scanf("%d", &version);


	//VERSION 0 - sequential

	if(version==0){

	{

	int i,j,k;

	/* constants */
	int ione = 1;
	double done = 1.0;
	double dmone = -1.0;


	/* for each row block */
	int i_block;
	for(i_block=0; i_block < it_matrix->Nrblocks; i_block++){

	  /* if block is active - necessary to allow for inverting diagonal small block */
	  if(it_matrix->ilu_row_ptr[i_block]<it_matrix->ilu_row_ptr[i_block+1]){

/*kbw
	printf("\nin row block i_block %d, sm_blocks %d-%d, diag %d\n",
		   i_block, it_matrix->ilu_row_ptr[i_block],
		   it_matrix->ilu_row_ptr[i_block+1],
		   it_matrix->ilu_diag_ptr[i_block]);
/*kew*/


	int sm_block;
	// for all k<i - all row blocks above i for which small blocks a_ik exist
	// if a_ik is zero than  we do not have to subtract i-th row block
	// from k-th row block (to zero a_ik)
	for(sm_block = it_matrix->ilu_row_ptr[i_block];
		sm_block < it_matrix->ilu_diag_ptr[i_block];
		sm_block++){

	  int k_block = it_matrix->ilu_col_ind[sm_block];

	  /* compute product a_ik = a_ik * a_kk^-1 */
	  // new a_ik is used in calculations and left in lower part of LU
	  // to be used again in forward reduction

	  // pointer to a_ik -> sm_block*BLS*BLS
	  int a_ik_p = sm_block*BLS*BLS;
	  // pointer to a_kk -> it_matrix->ilu_diag_ptr[k_block]*BLS*BLS
	  int a_kk_p = it_matrix->ilu_diag_ptr[k_block]*BLS*BLS;
	  // a_kk holds original a_kk inverted (in last step in LU decomposition)

/*kbw
		if(i_block==6 && k_block==5)
	  //if(sm_block==55)
	  //if(i_block==352)
	  {
		int i;
		printf("\ni_block %d, k_block %d, sm_block %d, Dia_kk before solution:\n",
		   i_block, k_block, it_matrix->ilu_diag_ptr[k_block] );
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_kk_p+i]);
		}
		printf("\n");
		printf("i_block %d, k_block %d, sm_block %d, Aux_ik before solution:\n",
		   i_block, k_block, sm_block);
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_ik_p+i]);
		}
		printf("\n");
	  }
/*kew*/

	  {
	  // MUST BE VECTORIZED!!!
	  double a_work[BLS*BLS];
	  int i,j,k;
	  for(i=0;i<BLS;i++){
		for(j=0;j<BLS;j++){
		  double daux = 0.0;
		  for(k=0;k<BLS;k++){

		daux += it_matrix->ilu_val[a_ik_p+i+BLS*k] *
		  it_matrix->ilu_val[a_kk_p+k+j*BLS];

		  }

		  a_work[i+j*BLS]=daux;

		}
	  }
	  for(i=0;i<BLS*BLS;i++) it_matrix->ilu_val[a_ik_p+i] = a_work[i];


#ifdef TIME_TEST_ILU
	nr_oper += 2.0*BLS*BLS*BLS; // dgemm
	nr_access += 5.0*BLS*BLS; // optimal version - each matrix is retrieved only once
#endif

	  } // end single thread region with barrier at the end
/*kbw
		if(i_block==6 && k_block==5)
		  //if(sm_block==55)
	  //if(i_block==352)
	  {
		int i;
		printf("i_block %d, k_block %d, sm_block %d, Aux_ik after solution:\n",
		   i_block, k_block, sm_block);
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_ik_p+i]);
		}
		printf("\n");
		getchar(); getchar();
	  }
/*kew*/


	  // now k becomes an index for the k-th row block and i-th row block
	  // is subtracted from the k-th row block (with zero blocks ommitted!!!)

	  // we have to advance the pointer in k-th row until j-th column is found
	  // since j_block > k_block we start from the diagonal in k-th row
	  int smbl_kj = it_matrix->ilu_diag_ptr[k_block];

	  /* loop over all neighbors of i_block, with indices greater than k */
	  int j_smbl;
	  for(j_smbl = sm_block+1; j_smbl < it_matrix->ilu_row_ptr[i_block+1]; j_smbl++){

		int j_block = it_matrix->ilu_col_ind[j_smbl];

		// pointer to a_ik -> sm_block*BLS*BLS
		// a_ik now holds a_ik * a_kk^-1
		int a_ik_p = sm_block*BLS*BLS;

		// pointer to a_ij
		int a_ij_p = j_smbl*BLS*BLS;


/*kbw
			smbl_kj = it_matrix->ilu_diag_ptr[k_block];
		printf("\ntest starting position %d (%d)\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj]);

		// each thread advance its own counter for kj blocks
		while(it_matrix->ilu_col_ind[smbl_kj] < j_block  //) smbl_kj++;
		  &&  smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;
		//smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;

		printf("test finished position %d (%d), j_block %d, end row %d\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj], j_block,
		   it_matrix->ilu_row_ptr[k_block+1]);
		int test_sm=smbl_kj;

/*kew*/
		// we have to advance the pointer in k-th row until j-th column is found
		// since j_block > k_block we start from the diagonal in k-th row
		//smbl_kj = it_matrix->ilu_diag_ptr[k_block];
		// instead of resetting for each new j_block we can just go back
		// one block from where we finished for the previous j_block
		if(smbl_kj<2) smbl_kj--;

/*kbw
		printf("\nstarting position %d (%d)\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj]);
/*kew*/

		// each thread advance its own counter for kj blocks
		while(it_matrix->ilu_col_ind[smbl_kj] < j_block  //) smbl_kj++;
		  &&  smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;

/*kbw
		printf("finished position %d (%d), j_block %d, end row %d\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj], j_block,
		   it_matrix->ilu_row_ptr[k_block+1]);
		if(smbl_kj!=test_sm) {getchar();getchar();}
/*kew*/

		// if a_kj != 0 ; i.e. small block exist for index kj
		if(it_matrix->ilu_col_ind[smbl_kj]==j_block &&
		   smbl_kj < it_matrix->ilu_row_ptr[k_block+1]){

		  // pointer to a_kj
		  int a_kj_p = smbl_kj*BLS*BLS;

/*kbw
	  if(j_smbl==55)
		{
		  int i;
		  printf("\ni %d, j %d (small block %d, pointer %d), a_ij before solution:\n",
			 i_block, j_block, j_smbl, a_ij_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ij_p+i]);
		  }
		  printf("\n");
		  printf("i %d, k %d (small block %d, pointer %d), a_ik before solution:\n",
			 i_block, k_block, sm_block, a_ik_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ik_p+i]);
		  }
		  printf("\n");
		  printf("k %d, j %d (small block %d, pointer %d), a_kj before solution:\n",
			 k_block, j_block, smbl_kj, a_kj_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_kj_p+i]);
		  }
		  printf("\n");
		 getchar();
		}
/*kew*/

/* compute a_ij = a_ij - a_ik*a_kj */

// LAPACK VERSION
		  int bls = BLS;
		  double dmone = -1.0;
		  double done = 1.0;
		  char c = 'N';
		  dgemm_(&c, &c, &bls, &bls, &bls, &dmone,
				 &it_matrix->ilu_val[a_ik_p], &bls,
				 &it_matrix->ilu_val[a_kj_p], &bls, &done,
				 &it_matrix->ilu_val[a_ij_p], &bls);


		  // j,k,i VERSION

/* 	      // MUST BE VECTORIZED!!! */
/* 	      int i,j,k; */
/* 	      //double daux[BLS*BLS]; */
/* 	      //double eaux[BLS*BLS]; */
/* 	      //double faux[BLS*BLS]; */
/* 	      //for(i=0;i<BLS*BLS;i++){daux[i] = 0.0;} */
/* 	      //for(i=0;i<BLS*BLS;i++){eaux[i] = it_matrix->ilu_val[a_ik_p+i];} */
/* 	      //for(i=0;i<BLS*BLS;i++){faux[i] = it_matrix->ilu_val[a_kj_p+i];} */
/* 	      for(j=0;j<BLS;j++){ */
/* 		for(k=0;k<BLS;k++){ */
/* 		  double faux = it_matrix->ilu_val[a_kj_p+k+j*BLS]; */
/* 		  //#pragma vector always */
/* #pragma ivdep */
/* 		  for(i=0;i<BLS;i++){ */
/* 		    it_matrix->ilu_val[a_ij_p+i+j*BLS] -= */
/* 		    //daux[i+j*BLS] +=  */
/* 		      it_matrix->ilu_val[a_ik_p+i+BLS*k] * faux; */
/* 		      //eaux[i+BLS*k] * faux; */
/* 		      //it_matrix->ilu_val[a_kj_p+k+j*BLS]; */
/* 		      //faux[k+j*BLS]; */

/* 		  } */


/* 		} */
/* 	      } */

		  //for(i=0;i<BLS*BLS;i++){
		  //it_matrix->ilu_val[a_ij_p+i] -= daux[i];
		  //}


		  // i,j,k VERSION

		  /* // MUST BE VECTORIZED!!! */
		  /* for(i=0;i<BLS;i++){ */
		  /* 	for(j=0;j<BLS;j++){ */
		  /* 	  double daux = 0.0; */
		  /* 	  for(k=0;k<BLS;k++){ */
		  /* 	    daux += it_matrix->ilu_val[a_ik_p+i+BLS*k] * */
		  /* 	      it_matrix->ilu_val[a_kj_p+k+j*BLS]; */

		  /* 	  } */
		  /* 	  it_matrix->ilu_val[a_ij_p+i+j*BLS] -= daux; */

		  /* 	} */
		  /* } */

#ifdef TIME_TEST_ILU
	nr_oper += 2.0*BLS*BLS*BLS; // dgemm
	nr_access += 4.0*BLS*BLS; // optimal version - each matrix is retrieved only once
							  // for reading and a_ij once for writing
#endif



/*kbw
	  if(j_smbl==55)
		{
		  printf("\ni %d, j %d (small block %d, pointer %d), a_ij after solution:\n",
			 i_block, j_block, j_smbl, a_ij_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ij_p+i]);
		  }
		  printf("\n");
		  getchar();
		}
/*kew*/


		} // end if non-zero a_kj found for a_ij

	  } // end parallel loop over small blocks in k-th row block


	} // end loop over all row blocks (k_blocks) above i-th block

/* invert Dia for i_block */
	// pointer to a_ii -> it_matrix->ilu_diag_ptr[i_block]*BLS*BLS
	int a_ii_p = it_matrix->ilu_diag_ptr[i_block]*BLS*BLS;

/*kbw
	  if(i_block==55)
	  {
	int i;
	printf("\ni_block %d, Dia before inverting:\n",i_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");
	  }
/*kew*/

	{
	double a_work[(BLS+1)*BLS];
	int i;
	for(i=0;i<(BLS+1)*BLS;i++) a_work[i]=0.0;
	int ips_work[BLS*BLS];
	int bls = BLS;
/*kb*/
	int iaux;
	dgetrf_(&bls, &bls, &it_matrix->ilu_val[a_ii_p],
		&bls, ips_work, &iaux);

	int lwork = bls*bls;
	dgetri_(&bls, &it_matrix->ilu_val[a_ii_p],
		&bls, ips_work , a_work, &lwork, &iaux);
/*kb*/

#ifdef TIME_TEST_ILU
	nr_oper += (2.0/3.0+4.0/3.0)*BLS*BLS*BLS; // dgetrf+dgetri
	nr_access += 6.0*BLS*BLS; // optimal version - a_ii accessed twice in dgetrf
								  // and once in dgetri, a_work twice in dgetri (+init)
#endif

/*kbw
	  //if(i_block>235)
	  if(i_block==55)
	  {
	int i;
	printf("i_block %d, Dia after inverting:\n",i_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");

	//getchar();
	  }
/*kew*/

	} // end single thread region with implicit barrier at the end

	  } // end if not ghost block

	} // end outer loop over row blocks - i_block

	} // end parallel region

	}
	else if(version==1){

	//VERSION I - PARALLEL IN EACH ROW - NOT EFFICIENT

#ifdef TIME_TEST_ILU
#pragma omp parallel default(none) firstprivate(it_matrix)  shared(nr_oper, nr_access)
#else
#pragma omp parallel default(none) firstprivate(it_matrix)
#endif
	{

	int i,j,k;

	/* constants */
	int ione = 1;
	double done = 1.0;
	double dmone = -1.0;


	/* for each row block */
	int i_block;
	for(i_block=0; i_block < it_matrix->Nrblocks; i_block++){

	  /* if block is active - necessary to allow for inverting diagonal small block */
	  if(it_matrix->ilu_row_ptr[i_block]<it_matrix->ilu_row_ptr[i_block+1]){

/*kbw
	printf("\nin row block i_block %d, sm_blocks %d-%d, diag %d\n",
		   i_block, it_matrix->ilu_row_ptr[i_block],
		   it_matrix->ilu_row_ptr[i_block+1],
		   it_matrix->ilu_diag_ptr[i_block]);
/*kew*/


	int sm_block;
	// for all k<i - all row blocks above i for which small blocks a_ik exist
	// if a_ik is zero than  we do not have to subtract i-th row block
	// from k-th row block (to zero a_ik)
	for(sm_block = it_matrix->ilu_row_ptr[i_block];
		sm_block < it_matrix->ilu_diag_ptr[i_block];
		sm_block++){

	  int k_block = it_matrix->ilu_col_ind[sm_block];

	  /* compute product a_ik = a_ik * a_kk^-1 */
	  // new a_ik is used in calculations and left in lower part of LU
	  // to be used again in forward reduction

	  // pointer to a_ik -> sm_block*BLS*BLS
	  int a_ik_p = sm_block*BLS*BLS;
	  // pointer to a_kk -> it_matrix->ilu_diag_ptr[k_block]*BLS*BLS
	  int a_kk_p = it_matrix->ilu_diag_ptr[k_block]*BLS*BLS;
	  // a_kk holds original a_kk inverted (in last step in LU decomposition)

/*kbw
		if(i_block==6 && k_block==5)
	  //if(sm_block==55)
	  //if(i_block==352)
	  {
		int i;
		printf("\ni_block %d, k_block %d, sm_block %d, Dia_kk before solution:\n",
		   i_block, k_block, it_matrix->ilu_diag_ptr[k_block] );
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_kk_p+i]);
		}
		printf("\n");
		printf("i_block %d, k_block %d, sm_block %d, Aux_ik before solution:\n",
		   i_block, k_block, sm_block);
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_ik_p+i]);
		}
		printf("\n");
	  }
/*kew*/

#pragma omp single
	  {
	  // MUST BE VECTORIZED!!!
	  double a_work[BLS*BLS];
	  int i,j,k;
	  for(i=0;i<BLS;i++){
		for(j=0;j<BLS;j++){
		  double daux = 0.0;
		  for(k=0;k<BLS;k++){

		daux += it_matrix->ilu_val[a_ik_p+i+BLS*k] *
		  it_matrix->ilu_val[a_kk_p+k+j*BLS];

		  }

		  a_work[i+j*BLS]=daux;

		}
	  }
	  for(i=0;i<BLS*BLS;i++) it_matrix->ilu_val[a_ik_p+i] = a_work[i];


#ifdef TIME_TEST_ILU
#pragma omp atomic
	nr_oper += 2.0*BLS*BLS*BLS; // dgemm
#pragma omp atomic
	nr_access += 5.0*BLS*BLS; // optimal version - each matrix is retrieved only once
#endif

	  } // end single thread region with barrier at the end
/*kbw
		if(i_block==6 && k_block==5)
		  //if(sm_block==55)
	  //if(i_block==352)
	  {
		int i;
		printf("i_block %d, k_block %d, sm_block %d, Aux_ik after solution:\n",
		   i_block, k_block, sm_block);
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_ik_p+i]);
		}
		printf("\n");
		getchar(); getchar();
	  }
/*kew*/


	  // now k becomes an index for the k-th row block and i-th row block
	  // is subtracted from the k-th row block (with zero blocks ommitted!!!)

	  // we have to advance the pointer in k-th row until j-th column is found
	  // since j_block > k_block we start from the diagonal in k-th row
	  int smbl_kj = it_matrix->ilu_diag_ptr[k_block];

	  /* loop over all neighbors of i_block, with indices greater than k */
	  int j_smbl;
#pragma omp for
	  for(j_smbl = sm_block+1; j_smbl < it_matrix->ilu_row_ptr[i_block+1]; j_smbl++){

		int j_block = it_matrix->ilu_col_ind[j_smbl];

		// pointer to a_ik -> sm_block*BLS*BLS
		// a_ik now holds a_ik * a_kk^-1
		int a_ik_p = sm_block*BLS*BLS;

		// pointer to a_ij
		int a_ij_p = j_smbl*BLS*BLS;


/*kbw
			smbl_kj = it_matrix->ilu_diag_ptr[k_block];
		printf("\ntest starting position %d (%d)\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj]);

		// each thread advance its own counter for kj blocks
		while(it_matrix->ilu_col_ind[smbl_kj] < j_block  //) smbl_kj++;
		  &&  smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;
		//smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;

		printf("test finished position %d (%d), j_block %d, end row %d\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj], j_block,
		   it_matrix->ilu_row_ptr[k_block+1]);
		int test_sm=smbl_kj;

/*kew*/
		// we have to advance the pointer in k-th row until j-th column is found
		// since j_block > k_block we start from the diagonal in k-th row
		//smbl_kj = it_matrix->ilu_diag_ptr[k_block];
		// instead of resetting for each new j_block we can just go back
		// one block from where we finished for the previous j_block
		if(smbl_kj<2) smbl_kj--;

/*kbw
		printf("\nstarting position %d (%d)\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj]);
/*kew*/

		// each thread advance its own counter for kj blocks
		while(it_matrix->ilu_col_ind[smbl_kj] < j_block  //) smbl_kj++;
		  &&  smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;

/*kbw
		printf("finished position %d (%d), j_block %d, end row %d\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj], j_block,
		   it_matrix->ilu_row_ptr[k_block+1]);
		if(smbl_kj!=test_sm) {getchar();getchar();}
/*kew*/

		// if a_kj != 0 ; i.e. small block exist for index kj
		if(it_matrix->ilu_col_ind[smbl_kj]==j_block &&
		   smbl_kj < it_matrix->ilu_row_ptr[k_block+1]){

		  // pointer to a_kj
		  int a_kj_p = smbl_kj*BLS*BLS;

/*kbw
	  if(j_smbl==55)
		{
		  int i;
		  printf("\ni %d, j %d (small block %d, pointer %d), a_ij before solution:\n",
			 i_block, j_block, j_smbl, a_ij_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ij_p+i]);
		  }
		  printf("\n");
		  printf("i %d, k %d (small block %d, pointer %d), a_ik before solution:\n",
			 i_block, k_block, sm_block, a_ik_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ik_p+i]);
		  }
		  printf("\n");
		  printf("k %d, j %d (small block %d, pointer %d), a_kj before solution:\n",
			 k_block, j_block, smbl_kj, a_kj_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_kj_p+i]);
		  }
		  printf("\n");
		 getchar();
		}
/*kew*/

/* compute a_ij = a_ij - a_ik*a_kj */

// LAPACK VERSION
		  int bls = BLS;
		  double dmone = -1.0;
		  double done = 1.0;
		  char c = 'N';
		  dgemm_(&c, &c, &bls, &bls, &bls, &dmone,
				 &it_matrix->ilu_val[a_ik_p], &bls,
				 &it_matrix->ilu_val[a_kj_p], &bls, &done,
				 &it_matrix->ilu_val[a_ij_p], &bls);


		  // j,k,i VERSION

/* 	      // MUST BE VECTORIZED!!! */
/* 	      int i,j,k; */
/* 	      //double daux[BLS*BLS]; */
/* 	      //double eaux[BLS*BLS]; */
/* 	      //double faux[BLS*BLS]; */
/* 	      //for(i=0;i<BLS*BLS;i++){daux[i] = 0.0;} */
/* 	      //for(i=0;i<BLS*BLS;i++){eaux[i] = it_matrix->ilu_val[a_ik_p+i];} */
/* 	      //for(i=0;i<BLS*BLS;i++){faux[i] = it_matrix->ilu_val[a_kj_p+i];} */
/* 	      for(j=0;j<BLS;j++){ */
/* 		for(k=0;k<BLS;k++){ */
/* 		  double faux = it_matrix->ilu_val[a_kj_p+k+j*BLS]; */
/* 		  //#pragma vector always */
/* #pragma ivdep */
/* 		  for(i=0;i<BLS;i++){ */
/* 		    it_matrix->ilu_val[a_ij_p+i+j*BLS] -= */
/* 		    //daux[i+j*BLS] +=  */
/* 		      it_matrix->ilu_val[a_ik_p+i+BLS*k] * faux; */
/* 		      //eaux[i+BLS*k] * faux; */
/* 		      //it_matrix->ilu_val[a_kj_p+k+j*BLS]; */
/* 		      //faux[k+j*BLS]; */

/* 		  } */


/* 		} */
/* 	      } */

		  //for(i=0;i<BLS*BLS;i++){
		  //it_matrix->ilu_val[a_ij_p+i] -= daux[i];
		  //}


		  // i,j,k VERSION

		  /* // MUST BE VECTORIZED!!! */
		  /* for(i=0;i<BLS;i++){ */
		  /* 	for(j=0;j<BLS;j++){ */
		  /* 	  double daux = 0.0; */
		  /* 	  for(k=0;k<BLS;k++){ */
		  /* 	    daux += it_matrix->ilu_val[a_ik_p+i+BLS*k] * */
		  /* 	      it_matrix->ilu_val[a_kj_p+k+j*BLS]; */

		  /* 	  } */
		  /* 	  it_matrix->ilu_val[a_ij_p+i+j*BLS] -= daux; */

		  /* 	} */
		  /* } */

#ifdef TIME_TEST_ILU
#pragma omp atomic
	nr_oper += 2.0*BLS*BLS*BLS; // dgemm
#pragma omp atomic
	nr_access += 4.0*BLS*BLS; // optimal version - each matrix is retrieved only once
							  // for reading and a_ij once for writing
#endif



/*kbw
	  if(j_smbl==55)
		{
		  printf("\ni %d, j %d (small block %d, pointer %d), a_ij after solution:\n",
			 i_block, j_block, j_smbl, a_ij_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ij_p+i]);
		  }
		  printf("\n");
		  getchar();
		}
/*kew*/


		} // end if non-zero a_kj found for a_ij

	  } // end parallel loop over small blocks in k-th row block


	} // end loop over all row blocks (k_blocks) above i-th block

/* invert Dia for i_block */
	// pointer to a_ii -> it_matrix->ilu_diag_ptr[i_block]*BLS*BLS
	int a_ii_p = it_matrix->ilu_diag_ptr[i_block]*BLS*BLS;

/*kbw
	  if(i_block==55)
	  {
	int i;
	printf("\ni_block %d, Dia before inverting:\n",i_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");
	  }
/*kew*/

#pragma omp single
	{
	double a_work[(BLS+1)*BLS];
	int i;
	for(i=0;i<(BLS+1)*BLS;i++) a_work[i]=0.0;
	int ips_work[BLS*BLS];
	int bls = BLS;
/*kb*/
	int iaux;
	dgetrf_(&bls, &bls, &it_matrix->ilu_val[a_ii_p],
		&bls, ips_work, &iaux);

	int lwork = bls*bls;
	dgetri_(&bls, &it_matrix->ilu_val[a_ii_p],
		&bls, ips_work , a_work, &lwork, &iaux);
/*kb*/

#ifdef TIME_TEST_ILU
#pragma omp atomic
	nr_oper += (2.0/3.0+4.0/3.0)*BLS*BLS*BLS; // dgetrf+dgetri
#pragma omp atomic
	nr_access += 6.0*BLS*BLS; // optimal version - a_ii accessed twice in dgetrf
								  // and once in dgetri, a_work twice in dgetri (+init)
#endif

/*kbw
	  //if(i_block>235)
	  if(i_block==55)
	  {
	int i;
	printf("i_block %d, Dia after inverting:\n",i_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");

	//getchar();
	  }
/*kew*/

	} // end single thread region with implicit barrier at the end

	  } // end if not ghost block

	} // end outer loop over row blocks - i_block

	} // end parallel region

	}
	else if(version==2){

// VERSION II of ILU - PARALLEL FOR MULTIPLE ROWS

		char str4[] = "original row %d, original column %d, block %d\n";
		char str5[] = "not found symmetric block in row %d\n";
		char str6[] = "error 777111222333 testing ilu bcrs matrices\n";
#ifdef TIME_TEST_ILU
#pragma omp parallel default(none) firstprivate(it_matrix)  shared(nr_oper, nr_access,str4,str5,str6)
#else
#pragma omp parallel default(none) firstprivate(it_matrix) shared(str4,str5,str6)
#endif
	{

	  int i,j,k;

	  /* constants */
	  int ione = 1;
	  double done = 1.0;
	  double dmone = -1.0;


	  /* for each row block */
	  int i_block;
	  for(i_block=0; i_block < it_matrix->Nrblocks; i_block++){

		/* if block is active - necessary to allow for inverting diagonal small block */
		if(it_matrix->ilu_row_ptr[i_block]<it_matrix->ilu_row_ptr[i_block+1]){

/*kbw
  printf("\nin row block i_block %d, sm_blocks %d-%d, diag %d\n",
  i_block, it_matrix->ilu_row_ptr[i_block],
  it_matrix->ilu_row_ptr[i_block+1],
  it_matrix->ilu_diag_ptr[i_block]);
/*kew*/


// 1. Invert diagonal block and put to temporary storage - private to each thread!!!

		  /* invert Dia for i_block */
		  // pointer to a_ii -> it_matrix->ilu_diag_ptr[i_block]*BLS*BLS
		  int a_ii_p = it_matrix->ilu_diag_ptr[i_block]*BLS*BLS;

/*kbw
	  if(i_block==55)
	  //if(ib_dia==352)
	  {
	int i;
	printf("\ni_block %d, Dia before inverting:\n",i_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");
	  }
/*kew*/

	  double a_work[BLS*BLS];
	  double a_work_1[(BLS+1)*BLS];
	  int i;
	  for(i=0;i<(BLS+1)*BLS;i++) a_work_1[i]=0.0;
	  for(i=0;i<BLS*BLS;i++) a_work[i]=it_matrix->ilu_val[a_ii_p+i];
	  int ips_work[BLS*BLS];
	  int bls = BLS;
	  int iaux;
	  dgetrf_(&bls, &bls, a_work,
			  &bls, ips_work, &iaux);

	  int lwork = bls*bls;
	  dgetri_(&bls, a_work,
			  &bls, ips_work , a_work_1, &lwork, &iaux);


#ifdef TIME_TEST_ILU
#pragma omp atomic
	nr_oper += (2.0/3.0+4.0/3.0)*BLS*BLS*BLS; // dgetrf+dgetri
#pragma omp atomic
	nr_access += 7.0*BLS*BLS; // optimal version - a_ii accessed twice in dgetrf
								  // and once in dgetri, a_work twice in dgetri
#endif


	  // for all k>i - all rows below i-th row (for which small blocks a_ik exist)
	  // if a_ik is zero than  we do not have to subtract i-th row block
	  // from k-th row block (to zero a_ik)
	  // structural symmetry of stiffness matrix is used here
	  int sm_block;

	  // AUTOMATIC VERSION WITH NOWAIT CLAUSE AND EXPLICIT BARRIER AT THE END !!!
	  //#pragma omp for schedule(static,1) nowait
#pragma omp for schedule(dynamic,1) nowait
	  for(sm_block = it_matrix->ilu_diag_ptr[i_block]+1;
			  sm_block < it_matrix->ilu_row_ptr[i_block+1];
			  sm_block++){

		// MANUAL VERSION WITH EXPLICIT BARRIER AT THE END !!!
	  /* int my_id = omp_get_thread_num(); */
	  /* int nr_threads = omp_get_num_threads();	     */
	  /* for(sm_block = it_matrix->ilu_diag_ptr[i_block]+1+my_id; */
		  /*     sm_block < it_matrix->ilu_row_ptr[i_block+1];   */
		  /*     sm_block+=nr_threads){ */

			// find column index for k_block
			int k_block = it_matrix->ilu_col_ind[sm_block];
		//pointer to a_ki-th block
		int a_ki_p;
		int smbl_ki=-1;

		// for all small blocks in k-th row
		int smbl_k;
		for(smbl_k = it_matrix->ilu_row_ptr[k_block];
				smbl_k < it_matrix->ilu_row_ptr[k_block+1];
				smbl_k++ ){

			  // find ki-th block
		  if(it_matrix->ilu_col_ind[smbl_k] == i_block){

				// pointer to a_ki -> smbl_k*BLS*BLS
				a_ki_p = smbl_k*BLS*BLS;
		smbl_ki = smbl_k;

/*kbw
		printf("for original row %d, original column %d, block %d\n",
			   i_block, it_matrix->ilu_col_ind[sm_block], sm_block);
		printf("found symmetric row %d, symmetric column %d, block %d\n",
			   k_block, it_matrix->ilu_col_ind[smbl_ki], smbl_ki);

/*kew*/

		break;
			  }

			}

		if(smbl_ki == -1){

		  printf(str4,
			 i_block, it_matrix->ilu_col_ind[sm_block], sm_block);
		  printf(str5,
			 k_block);

		  printf(str6); exit(-1);

		}

		/* compute product a_ki = a_ki * a_ii^-1 */
		// new a_ki is used in calculations and left in lower part of LU
		// to be used again in forward reduction


/*kbw
		if(k_block==6 && i_block==5)
		//if(smbl_ki==55)
	  //if(i_block==352)
	  {
		int i;
		printf("\ni_block %d, sm_block %d, Dia_ii before solution:\n",
	i_block, it_matrix->ilu_diag_ptr[i_block]);
		for(i=0;i<BLS*BLS;i++){
	  //printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
		  printf("%25.15lf",a_work[i]);
		}
		printf("\n");
		printf("k_block %d, i_block %d, sm_block %d, Aux before solution:\n",
	k_block, i_block, smbl_ki);
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_ki_p+i]);
		}
		printf("\n");
	  }
/*kew*/

		// MUST BE VECTORIZED!!!
		double a_work_2[BLS*BLS];
		for(i=0;i<BLS;i++){
			  for(j=0;j<BLS;j++){
				double daux = 0.0;
		for(k=0;k<BLS;k++){

				  daux += it_matrix->ilu_val[a_ki_p+i+BLS*k] *
				  a_work[k+j*BLS];
		  //		  it_matrix->ilu_val[a_ii_p+k+j*BLS];

				}

		a_work_2[i+j*BLS]=daux;

			  }
			}
		for(i=0;i<BLS*BLS;i++) it_matrix->ilu_val[a_ki_p+i] = a_work_2[i];

#ifdef TIME_TEST_ILU
#pragma omp atomic
	nr_oper += 2.0*BLS*BLS*BLS; // dgemm
#pragma omp atomic
	nr_access += 5.0*BLS*BLS; // optimal version - each matrix is retrieved only once
#endif

/*kbw
		if(k_block==6 && i_block==5)
		  //if(smbl_ki==55)
	  {
		int i;
		printf("k_block %d, sm_block %d, Aux after solution:\n",k_block,smbl_ki);
		for(i=0;i<BLS*BLS;i++){
		  printf("%25.15lf",it_matrix->ilu_val[a_ki_p+i]);
		}
		printf("\n");
		getchar(); getchar();
	  }
/*kew*/


		// i-th row block (from (i+1)-th column)
		// is subtracted from the k-th row block (with zero blocks ommitted!!!)

		// starting small block in i-th row is
		//int smbl_ij = it_matrix->ilu_diag_ptr[i_block]+1;

		// starting small block in k-th row is:

		// safe choice - start from the beginning
		//int smbl_kj_start = it_matrix->ilu_row_ptr[k_block];
		// beter choice - just after ki-th block
		int smbl_kj_start = smbl_ki+1;

		/* loop over all neighbors of i_block, with indices greater than i */
		int smbl_kj = smbl_kj_start;
		int smbl_ij;
		for(smbl_ij = it_matrix->ilu_diag_ptr[i_block]+1;
				smbl_ij < it_matrix->ilu_row_ptr[i_block+1];
		smbl_ij++){

			  int j_block = it_matrix->ilu_col_ind[smbl_ij];

		  // pointer to a_ij
		  int a_ij_p = smbl_ij*BLS*BLS;


/*kbw
			smbl_kj = it_matrix->ilu_diag_ptr[k_block];
		printf("\ntest starting position %d (%d)\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj]);

		// each thread advance its own counter for kj blocks
		while(it_matrix->ilu_col_ind[smbl_kj] < j_block  //) smbl_kj++;
		  &&  smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;
		//smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;

		printf("test finished position %d (%d), j_block %d, end row %d\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj], j_block,
		   it_matrix->ilu_row_ptr[k_block+1]);
		int test_sm=smbl_kj;

/*kew*/
		  // we have to advance the pointer in k-th row until j-th column is found
		  // int smbl_kj = smbl_kj_start;
		  // one block from where we finished for the previous j_block
		  if(smbl_kj>smbl_kj_start) smbl_kj--;

/*kbw
		printf("\nstarting position %d (%d)\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj]);
/*kew*/
		  int jfound = 0;
		  // each thread advance its own counter for kj blocks
		  while(it_matrix->ilu_col_ind[smbl_kj] < j_block  //) smbl_kj++;
				&&  smbl_kj < it_matrix->ilu_row_ptr[k_block+1]) smbl_kj++;

/*kbw
		printf("finished position %d (%d), j_block %d, end row %d\n",
		   smbl_kj, it_matrix->ilu_col_ind[smbl_kj], j_block,
		   it_matrix->ilu_row_ptr[k_block+1]);
		if(smbl_kj!=test_sm) {getchar();getchar();}
/*kew*/

		  // if a_kj != 0 ; i.e. small block exist for index kj
		  if(it_matrix->ilu_col_ind[smbl_kj]==j_block &&
			 smbl_kj < it_matrix->ilu_row_ptr[k_block+1]){

		jfound = 1;

/*kbw
		printf("for i_row %d, column %d, small_block %d\n",
			   i_block, j_block, smbl_ij);
		printf("\tfound k_row %d, column %d, small_block %d\n",
			   k_block, it_matrix->ilu_col_ind[smbl_kj], smbl_kj);
		//getchar();getchar();
/*kew*/

				// pointer to a_kj
				int a_kj_p = smbl_kj*BLS*BLS;

/*kbw
	  if(smbl_kj==55)
	{
		  printf("k %d, j %d (small block %d, pointer %d), a_kj before solution:\n",
			 k_block, j_block, smbl_kj, a_kj_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_kj_p+i]);
		  }
		  printf("\n");
		  printf("k %d, i %d (small block %d, pointer %d), a_ki before solution:\n",
			 k_block, i_block, smbl_ki, a_ki_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ki_p+i]);
		  }
		  printf("\n");
		  printf("\ni %d, j %d (small block %d, pointer %d), a_ij before solution:\n",
			 i_block, j_block, smbl_ij, a_ij_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_ij_p+i]);
		  }
		  printf("\n");
		  getchar();
	}
/*kew*/

/* compute a_kj = a_kj - a_ki*a_ij */

// LAPACK VERSION
		int bls = BLS;
		double dmone = -1.0;
		double done = 1.0;
		char c = 'N';
		dgemm_(&c, &c, &bls, &bls, &bls, &dmone,
				   &it_matrix->ilu_val[a_ki_p], &bls,
				   &it_matrix->ilu_val[a_ij_p], &bls, &done,
				   &it_matrix->ilu_val[a_kj_p], &bls);


/* compute a_ij = a_ij - a_ik*a_kj */
// LAPACK VERSION
		  /* int bls = BLS; */
		  /* double dmone = -1.0; */
		  /* double done = 1.0; */
		  /* char c = 'N'; */
		  /* dgemm_(&c, &c, &bls, &bls, &bls, &dmone,  */
		  /* 	     &it_matrix->ilu_val[a_ik_p], &bls, */
		  /* 	     &it_matrix->ilu_val[a_kj_p], &bls, &done, */
		  /* 	     &it_matrix->ilu_val[a_ij_p], &bls); */


		  // j,k,i VERSION

/* 	      // MUST BE VECTORIZED!!! */
/* 	      int i,j,k; */
/* 	      //double daux[BLS*BLS]; */
/* 	      //double eaux[BLS*BLS]; */
/* 	      //double faux[BLS*BLS]; */
/* 	      //for(i=0;i<BLS*BLS;i++){daux[i] = 0.0;} */
/* 	      //for(i=0;i<BLS*BLS;i++){eaux[i] = it_matrix->ilu_val[a_ik_p+i];} */
/* 	      //for(i=0;i<BLS*BLS;i++){faux[i] = it_matrix->ilu_val[a_kj_p+i];} */
/* 	      for(j=0;j<BLS;j++){ */
/* 		for(k=0;k<BLS;k++){ */
/* 		  double faux = it_matrix->ilu_val[a_kj_p+k+j*BLS]; */
/* 		  //#pragma vector always */
/* #pragma ivdep */
/* 		  for(i=0;i<BLS;i++){ */
/* 		    it_matrix->ilu_val[a_ij_p+i+j*BLS] -= */
/* 		    //daux[i+j*BLS] +=  */
/* 		      it_matrix->ilu_val[a_ik_p+i+BLS*k] * faux; */
/* 		      //eaux[i+BLS*k] * faux; */
/* 		      //it_matrix->ilu_val[a_kj_p+k+j*BLS]; */
/* 		      //faux[k+j*BLS]; */

/* 		  } */


/* 		} */
/* 	      } */

		  //for(i=0;i<BLS*BLS;i++){
		  //it_matrix->ilu_val[a_ij_p+i] -= daux[i];
		  //}


		  // i,j,k VERSION

		  /* // MUST BE VECTORIZED!!! */
		  /* int i,j,k; */
		  /* for(i=0;i<BLS;i++){ */
		  /* 	for(j=0;j<BLS;j++){ */
		  /* 	  double daux = 0.0; */
		  /* 	  for(k=0;k<BLS;k++){ */
		  /* 	    daux += it_matrix->ilu_val[a_ik_p+i+BLS*k] * */
		  /* 	      it_matrix->ilu_val[a_kj_p+k+j*BLS]; */

		  /* 	  } */
		  /* 	  it_matrix->ilu_val[a_ij_p+i+j*BLS] -= daux; */

		  /* 	} */
		  /* } */


#ifdef TIME_TEST_ILU
#pragma omp atomic
	nr_oper += 2.0*BLS*BLS*BLS; // dgemm
#pragma omp atomic
	nr_access += 4.0*BLS*BLS; // optimal version - each matrix is retrieved only once
							  // for reading and a_jk once for writing
#endif


/*kbw
	  if(smbl_kj==55)
	{
		  printf("\nk %d, j %d (small block %d, pointer %d), a_ij after solution:\n",
			 k_block, j_block, smbl_kj, a_kj_p);
		  for(i=0;i<BLS*BLS;i++){
		printf("%25.15lf", it_matrix->ilu_val[a_kj_p+i]);
		  }
		  printf("\n");
		  getchar();
	}
/*kew*/


			  } // end if non-zero a_kj found for a_ij
/*kbw
		  else {
		printf("for i_row %d, column %d, small_block %d\n",
			   i_block, j_block, smbl_ij);
		printf("\t\tnot found small block in k_row %d\n",
			   k_block);
		//getchar();getchar();
		  }
/*kew*/

		} // end loop over small blocks in k-th row block



		  } // end loop over all row blocks (k_blocks) below i-th block

	  // parallel for has nowait clause - hence we put implicit barrier here
#pragma omp barrier

	  // one thread has to rewrite a_work to a_ii^-1
#pragma omp single nowait
	  {
		// a_ii holds original a_ii inverted
		for(i=0; i<BLS*BLS; i++){
			  it_matrix->ilu_val[a_ii_p+i] = a_work[i];
			}

#ifdef TIME_TEST_ILU
#pragma omp atomic
	nr_access += 2.0*BLS*BLS;
#endif

/*kbw
	  if(i_block==55)
	  //if(i_block>235)
	  {
	int i;
	printf("i_block %d, Dia after inverting:\n",i_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");

	//getchar();
	  }
/*kew*/
		  } // end single thread task for rewriting a_ii^-1



		} // end if not ghost block

	  } // end outer loop over row blocks - i_block

//#pragma omp barrier

	} // end parallel region

	}


/*kbw
	  printf("AFTER ILU DECOMPOSITION:\n");
	  int ilu_row_gl;
	  for(ilu_row_gl=0; ilu_row_gl < it_matrix->Nrblocks; ilu_row_gl++){

		int sm_block;
		for(sm_block = it_matrix->ilu_row_ptr[ ilu_row_gl ];
			sm_block < it_matrix->ilu_row_ptr[ ilu_row_gl+1 ];
			sm_block++){

		  int ilu_col_gl = it_matrix->ilu_col_ind[ sm_block ];

		  int ilu_val_ind = sm_block*BLS*BLS;

		  printf("row %d, column %d, small_block %d:\n",
			 ilu_row_gl, ilu_col_gl, sm_block);

		  int i;
		  for(i=0;i<BLS*BLS;i++){
			printf("%25.15lf", it_matrix->ilu_val[ilu_val_ind+i]);
		  }
		  printf("\n");
	  }
	  getchar();
	}
/*kew*/


  }
  else{

	printf("not implemented preconditioner %d in BCRS module\n",
	   it_matrix->Precon);
	exit(-1);

  }

#ifdef TIME_TEST_ILU
  exec_time = time_clock() - exec_time;
  printf("Filling preconditioner: nr_oper %lf, nr_access %lf\n",
		  nr_oper, nr_access);
  printf("\t\t time %lf -> %lf GFlops, %lf GB/s\n",
	 exec_time, nr_oper/1.e9/exec_time, nr_access*sizeof(double)/1.e9/exec_time);
#endif

  return(0);

}



/*---------------------------------------------------------
  lar_free_preconditioner_bcrs - to free space for a block structure
---------------------------------------------------------*/
int lar_free_preconditioner_bcrs(
  int Matrix_id   /* in: matrix ID */
  )
{

  //printf("\n\n BCRS-lar_free_preconditioner \n\n");

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];

  if(it_matrix->Precon==BLOCK_JACOBI || it_matrix->Precon==BLOCK_GS){

	free(it_matrix->block_diag_ptr);
	free(it_matrix->block_diag_precon);
	free(it_matrix->block_diag_ips);

  } // end if BLOCK_JACOBI or BLOCK_GS preconditioner
  else if( it_matrix->Precon==MULTI_ILU || it_matrix->Precon==BLOCK_ILU ){

	free(it_matrix->ilu_diag_ptr);
	free(it_matrix->ilu_val);
	free(it_matrix->ilu_row_ptr);
	free(it_matrix->ilu_col_ind);

  } // end if MULTI_ILU or BLOCK_ILU preconditioner
  else{

	printf("not implemented preconditioner %d in BCRS module\n",
	   it_matrix->Precon);

  }

  return(0);
}

/*---------------------------------------------------------
  lar_free_SM_and_LV_bcrs - to free space for a block structure
---------------------------------------------------------*/
int lar_free_SM_and_LV_bcrs(
  int Matrix_id   /* in: matrix ID */
				)
{

  // to keep matrix management as simple as possible it is required that
  // the matrices are freed in the reverse order wrt creating (as in heap)
  if(Matrix_id != itv_nr_bcrs_matrices-1){
	printf("Matrix destroyed %d (from %d) is not the last in lar_free_SM_and_LV_bcrs!!!\n",
	   Matrix_id, itv_nr_bcrs_matrices);
	exit(-1);
  }


  //printf("\n\n BCRS-free SM \n\n");

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];

  free(it_matrix->block_val);
  free(it_matrix->block_row_ptr);
  free(it_matrix->block_col_ind);
  free(it_matrix->rhs);

  itv_nr_bcrs_matrices--;

  //printf("\n\n BCRS-free-end \n\n");

  return(0);
}

/*---------------------------------------------------------
lar_compute_residual_bcrs - to compute the residual of the system of equations,
	v = ( b - Ax ) (used also to compute the product v = -Ax)
---------------------------------------------------------*/
void lar_compute_residual_bcrs (
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

  int i;

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];
  int nrbl_inter=it_matrix->Nrbl_internal;

/*kbw
  printf("X on entrance\n" );
  for(i=0;i<40;i++) printf("%20.15lf",X[i]);
  printf("\n" );
  getchar();
  getchar();
/*kew*/

  //printf("\n\n BCRS-lar_compute_residual %d %d\n\n",Ndof, it_matrix->Nrdofgl);//getchar();

  assert(Ndof == it_matrix->Nrdofgl);


  //#pragma omp  for
#pragma omp parallel for default(none) firstprivate(Ndof, V)
  for(i=0;i<Ndof;i++) V[i]=0.0;

  if(!Ini_zero){ // if X non-zero

	int block_row_gl;
	//#pragma omp for
#pragma omp parallel for default(none) firstprivate(Ndof, V, X, B, it_matrix, Ini_zero, Use_rhs, nrbl_inter)
	for(block_row_gl=0; block_row_gl < nrbl_inter; block_row_gl++){

	  int sm_block;
	  for(sm_block = it_matrix->block_row_ptr[ block_row_gl ];
	  sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];
	  sm_block++){

	int block_col_gl = it_matrix->block_col_ind[ sm_block ];

	int block_val_ind = sm_block*BLS*BLS;

/*kbw
	int nrdofbl=BLS;
	printf("compres for blocks %d, %d before: X - \n", block_row_gl, block_col_gl);
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",X[ block_col_gl*BLS+i]);
	printf("\n");
	printf("compres before: V - \n");
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",V[ block_row_gl*BLS+i]);
	printf("\n");
	printf("compres before: SM - \n");
	for(i=0;i<nrdofbl*nrdofbl;i++) printf("%20.15lf",it_matrix->block_val[ block_val_ind+i]);
	printf("\n");
/*kew*/

	int k,l;

	for(l=0; l<BLS; l++){
	  for(k=0; k<BLS; k++){
		// small blocks are stored by columns
		V[ block_row_gl*BLS + k] -=
		  it_matrix->block_val[ block_val_ind + k + l*BLS] * X[ block_col_gl*BLS + l];

		// small blocks are stored by rows
	  //it_matrix->block_val[ block_val_ind + l + k*BLS] * X[ block_col_gl*BLS + l];
	  }
	}
/*kbw
	printf("compres after: V - \n");
	for(i=0;i<BLS;i++) printf("%20.15lf",V[block_row_gl*BLS+i]);
	printf("\n");
	getchar();
	getchar();
/*kew*/

	  } // end loop over horizontal blocks
	} // end loop over row blocks

/*kbw
  printf("X after mat-vec\n" );
  for(i=0;i<40;i++) printf("%20.15lf",V[i]);
  printf("\n" );
  getchar();
  getchar();
/*kew*/

  } // end if X non-zero

  if(Use_rhs==1){

/*kbw
	printf("compres RHS or B on entrance\n");
	if(B==NULL){
	  for(i=0;i<40;i++) printf("R%19.15lf",it_matrix->rhs[i]);
	}else{
	  for(i=0;i<40;i++) printf("B%19.15lf",B[i]);
	}
	printf("\n");
	getchar();
	getchar();
/*kew*/

	if(B==NULL){

	  //#pragma omp  for
#pragma omp parallel for default(none) firstprivate(Ndof, V, it_matrix)
	  for(i=0;i<Ndof;i++) V[i] += it_matrix->rhs[i];

	  //daxpy_(&Ndof, &done, it_matrix->rhs, &ione, V, &ione);

	}
	else{

	  //#pragma omp for
#pragma omp parallel for default(none) firstprivate(Ndof, V, B, it_matrix)
	  for(i=0;i<Ndof;i++) V[i] += B[i];

	  //daxpy_(&Ndof, &done, B, &ione, V, &ione);

	}
  }

/*kbw
  printf("V on leaving\n" );
  for(i=0;i<40;i++) printf("%20.15lf",V[i]);
  printf("\n" );
  getchar();
  getchar();
/*kew*/

/*     printf("compres: subsequent rows - \n"); */
/*     for(i=0;i<Ndof;i++) printf("%20.15lf",V[i]); */
/*     printf("\n"); */
/* getchar(); */

  return;
}


/*---------------------------------------------------------
lar_perform_BJ_or_GS_iterations_bcrs - to perform one iteration of block Gauss-Seidel
	algorithm (block Jacobi switched on by it_matrix->Precon==BLOCK_JACOBI)
		v_out = v_in + M^-1 * ( b - A * v_in )
---------------------------------------------------------*/
void lar_perform_BJ_or_GS_iterations_bcrs(
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

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];
  int nrbl_inter=it_matrix->Nrbl_internal;

  //int *block_diag_ips= (int*)malloc( (it_matrix->Nrblocks)*sizeof(int));

  int i_loop;		/* counter for loops for GS */
  int blsl=BLS,iaux;
  int ione = 1;
  double done = 1.;
  int i,ret,index;

  int block_row_gl, sm_block, block_col_gl,block_val_ind;

/*kbw
if(pcv_my_proc_id==1){
  printf("\n\n BCRS-lar_perform_BJ_or_GS_iterations Nr_prec=%d Ndof=%d Nrblocks=%d\n\n",Nr_prec,Ndof, it_matrix->Nrblocks);
  for(i=0;i<Ndof;i++) printf("(%.12lf)",V[i]);printf("\n\n");
  for(i=0;i<Ndof;i++) printf("w(%.12lf)",it_matrix->rhs[i]);printf("\n\n");
  getchar();
}
/**/


/*kbw
//if(pcv_my_proc_id==1){
#define SAVE_LOG_TO_FILE
#ifdef SAVE_LOG_TO_FILE
  char filename[256];
  const char *program_name = getenv("_");

  sprintf(filename,"BCRS_Solution.txt",program_name);

  static int file_bcrs_opened_by_program = 0;
  FILE *file = NULL;
  if(file_bcrs_opened_by_program == 0) {
	file = fopen(filename,"w+");
	file_bcrs_opened_by_program = 1;
  } else {
	file = fopen(filename,"a+");
  }

#else
  FILE *file = stdout;
#endif

  int row_ind,k,l;
  fprintf(file,"\n\n BCRS on entry: row, column, value\n");
  for(i=0;i<nrbl_inter;++i){
	for(l=0;l<blsl;l++){
	  for(row_ind=it_matrix->block_row_ptr[i];row_ind<it_matrix->block_row_ptr[i+1];row_ind++){
		for(k=0;k<blsl;k++){
		  fprintf(file,"%6d%6d%20.10lf\n",i*BLS+l,it_matrix->block_col_ind[ row_ind ]*BLS+k,it_matrix->block_val[row_ind*BLS*BLS+l+k*BLS]);
		}
	  }
	}
  }

  fprintf(file,"\n RHS on entry: index, value\n");
  if(B == NULL) {
	for(i=0;i<Ndof;++i){
	  fprintf(file,"%6d%20.10lf\n",i,it_matrix->rhs[i]);
	}
  } else {
	for(i=0;i<Ndof;++i){
	  fprintf(file,"%6d%20.10lf\n",i,B[i]);
	}
  }

  fprintf(file,"\n SOLUTION [before solve] on entry: index, value\n");
  for(i=0;i<Ndof;++i) {
	fprintf(file,"%6d%20.10lf\n",i,V[i]);
  }
#ifdef SAVE_LOG_TO_FILE
  fclose(file);
#endif
  // }

  // Exit - uncomment if needed for testing
  //exit(0);
/*kew*/

  assert(Ndof==it_matrix->Nrblocks*BLS);
  assert(Ndof==it_matrix->Nrdofgl);

  int nr_prec2;
  if(Nr_prec%2==0) nr_prec2=Nr_prec/2;
  else nr_prec2=Nr_prec/2+1;

  /* loop over loops */
  for(i_loop=0;i_loop<Nr_prec;i_loop++){

	  char N_string[] = "N";
	/* loop over all blocks - in standard order */
#pragma omp parallel default(none)					\
  firstprivate(V,B,it_matrix,i_loop,Use_rhs,Ini_zero,Ndof, blsl,ione,iaux,done,ret) \
  private(block_row_gl, sm_block, block_col_gl, block_val_ind, i) shared(N_string, nrbl_inter)
	{

	  // each thread gets its copy of input vector
	  double *vtemp_GS_BJ;		/* temporary vector*/
	  vtemp_GS_BJ = (double*)malloc(Ndof*sizeof(double));

#ifdef _OPENMP
	  int my_id = omp_get_thread_num();
	  int portion = (nrbl_inter*BLS)/omp_get_num_threads() + 1;
#else
	  int my_id = 0;
	  int portion = (nrbl_inter*BLS);
#endif

	  int my_first_x = (my_id*portion - it_matrix->Half_bandwidth - 10)*BLS;
	  if(my_first_x<0) my_first_x = 0;
	  int my_last_x = ((my_id+1)*portion + it_matrix->Half_bandwidth + 10)*BLS;
	  if(my_last_x>(nrbl_inter*BLS)) my_last_x=(nrbl_inter*BLS);
	  //for(i=my_first_x; i<my_last_x; i++) vtemp_GS_BJ[i]=V[i];
for(i=0; i<Ndof; i++)vtemp_GS_BJ[i]=V[i];

	  int my_first = my_id*portion;
	  if(my_first<0) my_first = 0;
	  int my_last = (my_id+1)*portion;
	  if(my_last>nrbl_inter) my_last=nrbl_inter;



	  //first forward loop

	  //OLD #pragma omp for
	  //OLD for (index=0; index < it_matrix->Nrblocks; index++){
	  for (block_row_gl = my_first; block_row_gl < my_last; block_row_gl++){

	//OLD if(i_loop%2==0)	block_row_gl = index;
	//OLD else block_row_gl = it_matrix->Nrblocks-1-index;

	double vloc[BLS];
	for(i=0;i<BLS;i++) vloc[i]=0.0;

	//Lx
	for(sm_block = it_matrix->block_row_ptr[block_row_gl];
		sm_block < it_matrix->block_diag_ptr[block_row_gl];
		sm_block++){

	  block_col_gl = it_matrix->block_col_ind[ sm_block ];

	  block_val_ind = sm_block*BLS*BLS;

	  assert(block_row_gl != block_col_gl);
	  int k,l;

	  if(it_matrix->Precon==BLOCK_JACOBI){

		for(l=0; l<BLS; l++){
		  for(k=0; k<BLS; k++){

		// small blocks stored by columns
		vloc[k] -= it_matrix->block_val[ block_val_ind + k + l*BLS] * V[ block_col_gl*BLS + l];
		//if( block_col_gl*BLS + l>1880)printf("(%d * %d) ",(block_val_ind + k*BLS + l),(block_col_gl*BLS + l));



/*kbw
	if(pcv_my_proc_id==1)
			{
//#pragma omp critical(printing)
		 {
			  printf("%6d %20.10lf %20.10lf\n",
				 block_col_gl*BLS + l,
				 it_matrix->block_val[ block_val_ind+k+l*BLS],
				 V[ block_col_gl*BLS + l]);
			  }
			  }
/**/


		  }
		}

	  }
	  else{

		for(l=0; l<BLS; l++){
		  for(k=0; k<BLS; k++){
		// small blocks stored by columns
		vloc[k] -= it_matrix->block_val[ block_val_ind+k+l*BLS] * vtemp_GS_BJ[ block_col_gl*BLS+l];
		//printf("(%d * %d, %20.10lf->%20.10lf) ",(block_val_ind+k*BLS+l),(block_col_gl*BLS + l),
		//vloc[k],vtemp_GS_BJ[ block_col_gl*BLS + l]);

		  }
		}

	  }
	}

	//Ux
	for(sm_block = it_matrix->block_diag_ptr[block_row_gl]+1;
		sm_block < it_matrix->block_row_ptr[block_row_gl+1];
		sm_block++){

	  block_col_gl = it_matrix->block_col_ind[ sm_block ];
	  block_val_ind = sm_block*BLS*BLS;

	  //assert(block_row_gl != block_col_gl);
	  int k,l;

	  for(k=0; k<BLS; k++){
		for(l=0; l<BLS; l++){
		  // small blocks are stored by columns
		  vloc[k] -= it_matrix->block_val[ block_val_ind + k + l*BLS] * V[ block_col_gl*BLS + l];
		  //printf("(%d * %d) ",(block_val_ind + k*BLS + l),(block_col_gl*BLS + l));
		//printf("%d k=%d = %20.10lf %20.10lf %20.10lf\n",block_col_gl*BLS + l,k,vloc[k],it_matrix->block_val[ block_val_ind + k + l*BLS],V[ block_col_gl*BLS + l]);
		//getchar();

/*kbw
	if(pcv_my_proc_id==1)
			{
//#pragma omp critical(printing)
		 {
			  printf("%6d %20.10lf %20.10lf\n",
				 block_col_gl*BLS + l,
				 it_matrix->block_val[ block_val_ind+k+l*BLS],
				 V[ block_col_gl*BLS + l]);
			  }
			  }
/**/

		}
	  }
	}

	// now RHS: b - Ax
	if(Use_rhs==1){

	  if(B==NULL){

		daxpy_(&blsl, &done, &it_matrix->rhs[block_row_gl*BLS], &ione, vloc, &ione);
		//vtemp_GS_BJ[block_row_gl]+=it_matrix->rhs[block_row_gl];
		// printf("RHS for subtraction\n");
		// for(i=0;i<BLS;i++) printf("%20.10lf", vloc[i]);
		// printf("\n");

	  }
	  else{

		daxpy_(&blsl, &done, &B[block_row_gl*BLS], &ione, vloc, &ione);
		//vtemp_GS_BJ[block_row_gl]+=B[block_row_gl];

	  }
	}

	// and finally apply preconditioner
	dgetrs_(N_string,&blsl,&ione,&it_matrix->block_diag_precon[block_row_gl*BLS*BLS],&blsl,
		&it_matrix->block_diag_ips[block_row_gl*BLS*BLS],vloc,&blsl,&iaux);


	//printf("\n (%6d,%20.10lf,%20.10lf) \n",block_row_gl,vloc[0],vloc[1]);

	// rewrite block to temporary vector
	int k;
	for(k=0; k<BLS; k++) vtemp_GS_BJ[block_row_gl*BLS+k]=vloc[k];
	//printf("\n (%6d,%20.10lf,%20.10lf) \n",(block_row_gl*BLS),vloc[0],vtemp_GS_BJ[block_row_gl]);

	  } // end loop over row blocks

#pragma omp barrier


	  //OLD#pragma omp for
	  //OLDfor (index=0; index < it_matrix->Nrblocks; index++){

	  for (block_row_gl = my_first; block_row_gl < my_last; block_row_gl++){

	//OLDif(i_loop%2==0)	block_row_gl = index;
	//OLDelse block_row_gl = it_matrix->Nrblocks-1-index;

	int k;
	for(k=0; k<BLS; k++)
	  V[block_row_gl*BLS+k]=vtemp_GS_BJ[block_row_gl*BLS+k];

	  }
	  free(vtemp_GS_BJ);

	  //for(i=0;i<Ndof;i++) V[i]=vtemp_GS_BJ[i];

	} // end parallel region

	if(Nr_prec%2==0){

	/* loop over all blocks - in reverse order*/
#pragma omp parallel default(none)					\
  firstprivate(V,B,it_matrix,i_loop,Use_rhs,Ini_zero,Ndof, blsl,ione,iaux,done,ret) \
  private(block_row_gl, sm_block, block_col_gl, block_val_ind, i) shared(N_string,nrbl_inter)
	 {

	  // each thread gets its copy of input vector
	  double *vtemp_GS_BJ;		/* temporary vector*/
	  vtemp_GS_BJ = (double*)malloc(Ndof*sizeof(double));

#ifdef _OPENMP
	  int my_id = omp_get_thread_num();
	  int portion = (nrbl_inter*BLS)/omp_get_num_threads() + 1;
#else
	  int my_id = 0;
	  int portion = (nrbl_inter*BLS);
#endif

	  int my_first_x = (my_id*portion - it_matrix->Half_bandwidth - 10)*BLS;
	  if(my_first_x<0) my_first_x = 0;
	  int my_last_x = ((my_id+1)*portion + it_matrix->Half_bandwidth + 10)*BLS;
	  if(my_last_x>(nrbl_inter*BLS)) my_last_x=(nrbl_inter*BLS);

	  //for(i=my_first_x; i<my_last_x; i++) vtemp_GS_BJ[i]=V[i];
	  for(i=0; i<Ndof; i++)vtemp_GS_BJ[i]=V[i];

	  int my_first = my_id*portion;
	  if(my_first<0) my_first = 0;
	  int my_last = (my_id+1)*portion;
	  if(my_last>nrbl_inter) my_last=nrbl_inter;

	  //now backward loop

	  //OLD #pragma omp for
	  //OLD for (index=0; index < it_matrix->Nrblocks; index++){
	  for (block_row_gl = my_last-1; block_row_gl >= my_first; block_row_gl--){

	//OLD if(i_loop%2==0)	block_row_gl = index;
	//OLD else block_row_gl = it_matrix->Nrblocks-1-index;

	double vloc[BLS];
	for(i=0;i<BLS;i++) vloc[i]=0.0;

	//Ux - now Upper part is first and split according to type of preconditioner
	for(sm_block = it_matrix->block_diag_ptr[block_row_gl]+1;
		sm_block < it_matrix->block_row_ptr[block_row_gl+1];
		sm_block++){

	  block_col_gl = it_matrix->block_col_ind[ sm_block ];
	  block_val_ind = sm_block*BLS*BLS;

	  assert(block_row_gl != block_col_gl);
	  int k,l;

	  if(it_matrix->Precon==BLOCK_JACOBI){

		for(l=0; l<BLS; l++){
		  for(k=0; k<BLS; k++){

		// small blocks stored by columns
		vloc[k] -= it_matrix->block_val[ block_val_ind + k + l*BLS] * V[ block_col_gl*BLS + l];
		//printf("(%d * %d) ",(block_val_ind + k*BLS + l),(block_col_gl*BLS + l));


/*kbw
	if(pcv_my_proc_id==1)
			{
//#pragma omp critical(printing)
		 {
			  printf("%6d %20.10lf %20.10lf\n",
				 block_col_gl*BLS + l,
				 it_matrix->block_val[ block_val_ind+k+l*BLS],
				 V[ block_col_gl*BLS + l]);
			  }
			  }
/**/


		  }
		}

	  }
	  else{

		for(l=0; l<BLS; l++){
		  for(k=0; k<BLS; k++){

		vloc[k] -= it_matrix->block_val[ block_val_ind+k+l*BLS] * vtemp_GS_BJ[ block_col_gl*BLS+l];
		//printf("(%d * %d, %20.10lf->%20.10lf) ",(block_val_ind+k*BLS+l),(block_col_gl*BLS + l),
		//vloc[k],vtemp_GS_BJ[ block_col_gl*BLS + l]);
		  }
		}

	  }

	}

	//Lx - lower part - now always with previous iteration vector
	for(sm_block = it_matrix->block_row_ptr[block_row_gl];
		sm_block < it_matrix->block_diag_ptr[block_row_gl];
		sm_block++){

	  block_col_gl = it_matrix->block_col_ind[ sm_block ];

	  block_val_ind = sm_block*BLS*BLS;

	  assert(block_row_gl != block_col_gl);
	  int k,l;

	  for(k=0; k<BLS; k++){
		for(l=0; l<BLS; l++){

		  vloc[k] -= it_matrix->block_val[ block_val_ind + k + l*BLS] * V[ block_col_gl*BLS + l];
		  //printf("(%d * %d) ",(block_val_ind + k*BLS + l),(block_col_gl*BLS + l));

/*kbw
	if(pcv_my_proc_id==1)
			{
//#pragma omp critical(printing)
		 {
			  printf("%6d %20.10lf %20.10lf\n",
				 block_col_gl*BLS + l,
				 it_matrix->block_val[ block_val_ind+k+l*BLS],
				 V[ block_col_gl*BLS + l]);
			  }
			  }
/**/

		}
	  }

	}

	// now RHS: b - Ax
	if(Use_rhs==1){

	  if(B==NULL){

		daxpy_(&blsl, &done, &it_matrix->rhs[block_row_gl*BLS], &ione, vloc, &ione);
		//vtemp_GS_BJ[block_row_gl]+=it_matrix->rhs[block_row_gl];
		// printf("RHS for subtraction\n");
		// for(i=0;i<BLS;i++) printf("%20.10lf", vloc[i]);
		// printf("\n");

	  }
	  else{

		daxpy_(&blsl, &done, &B[block_row_gl*BLS], &ione, vloc, &ione);
		//vtemp_GS_BJ[block_row_gl]+=B[block_row_gl];

	  }
	}

	// and finally apply preconditioner
	dgetrs_(N_string,&blsl,&ione,&it_matrix->block_diag_precon[block_row_gl*BLS*BLS],&blsl,
		&it_matrix->block_diag_ips[block_row_gl*BLS*BLS],vloc,&blsl,&iaux);

	//printf("\n (%6d,%20.10lf,%20.10lf) \n",block_row_gl,vloc[0],vloc[1]);

	// rewrite block to temporary vector
	int k;
	for(k=0; k<BLS; k++) vtemp_GS_BJ[block_row_gl*BLS+k]=vloc[k];
	//printf("\n (%6d,%20.10lf,%20.10lf) \n",(block_row_gl*BLS),vloc[0],vtemp_GS_BJ[block_row_gl]);

	  } // end loop over row blocks - in reverse order

#pragma omp barrier

	  //OLD#pragma omp for
	  //OLDfor (index=0; index < it_matrix->Nrblocks; index++){

	  for (block_row_gl = my_first; block_row_gl < my_last; block_row_gl++){

	//OLDif(i_loop%2==0)	block_row_gl = index;
	//OLDelse block_row_gl = it_matrix->Nrblocks-1-index;

	int k;
	for(k=0; k<BLS; k++)
	  V[block_row_gl*BLS+k]=vtemp_GS_BJ[block_row_gl*BLS+k];

	  }
	  free(vtemp_GS_BJ);

	  //for(i=0;i<Ndof;i++) V[i]=vtemp_GS_BJ[i];

	 } // end parallel region

	} // end if second loop over row blocks required - the one in reverse order

  } // end loop over preconditioner iterations

/*kbw
if(pcv_my_proc_id==1){
  printf("X on exit\n\n");
  for(i=0;i<Ndof;i++) printf("(%.12lf) ",V[i]);printf("\n\n");
}
/**/
  return;
}


/*---------------------------------------------------------
  lar_perform_rhsub_bcrs - to perform forward reduction and back-substitution for ILU
						   preconditioning: v_out = M^-1 * b
---------------------------------------------------------*/
void lar_perform_rhsub_bcrs(
  int Matrix_id,   /* in: matrix ID */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* out: vector of unknowns updated */
		/*      during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	)
{

#ifdef TIME_TEST_ILU
  double nr_oper = 0;
  double nr_access = 0;
  double exec_time = time_clock();

#endif

  //printf("\n\n BCRS-lar_perform_rhsub \n\n");

  itt_bcrs_matrices *it_matrix;
  it_matrix = &itv_bcrs_matrices[Matrix_id];

/* FORWARD REDUCTION */

  int version=0;
  //printf("Starting FW and BS - choose the version (current %d)\n", version);
  //scanf("%d", &version);


  //VERSION 0 - sequential (once was parallel in row)
  if(version==0){


	int i,j,k;

	/* loop over all row_blocks */
	int i_block;
	for(i_block=0; i_block < it_matrix->Nrblocks; i_block++){

	  /* if block is active */
	  if(it_matrix->ilu_row_ptr[i_block]<it_matrix->ilu_row_ptr[i_block+1]){

	/* initialize vloc */
	double vloc[BLS];
	/* substitute value of rhs to vloc */
	if(B==NULL){
	  int i;
	  for(i=0;i<BLS;i++) vloc[i]=it_matrix->rhs[i_block*BLS+i];
	}
	else{
	  int i;
	  for(i=0;i<BLS;i++) vloc[i]=B[i_block*BLS+i];
	}

	int sm_block;
	// for all j<i - lower subdiagonal entries
	for(sm_block = it_matrix->ilu_row_ptr[i_block];
		sm_block < it_matrix->ilu_diag_ptr[i_block];
		sm_block++){

	  int j_block = it_matrix->ilu_col_ind[sm_block];

	  // pointer to a_ij -> sm_block*BLS*BLS
	  int a_ij_p = sm_block*BLS*BLS;

	  // MUST BE VECTORIZED!!!
	  int i,j;
	  for(j=0;j<BLS;j++){
		for(i=0;i<BLS;i++){
		  vloc[i] -= it_matrix->ilu_val[a_ij_p+i+BLS*j] * V[j_block*BLS+j];
		}
	  }

/*kbw
	printf("V before updating with block %d (%d)\n", j_block, sm_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",V[j_block*BLS+i]);
	printf("\n");
	printf("i_block %d, j_block %d (sm_block %d), Aux before solution:\n",
		   i_block,j_block,sm_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ij_p+i]);
	}
	printf("\n");
	printf("vloc after updating with block %d (%d)\n", j_block, sm_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	//getchar();
/*kew*/

	}

	/* rewrite back to global vector */
	for(i=0;i<BLS;i++) V[i_block*BLS+i] = vloc[i];

/*kbw
	  printf("part of V %d-%d after updating with all small blocks\n",
		 i_block*BLS, i_block*BLS+BLS-1);
	  for(i=0;i<BLS;i++) printf("%25.15lf", V[i_block*BLS+i]);
	  //getchar();
/*kew*/

	  }

	}

/*kbw
  printf("V after FORWARD REDUCTION:\n");
  for(i=0;i<it_matrix->Nrblocks*BLS;i++) printf("%25.15lf",V[i]);
  printf("\n");
  getchar();getchar();
/*kew*/

/* BACK SUBSTITUTION */

	/* loop over all row_blocks - in reverse order */
	for(i_block=it_matrix->Nrblocks-1; i_block >= 0; i_block--){

	  double vloc[BLS];
	  int i;
	  for(i=0;i<BLS;i++) vloc[i] = V[i_block*BLS+i];

/*kbw
	printf("vloc initialized: i_block %d\n", i_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	getchar();
/*kew*/

	  int sm_block;
	  // for all j<i - lower subdiagonal entries
	  for(sm_block = it_matrix->ilu_diag_ptr[i_block] + 1;
	  sm_block < it_matrix->ilu_row_ptr[i_block+1];
	  sm_block++){

	int j_block = it_matrix->ilu_col_ind[sm_block];

	// pointer to a_ij -> sm_block*BLS*BLS
	int a_ij_p = sm_block*BLS*BLS;

	// MUST BE VECTORIZED!!!
	int i,j;
	for(j=0;j<BLS;j++){
	  for(i=0;i<BLS;i++){
		vloc[i] -= it_matrix->ilu_val[a_ij_p+i+BLS*j] * V[j_block*BLS+j];
	  }
	}

/*kbw
	printf("V before updating with block %d (%d)\n", j_block, sm_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",V[j_block*BLS+i]);
	printf("\n");
	printf("i_block %d, j_block %d (sm_block %d), Aux before solution:\n",
		   i_block,j_block,sm_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ij_p+i]);
	}
	printf("\n");
	printf("vloc after updating\n");
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	getchar();
/*kew*/

	  }


	  // update using A_ii

	  // pointer to a_ii -> it_matrix->ilu_diag_ptr[i_block]*BLS*BLS
	  int a_ii_p = it_matrix->ilu_diag_ptr[i_block]*BLS*BLS;

	  // MUST BE VECTORIZED!!!
	  for(i=0;i<BLS;i++) V[i_block*BLS+i] = 0.0;
	  for(j=0;j<BLS;j++){
	for(i=0;i<BLS;i++){
	  V[i_block*BLS+i] += it_matrix->ilu_val[a_ii_p+i+BLS*j] * vloc[j];
	}
	  }

/*kbw
	printf("vloc before solving with Dia\n");
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	printf("i_block %d (%d), Dia before solution:\n",
	   i_block, it_matrix->ilu_diag_ptr[i_block]);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");
	printf("V after updating\n");
	for(i=0;i<BLS;i++) printf("%25.15lf",V[i_block*BLS+i]);
	printf("\n");
	getchar();
/*kew*/


	} // end loop over row blocks in reverse order


  } else if(version==2){

	// APPROXIMATE SOLUTION - PARALLELIZED IN DOMAIN DECOMPOSITION MANNER


#pragma omp parallel default(none) firstprivate(it_matrix, V, B, Ndof)
  {

	int i,j,k;

	// each thread gets its copy of input vector
	double *vtemp_FR_BS;		/* temporary vector*/
	vtemp_FR_BS = (double*)malloc(Ndof*sizeof(double));

#ifdef _OPENMP
	int my_id = omp_get_thread_num();
	int portion = Ndof/omp_get_num_threads() + 1;
#else
	int my_id = 0;
	int portion = Ndof;
#endif

	int my_first_x = (my_id*portion - it_matrix->ILU_half_bandwidth - 10)*BLS;
	if(my_first_x<0) my_first_x = 0;
	int my_last_x = ((my_id+1)*portion + it_matrix->ILU_half_bandwidth + 10)*BLS;
	if(my_last_x>Ndof) my_last_x=Ndof;
	for(i=my_first_x; i<my_last_x; i++) vtemp_FR_BS[i]=V[i];

	int my_first = my_id*portion;
	if(my_first<0) my_first = 0;
	int my_last = (my_id+1)*portion;
	if(my_last>it_matrix->Nrblocks) my_last=it_matrix->Nrblocks;


	/* loop over all row_blocks */
	int i_block;
	/* for(i_block=0; i_block < it_matrix->Nrblocks; i_block++){ */
	for (i_block = my_first; i_block < my_last; i_block++){

	  /* if block is active */
	  if(it_matrix->ilu_row_ptr[i_block]<it_matrix->ilu_row_ptr[i_block+1]){

	/* initialize vloc */
	double vloc[BLS];

	/* substitute value of rhs to vloc */
	if(B==NULL){
	  int i;
	  for(i=0;i<BLS;i++) vloc[i]=it_matrix->rhs[i_block*BLS+i];
	}
	else{
	  int i;
	  for(i=0;i<BLS;i++) vloc[i]=B[i_block*BLS+i];
	}

	int sm_block;
	// for all j<i - lower subdiagonal entries
	for(sm_block = it_matrix->ilu_row_ptr[i_block];
		sm_block < it_matrix->ilu_diag_ptr[i_block];
		sm_block++){

	  int j_block = it_matrix->ilu_col_ind[sm_block];

//#ifdef DEBUG_LSM
	  if(j_block*BLS<my_first_x || j_block*BLS+BLS-1>=my_last_x){

			  printf("wrong my_first_x %d, my_last_x %d, ILU_half_bandwidth %d in FR\n",
				  my_first_x, my_last_x, it_matrix->ILU_half_bandwidth);
			  printf("i_block %d, j_block %d (*BLS %d), \n",
				  i_block, j_block, j_block*BLS);

			  exit(-1);
		  }

//#endif

	  // pointer to a_ij -> sm_block*BLS*BLS
	  int a_ij_p = sm_block*BLS*BLS;

	  // MUST BE VECTORIZED!!!
	  int i,j;
	  for(j=0;j<BLS;j++){
		for(i=0;i<BLS;i++){
		  vloc[i] -= it_matrix->ilu_val[a_ij_p+i+BLS*j] * vtemp_FR_BS[j_block*BLS+j];
		  //vloc[i] -= it_matrix->ilu_val[a_ij_p+i+BLS*j] * V[j_block*BLS+j];
		}
	  }

/*kbw
	printf("V before updating with block %d (%d)\n", j_block, sm_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",vtemp_FR_BS[j_block*BLS+i]);
	for(i=0;i<BLS;i++) printf("%25.15lf",V[j_block*BLS+i]);
	printf("\n");
	printf("i_block %d, j_block %d (sm_block %d), Aux before solution:\n",
		   i_block,j_block,sm_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ij_p+i]);
	}
	printf("\n");
	printf("vloc after updating with block %d (%d)\n", j_block, sm_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	//getchar();
/*kew*/

	}


	/* rewrite back to global vector */
	for(i=0;i<BLS;i++) vtemp_FR_BS[i_block*BLS+i] = vloc[i];

/*kbw
	  printf("part of V %d-%d after updating with all small blocks\n",
		 i_block*BLS, i_block*BLS+BLS-1);
	  for(i=0;i<BLS;i++) printf("%25.15lf", vtemp_FR_BS[i_block*BLS+i]);
	  for(i=0;i<BLS;i++) printf("%25.15lf", V[i_block*BLS+i]);
	  //getchar();
/*kew*/

	  } // end if not ghost block


	} // end loop over row blocks

#pragma omp barrier

	for (i_block = my_first; i_block < my_last; i_block++){

	  /* rewrite back to global vector */
	  for(i=0;i<BLS;i++) V[i_block*BLS+i] = vtemp_FR_BS[i_block*BLS+i];

/*kbw
	  printf("part of V %d-%d after updating with all small blocks\n",
		 i_block*BLS, i_block*BLS+BLS-1);
	  for(i=0;i<BLS;i++) printf("%25.15lf", V[i_block*BLS+i]);
	  //getchar();
/*kew*/
	}

/*kbw
#pragma omp single
	{
  printf("V after FORWARD REDUCTION:\n");
  for(i=0;i<it_matrix->Nrblocks*BLS;i++) printf("%25.15lf",V[i]);
  printf("\n");
  getchar();getchar();
	}
/*kew*/

#pragma omp barrier

/* BACK SUBSTITUTION */

	// each thread gets its copy of vector V after forward reduction ("exchange of DOFS")
	for(i=my_first_x; i<my_last_x; i++) vtemp_FR_BS[i]=V[i];

	/* loop over all row_blocks - in reverse order */
	//for(i_block=it_matrix->Nrblocks-1; i_block >= 0; i_block--){
	for (i_block = my_last-1; i_block >= my_first; i_block--){

	  double vloc[BLS];

	  for(i=0;i<BLS;i++) vloc[i] = vtemp_FR_BS[i_block*BLS+i];
	  //for(i=0;i<BLS;i++) vloc[i] = V[i_block*BLS+i];

/*kbw
	printf("vloc initialized: i_block %d\n", i_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	getchar();
/*kew*/


	int sm_block;
	// for all j<i - lower subdiagonal entries
	for(sm_block = it_matrix->ilu_diag_ptr[i_block] + 1;
	sm_block < it_matrix->ilu_row_ptr[i_block+1];
	sm_block++){

	  int j_block = it_matrix->ilu_col_ind[sm_block];

//#ifdef DEBUG_LSM
	  if(j_block*BLS<my_first_x || j_block*BLS+BLS-1>=my_last_x){
		printf("wrong my_first_x, my_last_x in BS\n"); exit(-1);
	  }
//#endif

	  // pointer to a_ij -> sm_block*BLS*BLS
	  int a_ij_p = sm_block*BLS*BLS;

	  // MUST BE VECTORIZED!!!
	  int i,j;
	  for(j=0;j<BLS;j++){
	for(i=0;i<BLS;i++){
	  vloc[i] -= it_matrix->ilu_val[a_ij_p+i+BLS*j] * vtemp_FR_BS[j_block*BLS+j];
	  //vloc[i] -= it_matrix->ilu_val[a_ij_p+i+BLS*j] * V[j_block*BLS+j];
	}
	  }


/*kbw
	printf("V before updating with block %d (%d)\n", j_block, sm_block);
	for(i=0;i<BLS;i++) printf("%25.15lf",vtemp_FR_BS[j_block*BLS+i]);
	printf("\n");
	printf("i_block %d, j_block %d (sm_block %d), Aux before solution:\n",
		   i_block,j_block,sm_block);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ij_p+i]);
	}
	printf("\n");
	printf("vloc after updating\n");
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	getchar();
/*kew*/

	}

	// pointer to a_ii -> it_matrix->ilu_diag_ptr[i_block]*BLS*BLS
	int a_ii_p = it_matrix->ilu_diag_ptr[i_block]*BLS*BLS;

	// MUST BE VECTORIZED!!!
	for(i=0;i<BLS;i++) vtemp_FR_BS[i_block*BLS+i] = 0.0;
	for(j=0;j<BLS;j++){
	  for(i=0;i<BLS;i++){
	vtemp_FR_BS[i_block*BLS+i] += it_matrix->ilu_val[a_ii_p+i+BLS*j] * vloc[j];
	//V[i_block*BLS+i] += it_matrix->ilu_val[a_ii_p+i+BLS*j] * vloc[j];
	  }
	}


/*kbw
	printf("vloc before solving with Dia\n");
	for(i=0;i<BLS;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	printf("i_block %d (%d), Dia before solution:\n",
	   i_block, it_matrix->ilu_diag_ptr[i_block]);
	for(i=0;i<BLS*BLS;i++){
	  printf("%25.15lf",it_matrix->ilu_val[a_ii_p+i]);
	}
	printf("\n");
	printf("V after updating\n");
	for(i=0;i<BLS;i++) printf("%25.15lf",vtemp_FR_BS[i_block*BLS+i]);
	printf("\n");
	getchar();
/*kew*/

  } // end loop over row blocks in reverse order

#pragma omp barrier

	for (i_block = my_first; i_block < my_last; i_block++){

	  /* rewrite back to global vector */
	  for(i=0;i<BLS;i++) V[i_block*BLS+i] = vtemp_FR_BS[i_block*BLS+i];

/*kbw
	  printf("part of V %d-%d after updating with all small blocks\n",
		 i_block*BLS, i_block*BLS+BLS-1);
	  for(i=0;i<BLS;i++) printf("%25.15lf", V[i_block*BLS+i]);
	  //getchar();
/*kew*/
	}

  } // end parallel region

  } // end version 2

#ifdef TIME_TEST_ILU
  nr_oper += it_matrix->ILU_nr_sm_blocks*2*BLS*BLS; // dgemm
  nr_access +=it_matrix->ILU_nr_sm_blocks*BLS*BLS + it_matrix->Nrblocks*BLS;
										  // optimal version - one access to V
  exec_time = time_clock() - exec_time;
  printf("FR and BS (version %d): nr_oper %lf, nr_access %lf\n",
	 version, nr_oper, nr_access);
  printf("\t\t time %lf -> %lf GFlops, %lf GB/s\n",
	 exec_time, nr_oper/1.e9/exec_time, nr_access*sizeof(double)/1.e9/exec_time);
#endif

/*kbw
  int i;
  printf("V after BACKWARD SUBSTITUTION:\n");
  for(i=0;i<it_matrix->Nrblocks*BLS;i++) printf("%25.15lf",V[i]);
  printf("\n");
  getchar();getchar();
/*kew*/


}


/*---------------------------------------------------------
  lar_get_SM_and_LV_crs_from_bcrs - convert bcrs matrix format to crs
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_crs_from_bcrs( // returns: flag - 0 if not needed delete crs matrices
  int Matrix_id, /* in: matrix ID */
  int offset, /* in: offset in crs_row and crs_col matrix */
  int** crs_row,	   /* out: matrix of rows in crs */
  int** crs_col,	   /* out: matrix of column in crs */
  double** crs_val,	   /* out: matrix of value in crs */
  double** rhs	   /* out: rhs vector */
			   )
{

  itt_bcrs_matrices *it_matrix = &itv_bcrs_matrices[Matrix_id];
  *rhs=it_matrix->rhs;

  // allocate new crs structure
  if(*crs_row==NULL){
	*crs_row = (int*)malloc( (it_matrix->Nrdofgl+1)*sizeof(int));
	if( *crs_row==NULL ) {
	  printf("Not enough memory for allocating crs structure for row_ptr vector\n");
	  exit(-1);
	}
  }
  if(*crs_col==NULL){
	*crs_col = (int*)malloc(it_matrix->Nr_sm_blocks*BLS*BLS*sizeof(int));
	if( *crs_col==NULL ) {
	  printf("Not enough memory for allocating crs structure for col_ind vector\n");
	  exit(-1);
	}
  }
  if(*crs_val==NULL){
	*crs_val= (double*)malloc(it_matrix->Nr_sm_blocks*BLS*BLS*sizeof(double));
	if( *crs_val==NULL ) {
	  printf("Not enough memory for allocating crs structure for val vector\n");
	  exit(-1);
	}
  }


  // converting from bcrs to crs
  int block_row_gl,k,sm_block,l,block_col_gl, block_val_ind;
  int ccolind=0, crowind=0;
  ccolind=0, crowind=0;
  for(block_row_gl=0; block_row_gl < it_matrix->Nrblocks; block_row_gl++){
	for(k=0;k<BLS;k++){
	  (*crs_row)[crowind]=ccolind+offset;
	  crowind++;
	  for(sm_block = it_matrix->block_row_ptr[ block_row_gl ];
		  sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];
		  sm_block++){

	block_col_gl = it_matrix->block_col_ind[ sm_block ]*BLS;

	block_val_ind = sm_block*BLS*BLS;

	for(l=0;l<BLS;l++){
	  assert(ccolind<it_matrix->Nr_sm_blocks*BLS*BLS);

	  (*crs_val)[ccolind]=it_matrix->block_val[ block_val_ind + k + l*BLS];
	  //assert(ccolind>nr_sm_blocks*BLS*BLS);
	  (*crs_col)[ccolind]=block_col_gl+l+offset;
	  ccolind++;
	}

	  }
	}
  }
  assert(crowind==it_matrix->Nrdofgl);
  (*crs_row)[crowind]=ccolind+offset;

  int info=1;

  return(info);

}

/*---------------------------------------------------------
  lar_get_SM_and_LV_coo_from_bcrs - convert bcrs matrix format to coo for external solver
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_coo_from_bcrs( // returns: flag - 0 if not needed delete coo matrices
  int Matrix_id, /* in: matrix ID */
  int* nnz,		/* out: Number of nonzero elements in matrix */
  int offset, /* in: offset in coo_row and coo_col matrix */
  int** coo_row,	   /* out: matrix of rows in coo */
  int** coo_col,	   /* out: matrix of column in coo */
  double** coo_val,	   /* out: matrix of value in coo */
  double** rhs	   /* out: rhs vector */
			   )
{
  itt_bcrs_matrices *it_matrix = &itv_bcrs_matrices[Matrix_id];
  int info=1;
  int i,j;
  *nnz=it_matrix->Nr_sm_blocks*BLS*BLS;
  *rhs=it_matrix->rhs;
  /*
  int BLSc=2;
  it_matrix->Nr_sm_blocks=5;
  it_matrix->Nrblocks=3;

 // it_matrix->block_row_ptr[4]={0,1,3,5};
  //it_matrix->block_col_ind[5]={0,1,2,0,2};
  //it_matrix->block_val[20]={1,0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};


  it_matrix->block_row_ptr[0]=0;
  it_matrix->block_row_ptr[1]=1;
  it_matrix->block_row_ptr[2]=3;
  it_matrix->block_row_ptr[3]=5;

  it_matrix->block_col_ind[0]=0;
  it_matrix->block_col_ind[1]=1;
  it_matrix->block_col_ind[2]=2;
  it_matrix->block_col_ind[3]=0;
  it_matrix->block_col_ind[4]=2;

  it_matrix->block_val[0]=1.0;
  it_matrix->block_val[1]=0.0;
  it_matrix->block_val[2]=2.0;
  it_matrix->block_val[3]=3.0;
  it_matrix->block_val[4]=4.0;
  it_matrix->block_val[5]=5.0;
  it_matrix->block_val[6]=6.0;
  it_matrix->block_val[7]=7.0;
  it_matrix->block_val[8]=8.0;
  it_matrix->block_val[9]=9.0;
  it_matrix->block_val[10]=10.0;
  it_matrix->block_val[11]=11.0;
  it_matrix->block_val[12]=12.0;
  it_matrix->block_val[13]=13.0;
  it_matrix->block_val[14]=14.0;
  it_matrix->block_val[15]=15.0;
  it_matrix->block_val[16]=16.0;
  it_matrix->block_val[17]=17.0;
  it_matrix->block_val[18]=18.0;
  it_matrix->block_val[19]=19.0;
  ////


  it_matrix->block_row_ptr[0]=0;
  it_matrix->block_row_ptr[1]=2;
  it_matrix->block_row_ptr[2]=3;
  it_matrix->block_row_ptr[3]=5;

  it_matrix->block_col_ind[0]=0;
  it_matrix->block_col_ind[1]=2;
  it_matrix->block_col_ind[2]=1;
  it_matrix->block_col_ind[3]=1;
  it_matrix->block_col_ind[4]=2;

  it_matrix->block_val[0]=1.0;
  it_matrix->block_val[1]=2.0;
  it_matrix->block_val[2]=0.0;
  it_matrix->block_val[3]=3.0;
  it_matrix->block_val[4]=12.0;
  it_matrix->block_val[5]=14.0;
  it_matrix->block_val[6]=13.0;
  it_matrix->block_val[7]=15.0;
  it_matrix->block_val[8]=4.0;
  it_matrix->block_val[9]=6.0;
  it_matrix->block_val[10]=5.0;
  it_matrix->block_val[11]=7.0;
  it_matrix->block_val[12]=8.0;
  it_matrix->block_val[13]=10.0;
  it_matrix->block_val[14]=9.0;
  it_matrix->block_val[15]=11.0;
  it_matrix->block_val[16]=16.0;
  it_matrix->block_val[17]=18.0;
  it_matrix->block_val[18]=17.0;
  it_matrix->block_val[19]=19.0;
  */

  //



  if(*coo_col==NULL){
	*coo_col = (int*)malloc(*nnz*sizeof(int));
	if( *coo_col==NULL ) {
	  printf("Not enough memory for allocating coo structure for col_ind vector\n");
	  exit(-1);
	}
  }
  if(*coo_row==NULL){
	*coo_row = (int*)malloc(*nnz*sizeof(int));
	if( *coo_row==NULL ) {
	  printf("Not enough memory for allocating coo structure for row_ptr vector\n");
	  exit(-1);
	}
  }
  if(*coo_val==NULL){
	*coo_val= (double*)malloc(*nnz*sizeof(double));
	if( *coo_val==NULL ) {
	  printf("Not enough memory for allocating coo structure for val vector\n");
	  exit(-1);
	}
  }

	// converting from bcrs to coo
  int block_row_gl,k,sm_block,l,block_col_gl, block_val_ind;
  int ccolind=0, crowind=0;
  ccolind=0, crowind=0;
  for(block_row_gl=0; block_row_gl < it_matrix->Nrblocks; block_row_gl++){
	for(k=0;k<BLS;k++){
	  //(*crs_row)[crowind]=ccolind+offset;
	  //crowind++;
	  for(sm_block = it_matrix->block_row_ptr[ block_row_gl ];
		  sm_block < it_matrix->block_row_ptr[ block_row_gl+1 ];
		  sm_block++){

	block_col_gl = it_matrix->block_col_ind[ sm_block ]*BLS;

	block_val_ind = sm_block*BLS*BLS;

	for(l=0;l<BLS;l++){
	  assert(ccolind<it_matrix->Nr_sm_blocks*BLS*BLS);

	  (*coo_val)[ccolind]=it_matrix->block_val[ block_val_ind + k + l*BLS];
	  //(*nnz)++;
	  (*coo_col)[ccolind]=block_col_gl+l+offset;
	  (*coo_row)[ccolind]=block_row_gl*BLS+k+offset;
	 // printf("(%d,%d,%d) ",ccolind,l,k);
	  ccolind++;
	}

	  }
	}
  }
  /*
  FILE *fp=fopen("mumps.txt", "w");
	k=1;
	for(i=0;i<*nnz;i++){
		if(k!=(*coo_row)[i])fprintf(fp,"\n");
		fprintf(fp,"(%d,%d,%lf) ",(*coo_row)[i]-1,(*coo_col)[i]-1,(*coo_val)[i]);
		k= (*coo_row)[i];

	}exit(0);

  */

  //assert(crowind==it_matrix->Nrdofgl);
  //(*crs_row)[crowind]=ccolind+offset;

	//
	//for(i=0;i<50;i++)printf("(%lf,%d,%d)\n",(*coo_val)[i],(*coo_row)[i],(*coo_col)[i]);//exit(0);

	//for(i=0;i<(*nnz);i++)printf("%d ",(*coo_col)[i]);exit(0);


	//printf("\n\n\n\n     nnz=%d\n\n\n\n",*nnz);

  return(info);
}

/************************ local utilities *****************/

/*---------------------------------------------------------
lar_util_chk_list_bcrs - to check whether Num is on the list List
	with length Ll (filled with numbers and LAC_END_OF_LIST at the end)
---------------------------------------------------------*/
int lar_util_chk_list_bcrs(	/* returns: */
			/* >=0 - position on the list */
					/* LAC_PUT_LIST_NOT_FOUND - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	)
{

  int i, il;

  for(i=0;i<Ll;i++){
	if((il=List[i])==LAC_END_OF_LIST) break;
	/* found on the list on i-th position - counting from 0 */
	if(Num==il) return(i);
  }
  /* not found on the list */
  return(LAC_PUT_LIST_NOT_FOUND);
}

/*---------------------------------------------------------
lar_util_put_list_bcrs - to put Num on the list List (with length Ll and filled with
					 meaningfull numbers and one or more LAC_END_OF_LIST at the end)
---------------------------------------------------------*/
int lar_util_put_list_bcrs( /* returns*/
		/*  >=0 - success, position on the list */
				/*  LAC_PUT_LIST_FOUND - not put - already on the list */
				/*  LAC_PUT_LIST_FULL - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
	)
{

  int i, il;

  for(i=0;i<Ll;i++){
	if((il=List[i])==LAC_END_OF_LIST) break;
	/* found on the list on i-th position - counting from 0 */
	if(Num==il) return(LAC_PUT_LIST_FOUND);
  }
  /* if list is full return error message */
  if(i==Ll) return(LAC_PUT_LIST_FULL);
  /* update the list and return*/
  List[i]=Num;
  return(i);
}

/*---------------------------------------------------------
  lar_sort_short - to sort an array by insertion (for small integer arrays)
---------------------------------------------------------*/
int lar_sort_short(
		   int* A, // array
		   int p,  //first index
		   int k   // last index
		   )
{
  int i,j;
  int t;

  for(i=p+1;i<=k;i++) {
	t=A[i];
	j=i-1;
	while( j>=p && A[j]>t ) {
	  A[j+1]=A[j];
	  j--;
	}
	A[j+1]=t;
  }
  return 0;
}
