/*************************************************************
File contains procedures:
 lar_allocate_SM_and_LV_block - to allocate space for stiffness matrix and load vector
  lar_initialize_SM_and_LV_block - to initialize stiffness matrix and/or load vector
  lar_get_storage_block - to compute storage of SM, LV and preconditioner
  lar_assemble_SM_and_LV_block - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  lar_allocate_preconditioner_block - to allocate space for preconditioner
  lar_fill_preconditioner_block - to fill preconditioner
  lar_free_preconditioner_block - to free space for a block structure
  lar_free_SM_and_LV_block - to free space for a block structure

lar_compute_residual_block - to compute the residual of the not preconditioned
	system of equations, v = ( b - Ax )
lar_perform_BJ_or_GS_iterations_block - to perform one iteration of block Gauss-Seidel
	or block Jacobi algorithm
	 v_out = v_in + M^-1 * ( b - A * v_in )
lar_perform_rhsub_block - to perform forward reduction and back-substitution for ILU
		   preconditioning

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

#include <modfem/uth_log.h>

// interface for linear algebra package
#include <modfem/ls_mkb/lah_intf.h>

/* internal information for the block module */
#include "./lah_block.h"

/* LAPACK and BLAS procedures */
#include <modfem/lin_alg_intf.h>


/* minimal number of dofs to make use of DGEMV reasonable */
#define MIN_DOF_DGEMV 10000
#define SMALL 1e-15 /* approximation of round-off error */
#define TOL 1e-9    /* default tolerance for different iterative processes */

#define Max_dof_block_dia 1500
#define MAX_BL_BL 20

const int LAC_MAX_BL_NGB = 200; /* maximal no of block's neighbors */
const int LAC_MAX_BL_DIA_NGB = 200;/* maximal number of dia block's neighbors */


/* GLOBAL VARIABLES */
int   itv_nr_matrices=0;        /* the number of solvers in the problem */
itt_matrices itv_matrices[LAC_MAX_MATRICES];        /* array of solvers */

/*---------------------------------------------------------
  lar_allocate_SM_and_LV_block - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
int lar_allocate_SM_and_LV_block( // returns: matrix index in itv_matrices array
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
	int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]<=0 - iblock is a ghost block with no rows associated
  int** L_offbl    /* in: list of lists of off diagonal blocks */
	)
{

/* local variables */
  itt_matrices *it_matrix;
  itt_blocks *block;	/* to simplify typing */

/* auxiliary variable */
  int i, iblock, bloff, nrblocks, nrdofbl, nroffbl;

/*++++++++++++++++ executable statements ++++++++++++++++*/

// create structures for a new matrix
  itv_nr_matrices++;

  if(itv_nr_matrices>=LAC_MAX_MATRICES){
	printf("Too much (%d) matrices requested in las_block_intf! (correct lah_block.h)\n",
	   itv_nr_matrices);
  }

  it_matrix = &itv_matrices[itv_nr_matrices-1];

  it_matrix->Max_SM_size = Max_SM_size;
  it_matrix->Nrblocks = Nrblocks;
  it_matrix->Nrdofgl = Nrdof_glob;

  /* allocate space for pointers to blocks */
  nrblocks = it_matrix->Nrblocks;
  it_matrix->Block = (itt_blocks **)malloc((nrblocks+1)*sizeof(itt_blocks *));
  if( it_matrix->Block==NULL ) {
	printf("Not enough memory for allocating block structure for stiffness matrix\n");
	exit(-1);
  }

  /* !!! offset 1 numbering of blocks in lad_block !!! */
  for(iblock=1;iblock<=nrblocks;iblock++){

	/* !!! offset 0 numbering of blocks in sid_mkb and lsd_mkb !!! */
	int index = iblock-1;

	/* allocate space for a new block */
	it_matrix->Block[iblock] = (itt_blocks *)malloc(sizeof(itt_blocks));
	if(it_matrix->Block[iblock]==NULL){
	  printf("it_create_blocks: no space for block %d\n", iblock);
	  exit(1);
	}

/*kbw
	printf("Allocating block %d, it_matrix %u, it_matrix->Block %u, it_matrix->Block[%d] %u\n",
		   iblock, it_matrix, it_matrix->Block, iblock, it_matrix->Block[iblock]);
/*kew*/

	block = it_matrix->Block[iblock];

	/* specify block data */
	block->Ndof = Nrdofbl[index];
	block->Posg = Posglob[index];
	if(Nroffbl[index]>0){
	  block->Lngb=lar_util_ivector(Nroffbl[index]+1,"Lngb in it_create_blocks");
	  block->Lngb[0] = Nroffbl[index];
	  for(i=1;i<=Nroffbl[index];i++){
	/* !!! offset 1 numbering of blocks in lad_block !!! */
	block->Lngb[i] = L_offbl[index][i-1] + 1;
	  }
	}
	else block->Lngb=NULL;
	// all ghost blocks are identified by block->Lngb==NULL

  }

  for(iblock=1;iblock<=nrblocks;iblock++){

	block = it_matrix->Block[iblock];

	nrdofbl = block->Ndof;


	/* for ordinary (not ghost) blocks */
	if(block->Lngb!=NULL){

	  /* allocate space for diagonal matrix, rhs vector and */
	  /* auxiliary ips vector for pivoting information */
	  block->Dia=lar_util_dvector(nrdofbl*nrdofbl,"Dia in it_create_blocks");
	  block->Rhs=lar_util_dvector(nrdofbl,"Rhs in it_create_blocks");
	  block->Ips=lar_util_ivector(nrdofbl,"ips it_create_blocks");

	  /* array of pointers to off-diagonal blocks */
	  block->Aux=(double **)malloc(block->Lngb[0]*sizeof(double *));

	  for(i=1;i<=block->Lngb[0];i++) {

	bloff=block->Lngb[i];
	nroffbl = it_matrix->Block[bloff]->Ndof;

	/* allocate space for a new off-diagonal matrix */
	block->Aux[i-1]=lar_util_dvector(nrdofbl*nroffbl,"Aux in it_create_blocks");

/*kbw
	printf("Allocated space for block %d and block %d (%d), size %dx%d\n",
		   iblock,bloff,i-1,nrdofbl,nroffbl);
/*kew*/

	  }
	}
	else{

	  block->Dia=NULL;
	  block->Rhs=NULL;
	  block->Ips=NULL;
	  block->Aux=NULL;

	}
  }



  return(itv_nr_matrices-1);
}

/*---------------------------------------------------------
  lar_initialize_SM_and_LV_block - to initialize stiffness matrix and/or load vector
---------------------------------------------------------*/
int lar_initialize_SM_and_LV_block(
  int Matrix_id,   /* in: matrix ID */
  int Scope    /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  )

{

/* local variables */
  itt_matrices *it_matrix;
  itt_blocks *block;	/* to simplify typing */

/* auxiliary variable */
  int i, iblock, ineig, bloff, nrblocks, nrdofbl, nroffbl;

/*++++++++++++++++ executable statements ++++++++++++++++*/


  it_matrix = &itv_matrices[Matrix_id];

  /* allocate space for pointers to blocks */
  nrblocks = it_matrix->Nrblocks;
  for(iblock=1;iblock<=nrblocks;iblock++){

	block = it_matrix->Block[iblock];

/*kbw
	printf("Initializing block %d, it_matrix %u, it_matrix->Block %u, it_matrix->Block[%d] %u\n",
		   iblock, it_matrix, it_matrix->Block, iblock, it_matrix->Block[iblock]);
/*kew*/

	/* for ordinary (not ghost) blocks */
	if(block->Lngb!=NULL){

	  nrdofbl = block->Ndof;

	  /* initialize blocks forming RHS vector */
	  lar_util_d_zero(block->Rhs,nrdofbl);

	  if(Scope==LAC_SCOPE_SM_AND_LV){

	/* initialize arrays forming stiffness matrix */
	lar_util_d_zero(block->Dia,nrdofbl*nrdofbl);

	for(ineig=1;ineig<=block->Lngb[0];ineig++) {

	  bloff=block->Lngb[ineig];
	  nroffbl = it_matrix->Block[bloff]->Ndof;

/*kbw
	  printf("Initializing block %d, Aux block %d, size %dx%d\n",
		 iblock, ineig-1, nrdofbl, nroffbl);
/*kew*/


	  for(i=0;i<nrdofbl*nroffbl;i++) {

		block->Aux[ineig-1][i]=0;

/*kbw
		if(iblock == 1 && ineig == 1 && i == 61){
		  printf("Initializing Aux block %d, entry %d\n", ineig-1, i);
		  getchar();
		  getchar();
		}
/*kew*/

	  }
	}
	  }
	}
  }

  return(0);
}

/*---------------------------------------------------------
  lar_get_storage_block - to compute storage of SM, LV and preconditioner
---------------------------------------------------------*/
double lar_get_storage_block( /* returns: storage in MB */
  int Matrix_id   /* in: matrix ID */
			   )
{

  itt_blocks *block;	      /* array of pointers to small blocks */
  itt_blocks_dia *block_dia; /* array of pointers to diagonal blocks */
  int i,j,k,iaux,jaux,kaux,iblock,nrdofbl,nroffbl;
  itt_matrices *it_matrix = &itv_matrices[Matrix_id];
  double storage = 0.0; double daux=0.0;

  for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){
	block = it_matrix->Block[iblock];
	if(block->Lngb!=NULL){
	  nrdofbl = block->Ndof;
	  /* integers: Ndof, Posg */
	  storage += 2.0*sizeof(int)/1024.0/1024.0;
	  /* pointers to arrays of integers: Lngb, Ips */
	  storage += 2.0*sizeof(int *)/1024.0/1024.0;
	  /* integer arrays: Lngb, Ips */
	  storage += (block->Lngb[0]+1.0+nrdofbl)*sizeof(int)/1024.0/1024.0;
	  /* pointers to arrays of doubles:Dia, Rhs, Aux, Aux[i] */
	  storage += (3.0+block->Lngb[0])*sizeof(double *)/1024.0/1024.0;
	  /* arrays of doubles: Rhs */
	  storage += nrdofbl*sizeof(double)/1024.0/1024.0;

	  /* arrays of doubles: Dia */
	  storage += nrdofbl*nrdofbl*sizeof(double)/1024.0/1024.0;
	  /* arrays of doubles: Aux */
	  for(j=0;j<block->Lngb[0];j++){
	nroffbl=it_matrix->Block[block->Lngb[j+1]]->Ndof;
	storage += nrdofbl*nroffbl*sizeof(double)/1024.0/1024.0;
	  }
	}
  }
  printf("Storage for %d elementary blocks %f MBytes\n",
	 it_matrix->Nrblocks,storage-daux);
  daux=storage;

  if(it_matrix->Precon != NO_PRECON){
	for(iblock=1;iblock<=it_matrix->Nrblocks_dia;iblock++){
	  block_dia=it_matrix->Block_dia[iblock];
	  iaux=block_dia->Lsmall[0];
	  jaux = block_dia->Lpos[0]; /* total number of dofs for block_dia */
	  if( it_matrix->Precon==MULTI_ILU ||it_matrix->Precon==BLOCK_ILU )
	kaux=block_dia->Lneig[0];
	  else kaux=0;
	  /* pointers to arrays: Lsmall, Lelem, Lneig, Lpos, Ips */
	  storage += 5.0*sizeof(int *)/1024.0/1024.0;
	  /* integer arrays: Lsmall, Lpos, Ips and Lneig*/
	  storage += (2.0*(iaux+1.0)+jaux+kaux+1.0)*sizeof(int)/1024.0/1024.0;
	  /* pointers to arrays: Dia and Aux */
	  storage += (2.0+kaux)*sizeof(double *)/1024.0/1024.0;
	  /* array of doubles: Dia */
	  storage += 1.0*jaux*jaux*sizeof(double)/1024.0/1024.0;
	  if( it_matrix->Precon==MULTI_ILU ||it_matrix->Precon==BLOCK_ILU ){
	/* arrays of doubles: Aux[i] */
	for(i=0;i<block_dia->Lneig[0];i++){
	  j = block_dia->Lneig[i+1];
	  k = it_matrix->Block_dia[j]->Lpos[0];
	  storage += jaux*k*sizeof(double)/1024.0/1024.0;
	}
	  }
	}

	printf("Storage for %d diagonal blocks %f MBytes\n",
	   it_matrix->Nrblocks_dia, storage-daux);
	daux=storage;

  }

  return(storage);
}


/*------------------------------------------------------------
  lar_assemble_SM_and_LV_block - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
!!!!!!!!!******************************************!!!!!!!!!!!!!!
lar_assemble_SM_and_LV allows for sending GHOST blocks on lists of blocks.
GHOST blocks  (indicated by block->Lngb==NULL) are blocks with associated block
structures, but without associated rows of GLOBAL stiffness matrix.
Such blocks are used for matrix-vector product with rows belonging to ordinary
(not ghost) blocks. In lar_assemble_SM_and_LV we nevertheless assume that provided
Stiff_mat is square, i.e. there are entries for all blocks. Procedures
calling lar_assemble_SM_and_LV must conform to this convention
!!!!!!!!!******************************************!!!!!!!!!!!!!!
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_block(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Scope,         /* in: indicator for the scope of computations: */
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
  itt_blocks *block_i, *block_j;	/* to simplify typing */

  /* auxiliary variables */
  int iblock, jblock, nrdof_i, nrdof_j, posbl_i, posbl_j;
  int posloc_i, posloc_j, nrdof, nrdof_glob, ibl;
  int i,j,k, iaux, jaux;

/*kbw
  static FILE *fp=NULL;
/*kew*/

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kbw
  if(fp==NULL) fp=fopen("stiff_mat","a");
/*kew*/


/*kbw
	printf("\n\nEntering lar_assemble_SM_and_LV_block!!!\n");
/*kew*/

  /* get pointer to the level structure */
  itt_matrices *it_matrix = &itv_matrices[Matrix_id];

  /* compute local stiffness matrix nrdof */
  nrdof=0;
  for(iblock=0;iblock<Nr_dof_bl;iblock++){
	nrdof += L_bl_nrdof[iblock];
  }

  /* loop over stiffness matrix blocks */
  posloc_i=0;
  for(iblock=0;iblock<Nr_dof_bl;iblock++){

/*kbw
	printf("\n\nassembling matrix %d: it_matrix %u, it_matrix->Block %u\n",
		   Matrix_id, it_matrix, it_matrix->Block);
/*kew*/

	// offset 1 numbering of blocks in lad_block !!!!!
	block_i =  it_matrix->Block[L_bl_id[iblock]+1];

	/* number of dofs for a block */
	nrdof_i = L_bl_nrdof[iblock];


/*kbw
	if(block_i->Posg >= 0){
	  printf("\nin group of columns %d (block_i - bl_id %d), posglob_i %d, nrdofs_i %d\n\n",
		 iblock,  L_bl_id[iblock]+1, block_i->Posg, L_bl_nrdof[iblock]);
	  getchar();
	}
/*kew*/

#ifdef DEBUG
	if(nrdof_i!=block_i->Ndof){
	  printf("Error 409 in block structure for SM in assemble!, i %d != %d\n",
		 nrdof_i, block_i->Ndof);
	  exit(-1);
	}
#endif

	posloc_j=0;
	for(jblock=0;jblock<Nr_dof_bl;jblock++){

	  // offset 1 numbering of blocks in lad_block !!!!!
	  block_j =  it_matrix->Block[L_bl_id[jblock]+1];

	  /* number of dofs for a block */
	  nrdof_j = L_bl_nrdof[jblock];


	  /* for ordinary (not ghost) blocks */
	  if(block_j->Lngb!=NULL){

#ifdef DEBUG
	if(nrdof_j!=block_j->Ndof){
	  printf("Error in block structure for SM in assemble!, j %d != %d\n",
		 nrdof_j, block_j->Ndof);
	  exit(-1);
	}
#endif

	iaux = 0;
	jaux = posloc_j + posloc_i*nrdof;

	if(block_i==block_j){

	  for(i=0;i<nrdof_i;i++){

/*kbw
		if(block_i->Posg >= 0){
		  printf("assembling Dia entries from local : %d - %d ",
			 jaux, jaux+nrdof_j-1 );
		  printf("to global by columns %d - %d: \n",
			 block_j->Posg + (block_i->Posg+i)*it_matrix->Nrdofgl,
			 block_j->Posg + (block_i->Posg+i)*it_matrix->Nrdofgl + nrdof_j-1);
		}
/*kew*/

		for(j=0;j<nrdof_j;j++){
/*kbw
		  if(block_i->Posg >= 0){
		printf("%20.15lf + %20.15lf",
			   Stiff_mat[jaux+j], block_i->Dia[iaux+j]);
		  }
/*kew*/
		  block_i->Dia[iaux+j] += Stiff_mat[jaux+j];

/*kbw
		  if(block_i->Posg >= 0){
		printf(" =  %20.15lf\n",
			   block_i->Dia[iaux+j]);
		  }
		  if((posloc_j+j)==21 &&
		 (block_j->Posg + (block_i->Posg+i)*it_matrix->Nrdofgl +j)==3301)
		   {printf("posloc 21\n"); getchar(); getchar();}
/*kew*/

		}

		jaux += nrdof;
		iaux += nrdof_j;
	  }


	}
	else{

/*kbw
	  if(block_i->Posg >= 0){
		printf("\nin group of rows %d (block_j - bl_id %d), posglob_j %d, nrdofs_j %d\n",
		   jblock, L_bl_id[jblock]+1, block_j->Posg, L_bl_nrdof[jblock]);
		getchar();
	  }
/*kew*/

	// offset 1 numbering of blocks in lad_block !!!!!
	  int lad_block_id = L_bl_id[iblock]+1;
	  ibl = lar_util_chk_list(lad_block_id,
				  &block_j->Lngb[1],block_j->Lngb[0]);

/*kbw
	  if(block_i->Posg >= 0){
		printf("block_i being neighbor no %d, assembling block_j->Aux[%d] - %d values:\n",
		   ibl-1, ibl-1, nrdof_j*nrdof_i );
	  }
/*kew*/

// due to ns_supg_ext we allow for supplying SM with associated lists of DOF blocks
// where some combinations of blocks are not considered neighbors!!!!!!
	  if(ibl>0){

	   for(i=0;i<nrdof_i;i++){

/*kbw
if(block_i->Posg >= 0){
	printf("\nassembling Aux[%d] entries from local : %d - %d ",
		   ibl-1, jaux, jaux+nrdof_j-1 );
	printf("to global %d - %d: \n",
		   block_j->Posg + (block_i->Posg+i) * it_matrix->Nrdofgl,
		   block_j->Posg + (block_i->Posg+i)*it_matrix->Nrdofgl + nrdof_j-1 );
}
/*kew*/

		for(j=0;j<nrdof_j;j++){

/*kbw
if(block_i->Posg >= 0){
		printf("%20.15lf + %20.15lf ",
			Stiff_mat[jaux+j], block_j->Aux[ibl-1][iaux+j]);
}
/*kew*/

		  block_j->Aux[ibl-1][iaux+j] += Stiff_mat[jaux+j];


/*kbw
if(block_i->Posg >= 0){
		printf("=  %20.15lf\n",
			block_j->Aux[ibl-1][iaux+j]);
}
		  if((posloc_j+j)==21 &&
		 (block_j->Posg + (block_i->Posg+i)*it_matrix->Nrdofgl +j)==3301)
		{ printf("posloc 21: block %d, Aux %d, entry %d\n",
			 L_bl_id[jblock]+1,ibl-1,iaux+j); getchar(); getchar();}
/*kew*/

		}

		jaux += nrdof;
		iaux += nrdof_j;
	  }

	 } // end if horizontal block (block_i) being neighbor

	} // end if AUX block (off-diagonal)

/*kbw
	if(block_i->Posg >= 0){
	  printf("\n");
	}
/*kew*/


	  } /* end if ordinary (not ghost) jblock */

	  posloc_j += nrdof_j;
	}

	/* for assembling RHSV loop over iblocks should be treated as loop over */
	/* groups of rows !!! */

	/* for ordinary (not ghost) blocks */
	if(block_i->Lngb!=NULL){

/*kbw
	  if(block_i->Posg >= 0){
	printf("\nin group of rows %d (%d), assembling Rhs entries \nfrom local : %d - %d ",
		   iblock, L_bl_id[iblock]+1, posloc_i, posloc_i+nrdof_i-1 );
	printf("to global %d - %d: \n",
		   block_i->Posg, block_i->Posg+nrdof_i-1 );
	  }
/*kew*/

	  for(i=0;i<nrdof_i;i++) {

/*kbw
	if(block_i->Posg >= 0){
	  printf("%20.15lf + %20.15lf",
		 Rhs_vect[posloc_i+i], block_i->Rhs[i]);
	}
/*kew*/

	/* assemble right hand side block's entries */
	block_i->Rhs[i] += Rhs_vect[posloc_i+i];

/*kbw
	if(block_i->Posg >= 0){
	  printf(" = %20.15lf\n",
		 block_i->Rhs[i]);
	}
/*kew*/

	  }

/*kbw
	  if(block_i->Posg >= 0){
	printf("\n");
	  }
/*kew*/


	} // end if not ghost block !!!

/*kbw
	if(block_i->Posg >= 0){
	  getchar();
	}
/*kew*/


	posloc_i += nrdof_i;
  }


  return(1);
}


/*---------------------------------------------------------
lar_allocate_preconditioner_block - to create preconditioner blocks corresponding
					   to small subdomains of neighboring elements
---------------------------------------------------------*/
int lar_allocate_preconditioner_block( /* returns:   >0 number of diagonal blocks */
						  /*	       <=0 - error */
  int Matrix_id,   /* in: matrix ID */
  int Precon,       /* in: type of preconditioner - see lah_block.h - line circa 45 */
  int Block_type  /* in: number of elementary DOF blocks in a solver block */
  )
{

  /* pointers to solver structure and stiffness matrix blocks' info */
  itt_matrices *it_matrix;

/* local variables */
  itt_blocks_dia *block_dia;/* to simplify typing */
  int iblock, ib_dia;	/* counters for blocks */
  int nrdof_dia;		/* number of dofs in a block */
  int nmfa;	        /* number of faces */

  int* tags;		/* to tag objects for which blocks were created */
  int nr_neighb, bl_neighb[100];/* number and list of block's neighbors */
  int **temp_list, i_temp_list;	/* temporary list for block_dia blocks */
  int nr_bl_block;	/* number of elementary blocks in a dia block */
  int neig_loc[LAC_MAX_BL_DIA_NGB]; /* temp array with dia block neighbors */
  int **big_block; /* temporary array with big block numbers for small blocks */


/* auxiliary variables */
  int i, j, k, l, iaux, jaux, kaux, laux, ibl, nrnos, nrneig, ineig;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*ok_kbw*/
  printf("Entering lar_allocate_preconditioner_block: Matrix_id %d, Block_type %d, Precon %d\n",
	 Matrix_id, Block_type, Precon);
  //getchar();
  //getchar();
/*kew*/

  it_matrix = &itv_matrices[Matrix_id];

  it_matrix->Block_type = Block_type;
  it_matrix->Precon = Precon;

  if(abs(it_matrix->Block_type)==1){

/* allocate space for pointers to blocks */
	it_matrix->Block_dia =
	  (itt_blocks_dia **)malloc((it_matrix->Nrblocks+1)
				   *sizeof(itt_blocks_dia *));

/* dia->Dia blocks are just inverted Dia blocks */
	it_matrix->Nrblocks_dia = 0;
	ib_dia = 0;
	for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){

	  /* for ordinary (not ghost) blocks */
	  if(it_matrix->Block[iblock]->Lngb!=NULL){

	it_matrix->Nrblocks_dia++;
	ib_dia++;

	/* allocate space for block_dia structure */
	it_matrix->Block_dia[ib_dia] =
	  (itt_blocks_dia *)malloc(sizeof(itt_blocks_dia));
	if(it_matrix->Block_dia[ib_dia]==NULL){
	  printf("it_create_blocks_dia: no space for block %d\n", ib_dia);
	  return(1);
	}

	/* number of small blocks (nodes) in dia is equal to 1*/
	nrnos=1;

	/* prepare diagonal block for assembling */
	block_dia = it_matrix->Block_dia[ib_dia];
	block_dia->Lsmall=lar_util_ivector(nrnos+1,"block_dia->Lsmall in create_blocks_dia");

	block_dia->Lelem=NULL;

	block_dia->Lpos=lar_util_ivector(nrnos+1,
				   "block_dia->Lpos in create_blocks_dia");

	/* number of small blocks */
	block_dia->Lsmall[0] = nrnos;

	/* put small block number */
	block_dia->Lsmall[1] = iblock;

	/* number of neighbors and their numbers */
	if( it_matrix->Precon==MULTI_ILU||it_matrix->Precon==BLOCK_ILU ){
	  iaux = it_matrix->Block[iblock]->Lngb[0];
	  nrneig=0;
	  for(i=1;i<=iaux;i++){
		jaux = it_matrix->Block[iblock]->Lngb[i];
		if(it_matrix->Block[jaux]->Lngb!=NULL) nrneig++;
	  }
	  block_dia->Lneig=lar_util_ivector(nrneig+1,"block_dia->Lneig in create_blocks_dia");
	  block_dia->Lneig[0] = nrneig;
	  nrneig=0;
	  for(i=1;i<=iaux;i++){
		jaux = it_matrix->Block[iblock]->Lngb[i];
		if(it_matrix->Block[jaux]->Lngb!=NULL) {
		  nrneig++;
		  block_dia->Lneig[nrneig]=it_matrix->Block[iblock]->Lngb[i];
		}
	  }
	}
	else block_dia->Lneig=NULL;

	/* indicate local position of first dof in dia subarray */
	block_dia->Lpos[1] = 0;

	/* total number of dofs for dia */
	nrdof_dia = it_matrix->Block[iblock]->Ndof;

	/* total number of dofs for dia */
	block_dia->Lpos[0] = nrdof_dia;

	/* allocate space for dia and Ips */
	block_dia->Dia=
	  lar_util_dvector(nrdof_dia*nrdof_dia,"Dia in create_blocks_dia");
	block_dia->Ips=
	  lar_util_ivector(nrdof_dia,"Ips in create_blocks_dia");

	  } /* if ordinary (not ghost) block */

	} /* for all diagonal blocks */


  }
  else if(abs(it_matrix->Block_type)>1){ /* for larger blocks */

/* subdomains composed of blocks and its neighbors */
	if(    abs(it_matrix->Block_type)==2 || abs(it_matrix->Block_type)==3
	|| abs(it_matrix->Block_type)==4 || abs(it_matrix->Block_type)==5
	|| abs(it_matrix->Block_type)==6 || abs(it_matrix->Block_type)==7 ){

	  temp_list = lar_util_imatrix(it_matrix->Nrblocks+1,LAC_MAX_BL_NGB+1,
				 "temp_list in create_block_dia");

/* for neighborhood blocks with different overlap */

/* tagged blocks are excluded from being internal block_dia blocks */
	  tags = lar_util_ivector(it_matrix->Nrblocks+1,"tags in create_blocks_dia");
	  for(ibl=1;ibl<=it_matrix->Nrblocks;ibl++) tags[ibl] = 0;

/* allocate space for diagonal blocks */
	  it_matrix->Nrblocks_dia = 0;

/* in a loop over all blocks */
	  for(ibl=1;ibl<=it_matrix->Nrblocks;ibl++){

/* if block ordinary (not ghost) and not tagged */
	if( it_matrix->Block[ibl]->Lngb != NULL && (tags[ibl]==0 ||
			 abs(it_matrix->Block_type)==2 ||
			 abs(it_matrix->Block_type)==3 )){

/* increase the number of diagonal blocks and tag a block */
	  it_matrix->Nrblocks_dia++;
	  tags[ibl]=it_matrix->Nrblocks_dia;

/* initialize the counter for block_dia blocks */
	  nr_bl_block = 0;

/* put block on temporary list of block_dia blocks  */
	  nr_bl_block++;
	  temp_list[it_matrix->Nrblocks_dia][nr_bl_block]=ibl;

/* loop over block's neighbors */
/* i.e. other blocks connected by entries in the global stiffness matrix */
	  for(iaux=0;iaux<it_matrix->Block[ibl]->Lngb[0];iaux++){

		ineig = it_matrix->Block[ibl]->Lngb[iaux+1];

/*kbw
printf("in neig %d (%d), LNGB %u, tag %d\n",
	   ineig,iaux,it_matrix->Block[ineig]->Lngb,tags[ineig]);
/*kew*/

/* if ordinary (not ghost) block and (not tagged yet or we include overlap) */
		if(  it_matrix->Block[ineig]->Lngb != NULL && (tags[ineig]==0 ||
			 abs(it_matrix->Block_type)==4 ||
			 abs(it_matrix->Block_type)==5 ) ){

/* put on temporary list of block_dia blocks  */
		  nr_bl_block++;
		  temp_list[it_matrix->Nrblocks_dia][nr_bl_block]=ineig;

/*kbw
printf("in neig %d (%d), put as block %d in block %d\n",
	ineig,iaux,nr_bl_block,it_matrix->Nrblocks_dia);
/*kew*/

/* tag */
		  if(abs(it_matrix->Block_type)==3 ||
		 abs(it_matrix->Block_type)==5 ||
		 abs(it_matrix->Block_type)==7 )
		tags[ineig] = it_matrix->Nrblocks_dia;


		} /* if neighbor included */

	  }  /* for all block's neighbors: iaux */

/* the total number of blocks in a diagonal block */
	  if(nr_bl_block>=LAC_MAX_BL_NGB){
		printf("Too much blocks in it_create_blocks_dia\n");
		printf("Change parameters for block creation or increase LAC_MAX_BL_NGB\n");
		return(-1);
	  }
	  temp_list[it_matrix->Nrblocks_dia][0] = nr_bl_block;

/*kbw
printf("Created block %d: %d blocks\n",
	it_matrix->Nrblocks_dia,nr_bl_block);
for(iaux=1;iaux<=nr_bl_block;iaux++)
	printf("%d   ",temp_list[it_matrix->Nrblocks_dia][iaux]);
printf("\n");
 getchar();
/*kew*/


	} /* if block ordinary and not tagged  */

	  } /* for all active blocks: ibl */

/*kbw
printf("Total number of blocks %d\n",it_matrix->Nrblocks_dia);
/*kew*/



/* allocate space for pointers to blocks */
	  it_matrix->Block_dia =
	(itt_blocks_dia **)malloc((it_matrix->Nrblocks_dia+1)
				 *sizeof(itt_blocks_dia *));

/* loop over diagonal blocks */
	  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

/* allocate space for block_dia structure */
	it_matrix->Block_dia[ib_dia] =
	  (itt_blocks_dia *)malloc(sizeof(itt_blocks_dia));
	if(it_matrix->Block_dia[ib_dia]==NULL){
		  printf("it_matrix gmres: no space for block %d\n", ib_dia);
		  exit(1);
	}

/* prepare diagonal block for assembling */
	block_dia = it_matrix->Block_dia[ib_dia];

/* renumbering: i_temp_list = it_matrix->Nrblocks_dia + 1 - ib_dia;  */
	i_temp_list = ib_dia;

 /* number of small blocks (blocks) in dia */
	iaux = temp_list[i_temp_list][0];
	block_dia->Lsmall=lar_util_ivector(iaux+1,"block_dia in create_blocks_dia");
	block_dia->Lpos=lar_util_ivector(iaux+1,"block_dia in create_blocks_dia");
	block_dia->Lelem=NULL;

/* number of small blocks */
	block_dia->Lsmall[0] = iaux;

/* compute the number of dofs for diagonal block */
	nrdof_dia = 0;

/* for each block in block_dia */
	for(ibl=1;ibl<=temp_list[i_temp_list][0];ibl++){

/* get small block number */
	  iblock=temp_list[i_temp_list][ibl];

/* put small block number */
	  block_dia->Lsmall[ibl] = iblock;

/* indicate local position of first dof in dia subarray */
	  block_dia->Lpos[ibl] = nrdof_dia;

/* update total number of dofs for dia */
	  nrdof_dia += it_matrix->Block[iblock]->Ndof;

/*kbw
printf("counting dofs of el %d (loc %d)\n",
	temp_list[i_temp_list][ibl],ibl);
printf("small block %d, nrdof %d, total nrdof %d\n",
	iblock,it_matrix->Block[iblock]->Ndof,nrdof_dia);
/*kew*/


	} /* for all  blocks in a block */


/* total number of dofs for dia */
	block_dia->Lpos[0] = nrdof_dia;

/* allocate space for dia and Ips */
	block_dia->Dia=
	  lar_util_dvector(nrdof_dia*nrdof_dia,"Dia in create_blocks_dia");
	block_dia->Ips=
	  lar_util_ivector(nrdof_dia,"Ips in create_blocks_dia");


/*kbw
  printf("In block %d, %d small blocks (blocks):\n",
	 ib_dia,iaux);
 for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){
   printf("%d (%d)  ",block_dia->Lsmall[iblock],
	  it_matrix->Mesh_elems[block_dia->Lsmall[iblock]]);
 }
 printf("\n");
/*kew*/

/*kbw
printf("allocated %d dofs in block %d\n",
	nrdof_dia*nrdof_dia, ib_dia);
getchar();
/*kew*/


	  } /* for all diagonal blocks: ib_dia */

	  for(i=0;i<=it_matrix->Nrblocks;i++) free(temp_list[i]);
	  free(temp_list);
	  free(tags);

	} /* if specific block type */


/* for block ILU factorization */
	if( it_matrix->Precon==MULTI_ILU ||  it_matrix->Precon==BLOCK_ILU ){
	  printf("large blocks not allowed for ILU preconditioning in it_bliter_blocks.c\n");
	  exit(-1);
	} /* end if block ILU preconditioner */
	else {
	  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){
	block_dia=it_matrix->Block_dia[ib_dia];
	block_dia->Lneig=NULL;
	  }
	}

  } /* if larger blocks */

/*kbw
  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){
	block_dia=it_matrix->Block_dia[ib_dia];
	printf("Big block %d, Neighbors:",ib_dia);
	for(i=1;i<=block_dia->Lneig[0];i++) printf("%d ",block_dia->Lneig[i]);
	printf("\n");
  }
/*kew*/

/* allocate space for arrays of doubles: Aux[i] */
  if( it_matrix->Precon==MULTI_ILU||it_matrix->Precon==BLOCK_ILU ){
	for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

	  block_dia = it_matrix->Block_dia[ib_dia];
	  nrdof_dia = block_dia->Lpos[0];

	  nrneig = block_dia->Lneig[0];
	  block_dia->Aux = (double **)malloc(nrneig*sizeof(double *));

	  for(i=0;i<nrneig;i++){
	j = block_dia->Lneig[i+1];

	kaux = it_matrix->Block_dia[j]->Lpos[0];

/*kbw
  if(ib_dia>0){
	printf("Big block %d, nr_dof %d neighbor %d, nr_dofs %d\n",
	   ib_dia,nrdof_dia,j,kaux);
  getchar();
  }
/*kew*/

	block_dia->Aux[i] =
	  lar_util_dvector(nrdof_dia*kaux,"Aux[] in create_blocks_dia");

	  }

	} /* for all diagonal blocks */
  } /* if ILU preconditioner */

/*kbw
  printf("Order of blocks:\n");
  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){
	block_dia=it_matrix->Block_dia[ib_dia];
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){
	  iaux=block_dia->Lsmall[iblock];
	  printf("block_dia %d block_small %d\n",ib_dia, iaux);
	  //printf("%d ",it_matrix->Mesh_elems[iaux]);
	}
 }
 printf("\n");
 getchar();
/*kew*/

  return(it_matrix->Nrblocks_dia);
}

/*---------------------------------------------------------
  lar_fill_preconditioner_block - to fill preconditioner (here by factorizing
							diagonal blocks or by block ILU(0) factorization)
---------------------------------------------------------*/
int lar_fill_preconditioner_block(
  int Matrix_id   /* in: matrix ID */
	)
{

  /* pointers to solver structure and stiffness matrix blocks' info */
  itt_matrices *it_matrix;

  int i,j,k;
  //double *a_work;

/*++++++++++++++++ executable statements ++++++++++++++++*/


  it_matrix = &itv_matrices[Matrix_id];

/*kbw
  printf("Entering lar_fill_preconditioner_block: Matrix_id %d, Block_type %d, Precon %d\n",
	 Matrix_id, itv_matrices[Matrix_id].Block_type, itv_matrices[Matrix_id].Precon);
  //getchar();
  //getchar();
/*kew*/

  /*kbw
  {

	itt_matrices *it_matrix;

	int iaux, jaux, kaux, iblock;

	it_matrix = &itv_matrices[Matrix_id];

	printf("before factorizing blocks_dia, Matrix_id %d\n", Matrix_id);
	for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){
	  jaux=it_matrix->Block[iblock]->Ndof;
	  printf("In block %d (%d rows, the first %d), %d blocks\n", iblock,
		 it_matrix->Block[iblock]->Ndof, it_matrix->Block[iblock]->Posg,
		 it_matrix->Block[iblock]->Lngb[0]+1);
	  if(it_matrix->Block[iblock]->Lngb!=NULL){

	printf("diagonal block (displayed by rows)\n");
	for(j=0;j<jaux;j++){
	  printf("diagonal block (displayed by rows): row %d\n", j);
	  for(k=0;k<jaux;k++){
		printf("%20.15lf", it_matrix->Block[iblock]->Dia[j+jaux*k]);
	  }
	  printf("\n");
	}

	iaux = it_matrix->Block[iblock]->Lngb[0];
	for(i=1;i<=iaux;i++){
	  kaux = it_matrix->Block[iblock]->Lngb[i];
	  printf("for block %d (from column %d) (entries displayed by rows!)\n",
		 kaux, it_matrix->Block[kaux]->Posg);
	  for(j=0;j<jaux;j++){
	  printf("for block %d (from column %d) (entries displayed by rows!): row %d\n",
		 kaux, it_matrix->Block[kaux]->Posg, j);

		for(k=0;k<jaux;k++){
		  printf("%20.15lf", it_matrix->Block[iblock]->Aux[i-1][j+jaux*k]);
		}
		printf("\n");
	  }
	}

	printf("RHS\n");
	for(j=0;j<jaux;j++){
	  printf("%20.15lf", it_matrix->Block[iblock]->Rhs[j]);
	}
	printf("\n");
	  }
	  getchar();
	}
  }
	/*kew*/

/*according to the value of the global variable it_matrix->Block_type */
  if(abs(itv_matrices[Matrix_id].Block_type)==0){

	int iblock;
	/* prepare block arrays */
#pragma omp parallel for default(none) firstprivate(Matrix_id, it_matrix) private(i,j,k)
	for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){


	  if(it_matrix->Block[iblock]!=NULL&&it_matrix->Block[iblock]->Lngb!=NULL) {

	itt_blocks *block;  /* to simplify typing */

		block = it_matrix->Block[iblock];
		int nrdofbl = block->Ndof;
		int nrngbl = block->Lngb[0];

	/* for standard iterations as preconditioner prepare inverted arrays */
		if(nrdofbl==1){
		  block->Rhs[0]=block->Rhs[0]/block->Dia[0];
		  for(i=0;i<nrngbl;i++)block->Aux[i][0]=block->Aux[i][0]/block->Dia[0];
		}
		else{
	  /* LU decompose block->Dia */

	  int iaux;
		  dgetrf_(&nrdofbl,&nrdofbl,block->Dia,&nrdofbl,block->Ips,&iaux);

	  /* resolve to find  block->Dia^(-1) * block->Rhs */
		  int kaux=1;

		  dgetrs_("N",&nrdofbl,&kaux,block->Dia,&nrdofbl,block->Ips,
		  block->Rhs,&nrdofbl,&iaux);

		  for(i=0;i<nrngbl;i++){

		/* find the number of dofs for the neighbor */
			int j = block->Lngb[i+1];
			kaux = it_matrix->Block[j]->Ndof;

		/* resolve to find  block->Dia^(-1) * block->Aux[i] */
			dgetrs_("N",&nrdofbl,&kaux,block->Dia,&nrdofbl,block->Ips,
			block->Aux[i],&nrdofbl,&iaux);

		  }
		}

	  } /* if block active */
	} /* for each block: iblock */

/*kbw
	  //if(it_matrix->Monitor>ERRORS)
	{

	  int iblock;

	  printf("\nAfter factorization of Dia blocks in block_dia:\n");
	  for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){
	//if(iblock==352)
	{

	  itt_blocks *block;
	  block = it_matrix->Block[iblock];
	  int nrdofbl = block->Ndof;
	  printf("Block %d, ndof %d, Posg %d\n   Neighbors:",
		 iblock,block->Ndof, block->Posg);
	  for(i=1;i<=block->Lngb[0];i++)
		printf("  %d",block->Lngb[i]);
	  printf("\n");
	  printf("Dia:\n");
	  for(i=0;i<nrdofbl*nrdofbl;i++){
		printf("%20.15lf",block->Dia[i]);
	  }
	  for(j=0;j<block->Lngb[0];j++){
		printf("\nAux[%d]:\n",j);
		for(i=0;i<nrdofbl*nrdofbl;i++){
		  printf("%20.15lf",block->Aux[j][i]);
		}
	  }
	  printf("\nRHS:\n");
	  for(i=0;i<nrdofbl;i++){
		printf("%20.15lf",block->Rhs[i]);
	  }
	  printf("\n\n");
	  getchar();
	}
	  }

	}
/*kew*/

  } /* end if Block_type=size == 0 (iterations with no dia blocks) */
  else { /* all other kinds of blocks */

	if(it_matrix->Precon!=MULTI_ILU&&it_matrix->Precon!=BLOCK_ILU){

	  int ib_dia;
/* for each diagonal block */
#pragma omp parallel for default(none) firstprivate(Matrix_id, it_matrix) private(i,j,k)
	  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

/* if block is active */
	if(it_matrix->Block_dia[ib_dia]!=NULL){

	  itt_blocks_dia *block_dia;/* to simplify typing */
/* to simplify typing */
	  block_dia=it_matrix->Block_dia[ib_dia];

/* assemble diagonal and suitable off diagonal submatrices from
   small blocks */
/* total number of dofs for DIA */
	  int nrdof_dia = block_dia->Lpos[0];

/* initialize Dia and ips matrices */
	  for(i=0;i<nrdof_dia*nrdof_dia;i++) block_dia->Dia[i]=0.0;
	  for(i=0;i<nrdof_dia;i++) block_dia->Ips[i]=0;

	  int ibs;
/* loop over small blocks */
	  for(ibs=1;ibs<=block_dia->Lsmall[0];ibs++){

/* get into small block */
		itt_blocks *block;  /* to simplify typing */
		block=it_matrix->Block[block_dia->Lsmall[ibs]];

/* number of dofs for small block */
		int ndof1 = block->Ndof;

/* position of first dof in diagonal submatrix */
		int posloc1 = block_dia->Lpos[ibs];

/* assemble diagonal subarray */
		for(i=0;i<ndof1;i++){
		  for(j=0;j<ndof1;j++){
		block_dia->Dia[posloc1+i+nrdof_dia*(posloc1+j)]
		  = block->Dia[i+ndof1*j];
		  }
		}

/* if more than one small block */
		if(block_dia->Lsmall[0]>1){

		  int iaux;
/* loop over small blocks */
		  for(iaux=1;iaux<=block_dia->Lsmall[0];iaux++){

/* if different from block ibs  */
		if(iaux!=ibs){

		  int ioff = lar_util_chk_list(block_dia->Lsmall[iaux],
					 &block->Lngb[1],block->Lngb[0]);

/* if there is an off-diagonal block */
		  if(ioff){

/* number of dofs for the off-diagonal submatrix */
			int ndof2 = it_matrix->Block[block_dia->Lsmall[iaux]]->Ndof;

/* position for the off-diagonal submatrix in dia */
			int posloc2 = block_dia->Lpos[iaux];

/* assemble off-diagonal subarray */
			for(i=0;i<ndof1;i++){
			  for(j=0;j<ndof2;j++){
			block_dia->Dia[posloc1+i+nrdof_dia*(posloc2+j)]
			  = block->Aux[ioff-1][i+ndof1*j];
			  }
			}

		  } /* if there is off-diagonal block */

		} /* if two different blocks: iaux != ibs */

		  } /* for all small blocks: iaux */

		} /* if more small blocks */

		} /* for all small blocks: ibs */

/* decompose Dia */

	  if(nrdof_dia==1) block_dia->Dia[0]= 1/block_dia->Dia[0];
	  else{

		int iaux;
/* LU decompose block->Dia */
		dgetrf_(&nrdof_dia,&nrdof_dia,block_dia->Dia,
			&nrdof_dia,block_dia->Ips,&iaux);

	  }

	} /* if block is active */

	  } /* for all blocks */

/*kbw
	  //if(it_matrix->Monitor>ERRORS)
	{

	  int iblock;

	printf("\nAfter factorization of Dia blocks in block_dia:\n");
	for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){
	  //if(iblock==352)
		{

	  itt_blocks *block;
	  block = it_matrix->Block[iblock];
	  int nrdofbl = block->Ndof;
	  printf("Block %d, ndof %d, Posg %d\n   Neighbors:",
		 iblock,block->Ndof, block->Posg);
	  for(i=1;i<=block->Lngb[0];i++)
		printf("  %d",block->Lngb[i]);
	  printf("\n");
	  printf("Dia:\n");
	  for(i=0;i<nrdofbl*nrdofbl;i++){
		printf("%20.15lf",block->Dia[i]);
	  }
	  for(j=0;j<block->Lngb[0];j++){
		printf("\nAux[%d]:\n",j);
		for(i=0;i<nrdofbl*nrdofbl;i++){
		  printf("%20.15lf",block->Aux[j][i]);
		}
	  }
	  printf("\nRHS:\n");
	  for(i=0;i<nrdofbl;i++){
		printf("%20.15lf",block->Rhs[i]);
	  }
	  printf("\n\n");
	  getchar();
	  }
	}

	for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

	  if(it_matrix->Block_dia[ib_dia]!=NULL)
	  //if(ib_dia==352)
		{

		itt_blocks_dia *block_dia;
		block_dia=it_matrix->Block_dia[ib_dia];

		int nrdofbl = block_dia->Lpos[0];
		int iaux=block_dia->Lsmall[0];
		printf("Block_dia %d\nSmall blocks:",ib_dia);
		for(i=1;i<=iaux;i++) printf("  %d",block_dia->Lsmall[i]);
		printf("\n");
		//int kaux=block_dia->Lneig[0];
		//printf("Neighbors:");
		//for(i=1;i<=kaux;i++) printf("  %d",block_dia->Lneig[i]);
		//printf("\n");
		printf("Dia:\n");
		for(i=0;i<nrdofbl*nrdofbl;i++){
		  printf("%20.15lf",block_dia->Dia[i]);
		}
		//printf("\n");
		//int nrneig = block_dia->Lneig[0];
		//int ineig;
		//for(ineig=0;ineig<block_dia->Lneig[0];ineig++){
		//  j = block_dia->Lneig[ineig+1];
		//  int kaux = it_matrix->Block_dia[j]->Lpos[0];
		//  printf("\nAux[%d]:\n",ineig);
		//  for(i=0;i<nrdofbl*kaux;i++)
		//	printf("%20.15lf",block_dia->Aux[ineig][i]);
		//}
		printf("\n\n");
		getchar();


	  }
	}
	  }
/*kew*/


	} /* end if standard iterations preconditioning */
	else if(it_matrix->Precon==MULTI_ILU||it_matrix->Precon==BLOCK_ILU){

	  /* prepare data for forward reduction and back substitution */
	  printf("Preparing data for iterations\n");

/*       /\* loop over all big diagonal blocks *\/ */
/* #pragma omp parallel for  private(block_dia,jb_dia,ioff) */
/*       for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){ */

/* 	/\* to simplify typing *\/ */
/* 	block_dia=it_matrix->Block_dia[ib_dia]; */

/* 	block_dia->Llowerneig = (int *) malloc (200*sizeof(int)); */
/* 	block_dia->Llowerneig[0] = 0; */


/* 	//for all neighbours */
/* 	int iaux_neig; */
/* 	for(iaux_neig=1;iaux_neig<=block_dia->Lneig[0];iaux_neig++){  */

/* 	  jb_dia = block_dia->Lneig[iaux_neig]; */

/* 	  ioff=iaux_neig; */

/* 	  if(jb_dia<ib_dia){ */

/* 	/\* /\\* for all factorized (lower subdiagonal) entries of A *\\/ *\/ */
/* 	/\* for(jb_dia=1;jb_dia<ib_dia;jb_dia++){ *\/ */

/* 	/\*   /\\* if two blocks are neighbors *\\/ *\/ */
/* 	/\*   ioff = lar_util_chk_list(jb_dia,&block_dia->Lneig[1],block_dia->Lneig[0]); *\/ */
/* 	/\*   if(ioff){ *\/ */

/* #ifdef DEBUG */
/* 	    if(block_dia->Lneig[ioff]!=jb_dia){ */
/* 	  printf("error in reading list of lower neighbors in decompose!\n"); */
/* 	  printf("ib_dia %d, jneig %d, jb_dia %d, ioff %d, ? %d\n", */
/*       ib_dia, block_dia->Llowerneig[0]+1, jb_dia, ioff, block_dia->Lneig[ioff]); */
/* 	      exit(1); */
/* 	    } */
/* 	    if(block_dia->Llowerneig[0]>=99){ */
/* 	      printf("increase the number of lower subdiagonal neighbors!\n"); */
/* 	      exit(1); */
/* 	    } */
/* #endif */

/* 	    block_dia->Llowerneig[0]++; */
/* 	    //block_dia->Llowerneig[2*block_dia->Llowerneig[0]-1] = jb_dia; */
/* 	    //block_dia->Llowerneig[2*block_dia->Llowerneig[0]] = ioff; */

/* 	    // we must sort an array */
/* 	    i=2*block_dia->Llowerneig[0]-1; */
/* 	    while(i>0 && jb_dia<block_dia->Llowerneig[i-2]){ */
/* 	      block_dia->Llowerneig[i] = block_dia->Llowerneig[i-2]; */
/* 	      block_dia->Llowerneig[i+1] = block_dia->Llowerneig[i-1]; */
/* 	      i -= 2; */
/* 	    } */
/* 	    block_dia->Llowerneig[i] = jb_dia; */
/* 	    block_dia->Llowerneig[i+1] = ioff; */


/* 	  } */
/* 	} */
/*       } */

	  int ib_dia;
	  /* loop over all big diagonal blocks */
#pragma omp parallel for default(none) firstprivate(Matrix_id, it_matrix) private(i,j,k)
	  //for(ib_dia=it_matrix->Nrblocks_dia;ib_dia>=1;ib_dia--){
	  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

/* if block not active */
	if(it_matrix->Block_dia[ib_dia]==NULL){
	  printf("Block %d not active in it_factor_block_dia\n",ib_dia);
	  exit(0);
	}

	itt_blocks_dia *block_dia;/* to simplify typing */
		/* to simplify typing */
	block_dia=it_matrix->Block_dia[ib_dia];

	block_dia->Lupperneig = (int *) malloc (200*sizeof(int));
	block_dia->Lupperneig[0] = 0;

	block_dia->Llowerneig = (int *) malloc (200*sizeof(int));
	block_dia->Llowerneig[0] = 0;

	//for all neighbours
	int iaux_neig;
	for(iaux_neig=1;iaux_neig<=block_dia->Lneig[0];iaux_neig++){

	  int jb_dia = block_dia->Lneig[iaux_neig];

	  int ioff=iaux_neig;

	  if(jb_dia>ib_dia){


		/* for all upper diagonal blocks */
		//for(jb_dia=ib_dia+1;jb_dia<=it_matrix->Nrblocks_dia;jb_dia++){

		/* if two blocks are neighbors */
		//ioff = lar_util_chk_list(jb_dia,&block_dia->Lneig[1],block_dia->Lneig[0]);

		//if(ioff){

#ifdef DEBUG
		if(block_dia->Lneig[ioff]!=jb_dia){
		  printf("error in reading list of upper neighbors in decompose!\n");
		  printf("ib_dia %d, jneig %d, jb_dia %d, ioff %d, ? %d\n",
			 ib_dia, block_dia->Lupperneig[0]+1, jb_dia, ioff,
			 block_dia->Lneig[ioff]);
		  exit(1);
		}
		if(block_dia->Lupperneig[0]>=99){
		  printf("increase the number of upper subdiagonal neighbors!\n");
		  exit(1);
		}
#endif

		block_dia->Lupperneig[0]++;
		//block_dia->Lupperneig[2*block_dia->Lupperneig[0]-1] = jb_dia;
		//block_dia->Lupperneig[2*block_dia->Lupperneig[0]] = ioff;

		// we must sort an array
		i=2*block_dia->Lupperneig[0]-1;
		while(i>2 && jb_dia<block_dia->Lupperneig[i-2]){
		  block_dia->Lupperneig[i] = block_dia->Lupperneig[i-2];
		  block_dia->Lupperneig[i+1] = block_dia->Lupperneig[i-1];
		  i -= 2;
		}
		block_dia->Lupperneig[i] = jb_dia;
		block_dia->Lupperneig[i+1] = ioff;

#ifdef DEBUG
		if(block_dia->Lneig[ioff]!=jb_dia){
		  printf("error in reading list of upper neighbors in decompose!\n");
		  printf("ib_dia %d, jneig %d, jb_dia %d, ioff %d, ? %d\n",
			 ib_dia, block_dia->Lupperneig[0]+1, jb_dia, ioff,
			 block_dia->Lneig[ioff]);
		  exit(1);
		}
		if(block_dia->Lupperneig[0]>=99){
		  printf("increase the number of lower subdiagonal neighbors!\n");
		  exit(1);
		}
/*kbw
		if(ib_dia <3){
		printf("ib_dia %d, iaux_neig %d, jb_dia %d, i %d\n",
			   ib_dia , iaux_neig, jb_dia, i);
		}
//*kew*/
#endif

	  }
	  else if(jb_dia<ib_dia){

#ifdef DEBUG
		if(block_dia->Lneig[ioff]!=jb_dia){
		  printf("error in reading list of lower neighbors in decompose!\n");
		  printf("ib_dia %d, jneig %d, jb_dia %d, ioff %d, ? %d\n",
			 ib_dia, block_dia->Llowerneig[0]+1, jb_dia, ioff,
			 block_dia->Lneig[ioff]);
		  exit(1);
		}
		if(block_dia->Llowerneig[0]>=99){
		  printf("increase the number of lower subdiagonal neighbors!\n");
		  exit(1);
		}
#endif

		block_dia->Llowerneig[0]++;
		//block_dia->Llowerneig[2*block_dia->Llowerneig[0]-1] = jb_dia;
		//block_dia->Llowerneig[2*block_dia->Llowerneig[0]] = ioff;
		// we must sort the array
		i=2*block_dia->Llowerneig[0]-1;
		while(i>2 && jb_dia<block_dia->Llowerneig[i-2]){
		  block_dia->Llowerneig[i] = block_dia->Llowerneig[i-2];
		  block_dia->Llowerneig[i+1] = block_dia->Llowerneig[i-1];
		  i -= 2;
		}
		block_dia->Llowerneig[i] = jb_dia;
		block_dia->Llowerneig[i+1] = ioff;

#ifdef DEBUG
		if(block_dia->Lneig[ioff]!=jb_dia){
		  printf("error in reading list of lower neighbors in decompose!\n");
		  printf("ib_dia %d, jneig %d, jb_dia %d, ioff %d, ? %d\n",
			 ib_dia, block_dia->Llowerneig[0]+1, jb_dia, ioff,
			 block_dia->Lneig[ioff]);
		  exit(1);
		}
		if(block_dia->Llowerneig[0]>=99){
		  printf("increase the number of lower subdiagonal neighbors!\n");
		  exit(1);
		}
/*kbw
		if(ib_dia <3){
		printf("ib_dia %d, iaux_neig %d, jb_dia %d, i %d\n",
			   ib_dia , iaux_neig, jb_dia, i);
		}
//*kew*/
#endif



	  }

	}
	  }



/* ASSEMBLE STIFFNESS MATRIX INTO BIG BLOCKS */
	  printf("Assembling preconditioner arrays\n");

/* for each diagonal block */
#pragma omp parallel for default(none) firstprivate(Matrix_id, it_matrix) private(i,j,k)
	  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

/* if block is active */
	if(it_matrix->Block_dia[ib_dia]!=NULL){

/* to simplify typing */
	  itt_blocks_dia *block_dia;/* to simplify typing */
	  block_dia=it_matrix->Block_dia[ib_dia];

/* total number of dofs for DIA */
	  int nrdof_dia = block_dia->Lpos[0];

/* assemble diagonal and suitable off diagonal submatrices from
   small blocks into Dia array for big block */
/* initialize Dia and ips matrices */
	  for(i=0;i<nrdof_dia*nrdof_dia;i++) block_dia->Dia[i]=0.0;
	  for(i=0;i<nrdof_dia;i++) block_dia->Ips[i]=0;

/* loop over small blocks */
	  int ibs;
	  for(ibs=1;ibs<=block_dia->Lsmall[0];ibs++){

/* get into small block */
		itt_blocks *block;  /* to simplify typing */
		block=it_matrix->Block[block_dia->Lsmall[ibs]];

/* number of dofs for small block */
		int ndof1 = block->Ndof;

/* position of first dof in diagonal submatrix */
		int posloc1 = block_dia->Lpos[ibs];

/* assemble diagonal subarray */
		for(i=0;i<ndof1;i++){
		  for(j=0;j<ndof1;j++){
		block_dia->Dia[posloc1+i+nrdof_dia*(posloc1+j)]
		  = block->Dia[i+ndof1*j];
		  }
		}

/* if more than one small block */

		if(block_dia->Lsmall[0]>1){

		  printf("Large blocks not working in ILU!!! Exiting\n");
		  exit(-1);

		} /* if more small blocks */

	  } /* for all small blocks: ibs */

/* assemble diagonal and suitable off diagonal submatrices from
   small blocks into Aux arrays for big block */

/* initialize Aux matrices */
/* loop over big block neighbors */
	  int nrneig = block_dia->Lneig[0];
	  int ineig;
	  for(ineig=0;ineig<block_dia->Lneig[0];ineig++){
		j = block_dia->Lneig[ineig+1];
		int kaux = it_matrix->Block_dia[j]->Lpos[0];
		for(i=0;i<nrdof_dia*kaux;i++) block_dia->Aux[ineig][i]=0.0;
	  }


/* loop over small blocks */
	  for(ibs=1;ibs<=block_dia->Lsmall[0];ibs++){

/* get into small block */
		itt_blocks *block;  /* to simplify typing */
		block=it_matrix->Block[block_dia->Lsmall[ibs]];

/* number of dofs for small block */
		int ndof1 = block->Ndof;

/* position of first dof in diagonal submatrix */
		int posloc1 = block_dia->Lpos[ibs];

/* loop over neighboring big blocks */
		int ineig;
		for(ineig=0;ineig<block_dia->Lneig[0];ineig++){
		  int j = block_dia->Lneig[ineig+1];
		  itt_blocks_dia *dia_neig;
		  dia_neig = it_matrix->Block_dia[j];

/* loop over small blocks from dia_neig */
		  int iaux;
		  for(iaux=1;iaux<=dia_neig->Lsmall[0];iaux++){
		int ismneig=dia_neig->Lsmall[iaux];

/* check whether two small blocks are neighbors (non-zero entries in A) */
		int ioff = lar_util_chk_list(ismneig, &block->Lngb[1],block->Lngb[0]);

/* if there is a common off-diagonal block */
		if(ioff){

/* number of dofs for the off-diagonal submatrix */
		  int ndof2 = it_matrix->Block[ismneig]->Ndof;

/* position for the off-diagonal submatrix in dia */
		  int posloc2 = dia_neig->Lpos[iaux];

/* assemble off-diagonal subarray */
		  for(i=0;i<ndof1;i++){
			for(j=0;j<ndof2;j++){
			  block_dia->Aux[ineig][posloc1+i+nrdof_dia*(posloc2+j)]
			= block->Aux[ioff-1][i+ndof1*j];
			}
		  }

		} /* if there is off-diagonal block */

		  } /* for all small blocks of a nieghbor: iaux */

		} /* for all big neighboring blocks */

	  } /* for all small blocks: ibs */

	} /* if block is active */

	  } /* for all blocks */


/*kbw
	  //if(it_matrix->Monitor>ERRORS)
	{

	  int iblock;

	for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){
	  //if(iblock==352)
		{

	  itt_blocks *block;
	  block = it_matrix->Block[iblock];
	  int nrdofbl = block->Ndof;
	  printf("Block %d, ndof %d, Posg %d\n   Neighbors:",
		 iblock,block->Ndof, block->Posg);
	  for(i=1;i<=block->Lngb[0];i++)
		printf("  %d",block->Lngb[i]);
	  printf("\n");
	  printf("Dia:\n");
	  for(i=0;i<nrdofbl*nrdofbl;i++){
		printf("%25.15lf",block->Dia[i]);
	  }
	  for(j=0;j<block->Lngb[0];j++){
		printf("\nAux[%d]:\n",j);
		for(i=0;i<nrdofbl*nrdofbl;i++){
		  printf("%25.15lf",block->Aux[j][i]);
		}
	  }
	  printf("\n\n");
	  getchar();
	  }
	}

	for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

	  if(it_matrix->Block_dia[ib_dia]!=NULL)
	  //if(ib_dia==352)
		{

		itt_blocks_dia *block_dia;
		block_dia=it_matrix->Block_dia[ib_dia];

		int nrdofbl = block_dia->Lpos[0];
		int iaux=block_dia->Lsmall[0];
		printf("Block_dia %d\nSmall blocks:",ib_dia);
		for(i=1;i<=iaux;i++) printf("  %d",block_dia->Lsmall[i]);
		printf("\n");
		int kaux=block_dia->Lneig[0];
		printf("Neighbors:");
		for(i=1;i<=kaux;i++) printf("  %d",block_dia->Lneig[i]);
		printf("\n");
		printf("Dia:\n");
		for(i=0;i<nrdofbl*nrdofbl;i++){
		  printf("%20.15lf",block_dia->Dia[i]);
		}
		printf("\n");
		int nrneig = block_dia->Lneig[0];
		int ineig;
		for(ineig=0;ineig<block_dia->Lneig[0];ineig++){
		  j = block_dia->Lneig[ineig+1];
		  int kaux = it_matrix->Block_dia[j]->Lpos[0];
		  printf("\nAux[%d]:\n",ineig);
		  for(i=0;i<nrdofbl*kaux;i++)
		printf("%20.15lf",block_dia->Aux[ineig][i]);
		}
		printf("\n\n");
		getchar();


	  }
	}
	  }
/*kew*/

/* PERFORM INCOMPLETE LU FACTORIZATION */
	  printf("Starting factorization\n");
	  /* constants */
	  int ione = 1;
	  double done = 1.0;
	  double dmone = -1.0;


/* for each block row (diagonal block) */
	  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

/*kbw
	printf("block %d (%lf\%)\n",ib_dia, (double)ib_dia*100/it_matrix->Nrblocks_dia);
/*kew*/


/* if block is active */
	if(it_matrix->Block_dia[ib_dia]!=NULL){

/* to simplify typing */
	  itt_blocks_dia *block_dia;
	  block_dia=it_matrix->Block_dia[ib_dia];

/* total number of dofs for DIA */
	  int nrdof_dia = block_dia->Lpos[0];

/* for all rows (blocks) already considered */
	  //for(kb_dia=1;kb_dia<ib_dia;kb_dia++){
	  //double *a_work = lar_util_dvector(lwork,"A_work in factor_dia");


/* if two blocks are neighbors */
	  //ioff = lar_util_chk_list(kb_dia,&block_dia->Lneig[1],block_dia->Lneig[0]);
	  //if(ioff){

	  int iaux_dia;
	  for(iaux_dia=1;iaux_dia<=block_dia->Llowerneig[0];iaux_dia++){

		int kb_dia = block_dia->Llowerneig[2*iaux_dia-1];
		int ioff = block_dia->Llowerneig[2*iaux_dia];

		{


#ifdef DEBUG
		  if(block_dia->Lneig[ioff]!=kb_dia){
		printf("error in reading list of neighbors in decompose! ib_dia %d, iaux_dia %d (<%d), ioff %d, Lneig[ioff] %d, kb_dia %d \n",
				ib_dia, iaux_dia, block_dia->Llowerneig[0], ioff, block_dia->Lneig[ioff], kb_dia);
		exit(1);
		  }
#endif

		  itt_blocks_dia *block_k;
		  block_k = it_matrix->Block_dia[kb_dia];

/* compute product a_ik = a_ik * a_kk^-1 */
		  int kaux = block_k->Lpos[0];

/*kbw
		  //if(ib_dia==352)
		{
		printf("\nkb_dia %d, Dia before solution:\n",kb_dia);
		for(i=0;i<kaux*kaux;i++){
		  printf("%25.15lf",block_k->Dia[i]);
		}
		printf("\n");
		printf("ib_dia %d, ioff %d, Aux before solution:\n",ib_dia,ioff);
		for(i=0;i<kaux*nrdof_dia;i++){
		  printf("%25.15lf",block_dia->Aux[ioff-1][i]);
		}
		printf("\n");
		  }
/*kew*/

		  double *a_work = lar_util_dvector(kaux*nrdof_dia,"A_work in factor_dia");

		  for(i=0;i<nrdof_dia;i++){
		for(j=0;j<kaux;j++){
		  double daux=0.0;
		  for(k=0;k<kaux;k++){
			daux +=
			  block_dia->Aux[ioff-1][i+nrdof_dia*k]*
			  block_k->Dia[k+kaux*j];
		  }
		  //block_dia->Aux[ioff-1][i+nrdof_dia*j]=daux;
		  a_work[i+nrdof_dia*j] = daux;
		}
		  }
		  for(i=0;i<nrdof_dia;i++){
			for(j=0;j<kaux;j++){
			  block_dia->Aux[ioff-1][i+nrdof_dia*j] = a_work[i+nrdof_dia*j];
			}
		  }

		  free(a_work);

/*kb!!!!!!!!!?????????????
		  dgetrs_("N",&kaux,&nrdof_dia,block_k->Dia,&kaux,block_k->Ips,
			  block_dia->Aux[ioff-1],&kaux,&iaux);
/*kb!!!!!!!!!?????????????*/

/*kbw
		  //if(ib_dia==352)
		  {
		printf("i %d (%d), k %d (%d), Aux after solution:\n",
			   ib_dia,nrdof_dia,kb_dia,kaux);
		for(i=0;i<kaux*nrdof_dia;i++){
		  printf("%25.15lf",block_dia->Aux[ioff-1][i]);
		}
		printf("\n");
		//getchar(); getchar();
		  }
/*kew*/


/* update Dia block */
		  int koff = lar_util_chk_list(ib_dia,&block_k->Lneig[1],
				 block_k->Lneig[0]);

/*kbw
		  //if(ib_dia==352)
		  {
		printf("i %d (%d), k %d (%d), block_k->Aux for multiplication\n",
			   ib_dia,nrdof_dia,kb_dia,kaux);
		for(i=0;i<kaux*nrdof_dia;i++){
		  printf("%25.15lf",block_k->Aux[koff-1][i]);
		}
		printf("\n");
		  }
/*kew*/
		  for(i=0;i<nrdof_dia;i++){
		for(j=0;j<nrdof_dia;j++){
		  double daux=0.0;
		  for(k=0;k<kaux;k++){
			daux +=
			  block_dia->Aux[ioff-1][i+nrdof_dia*k]*
			  block_k->Aux[koff-1][k+kaux*j];
		  }
		  block_dia->Dia[i+nrdof_dia*j] -= daux;
		}
		  }

/*kbw
		  //if(ib_dia==352)
		  {
		printf("ib_dia %d, Dia after update:\n",ib_dia);
		for(i=0;i<nrdof_dia*nrdof_dia;i++){
		  printf("%25.15lf",block_dia->Dia[i]);
		}
		printf("\n");
		//getchar();
		  }
/*kew*/

/* loop over all neighbors of dia block, with indices greater than k */
		  int nrneig = block_dia->Lneig[0];
		  int ineig;
#pragma omp parallel for default(none) private(i,j,k) \
		  firstprivate(block_dia, kb_dia, it_matrix, block_k, nrdof_dia, kaux, ioff)
		  for(ineig=0;ineig<block_dia->Lneig[0];ineig++){
		int jb_dia = block_dia->Lneig[ineig+1];
		if(jb_dia>kb_dia){

		  itt_blocks_dia *dia_neig;/* to simplify typing */
		  dia_neig = it_matrix->Block_dia[jb_dia];

/* if k-block and j-block are neighbors, i.e. Aux blocks (j,k) are non-zero */
		  int koff = lar_util_chk_list(jb_dia,&block_k->Lneig[1],
					 block_k->Lneig[0]);
		  if(koff){

			int jaux=dia_neig->Lpos[0];

/*kbw
  printf("\ni %d (%d), j %d (%d), a_ij before solution:\n",
	 ib_dia,nrdof_dia,jb_dia,jaux);
  for(i=0;i<jaux*nrdof_dia;i++){
	printf("%25.15lf",block_dia->Aux[ineig][i]);
  }
  printf("\n");
  printf("i %d (%d), k %d (%d), a_ik before solution:\n",
	 ib_dia,nrdof_dia,kb_dia,kaux);
  for(i=0;i<kaux*nrdof_dia;i++){
	printf("%25.15lf",block_dia->Aux[ioff-1][i]);
  }
  printf("\n");
  printf("k %d (%d), j %d (%d), a_kj before solution:\n",
	 kb_dia,kaux,jb_dia,jaux);
  for(i=0;i<kaux*jaux;i++){
	printf("%25.15lf",block_k->Aux[koff-1][i]);
  }
  printf("\n");
/*kew*/

/* compute a_ij = a_ij - a_ik*a_kj */

/*kb!!!!!!!!!?????????????
			dgemm_('N','N',&nrdof_dia,&jaux,&kaux,
			   &dmone,block_dia->Aux[ioff-1],&nrdof_dia,
			   block_k->Aux[koff-1],&kaux,
			   &done,block_dia->Aux[ineig],&nrdof_dia);
/*kb!!!!!!!!!?????????????*/

			for(i=0;i<nrdof_dia;i++){
			  for(j=0;j<jaux;j++){
			double daux=0.0;
			for(k=0;k<kaux;k++){
			  daux +=
				block_dia->Aux[ioff-1][i+nrdof_dia*k]*
				block_k->Aux[koff-1][k+kaux*j];
			}
			block_dia->Aux[ineig][i+nrdof_dia*j] -= daux;
			  }
			}

/*kbw
  printf("i %d, j %d, a_ij after solution:\n",ib_dia,jb_dia);
  for(i=0;i<jaux*nrdof_dia;i++){
	printf("%25.15lf",block_dia->Aux[ineig][i]);
  }
  printf("\n");
/*kew*/
/*kbw
  printf("i %d, k %d, a_ik after solution:\n",ib_dia,kb_dia);
  for(i=0;i<kaux*nrdof_dia;i++){
	printf("%20.15lf",block_dia->Aux[ioff-1][i]);
  }
  printf("\n");
  printf("k %d, j %d, a_kj after solution:\n",kb_dia,jb_dia);
  for(i=0;i<kaux*jaux;i++){
	printf("%20.15lf",block_k->Aux[koff-1][i]);
  }
  printf("\n");
  getchar();
/*kew*/

		  } /* end if non-zero entry */

		}

		  } /* end parallel for all neighbors of dia block */

		} /* end if non-zero entry */

		//free(a_work);

	  } /* end for all blocks already considered:  kb_dia (columns k)

/* decompose Dia */
	  if(nrdof_dia==1) block_dia->Dia[0]= 1/block_dia->Dia[0];
	  else{
/* LU decompose block->Dia */

/*kbw
		//if(ib_dia==352)
		  {
		  printf("\nib_dia %d, Dia before inverting:\n",ib_dia);
		  for(i=0;i<nrdof_dia*nrdof_dia;i++){
		printf("%25.15lf",block_dia->Dia[i]);
		  }
		  printf("\n");
		}
/*kew*/

		double *a_work = lar_util_dvector((nrdof_dia+1)*nrdof_dia,"A_work in factor_dia");
		for(i=0;i<(nrdof_dia+1)*nrdof_dia;i++) a_work[i]=0.0;

		/* lar_util_dgetrf(block_dia->Dia,nrdof_dia,block_dia->Ips); */

		/* int iaux; */
		/* for(iaux=0;iaux<nrdof_dia;iaux++){ */
		/*   for(i=(iaux+1)*nrdof_dia;i<(iaux+2)*nrdof_dia;i++){ */
		/* 	a_work[i]=0.0; */
		/*   } */
		/*   a_work[(iaux+1)*nrdof_dia+iaux]=1.0; */
		/*   lar_util_dgetrs(block_dia->Dia,nrdof_dia, */
		/* 		&a_work[(iaux+1)*nrdof_dia], */
		/* 		&a_work[iaux*nrdof_dia],block_dia->Ips); */
		/* } */

		/* for(i=0;i<nrdof_dia*nrdof_dia;i++){ */
		/*   block_dia->Dia[i]=a_work[i]; */
		/* } */

/*kb*/
		int iaux;
		dgetrf_(&nrdof_dia,&nrdof_dia,block_dia->Dia,
			&nrdof_dia,block_dia->Ips,&iaux);

		int lwork = nrdof_dia*nrdof_dia;
		dgetri_(&nrdof_dia,block_dia->Dia,
			&nrdof_dia,block_dia->Ips,a_work,&lwork,&iaux);
/*kb*/

		free(a_work);

/*kbw
		if(ib_dia>236)
		{

		  printf("ib_dia %d, Dia after inverting:\n",ib_dia);
		  for(i=0;i<nrdof_dia*nrdof_dia;i++){
		printf("%25.15lf",block_dia->Dia[i]);
		  }
		  printf("\n");
		  getchar();
		}
/*kew*/


	  }

	} /* if block is active */

	  } /* for all blocks */




	} /* end if ILU preconditioning */

  } /* end if block size > 0 (i.e. there is Block_dia structure) */

  return(0);

}


/*---------------------------------------------------------
  lar_free_preconditioner_block - to free space for a block structure
---------------------------------------------------------*/
int lar_free_preconditioner_block(
  int Matrix_id   /* in: matrix ID */
  )
{

/* local variables */
  itt_matrices *it_matrix;
  itt_blocks_dia *block_dia;  /* to simplify typing */
  int ib_dia, ineig, nrneig;  /* counters for blocks */

/*++++++++++++++++ executable statements ++++++++++++++++*/


  it_matrix = &itv_matrices[Matrix_id];

  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

	block_dia = it_matrix->Block_dia[ib_dia];

	if( it_matrix->Precon==MULTI_ILU||it_matrix->Precon==BLOCK_ILU ){
	  nrneig = block_dia->Lneig[0];
	  for(ineig=0;ineig<nrneig;ineig++){
	free(block_dia->Aux[ineig]);
	  }
	  free(block_dia->Aux);
	}

	free(block_dia->Dia);
	free(block_dia->Ips);

	free(block_dia->Lsmall);

	free(block_dia->Lpos);

	if( it_matrix->Precon==MULTI_ILU||it_matrix->Precon==BLOCK_ILU ){
	  free(block_dia->Lneig);
	  free(block_dia->Llowerneig);
	  free(block_dia->Lupperneig);
	}

	free(it_matrix->Block_dia[ib_dia]);

  } /* for all diagonal blocks */


  free(it_matrix->Block_dia);

  return(0);
}


/*---------------------------------------------------------
  lar_free_SM_and_LV_block - to free space for a block structure
---------------------------------------------------------*/
int lar_free_SM_and_LV_block(
  int Matrix_id   /* in: matrix ID */
  )
{

/* local variables */
  itt_matrices *it_matrix;
  itt_blocks *block;	/* to simplify typing */

/* auxiliary variable */
  int i, iblock, bloff, nrblocks, nrdofbl, nroffbl;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Matrix_id != itv_nr_matrices-1){

	mf_fatal_err("Requested deallocation of matrix with id %d - which is not the last (of %d matrices). ",
		Matrix_id, itv_nr_matrices);
  }

  it_matrix = &itv_matrices[Matrix_id];

  /* allocate space for pointers to blocks */
  nrblocks = it_matrix->Nrblocks;

  for(iblock=1;iblock<=nrblocks;iblock++){

	block = it_matrix->Block[iblock];

	if(block->Lngb!=NULL){

	  for(i=1;i<=block->Lngb[0];i++){
	free(block->Aux[i-1]);
	  }
	  free(block->Dia);
	  free(block->Aux);
	  free(block->Rhs);
	  free(block->Ips);

	  free(block->Lngb);

	}

	free(it_matrix->Block[iblock]);

  }

  free(it_matrix->Block);

  itv_nr_matrices--;


  return(0);
}


/*---------------------------------------------------------
lar_compute_residual_block - to compute the residual of the system of equations,
	v = ( b - Ax ) (used also to compute the product v = -Ax)
---------------------------------------------------------*/
void lar_compute_residual_block (
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
  double dmone = -1.;

/* auxiliary variables */
  itt_blocks *block;  	/* to simplify typing */
  int iblock;		/* counters for blocks */
  int ioffbl;		/* counter for off diagonal blocks */
  int nrdofbl;		/* numbers of dofs in a block */
  int ndofngb;		/* number of dofs for a neighbor */

  int i,j;
  int iaux, jaux, kaux, laux, iblaux;
  double daux;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  itt_matrices *it_matrix = &itv_matrices[Matrix_id];

/*kbw
  printf("compres X on entrance\n" );
  for(i=0;i<Ndof;i++) printf("%20.15lf",X[i]);
  printf("\n" );
  getchar();
  getchar();
/*kew*/

/* perform matrix vector multiplication using elementary blocks */
/* loop over all small blocks */
  for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){

/* to simplify typing */
	block = it_matrix->Block[iblock];

   /* if regular (not ghost) block */
   if(block->Lngb!=NULL){

/* number of degrees of freedom for the block */
	nrdofbl = block->Ndof;

/* position of the first iblock dof in the vector v */
	kaux = block->Posg;

/* initialize part of v */
	for(i=0;i<nrdofbl;i++) V[kaux+i]=0.;

	if(!Ini_zero){
/* if initial guess not zero */

/* loop over neighbors */
	  for(ioffbl=0;ioffbl<block->Lngb[0];ioffbl++){

/* global number of the subsequent neighbor */
	iaux = block->Lngb[ioffbl+1];

/* number of degrees of freedom for the neighbor */
	ndofngb = it_matrix->Block[iaux]->Ndof;

/* position of the first ingb dof in the global v vector */
	jaux = it_matrix->Block[iaux]->Posg;

/*kbw
	printf("compres for blocks %d, %d before: X - \n", iblock, iaux);
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",X[jaux+i]);
	printf("\n");
	printf("compres before: V - \n");
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",V[kaux+i]);
	printf("\n");
	printf("compres before: SM - \n");
	for(i=0;i<nrdofbl*nrdofbl;i++) printf("%20.15lf",block->Aux[ioffbl][i]);
	printf("\n");
/*kew*/


/* multiplication by the off-diagonal block */
	/* if(nrdofbl*ndofngb>MIN_DOF_DGEMV){ */
	/*   dgemv_("N", &nrdofbl, &ndofngb, &dmone,  */
	/* 	 block->Aux[ioffbl], &nrdofbl, &X[jaux], &ione,  */
	/* 	 &done, &V[kaux], &ione); */
	/* } */
	/* else { /\* if problem too small to use dgemv *\/ */
	  for(i=0;i<nrdofbl;i++){
		for(j=0;j<ndofngb;j++){
		  V[kaux+i] -= block->Aux[ioffbl][i+j*nrdofbl]*X[jaux+j];
		}
	  }
	  //	}
/*kbw
	printf("compres after: V - \n");
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",V[kaux+i]);
	printf("\n");
	getchar();
	getchar();
/*kew*/


	  }

/*kbw
	printf("compres for blocks %d, %d before: X - \n", iblock, iblock);
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",X[kaux+i]);
	printf("\n");
	printf("compres before: V - \n");
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",V[kaux+i]);
	printf("\n");
	printf("compres before: SM - \n");
	for(i=0;i<nrdofbl*nrdofbl;i++) printf("%20.15lf",block->Dia[i]);
	printf("\n");
/*kew*/

/* multiplication with diagonal block - result stored in vtemp*/
	  /* if(nrdofbl*nrdofbl>MIN_DOF_DGEMV){ */
	  /* 	dgemv_("N", &nrdofbl, &nrdofbl, &dmone,  */
	  /* 	       block->Dia, &nrdofbl, &X[kaux], &ione,  */
	  /* 	       &done, &V[kaux], &ione); */
	  /* } */
	  /* else { /\* if problem too small to use dgemv *\/ */
	for(i=0;i<nrdofbl;i++){
	  for(j=0;j<nrdofbl;j++){
		V[kaux+i] -= block->Dia[i+j*nrdofbl]*X[kaux+j];
	  }
	}
	//      }

	} /* end if not Ini_zero (not zero initial guess X) */

/*kbw
	printf("compres before RHS update: subsequent rows %d-%d - V \n", kaux, kaux+nrdofbl);
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",V[kaux+i]);
	printf("\n");
/*kew*/

/* update using the rhs vector */
	if(Use_rhs==1){
	  if(B==NULL){
	daxpy_(&nrdofbl, &done, block->Rhs, &ione,
		   &V[kaux], &ione);
	  }
	  else{
	daxpy_(&nrdofbl, &done, &B[kaux], &ione,
		   &V[kaux], &ione);
	  }
	}

/*kbw
	printf("compres: subsequent rows %d-%d - B and V \n", kaux, kaux+nrdofbl);
	for(i=0;i<nrdofbl;i++) printf("%20.15lf",V[kaux+i]);
	printf("\n");
	if(B==NULL){
	  for(i=0;i<nrdofbl;i++) printf("%20.15lf",block->Rhs[i]);
	  printf("\n");
	}else{
	  for(i=0;i<nrdofbl;i++) printf("%20.15lf",B[kaux+i]);
	  printf("\n");
	}
	getchar();
	getchar();
/*kew*/

   }   /* end if regular (not ghost) block */
  } /* end loop over blocks */

/*kbw
	printf("compres on leaving: V \n");
	for(i=0;i<Ndof;i++) printf("%20.15lf",V[i]);
	printf("\n");
getchar();
/*kew*/

  return;
}



/*---------------------------------------------------------
lar_perform_BJ_or_GS_iterations_block - to perform one iteration of block Gauss-Seidel
	algorithm (block Jacobi switched on by it_matrix->Precon==BLOCK_JACOBI)
		v_out = v_in + M^-1 * ( b - A * v_in )
---------------------------------------------------------*/
void lar_perform_BJ_or_GS_iterations_block(
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


/*++++++++++++++++ executable statements ++++++++++++++++*/

// the option for introducing reverse ordering in iterations - not used now
  if(Nr_prec<0) Nr_prec*= -1;

  itt_matrices *it_matrix = &itv_matrices[Matrix_id];

/*kbw
  printf(" lar_BJ_block (matrix_id %d, use_rhs %d, ini_zero %d, nr_prec %d, ndof %d)\n",
	 Matrix_id, Use_rhs, Ini_zero, Nr_prec, Ndof);
  printf("V on entrance (precon %d, block_type %d)\n",
	 it_matrix->Precon, it_matrix->Block_type);
  int i;
  for(i=0;i<40;i++) printf("%20.15lf",V[i]);
  for(i=Ndof-40;i<Ndof;i++) printf("%20.15lf",V[i]);
  //for(i=0;i<Ndof;i++) printf("%10.6lf",V[i]);
  printf("\n");
  if(Use_rhs==1 && B!=NULL){
	printf("B on entrance:\n");
	for(i=0;i<40;i++) printf("%20.15lf",B[i]);
	for(i=Ndof-40;i<Ndof;i++) printf("%20.15lf",B[i]);
	//for(i=0;i<Ndof;i++) printf("%10.6lf",B[i]);
	printf("\n");
  }
  getchar();
/*kew*/

/*kbw
  #define SAVE_LOG_TO_FILE
#ifdef SAVE_LOG_TO_FILE
  char filename[256];
  const char *program_name = getenv("_");


	  sprintf(filename,"Block_Solution.txt",program_name);

  FILE *file = fopen(filename,"a+");
#else
  FILE *file = stdout;
#endif

itt_blocks *blocki,*blockj;
itt_blocks_dia *block_dia;
  int row_ind;
  int ii,jj;
  int col_ind;
  int iiblock,iioffbl,iiaux;
  int posgi,nrdofbli,ndofngbi;
  int diag_flag;
  fprintf(file,"\n\n BLOCK on entry: row, column, value\n");



  for(iiblock=1;iiblock<=it_matrix->Nrblocks;++iiblock){  //loop over block row
	blocki = it_matrix->Block[iiblock];
	posgi=blocki->Posg;
	nrdofbli=blocki->Ndof;
	for(ii=0;ii<blocki->Ndof;++ii){        //loop over row in blocks
	  diag_flag=0;
	  for(iioffbl=0;iioffbl<blocki->Lngb[0];iioffbl++){ //loop over blocks (neighbor)


	  int iiaux = blocki->Lngb[iioffbl+1];
	  blockj=it_matrix->Block[iiaux];
	  ndofngbi = blockj->Ndof;
		int jjaux = blockj->Posg;

	  if(jjaux>posgi && diag_flag==0){      //if diagonal block isnt last
		diag_flag=1;
		for(jj=0;jj<ndofngbi;++jj){        //loop over columns in block
		  fprintf(file,"%6d%6d%20.10lf\n",posgi+ii,posgi+jj,blocki->Dia[ii+jj*ndofngbi]);
		}
	  }

		//offdiagonal blocks
		for(jj=0;jj<ndofngbi;++jj){        //loop over columns in block
		  fprintf(file,"%6d%6d%20.10lf\n",posgi+ii,jjaux+jj,blocki->Aux[iioffbl][ii+jj*nrdofbli]);
		}

	  if(diag_flag==0 && iioffbl+1==blocki->Lngb[0]){     //if diagonal block is last
		diag_flag=1;
		for(jj=0;jj<ndofngbi;++jj){        //loop over columns in block
		  fprintf(file,"%6d%6d%20.10lf\n",posgi+ii,posgi+jj,blocki->Dia[ii+jj*ndofngbi]);
		}
	  }

	  }
	}

  }
  for(ii=0;ii<Ndof;++ii)
	fprintf(file,"%20.10lf",V[ii]);


#ifdef SAVE_LOG_TO_FILE
  fclose(file);
#endif

/*kew*/

 if(it_matrix->Block_type==0){


  int i_loop;		/* counter for loops for GS */
  /* loop over loops */
  for(i_loop=0;i_loop<Nr_prec;i_loop++){

	int iblock;		/* counter for blocks - private as loop index */
	/* loop over all blocks */
	//#pragma omp parallel default(none) firstprivate(V,B,it_matrix,i_loop,Use_rhs,Ini_zero,Ndof)
	{

	  // each thread gets its copy of input vector
	  double *vtemp_GS_BJ;		/* temporary vector*/
	  {
	//vtemp_GS_BJ = lar_util_dvector(Ndof,"vtemp_GS in it_bliter\n");
	vtemp_GS_BJ = (double*)malloc(Ndof*sizeof(double));
	int i;
	for(i=0;i<Ndof;i++) vtemp_GS_BJ[i]=V[i];
	  }

	  //#pragma omp for
	  for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){

	/* constants */
	int ione = 1;
	double done = 1.;
	double dmone = -1.;

	/* local variables */
	itt_blocks *block;  	/* to simplify typing */
	int ioffbl;		/* counter for off diagonal blocks */
	int nrdofbl;		/* number of dofs in a block */
	int ndofngb;		/* number of dofs for a neighbor */

	/* auxiliary variables */
	int i,j;
	int iaux, jaux, kaux;


	/* to simplify typing */
	if(i_loop%2==0) block = it_matrix->Block[iblock];
	else  block = it_matrix->Block[it_matrix->Nrblocks+1-iblock];

	/* number of degrees of freedom for the block */
	nrdofbl = block->Ndof;

	/* position of the first iblock dof in the vector v */
	kaux = block->Posg;

/*kbw
	printf("block %d, thread %d, kaux %d, nrdof %d \n",
		   iblock, omp_get_thread_num(), kaux, nrdofbl);
/*kew*/

	/* initialize part of v */
	for(i=0;i<nrdofbl;i++) vtemp_GS_BJ[kaux+i]=0.;

	if(!Ini_zero||it_matrix->Precon==BLOCK_GS||i_loop>0){
	  /* if initial guess not zero */

	  /* loop over neighbors */
	  for(ioffbl=0;ioffbl<block->Lngb[0];ioffbl++){

		/* global number of the subsequent neighbor */
		iaux = block->Lngb[ioffbl+1];

		/* number of degrees of freedom for the neighbor */
		ndofngb = it_matrix->Block[iaux]->Ndof;

		/* position of the first ingb dof in the global v vector */
		jaux = it_matrix->Block[iaux]->Posg;


		/* solution of the local problem */
		/* if(nrdofbl*ndofngb>MIN_DOF_DGEMV){ */
		/*   if(it_matrix->Precon==BLOCK_JACOBI)  */
		/*     dgemv_("N", &nrdofbl, &ndofngb, &dmone,  */
		/* 	     block->Aux[ioffbl], &nrdofbl, &V[jaux], &ione,  */
		/* 	     &done, &vtemp[kaux], &ione); */
		/*   else dgemv_("N", &nrdofbl, &ndofngb, &dmone,  */
		/* 		block->Aux[ioffbl], &nrdofbl, &V[jaux], &ione,  */
		/* 		&done, &V[kaux], &ione); */
		/* } */
		/* else  */
		{ /* if problem too small to use dgemv */
		  if(it_matrix->Precon==BLOCK_JACOBI){
		for(i=0;i<nrdofbl;i++){
		  for(j=0;j<ndofngb;j++){
			vtemp_GS_BJ[kaux+i] -= block->Aux[ioffbl][i+j*nrdofbl]*V[jaux+j];
/*kbw
			  if(iblock<20){
			  printf("iblock %d, nrdofbl %d, kaux %d, ioffbl %d, iaux %d, jaux %d, i %d, j %d\n",
			  iblock, nrdofbl, kaux, ioffbl, iaux, jaux, i, j);
			  printf("20.15lf%20.15lf%20.15lf\n",
			  vtemp_GS_BJ[kaux+i], block->Aux[ioffbl][i+j*nrdofbl],
			  V[jaux+j]);
			  }
/*kew*/
		  }
		}
		  }
		  else {
		for(i=0;i<nrdofbl;i++){
		  for(j=0;j<ndofngb;j++){
			//V[kaux+i] -= block->Aux[ioffbl][i+j*nrdofbl]*V[jaux+j];
			// operations are performed on local copies
			vtemp_GS_BJ[kaux+i] -= block->Aux[ioffbl][i+j*nrdofbl]*vtemp_GS_BJ[jaux+j];
		  }
		}
		  }
		}
	  }
	}

	/* update using the rhs vector */
	if(Use_rhs==1){
	  if(B==NULL){
		if(it_matrix->Precon==BLOCK_JACOBI) {
		  daxpy_(&nrdofbl, &done, block->Rhs, &ione,
			 &vtemp_GS_BJ[kaux], &ione);
		}
		else {
		  //daxpy_(&nrdofbl, &done, block->Rhs, &ione,
		  //	   &V[kaux], &ione);
		  // operations are performed on local copies
		  daxpy_(&nrdofbl, &done, block->Rhs, &ione,
			 &vtemp_GS_BJ[kaux], &ione);
		}
	  }
	  else{
		if(it_matrix->Precon==BLOCK_JACOBI) {
		  daxpy_(&nrdofbl, &done, &B[kaux], &ione,
			 &vtemp_GS_BJ[kaux], &ione);
		}
		else {
		  //daxpy_(&nrdofbl, &done, &B[kaux], &ione,
		  //	   &V[kaux], &ione);
		  // operations are performed on local copies
		  daxpy_(&nrdofbl, &done, &B[kaux], &ione,
			 &vtemp_GS_BJ[kaux], &ione);
		}
	  }
	}

	  } /* loop over all blocks */

	// still in parallel region
	  //#pragma omp for
	  for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){

	/* constants */
	int ione = 1;
	double done = 1.;
	double dmone = -1.;

	/* local variables */
	itt_blocks *block;  	/* to simplify typing */
	int ioffbl;		/* counter for off diagonal blocks */
	int nrdofbl;		/* number of dofs in a block */
	int ndofngb;		/* number of dofs for a neighbor */

	/* auxiliary variables */
	int i,j;
	int iaux, jaux, kaux;


	/* to simplify typing */
	if(i_loop%2==0) block = it_matrix->Block[iblock];
	else  block = it_matrix->Block[it_matrix->Nrblocks+1-iblock];

	/* number of degrees of freedom for the block */
	nrdofbl = block->Ndof;

	/* position of the first iblock dof in the vector v */
	kaux = block->Posg;

	/* rewrite vtemp to V */
	for(i=0;i<nrdofbl;i++) V[kaux+i] = vtemp_GS_BJ[kaux+i];

/*kbw
	printf("block %d, thread %d, kaux %d, nrdof %d \n",
		   iblock, omp_get_thread_num(), kaux, nrdofbl);
/*kew*/


	  } // end loop over blocks


	  free(vtemp_GS_BJ);

	} // end parallel region


  } /* loop over preconditioner iterations */


 } // end if block_type == 0
 else{

  int i_loop;		/* counter for loops */
  for(i_loop=0;i_loop<Nr_prec;i_loop++){

	/* loop over all big diagonal blocks */
	int ib_dia;	/* counters for blocks */
	//#pragma omp parallel if(abs(it_matrix->Block_type)==1) default(none) shared(V,B,it_matrix)   firstprivate(i_loop,Use_rhs,Ini_zero,Ndof)
	{

	  // each thread gets its copy of input vector
	  double *vtemp_GS_BJ;		/* temporary vector*/
	  {
	vtemp_GS_BJ = lar_util_dvector(Ndof,"vtemp_GS_BJ in it_bliter\n");
	int i;
	for(i=0;i<Ndof;i++) vtemp_GS_BJ[i]=V[i];
	  }


	  //#pragma omp for
	  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

	/* constants */
	int ione = 1;
	double done = 1.;
	double dmone = -1.;

	/* local variables */
	itt_blocks *block;  	/* to simplify typing */
	itt_blocks_dia *block_dia;/* to simplify typing */
	int iblock;	/* counters for blocks */
	int ioffbl;		/* counter for off diagonal blocks */
	int nrdofbl, nrdof_dia;	/* numbers of dofs in a block */
	int ndofngb;		/* number of dofs for a neighbor */
	double vloc[Max_dof_block_dia];	/* vector for DIA dofs */

	/* auxiliary variables */
	int i,j;
	int iaux, jaux, kaux, laux, iblaux;

	if(i_loop%2==0) iblaux = ib_dia;
	else iblaux = it_matrix->Nrblocks_dia + 1 - ib_dia;

/*kbw
#pragma omp critical(printing)
	  {
	printf("ib_dia %d, iblaux %d, it_matrix->Precon %d, i_loop %d, Ini_zero %d, Use_rhs %d, Ndof %d\n",
		   ib_dia, iblaux, it_matrix->Precon, i_loop, Ini_zero, Use_rhs, Ndof);


	  }
/*kew*/

	//#pragma omp ordered
	if(it_matrix->Block_dia[iblaux]!=NULL){

	  /* to simplify typing */
	  block_dia=it_matrix->Block_dia[iblaux];

	  /* total number of dofs for DIA */
	  nrdof_dia = block_dia->Lpos[0];
	  if(nrdof_dia>Max_dof_block_dia){
		printf("Error - too much dofs in dia block, \n");
		printf("Increase Max_dof_block_dia in itb_bliter (it_matrix->Block.c)\n");
		exit(1);
	  }

	  /* initialize vloc */
	  for(i=0;i<nrdof_dia;i++) vloc[i]=0.;

	  /* first perform multiplication A*x */

	  /* loop over small blocks */
	  for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

		/* to simplify typing */
		block = it_matrix->Block[block_dia->Lsmall[iblock]];

		/* number of degrees of freedom for the block */
		nrdofbl = block->Ndof;

		/* position of the first iblock dof in the vector v */
		kaux = block->Posg;

		/* position in vloc */
		laux = block_dia->Lpos[iblock];

		if(!Ini_zero||it_matrix->Precon==BLOCK_GS||i_loop>0){
		  /* if initial guess not zero */

		  /* loop over neighbors */
		  for(ioffbl=1;ioffbl<=block->Lngb[0];ioffbl++){

		/* global number of the subsequent neighbor */
		iaux = block->Lngb[ioffbl];

		/* if not found on the list of DIA blocks */
		if(!lar_util_chk_list(iaux,&block_dia->Lsmall[1],
				block_dia->Lsmall[0])){

		  /* number of degrees of freedom for the neighbor */
		  ndofngb = it_matrix->Block[iaux]->Ndof;

		  /* position of the first ingb dof in the global v vector */
		  jaux = it_matrix->Block[iaux]->Posg;

		  /* multiplication by off-diagonal subarray */
		  /* if(nrdofbl*ndofngb>MIN_DOF_DGEMV){ */
		  /* if(nrdofbl*ndofngb>MIN_DOF_DGEMV){ */
		  /*   dgemv_("N", &nrdofbl, &ndofngb, &dmone,  */
		  /* 	 block->Aux[ioffbl-1], &nrdofbl, &V[jaux], &ione,  */
		  /* 	 &done, &vloc[laux], &ione); */
		  /* } */
		  /* else  */
//#pragma omp ordered
//#pragma omp critical(V_update)
		  { // if problem too small to use dgemv
			for(i=0;i<nrdofbl;i++){
			  for(j=0;j<ndofngb;j++){
			if(it_matrix->Precon==BLOCK_JACOBI) {
			  vloc[laux+i] -=block->Aux[ioffbl-1][i+j*nrdofbl]*V[jaux+j];
/*kbw
			  //if(iblaux<20)
			{
#pragma omp critical(printing)
			  {
			  printf("%6d%6d%6d%20.10lf%20.10lf%20.10lf\n",
				 iblaux, block_dia->Lsmall[iblock], ioffbl,
				 vloc[laux+i], block->Aux[ioffbl-1][i+j*nrdofbl],
				 V[jaux+j]);
			  }
			  }
/*kew*/
			}
			else{
			  // operations are performed on local copies of V stored in vtemp
			  vloc[laux+i] -=block->Aux[ioffbl-1][i+j*nrdofbl]*vtemp_GS_BJ[jaux+j];
/*kbw
			  //if(iblaux<20)
			{
#pragma omp critical(printing)
			  {
			  printf("%6d%6d%6d%20.10lf%20.10lf%20.10lf\n",
				 iblaux, block_dia->Lsmall[iblock], ioffbl,
				 vloc[laux+i], block->Aux[ioffbl-1][i+j*nrdofbl],
				 vtemp_GS_BJ[jaux+j]);
			  }
			  }
/*kew*/

			}
			  }
			}
		  }
/*kbw
	getchar();
/*kew*/

		} /* if not found on the list of DIA blocks */
		  } /* loop over neighbors */
		} /* if initial guess not zero */

/*kbw
	  if(iblaux<20)
//#pragma omp ordered
		{
#pragma omp critical(printing)
		  {
		  printf("V (vloc) after multiplication, i_loop %d, ib_dia %d, iblaux  %d\n",
		i_loop, ib_dia, iblaux);
	for(i=0;i<nrdof_dia;i++) printf("%20.10lf",vloc[i]);
	printf("\n");
	//getchar();
		  }
	  }
/*kew*/

/* update using the rhs vector */
		if(Use_rhs==1){
		  if(B==NULL){
		daxpy_(&nrdofbl, &done, block->Rhs, &ione,
			   &vloc[laux], &ione);
/*kbw
//	  if(iblaux<20)
{
#pragma omp critical(printing)
  {
	printf("RHS for subtraction\n");
	for(i=0;i<nrdofbl;i++) printf("%20.10lf", block->Rhs[i]);
	printf("\n");
	//getchar();
  }
	  }
/*kew*/

		  }
		  else{
		daxpy_(&nrdofbl, &done, &B[kaux], &ione,
			   &vloc[laux], &ione);
		  }


/*kbw
		  if(iblaux<20||iblaux>it_matrix->Nrblocks_dia-20)
{
#pragma omp critical(printing)
  {
	printf("V (vloc) after RHS subtraction, i_loop %d, ib_dia %d, iblaux  %d, Nrblocks_dia %d\n",
	   i_loop, ib_dia, iblaux, it_matrix->Nrblocks_dia);
	for(i=0;i<nrdof_dia;i++) printf("%10.6lf",vloc[i]);
	printf("\n");
	getchar();
  }
	  }
/*kew*/

		} /* if updated using the rhs vector */

	  } /* loop over small blocks */


	  /* resolve subdomain problem - solve for diagonal submatrix */
	  /*   v = block->Dia^(-1) * vloc */

	  if(nrdof_dia==1){
		vloc[0]*=block_dia->Dia[0];
	  }
	  else{
		dgetrs_("N",&nrdof_dia,&ione,block_dia->Dia,&nrdof_dia,
			block_dia->Ips,vloc,&nrdof_dia,&iaux);
	  }

	  /* rewrite back to global vector */
	  /* loop over small blocks */
	  for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

		/* to simplify typing */
		block = it_matrix->Block[block_dia->Lsmall[iblock]];

		/* number of degrees of freedom for the block */
		nrdofbl = block->Ndof;

		/* position of the first iblock dof in the vector v */
		kaux = block->Posg;

		/* position in vloc */
		laux = block_dia->Lpos[iblock];

		// updates are performed on local copies of V stored in vtemp
		for(i=0;i<nrdofbl;i++) vtemp_GS_BJ[kaux+i] = vloc[laux+i];
/*kbw
//	  if(iblaux<20)
{
#pragma omp critical(printing)
  {
  printf("V (vloc) after rewriting small block %d (posglob %d)\n", iblock, kaux);
	for(i=0;i<nrdof_dia;i++) printf("%20.10lf",vloc[i]);
	printf("\n");
  }
	//getchar();
	  }
/*kew*/

	  } /* loop over small blocks */

	} /* if block_dia is active */
	  } /* loop over big diagonal blocks */



	// still in parallel region
	  //#pragma omp for
	for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){


/*kbw
	  printf("block %d, thread %d\n", ib_dia, omp_get_thread_num());
/*kew*/

/* constants */
	  int ione = 1;
	  double done = 1.;
	  double dmone = -1.;

	  /* local variables */
	  itt_blocks *block;  	/* to simplify typing */
	  itt_blocks_dia *block_dia;/* to simplify typing */
	  int iblock;	/* counters for blocks */
	  int ioffbl;		/* counter for off diagonal blocks */
	  int nrdofbl, nrdof_dia;	/* numbers of dofs in a block */
	  int ndofngb;		/* number of dofs for a neighbor */

	  /* auxiliary variables */
	  int i,j;
	  int iaux, jaux, kaux, laux, iblaux;

	  if(i_loop%2==0) iblaux = ib_dia;
	  else iblaux = it_matrix->Nrblocks_dia + 1 - ib_dia;

/*kbw
#pragma omp critical(printing)
	  {
	printf("ib_dia %d, iblaux %d, it_matrix->Precon %d, i_loop %d, Ini_zero %d, Use_rhs %d, Ndof %d\n",
		   ib_dia, iblaux, it_matrix->Precon, i_loop, Ini_zero, Use_rhs, Ndof);


	  }
/*kew*/

	  //#pragma omp ordered
	  if(it_matrix->Block_dia[iblaux]!=NULL){

		/* to simplify typing */
		block_dia=it_matrix->Block_dia[iblaux];

		/* total number of dofs for DIA */
		nrdof_dia = block_dia->Lpos[0];
		if(nrdof_dia>Max_dof_block_dia){
		  printf("Error - too much dofs in dia block, \n");
		  printf("Increase Max_dof_block_dia in itb_bliter (it_matrix->Block.c)\n");
		  exit(1);
		}

		/* loop over small blocks */
		for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

		  /* to simplify typing */
		  block = it_matrix->Block[block_dia->Lsmall[iblock]];

		  /* number of degrees of freedom for the block */
		  nrdofbl = block->Ndof;

		  /* position of the first iblock dof in the vector v */
		  kaux = block->Posg;

		  // critical not necessary since multithreading is switched on only for BLOCK_SIZE=1
		  //#pragma omp critical(V_update)
		  //{
		//for(i=0;i<nrdofbl;i++) V[kaux+i] = 0.5*(V[kaux+i]+vtemp_GS_BJ[kaux+i]);
		  //}

		  for(i=0;i<nrdofbl;i++) V[kaux+i] = vtemp_GS_BJ[kaux+i];

/*kbw
	  printf("block %d, thread %d, kaux %d, nrdof %d \n",
		 ib_dia, omp_get_thread_num(), kaux, nrdofbl);
/*kew*/


		}

	  }

	} // end second loop over big blocks (GS only to rewrite from local thread copies of V)

	free(vtemp_GS_BJ);


	} // end parallel region

  } /* loop over preconditioner iterations */

 }

/*kbw
{
	printf("V on leaving\n");
	for(i=0;i<40;i++) printf("%20.15lf",V[i]);
	 for(i=Ndof-40;i<Ndof;i++) printf("%20.15lf",V[i]);
   //for(i=0;i<Ndof;i++) printf("%10.6lf",V[i]);
	printf("\n");
   getchar();
}
/*kew*/

 return;
}

/*---------------------------------------------------------
lar_perform_rhsub_block - to perform forward reduction and back-substitution for ILU
		   preconditioning: v_out = M^-1 * b
---------------------------------------------------------*/
void lar_perform_rhsub_block(
  int Matrix_id,   /* in: matrix ID */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* out: vector of unknowns updated */
		/*      during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	)
{

/* constants */
  int ione = 1;
  double done = 1.;
  double dmone = -1.;


/* local variables */
  itt_blocks *block, *blneig;  	/* to simplify typing */
  itt_blocks_dia *block_dia, *blneig_dia;/* to simplify typing */
  int iblock, ioffbl, ib_dia, jb_dia, ibl_sm, jbl_sm; /* counters for blocks */
  int nrdof1, nrdof2, nrdof_dia;	/* numbers of dofs in a block */
  int posloc1, posloc2, posglob1, posglob2; /* positions in dofs vectors */
  double vloc[Max_dof_block_dia];	/* vector for DIA dofs */
  double vloc1[Max_dof_block_dia];	/* vector for DIA dofs */

/* auxiliary variables */
  int i,j, jneig, iaux, ioff;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  itt_matrices *it_matrix = &itv_matrices[Matrix_id];

/* FORWARD REDUCTION */

/* loop over all big diagonal blocks */
  for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

/* to simplify typing */
	block_dia=it_matrix->Block_dia[ib_dia];

/* total number of dofs for DIA */
	nrdof_dia = block_dia->Lpos[0];
	if(nrdof_dia>Max_dof_block_dia){
	  printf("Error - too much dofs in dia block, \n");
	  printf("Increase Max_dof_block_dia in itb_bliter (it_matrix->Block.c)\n");
	  exit(1);
	}

/* initialize vloc */
	for(i=0;i<nrdof_dia;i++) vloc[i]=0.;

/* substitute value of rhs to vloc */

/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

	  ibl_sm = block_dia->Lsmall[iblock];

/* to simplify typing */
	  block = it_matrix->Block[ibl_sm];

/* number of degrees of freedom for the block */
	  nrdof1 = block->Ndof;

/* position of the first iblock dof in the vector v */
	  posglob1 = block->Posg;

/* position in vloc */
	  posloc1 = block_dia->Lpos[iblock];

/* update using the rhs vector */
	  if(B==NULL){
	daxpy_(&nrdof1, &done, block->Rhs, &ione,
		   &vloc[posloc1], &ione);
	  }
	  else{
	daxpy_(&nrdof1, &done, &B[posglob1], &ione,
		   &vloc[posloc1], &ione);
	  }

/* for all neighboring factorized (lower subdiagonal) entries of A */
	  for(jneig=1;jneig<=block_dia->Llowerneig[0];jneig++){

	jb_dia = block_dia->Llowerneig[2*jneig-1];
	ioff = block_dia->Llowerneig[2*jneig];

#ifdef DEBUG
	  if(block_dia->Lneig[ioff]!=jb_dia){
		printf("error in reading list of lower neighbors in rhsub!\n");
		printf("ib_dia %d, iblock %d, jneig %d, jb_dia %d, ioff %d, ? %d\n",
		   ib_dia, iblock, jneig, jb_dia, ioff, block_dia->Lneig[ioff]);
		exit(1);
	  }
#endif

	  blneig_dia = it_matrix->Block_dia[jb_dia];


/* loop over neighbor's small blocks */
	  for(ioffbl=1;ioffbl<=blneig_dia->Lsmall[0];ioffbl++){

/* global number of the subsequent neighbor */
		jbl_sm = blneig_dia->Lsmall[ioffbl];
		blneig = it_matrix->Block[jbl_sm];

/* number of degrees of freedom for the neighbor */
		  nrdof2 = blneig->Ndof;

/* position of the first ingb dof in the global v vector */
		  posglob2 = blneig->Posg;

/* position of the first ingb dof in the Aux array  */
		  posloc2 = blneig_dia->Lpos[ioffbl];

/* multiplication by off-diagonal subarray
		if(nrdofbl*ndofngb>MIN_DOF_DGEMV){
		  dgemv_("N", &nrdofbl, &ndofngb, &dmone,
			 block->Aux[ioffbl-1], &nrdofbl, &V[jaux], &ione,
			 &done, &vloc[laux], &ione);
		}
		else {
*/
		  for(i=0;i<nrdof1;i++){
		for(j=0;j<nrdof2;j++){
		  iaux = posloc1+i+(posloc2+j)*nrdof_dia;
		  vloc[posloc1+i] -=
			block_dia->Aux[ioff-1][iaux]*V[posglob2+j];
		}
		  }

/*		} */

/*kbw
		  printf("V before updating with block %d (%d)\n", jb_dia, ioff);
		  for(i=0;i<nrdof2;i++) printf("%25.15lf",V[posglob2+i]);
		  printf("\n");
		  printf("ib_dia %d, jb_dia %d, Aux before solution:\n",ib_dia,jb_dia);
		  for(i=0;i<blneig_dia->Lpos[0]*nrdof_dia;i++){
		printf("%25.15lf",block_dia->Aux[ioff-1][i]);
		  }
		  printf("\n");
		  printf("\nblock_dia %d small block %d after updating with small block %d\n",
			 ib_dia,ibl_sm,jbl_sm);
		  for(i=0;i<nrdof1;i++) printf("%25.15lf",vloc[i]);
		  printf("\n");
		  //getchar();
/*kew*/

	  } /* loop over small blocks of a neighboring big blocks */
	  } /* loop over neighboring big blocks */
	} /* loop over small blocks */


/*kbw
  printf("kb_dia %d, Dia before solution:\n",kb_dia);
  for(i=0;i<kaux*kaux;i++){
	printf("%20.15lf",block_k->Dia[i]);
  }
  printf("\n");
  printf("ib_dia %d, ioff %d, Aux before solution:\n",ib_dia,ioff);
  for(i=0;i<kaux*nrdof_dia;i++){
	printf("%20.15lf",block_dia->Aux[ioff-1][i]);
  }
  printf("\n");
/*kew*/



/* rewrite back to global vector */
/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

/* to simplify typing */
	  block = it_matrix->Block[block_dia->Lsmall[iblock]];

/* number of degrees of freedom for the block */
	  nrdof1 = block->Ndof;

/* position of the first iblock dof in the vector v */
	  posglob1 = block->Posg;

/* position in vloc */
	  posloc1 = block_dia->Lpos[iblock];

	  for(i=0;i<nrdof1;i++) V[posglob1+i] = vloc[posloc1+i];

/*kbw
	  printf("part of V %d-%d after updating with all small blocks\n",
		 posglob1, posglob1+nrdof1-1 );
	  for(i=0;i<nrdof1;i++) printf("%25.15lf", V[posglob1+i]);
	  //getchar();
/*kew*/

	} /* loop over small blocks */

  } /* loop over big diagonal blocks */

/*kbw
  printf("V after FORWARD REDUCTION:\n");
  for(i=0;i<Ndof;i++) printf("%25.15lf",V[i]);
  printf("\n");
  getchar();getchar();
/*kew*/


/* BACK SUBSTITUTION */
/* loop over all big diagonal blocks */
  for(ib_dia=it_matrix->Nrblocks_dia;ib_dia>=1;ib_dia--){

/* to simplify typing */
	block_dia=it_matrix->Block_dia[ib_dia];

/* total number of dofs for DIA */
	nrdof_dia = block_dia->Lpos[0];


/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

	  ibl_sm = block_dia->Lsmall[iblock];

/* to simplify typing */
	  block = it_matrix->Block[ibl_sm];

/* number of degrees of freedom for the block */
	  nrdof1 = block->Ndof;

/* position of the first iblock dof in the vector v */
	  posglob1 = block->Posg;

/* position in vloc */
	  posloc1 = block_dia->Lpos[iblock];

/* initialize vloc */
	  for(i=0;i<nrdof1;i++) vloc[posloc1+i] = V[posglob1+i] ;

/*kbw
	printf("\nvloc initialized: block_dia %d small block %d \n",
	   ib_dia,ibl_sm);
	for(i=0;i<nrdof_dia;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	//getchar();
/*kew*/


/* for all neighboring factorized (upper subdiagonal) entries of A */
	  for(jneig=1;jneig<=block_dia->Lupperneig[0];jneig++){

	jb_dia = block_dia->Lupperneig[2*jneig-1];
	ioff = block_dia->Lupperneig[2*jneig];

#ifdef DEBUG
	if(block_dia->Lneig[ioff]!=jb_dia){
	  printf("error in reading list of upper neighbors in rhsub!\n");
	  printf("ib_dia %d, iblock %d, jneig %d, jb_dia %d, ioff %d, ? %d\n",
		 ib_dia, iblock, jneig, jb_dia, ioff, block_dia->Lneig[ioff]);
	  exit(1);
	}
#endif

	  blneig_dia = it_matrix->Block_dia[jb_dia];

/* loop over neighbor's small blocks */
	  for(ioffbl=1;ioffbl<=blneig_dia->Lsmall[0];ioffbl++){

/* global number of the subsequent neighbor */
		jbl_sm = blneig_dia->Lsmall[ioffbl];
		blneig = it_matrix->Block[jbl_sm];

/* number of degrees of freedom for the neighbor */
		  nrdof2 = blneig->Ndof;

/* position of the first ingb dof in the global v vector */
		  posglob2 = blneig->Posg;

/* position of the first ingb dof in the Aux array  */
		  posloc2 = blneig_dia->Lpos[ioffbl];

/* multiplication by off-diagonal subarray
		if(nrdofbl*ndofngb>MIN_DOF_DGEMV){
		  dgemv_("N", &nrdofbl, &ndofngb, &dmone,
			 block->Aux[ioffbl-1], &nrdofbl, &V[jaux], &ione,
			 &done, &vloc[laux], &ione);
		}
		else {
*/
		  for(i=0;i<nrdof1;i++){
		for(j=0;j<nrdof2;j++){
		  iaux = posloc1+i+(posloc2+j)*nrdof_dia;
		  vloc[posloc1+i] -=
			block_dia->Aux[ioff-1][iaux]*V[posglob2+j];
		}
		  }

/*		} */

/*kbw
		  printf("V before updating with block %d (%d)\n", jb_dia, ioff);
		  for(i=0;i<nrdof2;i++) printf("%25.15lf",V[posglob2+i]);
		  printf("\n");
		  printf("ib_dia %d, jb_dia %d, Aux after solution:\n",ib_dia,jb_dia);
		  for(i=0;i<blneig_dia->Lpos[0]*nrdof_dia;i++){
		printf("%25.15lf",block_dia->Aux[ioff-1][i]);
		  }
		  printf("\n");
		  printf("\nblock_dia %d small block %d after updating with small block %d\n",
			 ib_dia,ibl_sm,jbl_sm);
		  for(i=0;i<nrdof_dia;i++) printf("%25.15lf",vloc[i]);
		  printf("\n");
		  //getchar();
/*kew*/

	  } /* loop over small blocks of a neighboring big blocks */
	  } /* loop over neighboring big blocks */
	} /* loop over small blocks */

/* resolve subdomain problem - solve for diagonal submatrix */
/*   v = block->Dia^(-1) * vloc */

	if(nrdof_dia==1){
	  vloc1[0]=block_dia->Dia[0]*vloc[0];
	}
	else{
		  for(i=0;i<nrdof_dia;i++){
		  vloc1[i]=0.0;
		for(j=0;j<nrdof_dia;j++){
			vloc1[i] +=
			  block_dia->Dia[i+nrdof_dia*j]*vloc[j];
		  }
		  }

	}

/*kbw
	printf("\nblock_dia %d small block %d before solving with Dia\n",
	   ib_dia,ibl_sm);
	for(i=0;i<nrdof_dia;i++) printf("%25.15lf",vloc[i]);
	printf("\n");
	printf("ib_dia %d, Dia before solution:\n",ib_dia);
	for(i=0;i<nrdof_dia*nrdof_dia;i++){
	  printf("%25.15lf",block_dia->Dia[i]);
	}
	printf("\n");
	printf("\nblock_dia %d small block %d after solving with Dia\n",
	   ib_dia,ibl_sm);
	for(i=0;i<nrdof_dia;i++) printf("%25.15lf",vloc1[i]);
	//getchar();
/*kew*/


/* rewrite back to global vector */
/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

/* to simplify typing */
	  block = it_matrix->Block[block_dia->Lsmall[iblock]];

/* number of degrees of freedom for the block */
	  nrdof1 = block->Ndof;

/* position of the first iblock dof in the vector v */
	  posglob1 = block->Posg;

/* position in vloc */
	  posloc1 = block_dia->Lpos[iblock];

	  for(i=0;i<nrdof1;i++) V[posglob1+i] = vloc1[posloc1+i];

	} /* loop over small blocks */

  } /* loop over big diagonal blocks */

/*kbw
  printf("V after BACKWARD SUBSTITUTION:\n");
  for(i=0;i<Ndof;i++) printf("%25.15lf",V[i]);
  printf("\n");
  getchar();getchar();
/*kew*/

  return;
}






// kept for sentimental reasons and bleak outlooks for the future - matrix free approach

/*---------------------------------------------------------
it_mfaiter - to perform matrix vector multiplication (possibly
	in matrix-free manner) and additive Schwarz approximate solve
---------------------------------------------------------*/
void it_mfaiter(
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

/* constants */
int ione = 1;
double done = 1.;
double dmone = -1.;

/* local variables */
itt_blocks *block;  	/* to simplify typing */
itt_blocks_dia *block_dia;/* to simplify typing */
int iblock, ib_dia;	/* counters for blocks */
int ioffbl;		/* counter for off diagonal blocks */
int nrdofbl, nrdof_dia;	/* numbers of dofs in a block */
int ndofngb;		/* number of dofs for a neighbor */
int i_loop;		/* counter for loops for GS */
double vloc[Max_dof_block_dia];	/* vector for DIA dofs */
double inveps;	/* inverse of small parameter for finite differences */

/* auxiliary variables */
int i,j;
int iaux, jaux, kaux, laux, iblaux, iel;
double daux, *vtemp, *vtemp1;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kb!!!!!! inveps = 1.0/it_matrix->Jaceps; */

itt_matrices *it_matrix = &itv_matrices[Matrix_id];

vtemp = lar_util_dvector(Ndof,"vtemp in itb_dditer\n");
vtemp1 = lar_util_dvector(Ndof,"vtemp in itb_dditer\n");

/* loop over loops */
for(i_loop=0;i_loop<Nr_prec;i_loop++){

/* first perform matrix vector multiplication using elementary blocks */
/* make global pointer it_matrix->u_MF to point to current vector v
it_matrix->u_MF=V; */

/* for each element active in a given mesh*/
  //for(iel=1;iel<=it_matrix->Mesh_elems[0];iel++){

/* compute element load vector for matrix free multiplication */
/* and assemble to block arrays Aux */
/*  pbr_assemble_matrix_free(iel); */

  //}

/* loop over all blocks */
for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){

/* to simplify typing */
  if(i_loop%2==0) block = it_matrix->Block[iblock];
  else if(i_loop%2==1) block = it_matrix->Block[it_matrix->Nrblocks+1-iblock];

/* number of degrees of freedom for the block */
  nrdofbl = block->Ndof;

/* position of the first iblock dof in the vector v */
  kaux = block->Posg;

/* product Ax in a matrix free manner */
  if(B==NULL){
	for(i=0;i<nrdofbl;i++){
	  vtemp[kaux+i] = inveps*(block->Aux[0][i]-block->Rhs[i]);
	}
  }
/* for matrix free formulation B may store modified RHS */
  else{
	for(i=0;i<nrdofbl;i++){
	  vtemp[kaux+i] = inveps*(block->Aux[0][i]-block->Rhs[i]);
	}
  }

/* update using the rhs vector */
  if(Use_rhs==1){
	if(B==NULL){
	  daxpy_(&nrdofbl, &done, block->Rhs, &ione,
		&vtemp[kaux], &ione);
	}
/* for matrix free formulation B may store modified RHS */
	else{
	  daxpy_(&nrdofbl, &done, block->Rhs, &ione,
		&vtemp[kaux], &ione);
	}
  }

} /* loop over all blocks */

/* clear Aux arrays for new iterations */
for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){
  for(iaux=0;iaux<it_matrix->Block[iblock]->Ndof;iaux++){
	it_matrix->Block[iblock]->Aux[0][iaux]=0;
  }
}

/* create vector for initial state before solving subdomain problems */
dcopy_(&Ndof, vtemp, &ione, vtemp1, &ione);

/* perform preconditioning step using big blocks */

/* loop over all big diagonal blocks */
for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

  if(i_loop%2==0) iblaux = ib_dia;
  else if(i_loop%2==1) iblaux = it_matrix->Nrblocks_dia + 1 - ib_dia;

  if(it_matrix->Block_dia[iblaux]!=NULL){

/* to simplify typing */
	block_dia=it_matrix->Block_dia[iblaux];

/* total number of dofs for DIA */
	nrdof_dia = block_dia->Lpos[0];
	if(nrdof_dia>Max_dof_block_dia){
	  printf("Error - too much dofs in dia block, \n");
	  printf("Increase Max_dof_block_dia in itb_bliter (it_matrix->Block.c)\n");
	  exit(1);
	}

/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

/* to simplify typing */
	  block = it_matrix->Block[block_dia->Lsmall[iblock]];

/* number of degrees of freedom for the block */
	  nrdofbl = block->Ndof;

/* position of the first iblock dof in the vector v */
	  kaux = block->Posg;

/* position in vloc */
	  laux = block_dia->Lpos[iblock];

/* rewrite vtemp to vloc */
	  for(i=0;i<nrdofbl;i++) vloc[laux+i] = vtemp1[kaux+i];

	} /* loop over small blocks */

/* resolve subdomain problem - solve for diagonal submatrix */
/*   vloc = block->Dia^(-1) * vloc */

	if(nrdof_dia==1){
	  vloc[0]*=block_dia->Dia[0];
	}
	else{
	  dgetrs_("N",&nrdof_dia,&ione,block_dia->Dia,&nrdof_dia,
		  block_dia->Ips,vloc,&nrdof_dia,&iaux);
	}

/* rewrite back to global vector */
/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

/* to simplify typing */
	  block = it_matrix->Block[block_dia->Lsmall[iblock]];

/* number of degrees of freedom for the block */
	  nrdofbl = block->Ndof;

/* position of the first iblock dof in the vector v */
	  kaux = block->Posg;

/* position in vloc */
	  laux = block_dia->Lpos[iblock];

/* we use the newest values in vloc to update vtemp */
	  for(i=0;i<nrdofbl;i++) vtemp[kaux+i] = vloc[laux+i];

	} /* loop over small blocks */

  } /* if block_dia is active */
} /* loop over big diagonal blocks */

/* add v to get -A_ii^-1 * SUM_j ( A_ij * x_j - b_i ) (j!=i) */
daux = 1;
daxpy_(&Ndof, &daux, vtemp, &ione, V, &ione);

} /* loop over standard iterations in preconditioning */

free(vtemp);
free(vtemp1);

return;
}

/*---------------------------------------------------------
it_mfmiter - to perform one iteration of block Gauss-Seidel
	algorithm (block Jacobi switched on by it_matrix->Precon==BLOCK_JACOBI)
	within Newton-Krylov-Schwarz algorithm
---------------------------------------------------------*/
void it_mfmiter(
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

/* constants */
int ione = 1;
double done = 1.;

/* local variables */
itt_blocks *block;  	/* to simplify typing */
itt_blocks_dia *block_dia;/* to simplify typing */
int iblock, ib_dia;	/* counters for blocks */
int nrdofbl, nrdof_dia;	/* numbers of dofs in a block */
int i_loop;		/* counter for loops for GS */
double vloc[Max_dof_block_dia];	/* vector for DIA dofs */
double inveps;	/* inverse of small parameter for finite differences */

/* auxiliary variables */
int i,iel;
int iaux, kaux, laux, iblaux;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* inveps = 1.0/it_Jaceps; */

itt_matrices *it_matrix = &itv_matrices[Matrix_id];

/* make global pointer it_u_MF to point to current vector v
it_u_MF=V; */

for(i_loop=0;i_loop<Nr_prec;i_loop++){

/* loop over all big diagonal blocks */
for(ib_dia=1;ib_dia<=it_matrix->Nrblocks_dia;ib_dia++){

  if(i_loop%2==0) iblaux = ib_dia;
  else if(i_loop%2==1) iblaux = it_matrix->Nrblocks_dia + 1 - ib_dia;

  if(it_matrix->Block_dia[iblaux]!=NULL){

/* to simplify typing */
	block_dia=it_matrix->Block_dia[iblaux];

/* total number of dofs for DIA */
	nrdof_dia = block_dia->Lpos[0];
	if(nrdof_dia>Max_dof_block_dia){
	  printf("Error - too much dofs in dia block, \n");
	  printf("Increase Max_dof_block_dia in itb_bliter (it_matrix->Block.c)\n");
	  exit(1);
	}

/* first perform multiplication A*x */

/* loop over all subdomain elements */
	for(iel=1;iel<=block_dia->Lelem[0];iel++){
/*      pbr_assemble_matrix_free(block_dia->Lelem[iel]); */
	}


/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

/* to simplify typing */
	  block = it_matrix->Block[block_dia->Lsmall[iblock]];

/* number of degrees of freedom for the block */
	  nrdofbl = block->Ndof;

/* position of the first iblock dof in the vector v */
	  kaux = block->Posg;

/* position in vloc */
	  laux = block_dia->Lpos[iblock];

/* product Ax in a matrix free manner */
	  if(B==NULL){
		for(i=0;i<nrdofbl;i++){
		  vloc[laux+i] = inveps*(block->Aux[0][i]-block->Rhs[i]);
		}
	  }
	  else{
/* for matrix free formulation B may store modified RHS */
		for(i=0;i<nrdofbl;i++){
		  vloc[laux+i] = inveps*(block->Aux[0][i]-block->Rhs[i]);
		}
	  }

/* update using the rhs vector */
	  if(Use_rhs==1){
		if(B==NULL){
		  daxpy_(&nrdofbl, &done, block->Rhs, &ione,
			&vloc[laux], &ione);
		}
		else{
/* for matrix free formulation B may store modified RHS */
		  daxpy_(&nrdofbl, &done, block->Rhs, &ione,
			&vloc[laux], &ione);
		}

	  } /* if updated using the rhs vector */

	} /* loop over small blocks */

/* resolve subdomain problem - solve for diagonal submatrix */
/*   v = block->Dia^(-1) * vloc */

	if(nrdof_dia==1){
	  vloc[0]*=block_dia->Dia[0];
	}
	else{
	  dgetrs_("N",&nrdof_dia,&ione,block_dia->Dia,&nrdof_dia,
		  block_dia->Ips,vloc,&nrdof_dia,&iaux);
	}

/* rewrite back to global vector */
/* loop over small blocks */
	for(iblock=1;iblock<=block_dia->Lsmall[0];iblock++){

/* to simplify typing */
	  block = it_matrix->Block[block_dia->Lsmall[iblock]];

/* number of degrees of freedom for the block */
	  nrdofbl = block->Ndof;

/* position of the first iblock dof in the vector v */
	  kaux = block->Posg;

/* position in vloc */
	  laux = block_dia->Lpos[iblock];

	  for(i=0;i<nrdofbl;i++) V[kaux+i] += vloc[laux+i];

	} /* loop over small blocks */

/* clean Aux arrays */
	for(iel=1;iel<=block_dia->Lelem[0];iel++){
/*      pbr_assemble_matrix_free(-block_dia->Lelem[iel]); */
	}

  } /* if block_dia is active */
} /* loop over big diagonal blocks */

} /* loop over preconditioner iterations */

return;
}

//////////////////////////////////////////////////////
int lar_block_print_matrix(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id          /* in: matrix ID */


){

  int nr_levels, ineg, nrdofgl, nrdofbl, max_nrdof;
  itt_blocks* block;         /* pointers to single blocks */
  itt_blocks_dia *block_dia; /* to simplify typing */
  //double *v_sol;
  int i,j,k, iaux, ient, iblock;
  double *row;


/*++++++++++++++++ executable statements ++++++++++++++++*/

  itt_matrices *it_matrix = &itv_matrices[Matrix_id];

  row=(double *)malloc(it_matrix->Nrblocks*sizeof(double));
  for(iblock=1;iblock<=it_matrix->Nrblocks;iblock++){
	block = it_matrix->Block[iblock];
	if(block->Lngb!=NULL){

		for(j=0;j<it_matrix->Nrblocks;j++)row[j]=0.0;


		for(j=1;j<=block->Lngb[0];j++){
		  ineg=block->Lngb[j]-1;
		  row[ineg] = block->Aux[j-1][0];

		}
		row[iblock-1]=block->Dia[0];
		for(j=0;j<it_matrix->Nrblocks;j++)
			printf("%15lf ",row[j]);
		printf("\n");
	 }

  }

return(0);

}

