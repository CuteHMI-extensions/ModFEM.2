/*******superlu_ext******************************************************
File contains procedures:
  lsr_mkb_direct_init - to create a new solver instance, read its control
			 parameters and initialize its data structure
  lsr_mkb_direct_solve - to solve a system of equations, given previously constructed
			 system matrix in crs

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

// lsd_mkb interface - including interface for direct solvers
#include <modfem/ls_mkb/lsh_mkb_intf.h>

/* internal information for the solver module */
#include "./lsh_mkb_superlu.h"

// API for linear algebra routines supporting MKB solver
#include <modfem/ls_mkb/lah_intf.h>
// SUPERLU interface needs access to internal data structures of SM
//#include "../lad_crs_generic/lah_crs_generic.h"
//#include "../lad_crs/lah_crs.h"

//#define THREADS 2
//#if THREADS == 1
#define MULTITHREADED_SUPERLU // should depend on standard MULTITHREADED option...
#ifdef MULTITHREADED_SUPERLU
  #include "./superlu_threads/slu_mt_ddefs.h"
  #include "./superlu_threads/slu_mt_util.h"
#else
  #include "./superlu_seq/slu_ddefs.h"
#endif

/* GLOBAL VARIABLES */
//extern int lsv_mkb_superlu_cur_solver_id;   /* ID of the current solver */
lst_mkb_superlu_solvers lsv_mkb_superlu_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */

/**-----------------------------------------------------------
   lsr_mkb_direct_init - to create a new solver instance, read its control
   parameters and initialize its data structure
------------------------------------------------------------*/
int lsr_mkb_direct_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  char* Filename,  /* in: name of the file with control parameters */
  int Monitoring_level /* Level of output, -1 for Filename */
  ){

/*kbw/
  printf("\nIn lsr_mkb_direct_init: Solver id %d, Filename %s, Monitor %d\n",
  Solver_id,  Filename, Monitoring_level);
/*kew*/

#ifdef PARALLEL
  //mf_log_warn("Using direct solver in distributed memory environment!");
  printf("MKB direct solvers not adapted to parallel (distributed memory) execution. Exiting!\n");
  exit(-1);
#endif


  /* initialize data structure
  lsv_mkb_superlu_solver[Solver_id].crs_row=NULL;
  lsv_mkb_superlu_solver[Solver_id].crs_col=NULL;
  lsv_mkb_superlu_solver[Solver_id].crs_val=NULL;
  lsv_mkb_superlu_solver[Solver_id].rhs=NULL;
  lsv_mkb_superlu_solver[Solver_id].offset=0;
  */

  printf("\nSuperLU\n");

  //lsv_mkb_superlu_cur_solver_id = Solver_id;
  lsv_mkb_superlu_solver[Solver_id].solver_id = Solver_id;
  lsv_mkb_superlu_solver[Solver_id].monitor = Monitoring_level;

  return(Solver_id);
}


/**-----------------------------------------------------------
   lsr_mkb_direct_create - to create a new solver instance, read its control
   parameters and initialize its data structure
------------------------------------------------------------*/
int lsr_mkb_direct_create( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  char* Filename,  /* in: name of the file with control parameters */
  int Monitoring_level /* Level of output, -1 for Filename */
  ){

/*kbw/
  printf("\nIn lsr_mkb_direct_init: Solver id %d, Filename %s, Monitor %d\n",
  Solver_id,  Filename, Monitoring_level);
/*kew*/

#ifdef PARALLEL
  //mf_log_warn("Using direct solver in distributed memory environment!");
  printf("MKB direct solvers not adapted to parallel (distributed memory) execution. Exiting!\n");
  exit(-1);
#endif


  /* initialize data structure */
  lsv_mkb_superlu_solver[Solver_id].crs_row=NULL;
  lsv_mkb_superlu_solver[Solver_id].crs_col=NULL;
  lsv_mkb_superlu_solver[Solver_id].crs_val=NULL;
  lsv_mkb_superlu_solver[Solver_id].rhs=NULL;
  lsv_mkb_superlu_solver[Solver_id].offset=0;
  lsv_mkb_superlu_solver[Solver_id].free_flag=0;



  //lsv_mkb_superlu_cur_solver_id = Solver_id;
  //lsv_mkb_superlu_solver[Solver_id].solver_id = Solver_id;
  //lsv_mkb_superlu_solver[Solver_id].monitor = Monitoring_level;

  return(Solver_id);
}





/**--------------------------------------------------------
  lsr_mkb_direct_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
---------------------------------------------------------*/
 int lsr_mkb_direct_solve( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
						/* <0 - error code */
	int Solver_id,  /* in: solver ID */
	int Comp_type,  /* in: indicator for the scope of computations: */
					/*   LSC_SOLVE - solve the system */
					/*   LSC_RESOLVE - resolve for the new right hand side */
	int Matrix_id,  /* in: matrix ID */
	int Ndof, 	/* in: 	the number of degrees of freedom */
	double* X, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* B,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int Monitor	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data, 3 - 2+iteration data */
	)
{

  lsv_mkb_superlu_solver[Solver_id].SM_and_LV_id = Matrix_id;
  //itt_crs_generic_matrices *it_matrix = &itv_crs_generic_matrices[Matrix_id];

  // if(Ndof!=it_matrix->Nrdofgl){
	// printf("Error passing argument Ndof %d (!= %d) in lsd_mkb_solve_direct. Exiting.\n",
	   // Ndof, it_matrix->Nrdofgl);
  // }

  int i,j,nnz,aux,kux;

  int offset;
  int *crs_row=NULL, *crs_col=NULL;
  double *crs_val=NULL, *rhs=NULL;

  SuperMatrix A, L, U, RHS;

  int *perm_r; /* row permutations from partial pivoting */
  int *perm_c; /* column permutation vector */


  /* Setup SUPERLU control parameters.*/
  int nrhs = 1, info; /* Number of right hand sides. */
  double ddum=0; /* Double dummy */
  int idum=0; /* Integer dummy. */
  int phase;

  int error = 0;

  int n = Ndof;
  //it_matrix->crs_row_ptr[n]=it_matrix->Nnz;

  //nnz=it_matrix->Nnz;

  offset=lsv_mkb_superlu_solver[Solver_id].offset;


  lsv_mkb_superlu_solver[Solver_id].free_flag = lar_get_SM_and_LV_crs( Matrix_id, offset,
	  &(lsv_mkb_superlu_solver[Solver_id].crs_row),
	  &(lsv_mkb_superlu_solver[Solver_id].crs_col),
	  &(lsv_mkb_superlu_solver[Solver_id].crs_val),
	  &(lsv_mkb_superlu_solver[Solver_id].rhs));

  crs_row=lsv_mkb_superlu_solver[Solver_id].crs_row;
  crs_col=lsv_mkb_superlu_solver[Solver_id].crs_col;
  crs_val=lsv_mkb_superlu_solver[Solver_id].crs_val;
  rhs=lsv_mkb_superlu_solver[Solver_id].rhs;

  nnz=crs_row[n];

  /* number of processors to use. */
  //int_t      nprocs=THREADS;
  int_t      nprocs;
  int_t      permc_spec;


  //if(getenv("OMP_NUM_THREADS")==NULL)nprocs=1;
  //else nprocs=atoi(getenv("OMP_NUM_THREADS"));
#ifdef MULTITHREADED
  nprocs=omp_get_max_threads();
#else
  nprocs = 1;
#endif

  /*kbw
	printf("nzval: ");
	for(i=0;i<it_matrix->Nnz;++i)printf("%lf  ",it_matrix->crs_val[i]);
	printf("\n\nKONIEC1 %lf\n",it_matrix->crs_val[it_matrix->Nnz-1]);

	for(i=0;i<=n;++i)printf("%d ",it_matrix->crs_row_ptr[i]);
	printf("\n\nKONIEC2 %d\n",it_matrix->crs_row_ptr[n]);

	for(i=0;i<it_matrix->Nnz;++i)printf("%d ",it_matrix->crs_col_ind[i]);
	printf("\n\nKONIEC3 %d\n",it_matrix->crs_col_ind[it_matrix->Nnz-1]);

	for(i=0;i<n;++i)printf("%lf ",it_matrix->rhs[i]);
	printf("\n\nKONIEC4 %lf\n",it_matrix->rhs[n-1]);
	/*kew*/

  dCreate_CompCol_Matrix(&A, n, n, nnz, crs_val, crs_col, crs_row, SLU_NR, SLU_D, SLU_GE);

dCreate_Dense_Matrix(&RHS, n, nrhs, rhs, n, SLU_DN, SLU_D, SLU_GE);

  //#if THREADS == 1
#ifndef MULTITHREADED_SUPERLU

  superlu_options_t options;
  SuperLUStat_t stat;

  if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = MMD_AT_PLUS_A;
  options.ReplaceTinyPivot=YES;

  options.Fact = DOFACT;
  options.Equil = YES;
  options.ColPerm = MMD_AT_PLUS_A;
  options.RowPerm = LargeDiag;
  options.ReplaceTinyPivot = YES;
  options.Trans = NOTRANS;
  options.SolveInitialized = NO;
  options.RefineInitialized = NO;
  options.PrintStat = NO;

  /* Initialize the statistics variables. */
  StatInit(&stat);
  /* Solve the linear system. */
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &RHS, &stat, &info);


  if(info==0){
	double *sol = (double*) ((DNformat*) RHS.Store)->nzval;
	//for(i=0;i<n;++i)printf("%lf ",sol[i]);
	// printf("przepisanie ");
	/* rewrite the solution */
	for (i=0;i<n;i++) {
	  X[i]=sol[i];

	}
  }else{printf("1234 SuperLU solver fails %d",info);}



  StatFree(&stat);


#else

  if ( !(perm_r = intMalloc(n)) ) SUPERLU_ABORT("Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(n)) ) SUPERLU_ABORT("Malloc fails for perm_c[].");


  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = 0: natural ordering
   *   permc_spec = 1: minimum degree ordering on structure of A'*A
   *   permc_spec = 2: minimum degree ordering on structure of A'+A
   *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
   */
  permc_spec = 1;
  get_perm_c(permc_spec, &A, perm_c);

  pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &RHS, &info);


  if(info==0){
	double *sol = (double*) ((DNformat*) RHS.Store)->nzval;
	//for(i=0;i<n;++i)printf("%lf ",sol[i]);

	/* rewrite the solution */
	for (i=0;i<n;i++) {
	  X[i]=sol[i];
	  //printf("X=%lf\n",X[i]);

	}
  }else{printf("18544 SuperLU_threads solver fails %d",info);}


#endif

  // /* De-allocate storage */
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_SuperMatrix_Store(&A);
  Destroy_SuperMatrix_Store(&RHS);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);


  return(1);

}


/**--------------------------------------------------------
  lsr_mkb_direct_free - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
int lsr_mkb_direct_free( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
			 )
{


  if(lsv_mkb_superlu_solver[Solver_id].free_flag != 0){
	free(lsv_mkb_superlu_solver[Solver_id].crs_row);
	free(lsv_mkb_superlu_solver[Solver_id].crs_col);
	free(lsv_mkb_superlu_solver[Solver_id].crs_val);
  }

  // currently nothing to do

  return(1);
}

/**--------------------------------------------------------
  lsr_mkb_direct_destroy - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
int lsr_mkb_direct_destroy( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
			 )
{

  // currently nothing to do

  return(1);
}


