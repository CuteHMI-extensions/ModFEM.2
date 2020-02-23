/*************************************************************
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
#include "./lsh_mkb_pardiso.h"

// API for linear algebra routines supporting MKB solver
#include <modfem/ls_mkb/lah_intf.h>
// PARDISO interface needs access to internal data structures of SM
//#include "../lad_crs_generic/lah_crs_generic.h"
//#include "../lad_crs/lah_crs.h"

/* GLOBAL VARIABLES */
lst_mkb_pardiso_solvers lsv_mkb_pardiso_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */

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

  // set default value for pardiso config
  lsv_mkb_pardiso_solver[Solver_id].pardiso_config.mtype=1;
  lsv_mkb_pardiso_solver[Solver_id].pardiso_config.maxfct=1;
  lsv_mkb_pardiso_solver[Solver_id].pardiso_config.mnum=1;
  if(Monitoring_level==0 || Monitoring_level==1){
	lsv_mkb_pardiso_solver[Solver_id].pardiso_config.msglvl=Monitoring_level;
  }
  else{
	lsv_mkb_pardiso_solver[Solver_id].pardiso_config.msglvl=0;
  }
  lsv_mkb_pardiso_solver[Solver_id].pardiso_config.iparm1=0;


  //lsv_mkb_pardiso_cur_solver_id = Solver_id;
  lsv_mkb_pardiso_solver[Solver_id].solver_id = Solver_id;
  lsv_mkb_pardiso_solver[Solver_id].monitor = Monitoring_level;
  lsv_mkb_pardiso_solver[Solver_id].isPardisoConfigInitializedext = false;

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
  lsv_mkb_pardiso_solver[Solver_id].crs_row=NULL;
  lsv_mkb_pardiso_solver[Solver_id].crs_col=NULL;
  lsv_mkb_pardiso_solver[Solver_id].crs_val=NULL;
  lsv_mkb_pardiso_solver[Solver_id].rhs=NULL;
  lsv_mkb_pardiso_solver[Solver_id].offset=0;
  lsv_mkb_pardiso_solver[Solver_id].free_flag=0;
  lsv_mkb_pardiso_solver[Solver_id].analysis=1;

  lsv_mkb_pardiso_solver[Solver_id].iparm=(int*) malloc(64 * sizeof(int));
  int i;
  for (i = 0; i < 64; ++i) {
	lsv_mkb_pardiso_solver[Solver_id].pt[i] = 0;
	lsv_mkb_pardiso_solver[Solver_id].iparm[i] = 0;
  }
  /* iparm[1] = 1 for disabled default pardiso values */
  lsv_mkb_pardiso_solver[Solver_id].iparm[0] = 1;
  /* iparm[35] = 1 for zero-based indexing; = 0 for one-based indexing */
  lsv_mkb_pardiso_solver[Solver_id].iparm[34] = 1;

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
	){


  int i,j;
  int offset;
  int *crs_row, *crs_col;
  double *crs_val, *rhs;

  lsv_mkb_pardiso_solver[Solver_id].SM_and_LV_id = Matrix_id;

  offset=lsv_mkb_pardiso_solver[Solver_id].offset;

  rhs=lsv_mkb_pardiso_solver[Solver_id].rhs;

  lsv_mkb_pardiso_solver[Solver_id].free_flag = lar_get_SM_and_LV_crs( Matrix_id, offset,
	  &(lsv_mkb_pardiso_solver[Solver_id].crs_row),
	  &(lsv_mkb_pardiso_solver[Solver_id].crs_col),
	  &(lsv_mkb_pardiso_solver[Solver_id].crs_val),&rhs);

  crs_row=lsv_mkb_pardiso_solver[Solver_id].crs_row;
  crs_col=lsv_mkb_pardiso_solver[Solver_id].crs_col;
  crs_val=lsv_mkb_pardiso_solver[Solver_id].crs_val;



  // if(Ndof!=it_matrix->Nrdofgl){
	// printf("Error passing argument Ndof %d (!= %d) in lsd_mkb_solve_direct. Exiting.\n",
	   // Ndof, it_matrix->Nrdofgl);
  // }

  // // OFFSET TO LAR???
  // for(i=0;i<it_matrix->Nnz;i++)it_matrix->crs_col_ind[i]=it_matrix->crs_col_ind[i]+1;
  // for(i=0;i<=it_matrix->Nrdofgl;i++)it_matrix->crs_row_ptr[i]=it_matrix->crs_row_ptr[i]+1;


  // if(it_matrix->crs_row_ptr[Ndof] != it_matrix->Nnz+1){
	// printf("Error 0458 (%d != %d) in lsd_mkb_solve_direct. Exiting.\n",
	   // it_matrix->crs_row_ptr[Ndof], it_matrix->Nnz+1);
  // }




  /* Setup PARDISO control parameters.*/
  int nrhs = 1; /* Number of right hand sides. */
  double ddum=0; /* Double dummy */
  int idum=0; /* Integer dummy. */
  int phase;

  int *pt=lsv_mkb_pardiso_solver[Solver_id].pt;
  int *iparm=lsv_mkb_pardiso_solver[Solver_id].iparm;

  int mtype = lsv_mkb_pardiso_solver[Solver_id].pardiso_config.mtype;
  int maxfct = lsv_mkb_pardiso_solver[Solver_id].pardiso_config.maxfct;
  int mnum = lsv_mkb_pardiso_solver[Solver_id].pardiso_config.mnum;
  int error = 0;

  // ?????????
  //int msglvl = lsv_mkb_pardiso_solver[Solver_id].pardiso_config.msglvl;
  //msglvl=lsv_mkb_pardiso_solver[Solver_id].monitor;
  int msglvl = 1;

	int n = Ndof;

  if(lsv_mkb_pardiso_solver[Solver_id].analysis){
	lsv_mkb_pardiso_solver[Solver_id].analysis=0;
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n,
	  crs_val, crs_row,
	  crs_col, &idum, &nrhs, iparm,
	  &msglvl, rhs, X, &error);
	lsv_mkb_pardiso_solver[Solver_id].Ndof=Ndof;

  }

  phase = 23;//33;

  /*kbw
#define SAVE_LOG_TO_FILE
#ifdef SAVE_LOG_TO_FILE
  char filename[256];
  const char *program_name = getenv("_");

  if(strstr(program_name,"MOD_FEM_ns_mixed_prism_std_quad_d") != NULL ||
	 strstr(program_name,"MOD_FEM_ns_mixed_prism_std_quad") != NULL)
	{
	  sprintf(filename,"NS_MIXED_CRS_PARDISO_Solution.txt",program_name);
	}
  else if(strstr(program_name,"MOD_FEM_ns_supg_prism_std_quad_d") != NULL ||
	  strstr(program_name,"MOD_FEM_ns_supg_prism_std_quad") != NULL ||
	  strstr(program_name,"MOD_FEM_ns_supg_prism_std_d") != NULL ||
	  strstr(program_name,"MOD_FEM_ns_supg_prism_std") != NULL)
	{
	  sprintf(filename,"NS_SUPG_CRS_PARDISO_Solution.txt",program_name);
	}
  else
	{
	  sprintf(filename,"CRS_PARDISO_Solution.txt",program_name);
	}
  static int file_paraview_opened_by_program = 0;
  FILE *file = NULL;
  if(file_paraview_opened_by_program == 0) {
	file = fopen(filename,"w+");
	file_paraview_opened_by_program = 1;
  } else {
	file = fopen(filename,"a+");
  }

#else
  FILE *file = stdout;
#endif

  int row_ind;
  fprintf(file,"\n\n CRS on entry: row, column, value\n");
  for(i=0;i<n;++i){
	for(row_ind=crs_row[i];row_ind<crs_row[i+1];row_ind++){
	  fprintf(file,"%6d%6d%20.10lf\n",i,crs_col[row_ind],crs_val[row_ind]);
	}
  }

  fprintf(file,"\n RHS on entry: index, value\n");
  for(i=0;i<n;++i){
	fprintf(file,"%6d%20.10lf\n",i,rhs[i]);
  }

  fprintf(file,"\n SOLUTION [before solve] on entry: index, value\n");
  for(i=0;i<n;++i){
	fprintf(file,"%6d%20.10lf\n",i,X[i]);
  }

#ifdef SAVE_LOG_TO_FILE
  fclose(file);
#endif

  // Exit - uncomment if needed for testing
  // exit(0);

  //old version:
	//printf("\n n=%d num_nnz=%d mnum=%d mtype=%d phase=%d maxfct=%d msglvl=%d\n\n",it_matrix->Nrdofgl,it_matrix->Nnz, mnum, mtype,phase,maxfct,msglvl );

	//for(i=0;i<it_matrix->Nnz;++i)printf("%lf ",it_matrix->crs_val[i]);
	//printf("\n\nKONIEC1 %.2lf\n",it_matrix->crs_val[it_matrix->Nnz-1]);

	//for(i=0;i<=n;++i)printf("%d ",it_matrix->crs_row_ptr[i]);
	//printf("\n\nKONIEC2 %d\n",it_matrix->crs_row_ptr[n]);

	//for(i=0;i<it_matrix->Nnz;++i)printf("%d ",it_matrix->crs_col_ind[i]);
	//printf("\n\nKONIEC3 %d\n",it_matrix->crs_col_ind[it_matrix->Nnz-1]);

	//for(i=0;i<n;++i)printf("%lf ",it_matrix->rhs[i]);
	//printf("\n\nKONIEC4 %lf\n",it_matrix->rhs[n-1]);
/*kew*/


  if(msglvl == 1)
	{
	  printf("Starting PARDISO phase %i \n", phase);
	}

  assert(sizeof(_INTEGER_t) == sizeof(mnum));
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n,
	  crs_val, crs_row,
	  crs_col, &idum, &nrhs, iparm,
	  &msglvl, rhs, X, &error);

  if(msglvl == 1)
	{
	  printf("Running PARDISO finished with result: %i (0 - means no error)\n", error);
	  printf("Number of threads used by PARDISO: %i \n", iparm[2]);
	}



  /*kbw
#ifdef SAVE_LOG_TO_FILE
	file = fopen(filename,"a+");
#else
	file = stdout;
#endif
	fprintf(file,"\n SOLUTION [after solve] on entry: index, value\n");
	for(i=0;i<n;i++) { fprintf(file,"%6d%20.10lf\n",i,X[i]); }

#ifdef SAVE_LOG_TO_FILE
	fclose(file);
	printf("Press ENTER to continue...\n");
	getchar(); getchar();
#endif
  /*kew*/

  //free(iparm);

  return(1);

}

/**--------------------------------------------------------
  lsr_mkb_direct_free - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
int lsr_mkb_direct_free( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
			 )
{

	/* Setup PARDISO control parameters.*/
  int nrhs = 1; /* Number of right hand sides. */
  double ddum=0; /* Double dummy */
  int idum=0; /* Integer dummy. */
  int phase;
  int error = 0;
  double* X=NULL;
  int i;



  int n=lsv_mkb_pardiso_solver[Solver_id].Ndof;

  //Release all internal memory for all matrices in pardiso
  phase = -1;
  PARDISO(lsv_mkb_pardiso_solver[Solver_id].pt, &(lsv_mkb_pardiso_solver[Solver_id].pardiso_config.maxfct), &(lsv_mkb_pardiso_solver[Solver_id].pardiso_config.mnum), &(lsv_mkb_pardiso_solver[Solver_id].pardiso_config.mtype), &phase, &n, lsv_mkb_pardiso_solver[Solver_id].crs_val,
  lsv_mkb_pardiso_solver[Solver_id].crs_row, lsv_mkb_pardiso_solver[Solver_id].crs_col, &idum, &nrhs, lsv_mkb_pardiso_solver[Solver_id].iparm,
  &(lsv_mkb_pardiso_solver[Solver_id].pardiso_config.msglvl), lsv_mkb_pardiso_solver[Solver_id].rhs, X, &error);


  free(lsv_mkb_pardiso_solver[Solver_id].iparm);

  if(lsv_mkb_pardiso_solver[Solver_id].free_flag != 0){
	//printf("\nPardiso free\n");
	free(lsv_mkb_pardiso_solver[Solver_id].crs_row);
	free(lsv_mkb_pardiso_solver[Solver_id].crs_col);
	free(lsv_mkb_pardiso_solver[Solver_id].crs_val);
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

  //printf("\nPardiso destroy\n");
  // currently nothing to do

  return(1);
}

