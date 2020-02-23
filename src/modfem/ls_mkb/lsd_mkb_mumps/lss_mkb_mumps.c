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
#include "./lsh_mkb_mumps.h"
//#include "./mumps/dmumps_c.h"
//#include "./mumps/mpi.h"

//#include "/usr/include/dmumps_c.h"
//#include "/usr/include/mumps_c_types.h"
//#include "dmumps_struc.h"


#include "dmumps_c.h"
//#include "mpi.h"

// API for linear algebra routines supporting MKB solver
#include <modfem/ls_mkb/lah_intf.h>

/* GLOBAL VARIABLES */
lst_mkb_mumps_solvers lsv_mkb_mumps_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */


#define DEBUG
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD 1987654

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
  /*
  // set default value for pardiso config
  lsv_mkb_mumps_solver[Solver_id].pardiso_config.mtype=1;
  lsv_mkb_mumps_solver[Solver_id].pardiso_config.maxfct=1;
  lsv_mkb_mumps_solver[Solver_id].pardiso_config.mnum=1;
  lsv_mkb_mumps_solver[Solver_id].pardiso_config.msglvl=1;
  lsv_mkb_mumps_solver[Solver_id].pardiso_config.iparm1=0;

 */

  printf("\nMUMBS %d\n",Solver_id);

  lsv_mkb_mumps_solver[Solver_id].solver_id = Solver_id;
  lsv_mkb_mumps_solver[Solver_id].monitor = Monitoring_level;


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
  lsv_mkb_mumps_solver[Solver_id].coo_row=NULL;
  lsv_mkb_mumps_solver[Solver_id].coo_col=NULL;
  lsv_mkb_mumps_solver[Solver_id].coo_val=NULL;
  lsv_mkb_mumps_solver[Solver_id].rhs=NULL;
  lsv_mkb_mumps_solver[Solver_id].offset=1;
  //lsv_mkb_mumps_solver[Solver_id].analysis=1;

 /*
  lsv_mkb_mumps_solver[Solver_id].iparm=(int*) malloc(64 * sizeof(int));
  int i;
  for (i = 0; i < 64; ++i) {
	lsv_mkb_mumps_solver[Solver_id].pt[i] = 0;
	lsv_mkb_mumps_solver[Solver_id].iparm[i] = 0;
  }
  */

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
  int *coo_row, *coo_col;
  double *coo_val, *rhs;
  //double* rhs_copy;

  DMUMPS_STRUC_C id;
  MUMPS_INT n = Ndof;
  MUMPS_INT nnz=0;

  //rhs_copy= (double *) malloc(n*sizeof(double));


  lsv_mkb_mumps_solver[Solver_id].SM_and_LV_id = Matrix_id;

  offset=lsv_mkb_mumps_solver[Solver_id].offset;

  //rhs=lsv_mkb_mumps_solver[Solver_id].rhs;

  lar_get_SM_and_LV_coo( Matrix_id, &nnz, offset,
	  &(lsv_mkb_mumps_solver[Solver_id].coo_row),
	  &(lsv_mkb_mumps_solver[Solver_id].coo_col),
	  &(lsv_mkb_mumps_solver[Solver_id].coo_val),
	  &(lsv_mkb_mumps_solver[Solver_id].rhs));

  coo_row=lsv_mkb_mumps_solver[Solver_id].coo_row;
  coo_col=lsv_mkb_mumps_solver[Solver_id].coo_col;
  coo_val=lsv_mkb_mumps_solver[Solver_id].coo_val;
  rhs=lsv_mkb_mumps_solver[Solver_id].rhs;


/*
  FILE *fp=fopen("mumps.txt", "w");
	int k=1;

	for(i=0;i<nnz;i++){
		if(k!=(coo_row)[i])fprintf(fp,"\n");
		fprintf(fp,"(%d,%d,%.12lf) ",(coo_row)[i]-1,(coo_col)[i]-1,(coo_val)[i]);
		k= (coo_row)[i];

	}fprintf(fp,"\n");

	for(i=0;i<Ndof;i++)fprintf(fp,"(%d %.12lf) ",i,rhs[i]);

	exit(0);
/**/



//for(i=0;i<n;i++)X[i]=rhs[i];


  #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

  id.par=1; id.comm_fortran=USE_COMM_WORLD;




// Matrix type
//  0 :A is unsymmetric.
//  1 :A is assumed to be symmetric positive definite so that pivots are taken from the diagonal
//  without numerical pivoting during the factorization.  With this option, non-positive definite
//  matrices that do not require pivoting can also be treated in certain cases (see remark below).
//  2 :A is general symmetric
id.sym=0;

  id.job=JOB_INIT;
  dmumps_c(&id);



	id.n = n; id.nz =nnz; id.irn=coo_col; id.jcn=coo_row;
  id.a = coo_val; id.rhs = rhs;


// MAtrix format
//  0 : assembled format
//  1 : elemental format
id.ICNTL(5)=0;





// Mumps debug options
  if(lsv_mkb_mumps_solver[Solver_id].monitor<=0){
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
  }
  else{
	id.ICNTL(1)=6; id.ICNTL(2)=6; id.ICNTL(3)=6; id.ICNTL(4)=0;
  }

#ifdef DEBUG
 //Outputs for debugging
 id.ICNTL(1)=6; id.ICNTL(2)=6; id.ICNTL(3)=6; id.ICNTL(4)=0;

 id.ICNTL(11)=2; //compute main statistics (norms, residuals, componentwise backward errors)
#endif




 // Call the MUMPS package.
 id.job=6;
 dmumps_c(&id);


	/* rewrite the solution */
	for (i=0;i<Ndof;i++) {
	  X[i]=rhs[i];
	}


 // Terminate instance
 id.job=JOB_END; dmumps_c(&id);

//for(i=0;i<50;++i)printf("r%d =%.12lf\n",i,X[i]);

  return(1);

}

/**--------------------------------------------------------
  lsr_mkb_direct_free - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
int lsr_mkb_direct_free( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
			 )
{


  //if(lsv_mkb_mumps_solver[Solver_id].offset!=0){
	free(lsv_mkb_mumps_solver[Solver_id].coo_row);
	free(lsv_mkb_mumps_solver[Solver_id].coo_col);
	free(lsv_mkb_mumps_solver[Solver_id].coo_val);
  //}
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

