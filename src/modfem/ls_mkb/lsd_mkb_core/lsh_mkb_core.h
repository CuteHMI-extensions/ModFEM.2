/* lsh_mkb_core.h - internal data of the multigrid iterative solver */
/*    or multigrid preconditioned Krylow methods solver (currently GMRES only) */
/*    based on standard iterative methods  (Jacobi, Gauss-Seidel,  */
/*    additive Schwarz, multiplicative Schwarz) and ILU(0) as solvers, */
/*    smoothers and preconditioners: definition of parameters, */
/*    data types, global variables and internal functions    */

// The core procedures are usually called by some adapter like e.g.
// the ones in lsd_mkb or lsd_ns_supg directories, that implement
// an interface known to FEM code - in particular in ModFEM
// the FEM part of code interacting with the solver is a solver
// interface module with the current practice of one-to-one
// correspondence between interface and adapter modules (e.g.
// sid_mkb for lsd_mkd)

/* PROCEDURES */

/* in file "lss_mkb_core.c":

- solver management routines
lsr_mkb_core_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
lsr_mkb_core_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
lsr_mkb_core_destroy - to delete a solver instance

- preconditioner management routines
lsr_mkb_core_create_precon - to create preconditioner
lsr_mkb_core_fill_precon - to prepare preconditioner by factorizing the stiffness
					matrix, either only diagonal blocks or block ILU(0)
lsr_mkb_core_destroy_precon - to free preconditioner data structure

- solver algorithms
lsr_mkb_core_comp_norm_rhs - to compute the norm of the preconditioned rhs
lsr_mkb_core_standard - to solve the problem by standard iterations
lsr_mkb_core_vcycle - to perform one V-cycle of multigrid
lsr_mkb_core_smooth - to perform one smoothing step using different algorithms
lsr_mkb_core_precon - to perform preconditioning using different algorithms
lsr_mkb_core_compreres - to compute the residual of the left preconditioned
	system of equations, v = M^-1 * ( b - Ax )
		( used also to compute the product v = M^-1 * Ax)
lsr_mkb_core_solve_coarse - to launch solver for the coarse grid problem
//////////////////// basic restarted GMRES solver procedure ////////////////////
lsr_mkb_core_gmres - to solve a system of linear equations Ax = b
				 using the left preconditioned GMRES method

- additional utilities:
lsr_mkb_core_util_dvector - to allocate a double vector: name[0..ncom-1]:
lsr_mkb_core_util_ivector - to allocate an integer vector: name[0..ncom-1]:
lsr_mkb_core_util_imatrix - to allocate an integer matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
lsr_mkb_core_util_dmatrix - to allocate a double matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
lsr_mkb_core_util_chk_list - to check whether a number is on the list
lsr_mkb_core_util_put_list - to put Num on the list List with length Ll
lsr_mkb_core_util_d_zero - to zero a double vector
lsr_mkb_core_util_i_zero - to zero an integer vector
lsr_mkb_core_util_sort - to heap-sort an array
lsr_mkb_core_util_dgetrf - quasi-LU decomposition of a matrix
lsr_mkb_core_util_dgetrs - to perform forward reduction and back substitution
	of the RHS vector for solving a system of linear equations

///////////////////////////// GEOMETRIC MULTIGRID UTILITY //////////////////////
  lsr_mkb_core_get_pdeg_coarse - to get enforced pdeg for the coarse mesh
---------------------------------------------------------*/


#ifndef _lsh_mkb_core_
#define _lsh_mkb_core_

#include <stdio.h>

// API for all mkb solvers with some constants...
#include <modfem/ls_mkb/lsh_mkb_intf.h>


// API for linear algebra routines supporting MKB solver
#include <modfem/ls_mkb/lah_intf.h>


/* Parameters */
// The module can handle several solvers, each solver can have several subsystems
// and each subsystem can have several levels
#define LSC_MAX_NUM_SUBS 10


/* BLAS and LAPACK names with or without underscore */
#ifdef WITHOUT_
#define ddot_ ddot
#define dnrm2_ dnrm2
#define dscal_ dscal
#define dcopy_ dcopy
#define daxpy_ daxpy
#define dgemv_ dgemv
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#define drot_ drot
#define drotg_ drotg
#define dtrsv_ dtrsv
#endif


/* GMRES type */
#define STANDARD_GMRES	0
#define MATRIX_FREE	1

/* preconditioners-smoothers !!! MUST BE THE SAME AS IN LA PACKAGE */
// definition in lah_intf.h !!!

/* Types of convergence measures */
#define REL_RES_INI	0
#define ABS_RES		1
#define REL_RES_RHS	2


/*** Data types ***/


/* definition of type lst_mkb_core_levels - data structure for mesh levels */
/*	and associated solvers */
typedef struct {

/* control variables */
  int Solver;	/* linear equations solver: */
		/*	0 - direct solver */
		/*	1 - GMRES */
		/*	2 - multi-level GMRES */
  int GMRES_type;/* type of GMRES */
		/*	0 - standard */
		/*	1 - matrix free */
  int Krylbas;	/* number of Krylov space vectors for GMRES */
  int Pdeg;	/* indicator for choosing degree of approximation for coarse */
				/* meshes (if >=0 indicates the actual degree of approximation */
  int Precon;  	/* preconditioning type */
		/*	0 - block Jacobi */
		/*	1 - block Gauss-Seidel */
		/*	2 - additive Schwarz (BJ with global product) */
  int ILU_k;   // number of neighborhood rings for ILU(k)
  int Nr_prec;	/* number of block iterations for preconditioning */
  int Block_type;/* preconditioner block types: number of nodes in a block */
		/*    or some other indicator application dependent */
  int Max_iter;	/* maximum number of GMRES iterations to be performed */
  int Conv_type; /* residual convergence criterium: */
		/* 	0 - relative to initial residual */
		/* 	1 - absolute residual */
		/* 	2 - relative to rhs */
  double Conv_meas;	/* allowable convergence measure */
  int Monitor; 	/* monitoring flag: */
		/*	0 - silent run, 1 - warning messages */
		/*	2 - 1+restart data 3 - iteration data */

  int Nrdofgl;		/* total number of degrees of freedom */
  int SM_and_LV_id; // identifier assigned by linear algebra package to a triple:
					// system matrix, rhs vector, preconditioner data structure
  int Storage_type; // type of storage (LSC_STORAGE_ constants defined in lsh_mkb_intf.h)
  int Nr_pre_smooth;   // number of pre-smooth iterations
  int Nr_post_smooth;  // number of post-smooth iterations

  int Nr_dof_blocks;  // USED BY NS_SUPG_EXT ONLY!!!
} lst_mkb_core_levels;


/* definition of lst_mkb_core_solvers - data type for multi-level iterative solver */
typedef struct {

  int nr_level;	           /* number of levels in multi-level GMRES */
  int cur_level;           /* current level number in multi-level GMRES */
  lst_mkb_core_levels level[LSC_MAX_NUM_LEV];   /* array of solver data structures */
						   /* corresponding to different levels */
} lst_mkb_core_subsystems;

/* definition of lst_mkb_core_solvers - data type for multi-level iterative solver */
typedef struct {

  int solver_id;           /* solver_id */
  int parallel;            /* parameter specifying sequential (LSC_SEQUENTIAL) */
						   /* or parallel (LSC_PARALLEL) execution */
  int nr_subsystems;
  int cur_subsystem;
  lst_mkb_core_subsystems subsystem[LSC_MAX_NUM_SUBS];
						   /* array of subsystems handled by a solver */
} lst_mkb_core_solvers;


/* GLOBAL VARIABLES */
extern int lsv_mkb_core_cur_solver_id;   /* ID of the current solver */
extern lst_mkb_core_solvers lsv_mkb_core_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */


/* PROCEDURES */

/* in file "lss_mkb_core.c":

- solver management routines
lsr_mkb_core_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
lsr_mkb_core_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
lsr_mkb_core_destroy - to delete a solver instance

*/

extern int lsr_mkb_core_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Solver_type, // type of solver (as defined in problem input file)
  int Parallel,      /* parameter specifying sequential (LSC_SEQUENTIAL) */
					 /* or parallel (LSC_PARALLEL) execution */
  int* Max_num_levels_p,  /* in: number of levels for multigrid: */
				  /*     1 - enforce single level solver */
				  /*     >1 - enforce the number of levels */
				  /* out: actual number of levels !!! */
  char* Filename,  /* in: name of the file with control parameters */
  int Max_iter, /* maximal number of iterations, -1 for values from Filename */
  int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
  double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */
  int Monitoring_level /* Level of output, -1 for Filename */
  );


extern int lsr_mkb_core_solve( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
						/* <0 - error code */
	int Solver_id,      /* in: solver ID */
	int Comp_type,  /* in: indicator for the scope of computations: */
					/*   LSC_SOLVE - solve the system */
					/*   LSC_RESOLVE - resolve for the new right hand side */
	int* L_matrix_id, /* in: list of identifiers of SM and LV in lad_... module */
	int* L_ndof, 	/* in: 	the number of degrees of freedom */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* B,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int* Nr_iter, 	/* in:	the maximum iterations to be performed */
			/* out:	actual number of iterations performed */
	double* Toler, 	/* in:	tolerance level for chosen measure */
			/* out:	the final value of convergence measure */
	int Monitor,	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data, 3 - 2+iteration data */
	double* Conv_rate /* out: convergence rate */
			  );

extern int lsr_mkb_core_destroy( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
  );

/* in file "lss_mkb_core.c":

- preconditioner management routines
lsr_mkb_core_create_precon - to create preconditioner
lsr_mkb_core_fill_precon - to prepare preconditioner by factorizing the stiffness
					matrix, either only diagonal blocks or block ILU(0)
lsr_mkb_core_destroy_precon - to free preconditioner data structure

*/

extern int lsr_mkb_core_create_precon( /* returns:   >0 number of diagonal blocks */
						  /*	       <=0 - error */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,           /* in: level ID */
  int SM_and_LV_id
			   );

extern int lsr_mkb_core_fill_precon(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,           /* in: level ID */
  int SM_and_LV_id
			  );

extern int lsr_mkb_core_destroy_precon(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,           /* in: level ID */
  int SM_and_LV_id
  );


///////////////////////////// Solver algorithms
/* in file "lss_mkb_core.c":

- solver algorithms
lsr_mkb_core_gmres - to solve a system of linear equations Ax = b
	using the left preconditioned GMRES method
lsr_mkb_core_comp_norm_rhs - to compute the norm of the preconditioned rhs
lsr_mkb_core_standard - to solve the problem by standard iterations
lsr_mkb_core_vcycle - to perform one V-cycle of multigrid
lsr_mkb_core_smooth - to perform one smoothing step using different algorithms
lsr_mkb_core_precon - to perform preconditioning using different algorithms
lsr_mkb_core_compreres - to compute the residual of the left preconditioned
	system of equations, v = M^-1 * ( b - Ax )
		( used also to compute the product v = M^-1 * Ax)
lsr_mkb_core_solve_coarse - to launch solver for the coarse grid problem
*/

extern int lsr_mkb_core_gmres(	/* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Ndof, 	/* in: 	the number of degrees of freedom */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* B,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int Krylbas, 	/* in: 	number of Krylov space vectors */
	int*  Iter, 	/* in:	the maximum GMRES iterations to be performed */
			/* out:	actual number of iterations performed */
	double* Resid, 	/* in:	tolerance level for residual */
			/* out:	the final value of residual */
	double* Rel_res, /* in:	tolerance level for relative residual */
			/* out:	the final value of relative residual */
	double* Rhs_res, /* in:	tolerance level for ratio rhs/residual */
			/* out:	the final value of the ratio rhs/residual */
	int Monitor,	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data, 3 - 2+iteration data */
	double* Pconvr	/* out: convergence rate */
	);

extern double lsr_mkb_core_comp_norm_rhs (
				/* returns: the norm of the rhs vector */
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id	/* in: index of the mesh (level) */
	);

extern int lsr_mkb_core_standard( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Ndof, 	/* in: 	the number of degrees of freedom */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* B,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int*  Iter, 	/* in:	the maximum iterations to be performed */
			/* out:	actual number of iterations performed */
	double* Toler, 	/* in:	tolerance level for max norm of update */
			/* out:	the final value of max norm of update */
	int Monitor,	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data, 3 - 2+iteration data */
	double* Pconvr	/* out: convergence rate */
	);

extern void lsr_mkb_core_vcycle (
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Use_rhs,	/* in: flag for considering RHS */
	int Ini_zero,	/* in: flag for zero initial guess */
	double* X, 	/* in/out: initial guess vector and solution */
		double* B	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	);

extern void lsr_mkb_vcycle_amg(
	int Solver_id,      /* in: solver ID */
	int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Use_rhs,	/* in: flag for considering RHS */
	int Ini_zero,	/* in: flag for zero initial guess */
	double* X, 	/* in/out: initial guess vector and solution */
	double* B,	/* in:  the rhs vector, if NULL take rhs */
	int nr_pre, /* in: nr of presmoothing iterations */
	int nr_post /* in: nr of postsmoothing iterations */
	);

extern void lsr_mkb_core_smooth(
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Use_rhs,	/* in: flag for considering RHS */
	int Ini_zero,	/* in: flag for zero initial guess */
	double* X, 	/* in/out: initial guess vector and solution */
		double* B	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	);

extern void lsr_mkb_core_multi_precon (
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	double* X, 	/* in/out: initial guess vector and solution */
		double* B	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	);

extern void lsr_mkb_core_precon(
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	double* X, 	/* in/out: initial guess vector and solution */
		double* B	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	);

extern void lsr_mkb_core_compreres (
		int Solver_id,      /* in: pointer to solver data structure to be passed
						   to data structure dependent routines */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Control,	/* in: indicator whether to compute residual (1)
				   or matrix-vector product (0)  */
	int Ini_zero,	/* in: flag for zero input vector X */
	int Ndof, 	/* in: number of unknowns (components of X and V) */
	double* X, 	/* in: input vector */
		double* B,	/* in: the rhs vector, if NULL take rhs */
			/*     from block data structure */
	double* V 	/* out: output vector, V = M^-1 * ( B - A*X ) */
	);


extern int lsr_mkb_core_solve_coarse( /* returns: 1 - success; <=0 - error code*/
	int Solver_id,      /* in: solver ID */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X,	/* in/out: initial guess and solution vector */
	double* B	/* in: rhs vector */
	);


/* additional utilities:
lsr_mkb_core_util_dvector - to allocate a double vector: name[0..ncom-1]:
lsr_mkb_core_util_ivector - to allocate an integer vector: name[0..ncom-1]:
lsr_mkb_core_util_imatrix - to allocate an integer matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
lsr_mkb_core_util_dmatrix - to allocate a double matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
lsr_mkb_core_util_chk_list - to check whether a number is on the list
lsr_mkb_core_util_put_list - to put Num on the list List with length Ll
lsr_mkb_core_util_d_zero - to zero a double vector
lsr_mkb_core_util_i_zero - to zero an integer vector
lsr_mkb_core_util_sort - to heap-sort an array
lsr_mkb_core_util_dgetrf - quasi-LU decomposition of a matrix
lsr_mkb_core_util_dgetrs - to perform forward reduction and back substitution
	of the RHS vector for solving a system of linear equations
*/

extern double *lsr_mkb_core_util_dvector(
	/* return: pointer to allocated vector */
	int Ncom,  	/* in: number of components */
	char Error_text[]/* in: error text to be printed */
	);

extern int *lsr_mkb_core_util_ivector(
	/* return: pointer to allocated vector */
	int Ncom, 	/* in: number of components */
	char Error_text[]/* in: error text to be printed */
	);

extern int **lsr_mkb_core_util_imatrix( /* returns: pointer to array of pointers to integers */
	int Nrow, 	/* in: number of rows */
	int Ncol, 	/* in: number of columns */
	char Error_text[]/* in: text to print in case of error */
	);

extern double **lsr_mkb_core_util_dmatrix( /* returns: pointer to array of pointers to doubles */
	int Nrow, 	/* in: number of rows */
	int Ncol, 	/* in: number of columns */
	char Error_text[]/* in: text to print in case of error */
	);

extern int lsr_mkb_core_util_chk_list(	/* returns: */
			/* >0 - position on the list */
					/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	);

extern int lsr_mkb_core_util_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
				/*   0 - put on the list */
				/*  -1 - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
	);

extern void lsr_mkb_core_util_d_zero(double *Vec, int Num);

extern void lsr_mkb_core_util_i_zero(int *Vec, int Num);

extern void lsr_mkb_core_util_sort(
   int    *Ind_array,    /* in/out: index array for sorting */
   double *Val_array     /* in: array of values used for sorting */
   );

extern void lsr_mkb_core_util_dgetrf(double* a, int m, int* ips);

extern void lsr_mkb_core_util_dgetrs(double* a, int m, double* b, double* x, int* ips);



/**--------------------------------------------------------
lsr_mkb_core_util_skip_rest_of_line - to allow for comments in input files
---------------------------------------------------------*/
extern void lsr_mkb_core_util_skip_rest_of_line(
			  FILE *Fp  /* in: input file */
			  );


///////////////////////////// GEOMETRIC MULTIGRID UTILITY //////////////////////

/*---------------------------------------------------------
  lsr_mkb_core_get_pdeg_coarse - to get enforced pdeg for the coarse mesh
---------------------------------------------------------*/
extern int lsr_mkb_core_get_pdeg_coarse( // returns: enforced pdeg for the coarse mesh
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id /* in: level number */
				  );



// kept for sentimental reasons and bleak outlooks for the future - matrix free approach
/* extern void lsr_mkb_core_core_mfaiter( */
/*  lst_mkb_core_levels* level,/\* in: pointer to current level data structure *\/ */
/* 	int Use_rhs,	/\* in: 0 - no rhs, 1 - with rhs *\/ */
/* 	int Ini_zero,	/\* in: flag for zero initial guess *\/  */
//      int Nr_prec,  /* in: number of preconditioner iterations */
/* 	int Ndof,	/\* in: number of unknowns (components of v*) *\/  */
/* 	double* V,	/\* in,out: vector of unknowns updated *\/ */
/* 			/\* during the loop over subdomains *\/ */
/*         double* B	/\* in:  the rhs vector, if NULL take rhs *\/ */
/* 			/\*      from block data structure *\/ */
/* 	); */

/* extern void lsr_mkb_core_core_mfmiter( */
/*  lst_mkb_core_levels* level,/\* in: pointer to current level data structure *\/ */
/* 	int Use_rhs,	/\* in: 0 - no rhs, 1 - with rhs *\/ */
/* 	int Ini_zero,	/\* in: flag for zero initial guess *\/  */
//      int Nr_prec,  /* in: number of preconditioner iterations */
/* 	int Ndof,	/\* in: number of unknowns (components of v*) *\/  */
/* 	double* V,	/\* in,out: vector of unknowns updated *\/ */
/* 			/\* during the loop over subdomains *\/ */
/*         double* B	/\* in:  the rhs vector, if NULL take rhs *\/ */
/* 			/\*      from block data structure *\/ */
/* 	); */



#endif
