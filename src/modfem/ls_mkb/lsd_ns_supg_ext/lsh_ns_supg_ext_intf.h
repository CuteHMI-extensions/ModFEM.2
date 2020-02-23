/*************************************************************
lsh_ns_supg_ext_intf.h - interface of the extensions of the mkb solver of 
   linear equations for ns_supg problem module in ModFEM. Similar to mkb, 
   the extension uses la... modules for storage of the stiffness matrix and basic
   linear algebra operations (see interface include/lah_intf.h ) and
   employs Block versions of standard iterative methods (Jacobi, Gauss-Seidel, 
   additive Schwarz, multiplicative Schwarz) and ILU(0) for preconditioning.
   The extension is switched on by providing special input file
   with the first line: NS_SUPG_SOLVER_DATA ,
   otherwise standard mkb operations are performed.

   The file contains the provided interface with the headers of routines
   called from the FEM code

   The actual solution phase (lsr_ns_supg_ext_solve procedure) has a parameter
   that controls whether the solver is a distributed memory parallel solver
   or a standard sequential solver. In the first case the solver calls
   the finite element external part to perform three operations:
       fem_vec_norm - to compute a global (inter-processor) vector norm
       fem_sc_prod - to compute a global (inter-processor) scalar product
       fem_exchange_dofs - to exchange the values of degrees of freedom
                      between subdomains (in overlapping Schwarz manner)
   For multigrid solution the solver calls the finite element part with
       fem_proj_sol_lev - to project solution between levels (grids)
   The required interface of the solver containing headers of the
   procedures above is defined in the file "lsh_ns_supg_fem_intf.h"

Contents (declarations of the following routines):
  lsr_ns_supg_ext_init - to create a new solver instance, read its control 
                         parameters and initialize its data structure
  lsr_ns_supg_ext_create_matrix - to allocate space for a global system matrix
  lsr_ns_supg_ext_create_precon - to create preconditioner blocks corresponding
                                  to small subdomains of neighboring elements
  lsr_ns_supg_ext_clear_matrix - to initialize data structure of system matrix
  lsr_ns_supg_ext_assemble_local_sm - to assemble entries to the global stiffness
                            matrix and the global load vector using the provided 
                            local stiffness matrix and load vector
  lsr_ns_supg_ext_show_matrix - to show matrix structure using graphic library
  lsr_ns_supg_ext_fill_precon - to prepare special preconditioner for ns_supg 
  lsr_ns_supg_ext_solve - to solve a system of equations, given previously 
                          constructed system matrix and preconditioner
  lsr_ns_supg_ext_free_matrix - to free space for solver data structure (system
                                matrix, preconditioner and possibly other)

History:
        08.2015 - Krzysztof Banas, initial version

*************************************************************/
	  

#ifndef _lsh_ns_supg_intf_
#define _lsh_ns_supg_intf_


#ifdef __cplusplus
extern "C" {
#endif


/**-----------------------------------------------------------
  lsr_ns_supg_ext_init - to create a new solver instance, read its control 
                         parameters and initialize its data structure
------------------------------------------------------------*/
extern int lsr_ns_supg_ext_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
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

/**--------------------------------------------------------
  lsr_ns_supg_ext_create_matrix - to allocate space for a global system matrix
---------------------------------------------------------*/
extern int lsr_ns_supg_ext_create_matrix( 
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id,    /* in: level ID */
  int Nrblocks,    /* in: number of DOF blocks */
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_sm_size, /* in: maximal size of the stiffness matrix */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
  );

/**--------------------------------------------------------
lsr_ns_supg_ext_create_precon - to create preconditioner 
---------------------------------------------------------*/
extern int lsr_ns_supg_ext_create_precon( /* returns:   >0 success code */
                          /*	       <=0 - error */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id           /* in: level ID */
  );

/**--------------------------------------------------------
  lsr_ns_supg_ext_clear_matrix - to the structure of system matrix
---------------------------------------------------------*/
extern int lsr_ns_supg_ext_clear_matrix( 
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id,    /* in: level ID */
  int Comp_type    /* in: indicator for the scope of computations: */
                   /*   LSC_SOLVE - solve the system */
                   /*   LSC_RESOLVE - resolve for the new rhs vector */
  );

/**-----------------------------------------------------------
  lsr_ns_supg_ext_assemble_local_sm - to assemble entries to the global stiffness
                            matrix and the global load vector using the provided
                            local stiffness matrix and load vector
------------------------------------------------------------*/
extern int lsr_ns_supg_ext_assemble_local_sm( 
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator for the scope of computations: */
                         /*   LSC_SOLVE - solve the system */
                         /*   LSC_RESOLVE - resolve for the new rhs vector */
  int Nr_dof_bl,         /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
  );

/**--------------------------------------------------------
lsr_ns_supg_ext_show_matrix - to show matrix structure using graphic library
---------------------------------------------------------*/
extern int lsr_ns_supg_ext_show_matrix( 
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id          /* in: level ID */

  );
  
  
/**--------------------------------------------------------
  lsr_ns_supg_ext_fill_precon - to prepare special ns_supg preconditioner 
---------------------------------------------------------*/
extern int lsr_ns_supg_ext_fill_precon(  
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id           /* in: level ID */
  );

/**--------------------------------------------------------
  lsr_ns_supg_ext_solve - to solve a system of equations, given previously 
             constructed system matrix and preconditioner
---------------------------------------------------------*/
extern int lsr_ns_supg_ext_solve( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
                        /* <0 - error code */
	int Solver_id,  /* in: solver ID */
	int Ndof, 	/* in: 	the number of degrees of freedom */
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

/**--------------------------------------------------------
  lsr_ns_supg_ext_free_matrix - to free space for solver data structure (system
                       matrix, preconditioner and possibly other)
---------------------------------------------------------*/
extern int lsr_ns_supg_ext_free_matrix(
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
  );


/**--------------------------------------------------------
lsr_ns_supg_ext_compreres - to compute the residual of the left preconditioned 
	system of equations, v = M^-1 * ( b - Ax )
        ( used also to compute the product v = M^-1 * Ax)
---------------------------------------------------------*/
extern void lsr_ns_supg_ext_compreres ( 
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

extern double lsr_ns_supg_ext_compres(
	/* returns: residual norm */
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	double* X, /* in: input vector */
	int Ndof /* in: input vector size */
	);

extern void lsr_ns_supg_ext_compres_vector(
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	int Ndof, 	/* in: number of unknowns (components of x) */
	double* X, 	/* in: input vector (may be NULL if Ini_zero==0) */
	double* B,	/* in:  the rhs vector, if NULL take rhs */
					/*      from block data structure ( B is not taken */
					/*      into account if Use_rhs!=1) */
	double* V 	/* out: v = b-Ax */
	);

#ifdef __cplusplus
}
#endif

#endif
