/************************************************************************
File sih_intf.h - interface with the linear solver module of the code.

Contains declarations of constants and interface routines:
  sir_module_introduce - to return the solver name
  sir_solve_lin_sys - to perform the five steps below in one call
  sir_init - to create a new solver instance and read its control parameters
  sir_create - to create and initialize solver data structure
  sir_solve - to solve the system for a given data
  sir_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
  sir_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM
  sir_free - to free memory for a given solver instance
  sir_destroy - to destroy a solver instance

  sir_direct_solve_lin_sys - to solve the system of linear equations 
          using a direct solver for the current
          setting of parameters (with initialization and finalization phases)

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _sih_intf_
#define _sih_intf_

#ifdef __cplusplus
extern "C" {
#endif 

/** @defgroup SIM Solver
 *
 *  @{
 */

/*** CONSTANTS ***/

/** Options for parallel execution */
#define SIC_SEQUENTIAL 0    
#define SIC_PARALLEL 1

/** Options for solution procedure */
#define SIC_SOLVE 0     /** solve the system */ 
#define SIC_RESOLVE 1   /** resolve for the new right hand side */ 

#define SIC_SOLVER_ASSEMBLE 0
#define SIC_PROBLEM_ASSEMBLE 1

/** Monitoring options */
#define SIC_PRINT_NOT 0     /** do not print anything */ 
#define SIC_PRINT_ERRORS 1  /** print error messages only */
#define SIC_PRINT_INFO 2    /** print most important information */
#define SIC_PRINT_ALLINFO 3 /** print all available information */

#define SIC_MAX_NUM_SOLV 10 /** maximal number of solvers */
#define SIC_MAX_DOF_PER_INT 70 /*maximal number of dof_ent (blocks) per int_ent*/
  // 15 - tetra, 21 - prism, 27 - brick, but all without constrained approximation
  // 40 - testing (working for linear approximation with constraints) (does not work for constraint quadratic approximation [too small table])...
  // 54 - working for linear approximation with constraints (does not work for constraint quadratic approximation [too small table])...
  // 70 - working
#define SIC_MAX_INT_PER_DOF 100 /*maximal number of int_ent per dof_ent (block)*/
#define SIC_MAX_DOF_STR_NGB 296  /** maximal number of block's neighbors */
  // 200 does not work for constraint quadratic approximation [too small list])
  // 296 - working

#define SIC_LIST_END_MARK -1 /** number marking the end of list - cannot be put */

/*** INTERFACE ROUTINES ***/

/**-----------------------------------------------------------
  sir_module_introduce - to return the solver name
------------------------------------------------------------*/
extern int sir_module_introduce(  
                    /** returns: >=0 - success code, <0 - error code */
  char* Solver_name /** out: the name of the solver */
  );

/**-----------------------------------------------------------
  sir_solve_lin_sys - to solve the system of linear equations for the current
          setting of parameters (with initialization and finalization phases)
------------------------------------------------------------*/
extern int sir_solve_lin_sys( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,    /** ID of the problem associated with the solver */
  int Solver_type, /* type of solver: 0 - direct, >0 - iterative (passed to mkb_core) */
  int Parallel,      /** parameter specifying sequential (SIC_SEQUENTIAL) */
                     /** or parallel (SIC_PARALLEL) execution */
  int Nr_levels, // maximal number of levels for multigrid
  char* Filename,  /** in: name of the file with control parameters */
  int Max_iter, /** maximal number of iterations, -1 for values from Filename */
  int Error_type, /** type of error norm (stopping criterion), -1 for Filename*/
  double Error_tolerance, /** value for stopping criterion, -1.0 for Filename */ 
  int Monitoring_level /** Level of output, -1 for Filename */
  );


/**-----------------------------------------------------------
  sir_init - to create a new iterative solver, read its control parameters
             and initialize its data structure
------------------------------------------------------------*/
extern int sir_init( /** returns: >0 - solver ID, <0 - error code */
  int Solver_type, /* type of solver: 0 - direct, >0 - iterative (passed to mkb_core) */
  int Parallel,      /** parameter specifying sequential (SIC_SEQUENTIAL) */
                     /** or parallel (SIC_PARALLEL) execution */
  int Nr_levels, // maximal number of levels for multigrid
  char* Filename,  /** in: name of the file with control parameters */
  int Max_iter, /** maximal number of iterations, -1 for values from Filename */
  int Error_type, /** type of error norm (stopping criterion), -1 for Filename*/
  double Error_tolerance, /** value for stopping criterion, -1.0 for Filename */ 
  int Monitoring_level /** Level of output, -1 for Filename */
  );


/**-----------------------------------------------------------
  sir_create - to create and initialize solver data structure
------------------------------------------------------------*/
extern int sir_create( /** returns: >=0 - success code, <0 - error code */
  int Solver_id,    /** in: solver identification */
  int Problem_id    /** ID of the problem associated with the solver */
  );


/**-----------------------------------------------------------
sir_solve - to solve the system for a given data
------------------------------------------------------------*/
extern int sir_solve(/** returns: >=0 - success code, <0 - error code */
  int Solver_id,     /** in: solver identification */
  int Comp_type,     /** in: indicator for the scope of computations: */
                     /**   SIC_SOLVE - solve the system */
                     /**   SIC_RESOLVE - resolve for the new right hand side */
  int Ini_guess,     /** in: indicator whether to set initial guess (1), */
                     /**     or to initialize it to zero (0) */
  int Monitor,       /** in: monitoring flag with options: */
                     /**   SIC_PRINT_NOT - do not print anything */ 
                     /**   SIC_PRINT_ERRORS - print error messages only */
                     /**   SIC_PRINT_INFO - print most important information */
                     /**   SIC_PRINT_ALLINFO - print all available information */
  int* Nr_iter,      /** in:	the maximum iterations to be performed */
                     /** out:	actual number of iterations performed */
  double* Conv_meas, /** in:	tolerance level for chosen measure */
		     /** out:	the final value of convergence measure */
  double *Conv_rate  /** out (optional): the total convergence rate */ 
  );


/**-----------------------------------------------------------
  sir_free - to free memory for stiffness and preconditioner matrices
------------------------------------------------------------*/
extern int sir_free(/** returns: >=0 - success code, <0 - error code */
  int Solver_id     /** in: solver identification */
  );


/**-----------------------------------------------------------
  sir destroy - to destroy a solver instance and make room for next solvers
------------------------------------------------------------*/
extern int sir_destroy(/** returns: >=0 - success code, <0 - error code */
  int Solver_id   /** in: solver identification */
  );



/**-----------------------------------------------------------
  sir_direct_solve_lin_sys - to solve the system of linear equations 
          using a direct solver for the current
          setting of parameters (with initialization and finalization phases)
------------------------------------------------------------*/
extern int sir_direct_solve_lin_sys( /*returns: >=0-success code, <0-error code*/
  int Problem_id,    /** ID of the problem associated with the solver */
  int Parallel,      /** parameter specifying sequential (SIC_SEQUENTIAL) */
                     /** or parallel (SIC_PARALLEL) execution */
  char* Filename  /** in: name of the file with control parameters */
  );


/**-----------------------------------------------------------
  sir_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
------------------------------------------------------------*/
int sir_assemble_local_stiff_mat_with_table( 
                         /** returns: >=0 - success code, <0 - error code */
  int Solver_id,         /** in: solver ID (used to identify the subproblem) */
  int Level_id,          /** in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  int* Local_to_global_int_ent, /* part of the global local_to_global table */
  double* Stiff_mat,     /** in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /** in: rhs vector */
  char* Rewr_dofs         /** in: flag to rewrite or sum up entries */
                         /**   'T' - true, rewrite entries when assembling */
                         /**   'F' - false, sum up entries when assembling */
);

/**-----------------------------------------------------------
  sir_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM
------------------------------------------------------------*/
int sir_assemble_local_stiff_mat( 
                         /** returns: >=0 - success code, <0 - error code */
  int Solver_id,         /** in: solver ID (used to identify the subproblem) */
  int Level_id,          /** in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_dof_ent_type,   /* in: list of dof blocks' IDs */
  int* L_dof_ent_id,     /* in: list of dof blocks' IDs */
  int* L_dof_ent_nrdofs, /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /** in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /** in: rhs vector */
  char* Rewr_dofs         /** in: flag to rewrite or sum up entries */
                         /**   'T' - true, rewrite entries when assembling */
                         /**   'F' - false, sum up entries when assembling */
);
  
/** @} */ // end of group


#ifdef __cplusplus
}
#endif 

#endif
