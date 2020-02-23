/*************************************************************
lsh_mkb_intf.h - interface of the set of solvers that provide implementations
   to linear solver routines required by the sid_mkb solver interface package.
   The set currently includes direct solvers such as SuperLU, PARDISO, ViennaCL,
   as well as iterative solvers, with multigrid preconditioned (geometrically
   or algebraically) Krylow method solver (currently GMRES only) as default.
   The solvers use their own data structures or lad_... modules for storage
   of the stiffness matrix. The implementations reside in different subdirectories
   of lsd_mkb (with lsd_mkb_core containing multigrid preconditioned GMRES).
   lad_... modules additionally provides solvers with a set of linear algebra
   operations (see interface lah_intf.h ) to support basic operations on data
   structures (operations related mainly to iterative solvers).
   The file contains the provided interface with the headers of routines
   called from the FEM code

   The solver can handle several instances of solver data, each with
   a different set of control parameters, system matrix, etc.
   Hence it can be used to solve coupled problems (also transient)
   with two or more finite element problems solved simultaneously

   The actual solution phase (lsr_mkb_solve procedure) has a parameter
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
   procedures above is defined in the file "lsh_mkb_fem_intf.h"

Contents (declarations of the following routines):

- solver management
  lsr_mkb_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
  lsr_mkb_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
  lsr_mkb_destroy - to destroy a particular instance of the solver

- SM and LV management (using lah_intf.h API)
  lsr_mkb_create_matrix - to allocate space for a global system matrix
  lsr_mkb_clear_matrix - to initialize data structure of system matrix
  lsr_mkb_fill_assembly_table_int_ent - to create a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
  lsr_mkb_assemble_local_stiff_mat_with_table - to assemble entries to the global
						   stiffness matrix and the global load vector using the
						   provided local stiffness matrix, load vector and
						   the proper part of the global assembly table
  lsr_mkb_assemble_local_sm - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  lsr_mkb_free_matrix - to free space for solver data structure (system
					   matrix, preconditioner and possibly other)

//////////////////// utility to create data structures on GPU ///////////////
  lsr_create_solver_structures_accel
//////////////////// utility to get CRS parameters and pointers ///////////////
  lsr_get_crs_data

- preconditioner management (for iterative solvers)
  lsr_mkb_create_precon - to create preconditioner data structures
  lsr_mkb_fill_precon - to prepare preconditioner by factorizing the stiffness
					matrix, either only diagonal blocks or block ILU(0)

- multigrid utility
  lsr_mkb_get_pdeg_coarse - to get enforced pdeg for the coarse mesh

// three procedures for direct solvers acting through mkb interface
// IMPLEMENTED not in lss_mkb_intf but in lsd_direct_... modules...
  lsr_mkb_direct_init - to create a new solver instance, read its control
						 parameters and initialize its data structure
  lsr_mkb_direct_create - to create a new solver instance, read its control
   parameters and initialize its data structure
  lsr_mkb_direct_solve - to solve a system of equations, given previously constructed
						 system matrix, preconditioner
  lsr_mkb_direct_free - to destroy a particular instance of the direct solver
  lsr_mkb_direct_destroy - to destroy a particular instance of the direct solver


History:
		08.2013 - Krzysztof Banas, initial version

*************************************************************/


#ifndef _lsh_mkb_intf_
#define _lsh_mkb_intf_


#ifdef __cplusplus
extern "C" {
#endif

// linear algebra supporting interface
#include "lah_intf.h"

/* Scope of calculations */
#define LSC_SOLVE    0
#define LSC_RESOLVE  1

#define LSC_SEQUENTIAL 0
#define LSC_PARALLEL   1

#define LSC_MAX_NUM_SOLV 10
#define LSC_MAX_NUM_LEV  20

/* Monitoring options */
#define  LSC_SILENT      0
#define  LSC_ERRORS      1
#define  LSC_INFO        2
#define  LSC_ALLINFO     3

/* storage options */
#define LSC_STORAGE_UNDEFINED LAC_STORAGE_UNDEFINED
#define LSC_STORAGE_BLOCK LAC_STORAGE_BLOCK
#define LSC_STORAGE_CRS   LAC_STORAGE_CRS
#define LSC_STORAGE_BCRS  LAC_STORAGE_BCRS
#define LSC_STORAGE_CRS_GENERIC LAC_STORAGE_CRS_GENERIC
#define LSC_STORAGE_PETSC  LAC_STORAGE_PETSC

/* solvers */
#define DIRECT       0
#define GMRES        1
#define MULTI_GMRES  2
#define STANDARD_IT 10
#define MULTI_GRID  20
#define MKB_CORE_NS_SUPG_SOLVER 99
#define MULTI_GRID_AMG  125

/* definition of type lst_mkb_levels - data structure for mesh levels */
/*	and associated solvers and matrices */
typedef struct {
// the only relevant right now
  int SM_and_LV_id; // identifier assigned by linear algebra package to a triple:
                    // system matrix, rhs vector, preconditioner data structure
  int nrdofgl;		/* total number of degrees of freedom */
// possibly known - but not used
  //int nr_dof_blocks;	/* total number of blocks of DOFs */
  //int block_size; // size of SM blocks (-1 - non-constant size of blocks)
  int storage_type; // indicator of storage type: CRS, BCRS, block
                    // LSC_STORAGE_ constants defined above
} lst_mkb_levels;

/* definition of lst_mkb_solvers - data type for multi-level iterative solver */
typedef struct {
  int solver_id;           /* solver_id */
  int solver_type;  /* linear equations solver: */
            /*	0 - direct solver */
            /*	1 - GMRES */
            /*	2 - multi-level GMRES */
                    /* 99 - MKB_CORE_NS_SUPG_SOLVER */
  int parallel;            /* parameter specifying sequential (LSC_SEQUENTIAL) */
                           /* or parallel (LSC_PARALLEL) execution */
  int nr_levels;           /* number of levels in multi-level mkb */
  lst_mkb_levels level[LSC_MAX_NUM_LEV];   /* array of solver data structures */
} lst_mkb_solvers;

/* GLOBAL VARIABLES */
extern int lsv_mkb_nr_solvers; /* the number of solvers handled by the module */
extern int lsv_mkb_cur_solver_id;   /* ID of the current solver */
extern lst_mkb_solvers lsv_mkb_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */

/*** Parameters ***/

/**-----------------------------------------------------------
  lsr_mkb_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
------------------------------------------------------------*/
extern int lsr_mkb_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_type, // type of solver (as defined in problem input file)
  int Parallel,      /* parameter specifying sequential (LSC_SEQUENTIAL) */
                     /* or parallel (LSC_PARALLEL) execution */
  int* Max_num_levels_p,  /* in: number of levels for multigrid: */
                  /*     1 - enforce single level solver */
                  /*     >1 - enforce the number of levels */
                  /* out: actual number of levels !!! */
  int* Storage_type, /* in: requested storage type (NOT IMPLEMENTED YET) */
                      /* out: actual storage type */
  char* Filename,  /* in: name of the file with control parameters */
  int Max_iter, /* maximal number of iterations, -1 for values from Filename */
  int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
  double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */
  int Monitoring_level /* Level of output, -1 for Filename */
  );

/**--------------------------------------------------------
  lsr_mkb_destroy - to destroy a particular instance of the solver
---------------------------------------------------------*/
extern int lsr_mkb_destroy(
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
  );


/**--------------------------------------------------------
  lsr_mkb_create_matrix - to allocate space for a global system matrix
---------------------------------------------------------*/
extern int lsr_mkb_create_matrix(
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id,    /* in: level ID */
  int Storage_type,/* in: enforced storage type; if LSC_STORAGE_UNDEFINED */
                   /*     storage type is taken from input files or */
                   /*     decided by software (usually based on Block_size) */
  int Nrblocks,    /* in: number of DOF blocks */
  int Nrdof_glob,  /* in: total number of DOFs */
  int Block_size,  /* in: size of SM blocks (-1 - non-constant size) */
  int Max_sm_size, /* in: maximal size of the stiffness matrix */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
  );

/**--------------------------------------------------------
  lsr_mkb_clear_matrix - to initialize block structure of system matrix
---------------------------------------------------------*/
extern int lsr_mkb_clear_matrix(
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id,    /* in: level ID */
  int Comp_type    /* in: indicator for the scope of computations: */
                   /*   LSC_SOLVE - solve the system */
                   /*   LSC_RESOLVE - resolve for the new rhs vector */
  );

/**-----------------------------------------------------------
  lsr_mkb_fill_assembly_table_int_ent - to create a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
extern int lsr_mkb_fill_assembly_table_int_ent(
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int *Assembly_table_int_ent /* part of the global assembly table */
  );

/**-----------------------------------------------------------
  lsr_mkb_assemble_local_stiff_mat_with_table - to assemble entries to the global
						   stiffness matrix and the global load vector using the
						   provided local stiffness matrix, load vector and
						   the proper part of the global assembly table
------------------------------------------------------------*/
extern int lsr_mkb_assemble_local_stiff_mat_with_table(
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator for the scope of computations: */
                         /*   LSC_SOLVE - solve the system */
                         /*   LSC_RESOLVE - resolve for the new rhs vector */
  int Nr_dof_bl,         /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  int* Local_to_global_int_ent, /* part of the global local_to_global table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
  );

/**-----------------------------------------------------------
  lsr_mkb_assemble_local_sm - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
------------------------------------------------------------*/
extern int lsr_mkb_assemble_local_sm(
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
  lsr_mkb_free_matrix - to free space for solver data structure (system
					   matrix, preconditioner and possibly other)
---------------------------------------------------------*/
extern int lsr_mkb_free_matrix(
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
  );


/**--------------------------------------------------------
lsr_mkb_create_precon - to create preconditioner data structures
---------------------------------------------------------*/
extern int lsr_mkb_create_precon( /* returns:   >0 number of diagonal blocks */
                          /*	       <=0 - error */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id           /* in: level ID */
  );

/**--------------------------------------------------------
  lsr_mkb_fill_precon - to prepare preconditioner by factorizing the stiffness matrix,
				   either only diagonal blocks or block ILU(0)
---------------------------------------------------------*/
extern int lsr_mkb_fill_precon(
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id           /* in: level ID */
  );


/**--------------------------------------------------------
  lsr_mkb_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
---------------------------------------------------------*/
extern int lsr_mkb_solve( /* returns: convergence indicator: */
            /* 1 - convergence */
            /* 0 - noconvergence */
                        /* <0 - error code */
    int Solver_id,  /* in: solver ID */
    int Comp_type,  /* in: indicator for the scope of computations: */
                    /*   LSC_SOLVE - solve the system */
                    /*   LSC_RESOLVE - resolve for the new right hand side */
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


// three procedures for direct solvers acting through mkb interface
// IMPLEMENTED not in lss_mkb_intf but in lsd_direct_... modules...

/**-----------------------------------------------------------
   lsr_mkb_direct_init - to create a new solver instance, read its control
   parameters and initialize its data structure
------------------------------------------------------------*/
extern int lsr_mkb_direct_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  char* Filename,  /* in: name of the file with control parameters */
  int Monitoring_level /* Level of output, -1 for Filename */
             );

  /**-----------------------------------------------------------
   lsr_mkb_direct_create - to create a new solver instance, read its control
   parameters and initialize its data structure
------------------------------------------------------------*/
extern int lsr_mkb_direct_create( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  char* Filename,  /* in: name of the file with control parameters */
  int Monitoring_level /* Level of output, -1 for Filename */
               );

/**--------------------------------------------------------
  lsr_mkb_direct_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
---------------------------------------------------------*/
extern int lsr_mkb_direct_solve( /* returns: convergence indicator: */
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
              );

  /**--------------------------------------------------------
  lsr_mkb_direct_free - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
int lsr_mkb_direct_free( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
             );

/**--------------------------------------------------------
  lsr_mkb_direct_destroy - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
extern int lsr_mkb_direct_destroy( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
             );


//////////////////// utility to create data structures on GPU ///////////////
extern int lsr_create_solver_structures_accel( /* returns: >0 - success code, <0 - error code */
  int Solver_id,  /* in: solver ID */
  int Level_id    /* in: level ID */
                      );


//////////////////// utility to get CRS parameters and pointers ///////////////
extern int lsr_get_crs_data( /* returns: >0 - success code, <0 - error code */
  int Solver_id,  /* in: solver ID */
  int Level_id,    /* in: level ID */
  int* Nrdof_glob_p, /* output pointers: nrdof, nnz etc. - self explanatory */
  int* Nnz_p,
  int** Crs_col_ind_p,
  int** Crs_row_ptr_p,
  double** Crs_val_p,
  double** Rhs_p
              );

//////////////////////////// Multigrid utility ///////////////

/**--------------------------------------------------------
  lsr_mkb_get_pdeg_coarse - to get enforced pdeg for the coarse mesh
---------------------------------------------------------*/
extern int lsr_mkb_get_pdeg_coarse( // returns: enforced pdeg for the coarse mesh
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id /* in: level number */
  );

/** the end of provided interface */


#ifdef __cplusplus
}
#endif

#endif
