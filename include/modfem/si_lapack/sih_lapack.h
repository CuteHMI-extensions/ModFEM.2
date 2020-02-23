/************************************************************************
File sih_lapack.h - internal information for the interface with the lapack
					linear solver for dense matrices

Contains:
  - constants
  - data types
  - global variables (for the whole module)
  - function headers:
  sir_direct_module_introduce - to return the direct solver name
  sir_direct_init - to create a new iterative solver, read its control
			 parameters and initialize its data structure
  sir_direct_create - to create and initialize solver data structure
  sir_direct_solve - to solve the system for a given data
  sir_direct_free - to free memory for stiffness and preconditioner matrices
  sir_direct_destroy - to destroy a solver instance and make room for next solvers

auxiliary:
  sir_direct_select_solver - to return the pointer to a given solver
  sir_direct_assemble_stiff_mat - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  sir_direct_init - to create a new iterative solver, read its control
			 parameters and initialize its data structure
  sir_direct_create - to create and initialize solver data structure
  sir_direct_solve - to solve the system for a given data
  sir_direct_free - to free memory for stiffness and preconditioner matrices
  sir_direct_destroy - to destroy a solver instance and make room for next solvers

------------------------------
History:
		02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#ifndef _sih_lapack_
#define _sih_lapack_

#include <modfem/sih_intf.h>

/*** DATA TYPES ***/

/* dof structure with data useful for creating flexible interfaces*/
/*  between FEM code and different solvers */
typedef struct {
	int dof_ent_type;  /* type of the associated FEM code (mesh) entity */
	int dof_ent_id;    /* ID of the associated FEM code (mesh) entity */
	int nr_int_ent; /* number of  integration entities providing SMs and LVs*/
	int l_int_ent_index[SIC_MAX_INT_PER_DOF];
	/* list of integration entities providing SMs and LVs*/
	int nrdofs;      /* number of DOFs */
	int posglob;    /* position in a global stiffness matrix */
	int nrneig;     /* number of neighboring DOF structures */
	int l_neig[SIC_MAX_DOF_STR_NGB];  /* list of neighboring DOF structures */
}  sit_direct_dof_struct;


/* solver data for dense LAPACK procedures */
typedef struct {

	int problem_id;   /* ID of the problem associated with the solver */
	int nr_int_ent;   /* number of integration entities - entities that */
	/* provide solver with stiffness matrices and load vectors*/
	int nr_dof_ent;   /* number of dof entities - mesh entities with which */
	/* degrees of freedom are associated */
	int nrdofs_glob;     /* the global number of degrees of freedom */
	int max_dofs_int_ent; /* maximal number of dofs per integration entity, i.e. */
	/* maximal size of the local stiffness matrix */

	/* arrays for assembling local stiffness matrices into global stiffness matrix*/
	int * l_int_ent_type; /*list of types of entities providing local SMs and LVs */
	int * l_int_ent_id; /* list of ID's of entities providing local SMs and LVs */

	sit_direct_dof_struct * l_dof_struct; /* list of dof structures with data useful for */
	/* creating flexible interfaces between FEM code*/
	/* and different solvers */

	/* for each possible type of dof entity - its corresponding dof structure */
	int * l_dof_vert_to_struct;
	int * l_dof_edge_to_struct;
	int * l_dof_face_to_struct;
	int * l_dof_elem_to_struct;
	/* dimensions of the above arrays */
	int max_dof_vert_id;
	int max_dof_edge_id;
	int max_dof_face_id;
	int max_dof_elem_id;

	double * stiff_mat; /* the global stiffness matrix */
	double * rhs_vect;  /* the global right hand side vector */
	double * sol_vect;  /* the global solution vector */
	int  *  aux_vect;   /* auxiliary vector for pivoting information */


} sit_direct_solver;


/*** GLOBAL VARIABLES (for the solver module only) ***/

extern int siv_direct_nr_solvers;     /* the number of solvers in the problem */
extern int siv_direct_cur_solver_id;  /* ID of the current solver */
extern sit_direct_solver siv_direct_solvers[SIC_MAX_NUM_SOLV]; /* array of solvers */


/*** INTERFACE ROUTINES ***/

/**-----------------------------------------------------------
  sir_direct_module_introduce - to return the direct solver name
------------------------------------------------------------*/
extern int sir_direct_module_introduce(
		/** returns: >=0 - success code, <0 - error code */
		char * Solver_name /** out: the name of the solver */
);

/**-----------------------------------------------------------
  sir_direct_init - to create a new iterative solver, read its control parameters
			 and initialize its data structure
------------------------------------------------------------*/
extern int sir_direct_init( /** returns: >0 - solver ID, <0 - error code */
		int Parallel,      /** parameter specifying sequential (SIC_SEQUENTIAL) */
		/** or parallel (SIC_PARALLEL) execution */
		char * Filename /** in: name of the file with control parameters */
);


/**-----------------------------------------------------------
  sir_direct_create - to create and initialize solver data structure
------------------------------------------------------------*/
extern int sir_direct_create( /** returns: >=0 - success code, <0 - error code */
		int Solver_id,    /** in: solver identification */
		int Problem_id    /** ID of the problem associated with the solver */
);


/**-----------------------------------------------------------
sir_direct_solve - to solve the system for a given data
------------------------------------------------------------*/
extern int sir_direct_solve(/** returns: >=0 - success code, <0 - error code */
		int Solver_id,     /** in: solver identification */
		int Comp_type,     /** in: indicator for the scope of computations: */
		/**   SIC_SOLVE - solve the system */
		/**   SIC_RESOLVE - resolve for the new right hand side */
		int Monitor       /** in: monitoring flag with options: */
		/**   SIC_PRINT_NOT - do not print anything */
		/**   SIC_PRINT_ERRORS - print error messages only */
		/**   SIC_PRINT_INFO - print most important information */
		/**   SIC_PRINT_ALLINFO - print all available information */
);


/**-----------------------------------------------------------
  sir_direct_free - to free memory for stiffness and preconditioner matrices
------------------------------------------------------------*/
extern int sir_direct_free(/** returns: >=0 - success code, <0 - error code */
		int Solver_id     /** in: solver identification */
);


/**-----------------------------------------------------------
  sir_direct_destroy - to destroy a solver instance and make room for next solvers
------------------------------------------------------------*/
extern int sir_direct_destroy(/** returns: >=0 - success code, <0 - error code */
		int Solver_id   /** in: solver identification */
);


/*** AUXILIARY LOCAL PROCEDURES ***/

/**-----------------------------------------------------------
  sir_direct_select_solver - to return the pointer to a given solver
------------------------------------------------------------*/
extern sit_direct_solver * sir_direct_select_solver( /*returns: pointer to solver*/
		int Solver_id    /* in: solver identification */
);


/**-----------------------------------------------------------
  sir_direct_assemble_stiff_mat - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
------------------------------------------------------------*/
extern int sir_direct_assemble_stiff_mat(
		/* returns: >=0 - success code, <0 - error code */
		int Solver_id,         /* in: solver ID (used to identify the subproblem) */
		int Comp_type,         /* in: indicator for the scope of computations: */
		/*   SIC_SOLVE - solve the system */
		/*   SIC_RESOLVE - resolve for the new rhs vector */
		int Nr_dof_bl,         /* in: number of global dof blocks */
		/*     associated with the local stiffness matrix */
		int * L_bl_nrdofs,      /* in: list of dof blocks' nr dofs */
		int * L_bl_posglob,    /* in: list of blocks' global positions */
		double * Stiff_mat,    /* in: stiffness matrix stored columnwise */
		double * Rhs_vect,     /* in: rhs vector */
		char * Rewr_dofs        /* in: flag to rewrite or sum up entries */
		/*   'T' - true, rewrite entries when assembling */
		/*   'F' - false, sum up entries when assembling */
);

/**-----------------------------------------------------------
  sir_direct_init - to create a new iterative solver, read its control parameters
			 and initialize its data structure
------------------------------------------------------------*/
extern int sir_direct_init( /** returns: >0 - solver ID, <0 - error code */
		int Parallel,      /** parameter specifying sequential (SIC_SEQUENTIAL) */
		/** or parallel (SIC_PARALLEL) execution */
		char * Filename /** in: name of the file with control parameters */
);


/**-----------------------------------------------------------
  sir_direct_create - to create and initialize solver data structure
------------------------------------------------------------*/
extern int sir_direct_create( /** returns: >=0 - success code, <0 - error code */
		int Solver_id,    /** in: solver identification */
		int Problem_id    /** ID of the problem associated with the solver */
);


/**-----------------------------------------------------------
sir_direct_solve - to solve the system for a given data
------------------------------------------------------------*/
extern int sir_direct_solve(/** returns: >=0 - success code, <0 - error code */
		int Solver_id,     /** in: solver identification */
		int Comp_type,     /** in: indicator for the scope of computations: */
		/**   SIC_SOLVE - solve the system */
		/**   SIC_RESOLVE - resolve for the new right hand side */
		int Monitor       /** in: monitoring flag with options: */
		/**   SIC_PRINT_NOT - do not print anything */
		/**   SIC_PRINT_ERRORS - print error messages only */
		/**   SIC_PRINT_INFO - print most important information */
		/**   SIC_PRINT_ALLINFO - print all available information */
);


/**-----------------------------------------------------------
  sir_direct_free - to free memory for stiffness and preconditioner matrices
------------------------------------------------------------*/
extern int sir_direct_free(/** returns: >=0 - success code, <0 - error code */
		int Solver_id     /** in: solver identification */
);


/**-----------------------------------------------------------
  sir_direct_destroy - to destroy a solver instance and make room for next solvers
------------------------------------------------------------*/
extern int sir_direct_destroy(/** returns: >=0 - success code, <0 - error code */
		int Solver_id   /** in: solver identification */
);

#endif

