/************************************************************************
File: pdh_intf.h - headers of routines providing the interface of the
                   problem dependent module that can be used by other modules,
                   especially solver interface module and the part
                   of approximation module performing numerical integration

Contains declarations of routines (some routines are just wrappers for 
                 approximation module routines - see example implementations):
  pdr_module_introduce - to return the problem module's name
  pdr_get_list_ent - to return to the solver module :
                          1. the list of integration entities - entities
                             for which stiffness matrices and load vectors are
                             provided by the FEM code
                          2. the list of DOF entities - entities with which  
                             there are dofs associated by the given approximation
  pdr_get_list_ent_coarse - the same as above but for COARSE level and
                            given the corresponding lists from the fine level
  pdr_get_max_num_grid_levels - for limiting nr_levels in multigrid 
                              based on mesh and field data
  pdr_create_assemble_stiff_mat - to create element stiffness matrices
                                   and assemble them to the global SM
  pdr_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
  pdr_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM

  pdr_comp_stiff_mat - to construct a stiffness matrix and a load vector for
                      some given mesh entity
  pdr_comp_stiff_mat_uncon - to construct UNCONSTRAINED stiffness matrix and 
                      a load vector for some given mesh entity
  pdr_err_indi - to compute an error indicator for an element
  pdr_read_sol_dofs - to read a vector of dofs associated with a given 
                   mesh entity from approximation field data structure
  pdr_write_sol_dofs - to write a vector of dofs associated with a given 
                   mesh entity to approximation field data structure
  pdr_L2_proj_sol - to project solution between elements of different generations
  pdr_renum_coeff - to return a coefficient being a basis for renumbering
  pdr_get_ent_pdeg - to return the degree of approximation index 
                      associated with a given mesh entity
  pdr_dof_ent_sons - to return a list of dof entity sons

  pdr_proj_sol_lev - to project solution between mesh levels
  pdr_vec_norm - to compute a norm of global vector in parallel
  pdr_sc_prod - to compute a scalar product of two global vectors
  pdr_create_exchange_tables - to create tables to exchange dofs 
  pdr_exchange_dofs - to exchange dofs between processors

Routines to manage coefficients of PDE's
  pdr_select_el_coeff - to select coefficients returned to approximation
                        routines for element integrals
  pdr_el_coeff - to return coefficients for element integrals

------------------------------  			
History:
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _pdh_intf_
#define _pdh_intf_

#ifdef __cplusplus 
extern "C" 
{
#endif

/** @defgroup PRD Problem Definition
 *
 *  @{
 */

#define PDC_PI 3.141592654
#define PDC_SQRT_PI 1.772453851
#define PDC_SQRT_3 1.732050808

#define PDC_MAXEQ 30 /** maximal number of equations - for allocating space */
#define PDC_MAXBCVAL 10 /** maximal number of BC of given type (Dirichlet etc.) */

#define PDC_INTERIOR -1 /** interior face */
#define PDC_BC_DIRI   1 /** Dirichlet boundary face */
#define PDC_BC_NEUM   2 /** Neumann boundary face */
#define PDC_BC_MIXED  3 /** Robin boundary face */

#define PDC_USE_CURRENT_PROBLEM_ID -1 // indicator passed as function argument to use
  // the ID of the current problem stored in pdv_...._current_problem_id 

/** types of integration and dof entities */
extern const int PDC_ELEMENT;
extern const int PDC_FACE   ;
extern const int PDC_EDGE   ;
extern const int PDC_VERTEX ;
extern const int PDC_MIXED_ELEMENT;
extern const int PDC_MIXED_FACE   ;
extern const int PDC_MIXED_EDGE   ;
extern const int PDC_MIXED_VERTEX ;

/** options for create_stiff_mat procedure */
extern const int PDC_NO_COMP  ; /** do not compute stiff matrix and rhs vector */
extern const int PDC_COMP_SM  ; /** compute entries to stiff matrix only */
extern const int PDC_COMP_RHS ; /** compute entries to rhs vector only */
extern const int PDC_COMP_BOTH; /** compute entries for sm and rhsv */
extern const int PDC_COMP_MM; /** compute entries for sm and rhsv */

// when solving for time step n (n=1,2,3,...) and iteration k (k=1,2,...),
// (i.e. finding U^n_k ); the parameters below are used as the arguments
// for apr_read/write_ent_dofs and apr_get_el_dofs
#define Current_solution_ID 1 //  U^(n)_(k) after solver
#define Previous_iteration_sol_ID 2 //  U^(n)_(k-1)
#define Previous_time_step_sol_ID 3 //  U^(n)_0 = U^(n-1)
// ! for k==1 could be U^(n-2) since U^(n-1) is in previous iteration U^(n)_0 !
#define Geometry_dofs_ID 999


  // types of adaptations
#define PDC_UNIFORM_DEREF -2 // uniform mesh refinement
#define PDC_UNIFORM_REF   -1 // uniform mesh refinement
#define PDC_NO_ADAPT    0 // adaptations based on the knowledge of exact solution
#define PDC_ADAPT_EXACT 1 // adaptations based on the knowledge of exact solution
#define PDC_ADAPT_ZZ    2 // adaptations based on Zienkiewicz-Zhu error estimate
#define PDC_ADAPT_EXPL  3 // adaptations based on explicit residual error estimate

//#define PDC_MAX_DOF_PER_INT 27 /*maximal number of dof_ent (blocks) per int_ent */
//#define PDC_MAX_DOF_PER_INT 40 // increased due to P2P2 mixed approximation
#define PDC_MAX_DOF_PER_INT 70 // increased to 70 for constraint quadratic approximation

  // TODO !!!!!!
// should be moved to local problem dependent include files
#include "uth_intf.h"
extern utt_patches *pdv_patches;


/*------------------------------------------------------------
  pdr_module_introduce - to return the problem module's name
------------------------------------------------------------*/
int pdr_module_introduce(
                  /* returns: >=0 - success code, <0 - error code */
  char* Problem_module_name /* out: the name of the problem module */
  );


/*********************************************************************
Procedures defining the interface between the FEM code (problem, approximation
and mesh modules) and the solver interface module - implementation of procedures
knows interfaces of all FEM code modules and the interface of the solver
interface module
*************************************************************************/

/**-----------------------------------------------------------
  pdr_get_list_ent - to return the list of integration entities - entities
                         for which stiffness matrices and load vectors are
                         provided by the FEM code to the solver module 
------------------------------------------------------------*/
extern int pdr_get_list_ent( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in:  problem (and solver) identification */
  int* Nr_int_ent,    /** out: number of integration entitites */
	/** GHOST ENTITIES HAVE NEGATIVE TYPE !!! */
  int** List_int_ent_type,/** out: list of types of integration entitites */
  int** List_int_ent_id,  /** out: list of IDs of integration entitites */
  int* Nr_dof_ent,    /** out: number of dof entities (entities with which there
		              are dofs associated by the given approximation) */
  int** List_dof_ent_type,/** out: list of types of integration entitites */
  int** List_dof_ent_id,  /** out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/** out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /** out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent/** out: maximal number of dofs per dof entity */
    );

/**-----------------------------------------------------------
  pdr_get_list_ent_coarse - to return the list of integration entities - entities
                         for which stiffness matrices and load vectors are
                         provided by the FEM code to the solver module,
                         and DOF entities - entities with which there are dofs 
                         associated by the given approximation for COARSE level
                         given the corresponding lists from the fine level 
------------------------------------------------------------*/
int pdr_get_list_ent_coarse( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,       /** in:  problem (and solver) identification */
  int Nr_int_ent_fine, /** in: number of integration entitites */
  int* List_int_ent_type_fine,/** in: list of types of integration entitites */
  int* List_int_ent_id_fine,  /** in: list of IDs of integration entitites */
  int Nr_dof_ent_fine, /** in: number of dof entities (entities with which there
		             are dofs associated by the given approximation) */
  int* List_dof_ent_type_fine,/** in: list of types of integration entitites */
  int* List_dof_ent_id_fine,  /** in: list of IDs of integration entitites */
  int* List_dof_ent_nrdofs_fine,/** in: list of no of dofs for 'dof' entity */
  int Nrdofs_glob_fine, /** in: global number of degrees of freedom (unknowns) */
  int Max_dof_per_ent_fine, /** in: maximal number of dofs per dof entity */
  int* Pdeg_coarse,  /** in: degree of approximation for coarse space */
  int* Nr_int_ent,    /** out: number of integration entitites */
  int** List_int_ent_type,/** out: list of types of integration entitites */
  int** List_int_ent_id,  /** out: list of IDs of integration entitites */
  int* Nr_dof_ent,    /** out: number of dof entities (entities with which there
		              are dofs associated by the given approximation) */
  int** List_dof_ent_type,/** out: list of types of integration entitites */
  int** List_dof_ent_id,  /** out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/** out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /** out: global number of degrees of freedom (unknowns) */
  int* Max_dof_per_ent/** out: maximal number of dofs per dof entity */
  );

/*------------------------------------------------------------
pdr_get_max_num_grid_levels - for limiting nr_levels in multigrid 
                              based on mesh and field data
------------------------------------------------------------*/
int pdr_get_max_num_grid_levels(
			   int Problem_id
			   );

/**-----------------------------------------------------------
   pdr_create_assemble_stiff_mat - to create element stiffness matrices
                                   and assemble them to the global SM
------------------------------------------------------------*/
extern int pdr_create_assemble_stiff_mat(
  int Problem_id, 
  int Level_id, 
  int Comp_type,         /** in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /** do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /** compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /** compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /** compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  int Nr_int_ent,
  int* L_int_ent_type,
  int* L_int_ent_id,
  int Nr_colors_elems, 
  int* L_color_index_elems,
  int Nr_colors_faces,
  int* L_color_index_faces,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int* Pos_first_dof_int_ent,
  int* Local_to_global,
  int Max_dofs_int_ent
);

/*------------------------------------------------------------
  pdr_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
------------------------------------------------------------*/
int pdr_assemble_local_stiff_mat_with_table( 
                         /* returns: >=0 - success code, <0 - error code */
  int Problem_id,        /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  int* Local_to_global_int_ent, /* part of the global local_to_global table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs        /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
					     );

/**-----------------------------------------------------------
  pdr_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM
------------------------------------------------------------*/
extern int pdr_assemble_local_stiff_mat( 
                         /** returns: >=0 - success code, <0 - error code */
  int Problem_id,         /** in: solver ID (used to identify the subproblem) */
  int Level_id,          /** in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof entities */
                         /*     associated with the local stiffness matrix */
  int* L_dof_ent_type,   /* in: list of dof entities' types */
  int* L_dof_ent_id,     /* in: list of dof entities' IDs */
  int* L_dof_ent_nrdofs, /* in: list of dof entities' numbers of dof */
  double* Stiff_mat,     /** in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /** in: rhs vector */
  char* Rewr_dofs         /** in: flag to rewrite or sum up entries */
                         /**   'T' - true, rewrite entries when assembling */
                         /**   'F' - false, sum up entries when assembling */
);

/**-----------------------------------------------------------
  pdr_comp_stiff_mat - to construct a stiffness matrix and a load vector for
                      some given mesh entity
------------------------------------------------------------*/
extern int pdr_comp_stiff_mat( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in: approximation field ID  */
  int Int_ent_type,    /** in: type of the integration entity */ 
  int Int_ent_id,    /** in: unique identifier of the integration entity */ 
  int Comp_sm,       /** in: indicator for the scope of computations: */
                     /**   PDC_NO_COMP  - do not compute anything */
                     /**   PDC_COMP_SM - compute entries to stiff matrix only */
                     /**   PDC_COMP_RHS - compute entries to rhs vector only */
                     /**   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg,        /** in: enforced degree of polynomial (if != NULL ) */
  int* Nr_dof_ent,   /** in: size of List_dof_ent_... arrays, */
                /** out: no of filled entries, i.e. number of mesh entities*/
                /** with which dofs and stiffness matrix blocks are associated */
  int* List_dof_ent_type, /** out: list of types for 'dof' entity */
  int* List_dof_ent_id,   /** out: list of IDs for 'dof' entity */
  int* List_dof_ent_nrdofs, /** out: list of no of dofs for 'dof' entity */
  int* Nrdofs_loc,        /** in(optional): size of Stiff_mat and Rhs_vect */
                /** out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /** out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect,       /** out(optional): rhs vector */
  char* Rewr_dofs        /** out(optional): flag to rewrite or sum up entries */
                         /**   'T' - true, rewrite entries when assembling */
                         /**   'F' - false, sum up entries when assembling */
  );

/**-----------------------------------------------------------
  pdr_comp_stiff_mat_uncon - to construct UNCONSTRAINED stiffness matrix and 
                      a load vector for some given mesh entity
------------------------------------------------------------*/
extern int pdr_comp_stiff_mat_uncon( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in: approximation field ID  */
  int Int_ent_type,    /** in: type of the integration entity */ 
  int Int_ent_id,    /** in: unique identifier of the integration entity */ 
  int Comp_sm,       /** in: indicator for the scope of computations: */
                     /**   PDC_NO_COMP  - do not compute anything */
                     /**   PDC_COMP_SM - compute entries to stiff matrix only */
                     /**   PDC_COMP_RHS - compute entries to rhs vector only */
                     /**   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg,        /** in: enforced degree of polynomial (if != NULL ) */
  int* Nrdofs_loc,        /** in(optional): size of Stiff_mat and Rhs_vect */
                /** out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /** out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect,       /** out(optional): rhs vector */
  char* Rewr_dofs        /** out(optional): flag to rewrite or sum up entries */
                         /**   'T' - true, rewrite entries when assembling */
                         /**   'F' - false, sum up entries when assembling */
  );

/**-----------------------------------------------------------
pdr_err_indi - to compute an error indicator for an element
------------------------------------------------------------*/
extern double pdr_err_indi (	/** returns error indicator for an element */
        int Problem_id,	/** in: data structure to be used  */
        int Mode,	/** in: mode of operation */
        int El		/** in: element number */
       // double *Sigm /*out*/
        );

/**-----------------------------------------------------------
  pdr_read_sol_dofs - to read a vector of dofs associated with a given 
                   mesh entity from approximation field data structure
------------------------------------------------------------*/
extern int pdr_read_sol_dofs(/** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in: solver ID (used to identify the subproblem) */
  int Dof_ent_type,
  int Dof_ent_id,
  int Nrdofs,
  double* Vect_dofs  /** in: dofs to be written */
  );

/**-----------------------------------------------------------
  pdr_write_sol_dofs - to write a vector of dofs associated with a given 
                   mesh entity to approximation field data structure
------------------------------------------------------------*/
extern int pdr_write_sol_dofs(/** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in: solver ID (used to identify the subproblem) */
  int Dof_ent_type,
  int Dof_ent_id,
  int Nrdofs,
  double* Vect_dofs  /** in: dofs to be written */
  );

/**--------------------------------------------------------
  pdr_L2_proj_sol - to project solution between elements of different generations
----------------------------------------------------------*/
extern int pdr_L2_proj_sol(
  int Problem_id, /** in: problem ID */
  int El,		/** in: element number */
  int *Pdeg_vec,	/** in: element degree of approximation */
  double* Dofs,	/** out: workspace for degress of freedom of El */
  /** 	NULL - write to  data structure */
  int* El_from,	/** in: list of elements to provide function */
  int *Pdeg_vec_from,	/** in: degree of polynomial for each El_from */
     /** if element Pdeg is a vector - its components must be suitably placed */
     /** in Pdeg_vec_from, according to conventions in approximation module */
  double* Dofs_from /** in: Dofs of El_from or...*/
  );

/**--------------------------------------------------------
pdr_renum_coeff - to return a coefficient being a basis for renumbering
----------------------------------------------------------*/
  
  extern int pdr_renum_coeff(
  int Problem_id, /** in: problem ID */
  int Ent_type,	/** in: type of mesh entity */
  int Ent_id,	/** in: mesh entity ID */
  double* Ren_coeff  /** out: renumbering coefficient */
  );


/**-----------------------------------------------------------
  pdr_get_ent_pdeg - to return the degree of approximation index 
                      associated with a given mesh entity
------------------------------------------------------------*/
extern int pdr_get_ent_pdeg( /** returns: >0 - approximation index,
                                          <0 - error code */
  int Problem_id,     /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id         /** in: mesh entity ID */
  );


/**--------------------------------------------------------
  pdr_dof_ent_sons - to return a list of dof entity sons
---------------------------------------------------------*/
extern int pdr_dof_ent_sons( /** returns: success >=0 or <0 - error code */
  int Problem_id,     /** in: problem ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id,        /** in: mesh entity ID */
  int *Ent_sons     /** out: list of dof entity sons */
               	     /** 	Ent_sons[0] - number of sons */
  );

/**--------------------------------------------------------
  pdr_proj_sol_lev - to L2 project solution dofs between mesh levels
---------------------------------------------------------*/
extern int pdr_proj_sol_lev( /** returns: >=0 - success; <0 - error code*/
  int Problem_id, /** in: problem ID */
  int Solver_id,        /** in: solver data structure to be used */
  int Ilev_from,    /** in: level number to project from */
  double* Vec_from, /** in: vector of values to project */
  int Ilev_to,      /** in: level number to project to */
  double* Vec_to    /** out: vector of projected values */
  );

/**--------------------------------------------------------
  pdr_vec_norm - to compute a norm of global vector (in parallel)
---------------------------------------------------------*/
extern double pdr_vec_norm( /** returns: L2 norm of global Vector */
  int Problem_id, /** in: problem ID */
  int Solver_id,        /** in: solver data structure to be used */
  int Level_id,         /** in: level number */
  int Nrdof,            /** in: number of vector components */
  double* Vector        /** in: local part of global Vector */
  );

/**--------------------------------------------------------
  pdr_sc_prod - to compute a scalar product of two global vectors 
---------------------------------------------------------*/
extern double pdr_sc_prod( /** retruns: scalar product of Vector1 and Vector2 */
  int Problem_id, /** in: problem ID */
  int Solver_id,        /** in: solver data structure to be used */
  int Level_id,         /** in: level number */
  int Nrdof,           /** in: number of vector components */
  double* Vector1,     /** in: local part of global Vector */
  double* Vector2      /** in: local part of global Vector */
  );

/**--------------------------------------------------------
  pdr_create_exchange_tables - to create tables to exchange dofs 
---------------------------------------------------------*/
extern int pdr_create_exchange_tables( 
                   /** returns: >=0 -success code, <0 -error code */
  int Problem_id, /** in: problem ID */
  int Solver_id,        /** in: solver data structure to be used */
  int Level_id,    /** in: level ID */
  int Nr_dof_ent,  /** in: number of DOF entities in the level */
  /** all four subsequent arrays are indexed by block IDs with 1(!!!) offset */
  int* L_dof_ent_type,/** in: list of DOF entities associated with DOF blocks */
  int* L_dof_ent_id,  /** in: list of DOF entities associated with DOF blocks */
  int* L_bl_nrdof,    /** in: list of nrdofs for each dof block */
  int* L_bl_posg,     /** in: list of positions within the global */
                      /**     vector of dofs for each dof block */
  int* L_elem_to_bl,  /** in: list of DOF blocks associated with DOF entities */
  int* L_face_to_bl,  /** in: list of DOF blocks associated with DOF entities */
  int* L_edge_to_bl,  /** in: list of DOF blocks associated with DOF entities */
  int* L_vert_to_bl  /** in: list of DOF blocks associated with DOF entities */
  );

/**--------------------------------------------------------
  pdr_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
extern int pdr_exchange_dofs(
  int Problem_id, /** in: problem ID */
  int Solver_id,        /** in: solver data structure to be used */
  int Level_id,         /** in: level number */
  double* Vec_dofs  /** in: vector of dofs to be exchanged */
);

/*********************************************************************
Routines to manage coefficients of PDE's
*************************************************************************/

/*!!!!!!!!!!! THE FUNCTION BELOW IS USED BY APPROXIMATION MODULE !!!!!!!!!!!*/
/** !!! DO NOT CHANGE ITS HEADER (OR CHANGE APPROXIMATION MODULE AS WELL) !!!*/
/**-----------------------------------------------------------
  pdr_select_el_coeff_vect - to select coefficients returned to approximation
                        routines for element integrals in weak formulation
           (the procedure indicates which terms are non-zero in weak form)
------------------------------------------------------------*/
int pdr_select_el_coeff_vect( // returns success indicator
  int Problem_id,
  int *Coeff_vect_ind	/** out: coefficient indicator */
  // input to pdr_select_el_coeff_vect:
  // 0 - perform substitutions in coeff_vect_ind only if coeff_vect_ind[0]==0
  // output from pdr_select_el_coeff_vect:
  // 0 - coeff_vect_ind[0]==1 - always
  // 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz, 
  // 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
  // 14 - tx, 15 - ty, 16 - tz, 17 - cval
  // 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval 
			      );


/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
/**-----------------------------------------------------------
  pdr_select_el_coeff - to select coefficients returned to approximation
                        routines for element integrals
------------------------------------------------------------*/
extern double* pdr_select_el_coeff( 
			  /** returns: pointer !=NULL to indicate selection */
  int Problem_id,
  double **Mval,	/** out: mass matrix coefficient */
  double **Axx,double **Axy,double **Axz, /** out:diffusion coefficients, e.g.*/
  double **Ayx,double **Ayy,double **Ayz, /** Axy denotes scalar or matrix */
  double **Azx,double **Azy,double **Azz, /** related to terms with dv/dx*du/dy */
  /** second order derivatives in weak formulation (scalar for scalar problems */
  /** matrix for vector problems) */
  /** WARNING: if axy==NULL only diagonal (axx, ayy, azz) terms are considered */
  /** in apr_num_int_el */
  /** OPTIONS: */
  /** azz!=NULL, axy!=NULL - all a.. matrices must be specified */
  /** azz!=NULL, axy==NULL - axx, ayy, azz matrices must be specified */
  /** azz==NULL - axx, axy, ayx, ayy matrices must be specified */
  double **Bx,double **By,double **Bz,	/** out: convection coefficients */
  /** Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
  double **Tx,double **Ty,double **Tz,	/** out: convection coefficients */
  /** Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
  double **Cval,/** out: reaction coefficients - for terms without derivatives */
  /**  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
  double **Lval,/** out: rhs coefficient for time term, Lval denotes scalar */
  /** or matrix corresponding to time derivative - similar as mass matrix but  */
  /** with known solution at the previous time step (usually denoted by u_n) */
  double **Qx,/** out: rhs coefficients for terms with derivatives */
  double **Qy,/** Qy denotes scalar or matrix corresponding to terms with dv/dy */
  double **Qz,/** derivatives in weak formulation */
  double **Sval	/** out: rhs coefficients without derivatives (source terms) */
  );

/*!!!!!!!!!!! THE FUNCTION BELOW IS USED BY APPROXIMATION MODULE !!!!!!!!!!!*/
/** !!! DO NOT CHANGE ITS HEADER (OR CHANGE APPROXIMATION MODULE AS WELL) !!!*/
/**--------------------------------------------------------
  pdr_el_coeff - to return coefficients for element integrals
----------------------------------------------------------*/
extern int pdr_el_coeff(
  int Problem_id,
  int Elem,	/** in: element number */
  int Mat_num,	/** in: material number */
  double Hsize,	/** in: size of an element */
  int Pdeg,	/** in: local degree of polynomial */
  double *X_loc,      /** in: local coordinates of point within element */
  double *Base_phi,   /** in: basis functions */
  double *Base_dphix, /** in: x-derivatives of basis functions */
  double *Base_dphiy, /** in: y-derivatives of basis functions */
  double *Base_dphiz, /** in: z-derivatives of basis functions */
  double *Xcoor,	/** in: global coordinates of a point */
  double *Uk_val, 	/** in: computed solution from previous iteration */
  double *Uk_x, 	/** in: gradient of computed solution Uk_val */
  double *Uk_y,   	/** in: gradient of computed solution Uk_val */
  double *Uk_z,   	/** in: gradient of computed solution Uk_val */
  double *Un_val, 	/** in: computed solution from previous time step */
  double *Un_x, 	/** in: gradient of computed solution Un_val */
  double *Un_y,   	/** in: gradient of computed solution Un_val */
  double *Un_z,   	/** in: gradient of computed solution Un_val */
  double *Mval,	/** out: mass matrix coefficient */
  double *Axx, double *Axy, double *Axz,  /** out:diffusion coefficients */
  double *Ayx, double *Ayy, double *Ayz,  /** e.g. Axy denotes scalar or matrix */
  double *Azx, double *Azy, double *Azz,  /** related to terms with dv/dx*du/dy */
  /** second order derivatives in weak formulation (scalar for scalar problems */
  /** matrix for vector problems) */
  /** WARNING: if axy==NULL only diagonal (axx, ayy, azz) terms are considered */
  /** in apr_num_int_el */
  /** OPTIONS: */
  /** azz!=NULL, axy!=NULL - all a.. matrices must be specified */
  /** azz!=NULL, axy==NULL - axx, ayy, azz matrices must be specified */
  /** azz==NULL - axx, axy, ayx, ayy matrices must be specified */
  double *Bx, double *By, double *Bz,	/** out: convection coefficients */
  /** Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
  double *Tx, double *Ty, double *Tz,	/** out: convection coefficients */
  /** Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
  double *Cval,	/** out: reaction coefficients - for terms without derivatives */
  /**  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
  double *Lval,	/** out: rhs coefficient for time term, Lval denotes scalar */
  /** or matrix corresponding to time derivative - similar as mass matrix but  */
  /** with known solution at the previous time step (usually denoted by u_n) */
  double *Qx, /** out: rhs coefficients for terms with derivatives */
  double *Qy, /** Qy denotes scalar or matrix corresponding to terms with dv/dy */
  double *Qz, /** derivatives in weak formulation */
  double *Sval	/** out: rhs coefficients without derivatives (source terms) */
);


/** @} */ // end of group


#ifdef __cplusplus
}
#endif

#endif
