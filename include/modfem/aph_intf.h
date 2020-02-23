/************************************************************************
File aph_intf.h - interface with approximation modules of the code

Contains declarations of constants and interface routines:
  apr_module_introduce - to return the approximation method's name
  apr_get_mesh_id - to return the ID of the associated mesh
  apr_get_mixed_second_field_id - to return ID of second field in mixed approximation
  apr_get_nreq - to return the number of components in solution vector
  apr_get_mixed_nreq - to return the number of components in solution vector for mixed field
  apr_get_nr_sol - to return the number of solution vectors stored
  apr_get_base_type - to return the type of basis functions
  apr_get_ent_pdeg - to return the degree of approximation symbol
					  associated with a given mesh entity
  apr_set_ent_pdeg - to set the degree of approximation index
					  associated with a given mesh entity
  apr_get_ent_numshap - to return the number of shape functions (vector DOFs)
					associated with a given mesh entity
  apr_get_ent_nrdofs - to return the number of dofs associated with
					  a given mesh entity
  apr_get_el_pdeg - to return the degree of approximation vector
					  associated with a given element
  apr_set_el_pdeg - to set the degree of approximation vector
					  associated with a given element
  apr_get_mixed_pdeg - to get approximation pdeg associated with 
                       first or second field in mixed approximation
  apr_get_el_geo_pdeg - to return the degree of approximation
                                          associated with geometric transformation of element
  apr_set_el_geo_pdeg - to set the degree of approximation
                                          associated with geometric transformation of element
  apr_move_node - to move real geometric node by vector, node is related to vertex / edge / face or element
  apr_node_coor - to return node cordinates related to vertices / edges or faces
  apr_get_el_pdeg_numshap - to return the number of shape functions
					(scalar DOFs) for an element given its
				degree of approximation symbol or vector pdeg
  apr_get_el_dofs - to return the number and the list of element's degrees
			of freedom (corresponding to standard shape functions)
  apr_get_el_geo_dofs - to get geometric coordinates for element
  apr_get_nr_dof_ents - to return a global number of dof_entities
  apr_get_nrdofs_glob - to return a global dimension of the problem
  apr_read_ent_dofs - to read a vector of dofs associated with a given
			  mesh entity from approximation field data structure
  apr_write_ent_dofs - to write a vector of dofs associated with a given
			   mesh entity to approximation field data structure
  apr_create_ent_dofs - to write a vector of dofs associated with a given
			   mesh entity to approximation field data structure
  apr_free_field - to free approximation field data structure

Constrained approximation handling:
  apr_create_constr_data - to create and fill constraint arrays
				   (assuming approximation and mesh data agree)
  apr_get_el_constr_data - to return the number and the list of element's real
			   nodes, with the corresponding constraint coefficients
  apr_get_constr_data - to return constraints data for a node (dof entity)

Numerical integration support
  apr_prepare_integration_parameters - used e.g. by apr_num_int_el
  apr_num_int_el - to perform numerical integration for an element
  apr_get_stiff_mat_data - to return data on dof entities for an element and
		  to compute or rewrite element's stiffness matrix and RHSV

Adaptivity supporting functions
  apr_limit_deref - to return whether derefinement is allowed or not
  apr_limit_ref - to return whether refinement is allowed or not
	the routines does not use approximation data structures but the result
	depends on the approximation method
  apr_refine - to refine an element or the whole mesh checking mesh irregularity
  apr_derefine - to derefine an element or the whole mesh with irregularity check

Utilities:
  apr_shape_fun_3D - to compute values of shape functions and their
	local derivatives at a point within the master 3D element
  apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
  apr_set_quadr_3D - to prepare quadrature data for a given element
  apr_set_quadr_2D - to prepare quadrature data for a given face
  apr_set_quadr_2D_penalty - to prepare penalty quadrature data for a given face
  apr_L2_proj - to L2 project a function onto an element
  apr_get_profile - to return profile of the solution along a specified line


Global (O(N)) functions - very fast, but can be (should be?) parallelized...
  apr_init_field - to initiate new approximation field and read its data
  apr_write_field - to dump-out field data in the standard MOD_FEM format
  apr_check_field - to check approximation field data structure
  apr_set_ini_con - to set an initial condition
  apr_proj_dof_ref - to rewrite dofs after modifying the mesh
  apr_rewr_sol - to rewrite solution from one vector to another
  apr_sol_xglob - to return the solution at a point with global
	coordinates specified

apr_create_constr_data - O(N) also!

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#ifndef _aph_intf_
#define _aph_intf_

#include <stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

/** @defgroup APR Approximation
 *
 *  @{
 */


/*** CONSTANTS ***/
#define APC_TRUE 1
#define APC_FALSE 0

/** maximal degree of polynomial for an element */
#define APC_MAXELP_COMP 5
#define APC_MAXELP_TENS 505 // 6*6*7/2 = 126 shape functions (for DG only scalar)

/** maximal number of vector dofs (shape functions) for an element */
#define APC_MAXELVD 130 // was: 300, for degree 707 

/** maximal number of scalar dofs (nreq*num_shap) e.g. 4*15=60 for quadratic prisms and ns_supg */
#define APC_MAXELSD 300 // was: 300, for degree 707
#define APC_DOUBLE_MAXELSD 600 // for satisfy plast_flow module in some situation
  
/** maximal number of geo dofs */
#define APC_MAX_GEO_DOFS 18 // maximum is 18 dofs

/*** CONSTANTS ***/
#define APC_CUR_FIELD_ID  0 /** indicator for the current active field */

/** degree of approximation symbol indicating not allocated dofs */
#define APC_NO_DOFS  -1

/** Basis functions types */
#define APC_BASE_TENSOR_DG  1
#define APC_BASE_COMPLETE_DG  2
#define APC_BASE_PRISM_STD  3   /** for linear prismatic elements */
#define APC_BASE_PRISM_QUAD_STD 32 /** for quadratic prismatic elements [element type 3 / approximation degree 2] */
#define APC_BASE_TETRA_STD  4  /** for linear tetrahedral elements */
#define APC_BASE_TETRA_QUAD_STD  42  /** for quadratic tetrahedral elements [element type 4 / approximation degree 2] */

/** Initialization options */
#define APC_ZERO  0
#define APC_READ  1
#define APC_INIT  2

/** Identifiers of mesh entities */
// SHOULD BE CHOSEN AS SUBSEQUENT INTEGERS SINCE USED FOR ARRAY INDEXING !!!
#define APC_ELEMENT  4 /** ELEMENT type ID */ 
#define APC_FACE     3 /** FACE type ID */ 
#define APC_EDGE     2 /** EDGE type ID */ 
#define APC_VERTEX   1 /** VERTEX type ID */ 
  
#define APC_MIXED_ELEMENT  8 /** Mixed approximation ELEMENT type ID in second field */
#define APC_MIXED_FACE     7 /** Mixed approximation FACE type ID in second field */
#define APC_MIXED_EDGE     6 /** Mixed approximation EDGE type ID in second field */
#define APC_MIXED_VERTEX   5 /** Mixed approximation VERTEX type ID in second field */

// NO SECOND FIELD
#define APC_MIXED_NO_SEC_FIELD -1 /** No second field ID */

/** Options for assembling - should not be changed !!!*/
// Exact values are used in apr_get_stiff_mat_data
#define APC_NO_COMP    0 /** do not compute stiff matrix and rhs vector */
#define APC_COMP_SM    1 /** compute entries to stiff matrix only */
#define APC_COMP_RHS   2 /** compute entries to rhs vector only */
#define APC_COMP_BOTH  3 /** compute entries for sm and rhsv */
#define APC_COMP_MM  10 /** compute entries for sm and rhsv */
#define APC_REWR_SM    4 /** rewrite entries to stiff matrix only */
#define APC_REWR_RHS   5 /** rewrite entries to rhs vector only */
#define APC_REWR_BOTH  6 /** rewrite entries for sm and rhsv */
#define APC_REWR_MM  13 /** rewrite entries for sm and rhsv */

/** symbols for defining whether refinements/derefinements should be allowed */
/** e.g. due to mesh regularity constraints */
#define APC_REF_ALLOWED 1
#define APC_REF_DENIED  2
#define APC_DEREF_ALLOWED 3
#define APC_DEREF_DENIED  4

/** approximation field type */
#define APC_LINEAR_APPROXIMATION_PDEG 1 /** linear approximation PDEG */
#define APC_QUADRATIC_HIERACHICAL_APPROXIMATION_PDEG 15 /** quadratic hierachical approximation PDEG */
#define APC_QUADRATIC_APPROXIMATION_PDEG 2 /** quadratic Lagrange'a approximation PDEG */
#define APC_MIXED_P1_P1_APPROXIMATION_PDEG 11 /** mixed approximation PDEG (pure mixed only for testing), primary field LINEAR ; second field LINEAR approximation */
#define APC_MIXED_P2_P1_APPROXIMATION_PDEG 21 /** mixed approximation PDEG, primary field QUADRATIC ; second field LINEAR approximation */
#define APC_MIXED_P2_P2_APPROXIMATION_PDEG 22 /** mixed approximation PDEG (pure mixed only for testing), primary field QUADRATIC ; second field QUADRATIC approximation */

/** unified approximation pdeg [geo + base] */
#define APC_P1_NONE_G1 101 /** PDEG values used in problem config files ; LINEAR base approximation and LINEAR approximation of geometry */
#define APC_P15_NONE_G1 151 /** PDEG values used in problem config files ; QUADRATIC HIERARCHICAL base approximation and LINEAR approximation of geometry */
#define APC_P2_NONE_G1 201 /** PDEG values used in problem config files ; QUADRATIC LAGRANGE'A base approximation and LINEAR approximation of geometry */
#define APC_P2_NONE_G2 202 /** PDEG values used in problem config files ; QUADRATIC LAGRANGE'A base approximation and QUADRATIC approximation of geometry */
#define APC_P1_P1_G1 111 /** PDEG values used in problem config files ; MIXED (pure mixed only for testing) LINEAR base approximation for PRIMARY and SECOND fields, LINEAR approximation of geometry  */
#define APC_P2_P2_G2 222 /** PDEG values used in problem config files ; MIXED (pure mixed only for testing) QUADRATIC base approximation for PRIMARY and SECOND fields, QUADRATIC approximation of geometry */
#define APC_P2_P1_G2 212 /** PDEG values used in problem config files ; MIXED QUADRATIC base approximation for PRIMARY field and LINEAR base approximation for SECOND field, QUADRATIC approximation of geometry */
#define APC_P2_P1_G1 211 /** PDEG values used in problem config files ; MIXED QUADRATIC base approximation for PRIMARY field and LINEAR base approximation for SECOND field, LINEAR approximation of geometry */
  
/** aproximation number of shape function for element type */
#define APC_TETRA_NUM_SHAP_LINEAR_APROXIMATION 4 /** Number of shape functions for LINEAR tetrahedral element */
#define APC_TETRA_NUM_SHAP_QUADRATIC_APROXIMATION 10 /** Number of shape functions for QUADRATIC tetrahedral element */
  
#define APC_PRISM_NUM_SHAP_LINEAR_APROXIMATION 6 /** Number of shape functions for LINEAR prismatic element */
#define APC_PRISM_NUM_SHAP_QUADRATIC_APROXIMATION 18 /** Number of shape functions for QUADRATIC LAGRANGE'a prismatic element */
#define APC_PRISM_NUM_SHAP_HIERARCHICAL_QUADRATIC_APROXIMATION 15 /** Number of shape functions for QUADRATIC HIERARCHICAL prismatic element */

#define APC_MAX_CONSTR_EL_LOC 10 /** Maximal size of local constraint vector for any dof */
  
/// Close proximity precision
/// /* how close to the boundary is on the boundary */
#define APC_CLOSE 1.0e-6
#define APC_SMALL 1.0e-9

/*** FUNCTION DECLARATIONS - headers for external functions ***/

/**-----------------------------------------------------------
  apr_module_introduce - to return the approximation method's name
------------------------------------------------------------*/
extern int apr_module_introduce(
		    /** returns: >=0 - success code, <0 - error code */
  char* Approx_name /** out: the name of the approximation method */
  );

/**-----------------------------------------------------------
  apr_init_field - to initiate new approximation field and
				   read its control parameters
------------------------------------------------------------*/
extern int apr_init_field(  /** returns: >0 - field ID, <0 - error code */
  char Field_type, /** options: s - standard */ 
		   /**          c - discontinuous with complete basis */
                   /**          t - discontinuous with tensor product basis */
  int Control,	 /** in: control variable: */
		 /**      APC_INIT - to initialize the field to zero */
		 /**      APC_READ - to read field values from the file */
		 /**      APC_INIT - to initialize the field using function */
		 /**                 provided by the problem dependent module */
  int Mesh_id,	 /** in: ID of the corresponding mesh */
  int Nreq,	 /** in: number of equations - solution vector components */
  int Nr_sol,	 /** in: number of solution vectors for each dof entity */
  int Pdeg_in,   /** in: degree of approximating polynomial and geometry in format */
                 /**     [psg] p - primary field degree / s - secondary field degree for mixe default 0 / g - geometry */
  //int GeoPdeg_in,   /* in: degree of approximating polynomial for geometry */
  char *Filename, /** in: name of the file to read approximation data */
  double (*Fun_p)(int, double*, int) /** pointer to function that provides */
		 /** problem dependent initial condition data */
  );

/**--------------------------------------------------------
  apr_write_field - to dump-out field data in the standard MOD_FEM format
---------------------------------------------------------*/
extern int apr_write_field( /** returns: >=0 - success code, <0 - error code */
  int Field_id,    /** in: field ID */
  int Nreq,        /** in: number of equations (scalar dofs) */
  int Select,      /** in: parameter to select written vectors */
  double Accuracy, /** in: parameter specyfying accuracy */
		   /** (put 0 for full accuracy in "%g" format) */  
  char *Filename   /** in: name of the file to write field data */
  );

/**-----------------------------------------------------------
  apr_check_field - to check approximation field data structure
------------------------------------------------------------*/
extern int apr_check_field(
  int Field_id    /** in: approximation field ID  */
  );

/**-----------------------------------------------------------
  apr_get_mesh_id - to return the ID of the associated mesh
------------------------------------------------------------*/
extern int apr_get_mesh_id( /** returns: >0 - ID of the associated mesh,
					<0 - error code */
  int Field_id     /** in: approximation field ID  */
  );

/**-----------------------------------------------------------
apr_get_mixed_second_field_id - to return ID of second field in mixed approximation
------------------------------------------------------------*/
extern int apr_get_mixed_second_field_id( /** returns: >0 - ID of the associated second field id in mixed,
					      <0 - error code */
  int Field_id     /** in: approximation primary field ID  */
  );

/**-----------------------------------------------------------
  apr_get_nreq - to return the number of components in solution vector
------------------------------------------------------------*/
int apr_get_nreq( /** returns: >0 - number of solution components,
			      <0 - error code */
  int Field_id     /** in: approximation field ID  */
  );

/*------------------------------------------------------------
  apr_get_mixed_nreq - to return the number of components in solution vector for mixed field
------------------------------------------------------------*/
int apr_get_mixed_nreq( /* returns: >0 - number of solution components,
                                        <0 - error code */
  int Field_id     /* in: approximation field ID  */
  );

/**-----------------------------------------------------------
  apr_get_nr_sol - to return the number of solution vectors stored
------------------------------------------------------------*/
extern int apr_get_nr_sol(  /** returns: >0 - number of solution vectors,
					<0 - error code */
  int Field_id     /** in: approximation field ID  */
  );


/**--------------------------------------------------------
  apr_get_el_pdeg_numshap - to return the number of shape functions
	(scalar DOFs) for an element given its
	degree of approximation symbol or vector pdeg
---------------------------------------------------------*/
extern int apr_get_el_pdeg_numshap(
		 /** returns: >=0 - success code, <0 - error code*/
  int Field_id,  /** in: field ID */
  int El_id,        /** in: element ID */
  int *Pdeg_vec       /** in: degree of approximation symbol or vector */
  );


/*---------------------------------------------------------
apr_get_el_pdeg_numshap_related_with_dof_type - return number of shape functions retated with 
                                                  dof type 
---------------------------------------------------------*/
extern int apr_get_el_pdeg_numshap_related_with_dof_type(
  		 /* returns: >=0 - success code, <0 - error code*/
  int Field_id,  /* in: field ID */
  int El_id,     /* in: element ID */
  int Dof_type,  /* in: dof type */
  int *Pdeg      /* in: degree of approximation symbol or vector */
);


/**--------------------------------------------------------
  apr_get_base_type - to return the type of basis functions for an element
REMARK: type of basis functions differentiates element types as well 
   examples for discontinuous Galerkin approximation (from include/aph_intf.h):
          #define APC_BASE_TENSOR_DG  1
          #define APC_BASE_COMPLETE_DG  2  
   examples for standard linear approximation (from include/aph_intf.h): 
          #define APC_BASE_PRISM_STD  3   // for linear prismatic elements
          #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements
---------------------------------------------------------*/
extern int apr_get_base_type(/** returns: >0 - type of basis functions,
					 <0 - error code */
  int Field_id,  /** in: field ID */
  int El_id /** in: element ID */
  );

/**-----------------------------------------------------------
  apr_get_ent_pdeg - to return the degree of approximation symbol
					  associated with a given mesh entity
------------------------------------------------------------*/
extern int apr_get_ent_pdeg( /** returns: >0 - approximation index for P1 and P2,
                                           0 - dof entity inactive (constrained) or aproximation index in P0 
                                          <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id         /** in: mesh entity ID */
  );

/**-----------------------------------------------------------
  apr_set_ent_pdeg - to set the degree of approximation index
			  associated with a given mesh entity
------------------------------------------------------------*/
extern int apr_set_ent_pdeg( /** returns: >0 - success code,
					 <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id,         /** in: mesh entity ID */
  int Pdeg          /** in: degree of approximation */
  );

/**-----------------------------------------------------------
  apr_get_el_pdeg - to return the degree of approximation vector
					  associated with a given element
------------------------------------------------------------*/
extern int apr_get_el_pdeg( /** returns: >0 - success code
					<0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int El_id,         /** in: element ID */
  int *Pdeg_vec       /** out: degree of approximation symbol or vector */
  );

/**-----------------------------------------------------------
  apr_set_el_pdeg - to set the degree of approximation vector
					  associated with a given element
------------------------------------------------------------*/
extern int apr_set_el_pdeg( /** returns: >0 - success code,
					<0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int El_id,         /** in: element ID */
  int *Pdeg       /** out: degree of approximation symbol */
  );

/**-----------------------------------------------------------
  apr_get_mixed_pdeg - to get approximation pdeg associated with 
                       first or second field in mixed approximation
------------------------------------------------------------*/
extern int apr_get_mixed_pdeg( /** returns: >0 - success code,
					<0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int *Pdeg_vec       /** out: degree of approximation symbol */
  );

/**-----------------------------------------------------------
  apr_get_el_geo_pdeg - to return the degree of approximation
                                          associated with geometric transformation of element
------------------------------------------------------------*/
extern int apr_get_el_geo_pdeg( /** returns: >0 - success code,
					<0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int El_id,         /** in: element ID */
  int *Pdeg_vec       /** out: degree of approximation symbol */
  );
  
/**-----------------------------------------------------------
  apr_set_el_geo_pdeg - to set the degree of approximation
                                          associated with geometric transformation of element
------------------------------------------------------------*/
extern int apr_set_el_geo_pdeg( /** returns: >0 - success code,
					<0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int El_id,         /** in: element ID */
  int *Pdeg_vec       /** in: degree of approximation symbol or vector */
  );

/**-----------------------------------------------------------
 apr_move_node - to move real geometric node by vector, 
                              node is related to vertex / edge or face
------------------------------------------------------------*/
extern int apr_move_node( /** returns: >0 - success code,
			               <0 - error code */
  int Field_id, /* in: field ID */
  int Node_id, /* in: node ID */
  int Node_type, /* in: node type - APC_VERTEX / APC_EDGE / APC_FACE */
  double *Shift_vec /* in: shift vector */
);
  
/**-----------------------------------------------------------
  apr_node_coor - to return node cordinates related to vertices / edges or faces
------------------------------------------------------------*/
extern int apr_node_coor( /** returns: >0 - success code,
			               <0 - error code */
 int Field_id, /* in: field ID */
 int Node_id, /* in: node ID */
 int Node_type, /* in: node type - APC_VERTEX / APC_EDGE / APC_FACE */
 double *coord /* out: cordinates of node */
);

/**-----------------------------------------------------------
  apr_get_el_dofs - to return the list of standard dofs for a given element
			each dof corresponds to a standard shape function
------------------------------------------------------------*/
extern int apr_get_el_dofs( /** returns: >0 - the number of dofs
					<0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int El_id,        /** in: element ID */
  int Vect_id,       /** in: vector ID in case of multiple solution vectors */
  double *El_dofs_std     /** out: the list of values of element dofs */
  );

/**-----------------------------------------------------------
apr_get_el_geo_dofs - to get geometric coordinates for a given element
------------------------------------------------------------*/
extern int apr_get_el_geo_dofs( /** returns: >0 - the number of dofs
				               <0 - error code */
  int Field_id, /** in: approximation field ID  */
  int El_id,    /** in: element ID */
  int *GeoDofs_id,    /** out: list of dofs in element */
		/*	 (Dofs[0] - number of dofs) */
  int *GeoDofs_type, /** out: type of dofs [APC_VERTEX / APC_EDGE / APC_FACE] */
                  /*     (Dofs_type[0] - geo_pdeg */
  double *GeoDofs /* out: coordinates of element dofs */
  );

/**-----------------------------------------------------------
  apr_get_ent_numshap - to return the number of shape functions (vector DOFs)
			associated with a given mesh entity
------------------------------------------------------------*/
extern int apr_get_ent_numshap( /** returns: >0 - the number of dofs,
					    <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id         /** in: mesh entity ID */
  );

/**-----------------------------------------------------------
  apr_get_ent_nrdofs - to return the number of dofs associated with
					  a given mesh entity
------------------------------------------------------------*/
extern int apr_get_ent_nrdofs( /** returns: >0 - the number of dofs,
					   <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id         /** in: mesh entity ID */
  );

/**--------------------------------------------------------
  apr_get_nrdofs_glob - to return a global dimension of the problem
---------------------------------------------------------*/
extern int apr_get_nrdofs_glob(	/** returns: global dimension of the problem */
  int Field_id    /** in: field ID */
  );

/**-----------------------------------------------------------
  apr_read_ent_dofs - to read a vector of dofs associated with a given
			  mesh entity from approximation field data structure
------------------------------------------------------------*/
extern int apr_read_ent_dofs(/** returns: >=0 - success code, <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id,        /** in: mesh entity ID */
  int Ent_nrdofs,     /** in: number of dofs associated with the entity */
  int Vect_id,       /** in: vector ID in case of multiple solution vectors */
  double* Vect_dofs  /** out: dofs read from data structure */
  );


/**-----------------------------------------------------------
  apr_write_ent_dofs - to write a vector of dofs associated with a given
			   mesh entity to approximation field data structure
------------------------------------------------------------*/
extern int apr_write_ent_dofs(/** returns: >=0 - success code, <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id,        /** in: mesh entity ID */
  int Ent_nrdofs,     /** in: number of dofs associated with the entity */
  int Vect_id,       /** in: vector ID in case of multiple solution vectors */
  double* Vect_dofs);


/**-----------------------------------------------------------
  apr_create_ent_dofs - to write a vector of dofs associated with a given
			   mesh entity to approximation field data structure
------------------------------------------------------------*/
extern int apr_create_ent_dofs(
			 /** returns: >=0 - success code, <0 - error code */
  int Field_id,      /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id,        /** in: mesh entity ID */
  int Ent_nrdofs,     /** in: number of dofs associated with the entity */
  int Vect_id,       /** in: vector ID in case of multiple solution vectors */
  double* Vect_dofs  /** in: dofs to be written */
  );

/**-----------------------------------------------------------
  apr_set_ini_con - to set an initial condition
------------------------------------------------------------*/
extern int apr_set_ini_con(/** returns: >=0 - success code, <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  double (*Fun_p)(int, double*, int) /** pointer to function that provides */
		 /** problem dependent initial condition data */
  );


/**-----------------------------------------------------------
  apr_prepare_integration_parameters - used e.g. by apr_num_int_el
------------------------------------------------------------*/
extern int apr_prepare_integration_parameters( 
  int Field_id,     /** in: approximation field ID  */
  int El_id,        /** in: element ID */
  int *Geo_order,   /** out: geometrical order of approximation */
  int *Num_geo_dofs,/** out: number of geometrical degrees of freedom */
  double* Geo_dofs, /** out: geometrical degrees of freedom */
  int *El_mate,     /** out: material index (in materials database) */
  int *Base,        /** out: type of basis functions */
  int *Pdeg_vec,    /** out: degree of approximation symbol or vector */
  int *Num_shap,    /** out: number of shape functions */
  int *Nreq,        /** out: number of components for unknowns */
  int *Num_dofs     /** out: number of scalar dofs (num_dofs = nreq*num_shap) */
  );


/**-----------------------------------------------------------
  apr_num_int_el - to perform numerical integration for an element
------------------------------------------------------------*/
extern int apr_num_int_el(
  int Problem_id,
  int Field_id,    /** in: approximation field ID  */
  int El_id,       /** in: unique identifier of the element */ 
  int Comp_sm,     /** in: indicator for the scope of computations: */
                   /**   APC_NO_COMP  - do not compute anything */
                   /**   APC_COMP_SM - compute entries to stiff matrix only */
                   /**   APC_COMP_RHS - compute entries to rhs vector only */
                   /**   APC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,        /** in: enforced degree of polynomial (if !=NULL ) */
  double *Sol_dofs_k,   /** in: solution dofs from previous iteration */
                        /**     (for nonlinear problems) */
  double *Sol_dofs_n,   /** in: solution dofs from previous time step */ 
			/**     (for nonlinear problems) */
  int *diagonal,   /** array of indicators whether coefficient matrices are diagonal */
  double *Stiff_mat,	/** out: stiffness matrix stored columnwise */
  double *Rhs_vect	/** out: rhs vector */
  );



/**-----------------------------------------------------------
  apr_get_stiff_mat_data - to return data on dof entities for an element and
		  to compute or rewrite element's stiffness matrix and RHSV
------------------------------------------------------------*/
extern int apr_get_stiff_mat_data(
  int Field_id,  /** in: approximation field ID  */
  int El_id,     /** in: unique identifier of the element */
  int Comp_sm,   /** in: indicator for the scope of computations: */
		 /**   APC_NO_COMP  - do not compute anything */
		 /**   APC_COMP_SM - compute entries to stiff matrix only */
		 /**   APC_COMP_RHS - compute entries to rhs vector only */
		 /**   APC_COMP_BOTH - compute entries for sm and rhsv */
		 /**   APC_REWR_SM - rewrite only entries to stiff matrix only */
		 /**   APC_REWR_RHS - rewrite only entries to rhs vector only */
		 /**   APC_REWR_BOTH - rewrite only entries for sm and rhsv */
  char Transpose,/** in: perform transposition while rewriting */
				  /**     'y' or 'Y' - yes, otherwise - no */
  int Pdeg_in,    /** in: enforced degree of polynomial (if > 0 ) */
  int Nreq_in,    /** in: enforced nreq (if > 0 ) */
  int* Nr_dof_ent,/** in: size of arrays, */
		/** out: no of filled entries, i.e. number of mesh entities*/
		/** with which dofs and stiffness matrix blocks are associated */
  int* List_dof_ent_type, /** out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_id,   /** out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_nrdofs,/** out: list of no of dofs for 'dof' entity */
  int* Nrdofs_loc,/** in(optional): size of Stiff_mat and Rhs_vect */
		/** out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /** out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect        /** out(optional): rhs vector */
				 );

/**-----------------------------------------------------------
  apr_proj_dof_ref - to rewrite dofs after modifying the mesh
------------------------------------------------------------*/
extern int apr_proj_dof_ref(
  int Field_id,    /** in: approximation field ID  */
  int El,	   /** in:  >0 - rewrite after one [de]refinement of el */
		   /**     <=0 - rewrite after massive [de]refinements */
  int Max_elem_id, /** in: maximal element (face, etc.) id */
  int Max_face_id, /**     before and after refinements */
  int Max_edge_id,
  int Max_vert_id
  );

/**-----------------------------------------------------------
  apr_rewr_sol - to rewrite solution from one vector to another
------------------------------------------------------------*/
extern int apr_rewr_sol(
  int Field_id,      /** in: data structure to be used  */
  int Sol_from,      /** in: ID of vector to read solution from */
  int Sol_to         /** in: ID of vector to write solution to */
  );

/**-----------------------------------------------------------
  apr_free_field - to free approximation field data structure
------------------------------------------------------------*/
extern int apr_free_field(
  int Field_id    /** in: approximation field ID  */
  );


/**-----------------------------------------------------------
  apr_limit_deref - to return whether derefinement is allowed or not
	the routine does not use approximation data structures but the result
	depends on the approximation method
------------------------------------------------------------*/
extern int apr_limit_deref(
  int Field_id,    /** in: approximation field ID  */
  int El_id      /** in: unique identifier of the element */
  );

/**-----------------------------------------------------------
  apr_limit_ref - to return whether refinement is allowed or not
	the routine does not use approximation data structures but the result
	depends on the approximation method
------------------------------------------------------------*/
extern int apr_limit_ref(
  int Field_id,    /** in: approximation field ID  */
  int El_id      /** in: unique identifier of the element */
  );

/**-----------------------------------------------------------
  apr_refine - to refine an element or the whole mesh checking mesh irregularity
------------------------------------------------------------*/
int apr_refine( /** returns: >=0 - success code, <0 - error code */
  int Field_id,  /** in: field ID */
  int El_id   /** in: element ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
		);

/**-----------------------------------------------------------
  apr_derefine - to derefine an element or the whole mesh with irregularity check
------------------------------------------------------------*/
int apr_derefine( /** returns: >=0 - success code, <0 - error code */
  int Field_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El_id  /** in: element ID or -2 (MMC_DO_UNI_DEREF) for uniform derefinement */
		  );

/****************** FUNCTIONS IN APS_...._UTIL.C ******************/

/**--------------------------------------------------------
apr_shape_fun_3D - to compute values of shape functions and their
	local derivatives at a point within the master 3D element
----------------------------------------------------------*/
int apr_shape_fun_3D( /** returns: the number of shape functions (<=0 - failure) */
	int Base_type,	   /** in: type of basis functions: */
	/** REMARK: type of basis functions differentiates element types as well */
	/**  #define APC_BASE_TENSOR_DG 1 - tensor product for prismatic elements and DG*/
	/**  #define APC_BASE_COMPLETE_DG 2 - complete polynomials for prismatic elements and DG*/
        /**  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements */
        /**  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements  */
	int Pdeg, 	   /** in: degree of polynomial - can be either */
			   /*	a single number, for isotropic p, */
			   /*	or a combination pdegy*10+pdegx */
	double *Eta,	   /** in: local coord of the considered point */
	double *Base_phi,  /** out: basis functions */
	double *Base_dphix,/** out: x derivative of basis functions */
	double *Base_dphiy,/** out: y derivative of basis functions */
	double *Base_dphiz /** out: z derivative of basis functions */
		      );


#define APC_ELEM_CALC_SHAPE_FUNC_N_VALUES 1
#define APC_ELEM_CALC_DERIVATIVES_AND_JACOBIAN 2
#define APC_ELEM_CALC_ON_FACE(face_number) 3+face_number
/**-----------------------------------------------------------------
apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
-------------------------------------------------------------------*/
extern double apr_elem_calc_3D(
	/** returns: Jacobian determinant at a point, either for */
	/** 	volume integration if Vec_norm==NULL,  */
	/** 	or for surface integration otherwise */
	int Control,	/** in: control parameter (what to compute): */
            /*	APC_ELEM_CALC_SHAPE_FUNC_N_VALUES 1  - shape functions and values */
            /*	APC_ELEM_CALC_DERIVATIVES_AND_JACOBIAN 2  - derivatives and jacobian */
            /** APC_ELEM_CALC_ON_FACE(face_number) >2 - computations on the  */
	    /*	                                       (Control-2)-th element's face */
	int Nreq,	    /** in: number of equations */
	int *Pdeg_vec,	    /** in: element degree of polynomial */
	int Base_type,	    /** in: type of basis functions: */
/**          #define APC_BASE_TENSOR_DG  1 */
/**          #define APC_BASE_COMPLETE_DG  2 */
/**          #define APC_BASE_PRISM_STD   3   // for prismatic elements */
/**          #define APC_BASE_PRISM_QUAD_STD 32 // for quadratic prismatic elements [element type 3 / approximation degree 2] */
/**          #define APC_BASE_TETRA_STD  4  // for tetrahedral elements */ 
/**          #define APC_BASE_TETRA_QUAD_STD  42  // for quadratic tetrahedral elements [element type 4 / approximation degree 2] */
	double *Eta,    /** in: local coordinates of the input point */
	double *Node_coor,  /** in: array of coordinates of vertices of element*/
	double *Sol_dofs,   /** in: array of element' dofs */
	double *Base_phi,   /** out: basis functions */
	double *Base_dphix, /** out: x-derivatives of basis functions */
	double *Base_dphiy, /** out: y-derivatives of basis functions */
	double *Base_dphiz, /** out: z-derivatives of basis functions */
	double *Xcoor,	    /** out: global coordinates of the point*/
	double *Sol,        /** out: solution at the point */
	double *Dsolx,      /** out: derivatives of solution at the point */
	double *Dsoly,      /** out: derivatives of solution at the point */
	double *Dsolz,      /** out: derivatives of solution at the point */
	double *Vec_nor     /** out: outward unit vector normal to the face */
	);

/**--------------------------------------------------------
apr_set_quadr_3D - to prepare quadrature data for a given element
---------------------------------------------------------*/
extern int apr_set_quadr_3D(
	int Base_type,	   /** in: type of basis functions: */
			   /** 	1 (APC_TENSOR) - tensor product */
			   /** 	2 (APC_COMPLETE) - complete polynomials */
	int *Pdeg_vec,	/** in: element degree of polynomial */
	int *Ngauss,	/** out: number of gaussian points */
	double *Xg,	/** out: coordinates of gaussian points */
	double *Wg	/** out: weights associated with points */
	);

/**--------------------------------------------------------
apr_set_quadr_2D - to prepare quadrature data for a given face
---------------------------------------------------------*/
extern int apr_set_quadr_2D(
	int Fa_type,	/** in: type of a face */
	int Base_type,	   /** in: type of basis functions: */
			   /** 	1 (APC_TENSOR) - tensor product */
			   /** 	2 (APC_COMPLETE) - complete polynomials */
	int *Pdeg_vec,	/** in: element degree of polynomial */
	int *Ngauss,	/** out: number of gaussian points */
	double *Xg,	/** out: coordinates of gaussian points */
	double *Wg	/** out: weights associated with points */
	);

/**--------------------------------------------------------
apr_set_quadr_2D_penalty - to prepare penalty quadrature data for a given face
---------------------------------------------------------*/
extern int apr_set_quadr_2D_penalty(
	int Fa_type,	/** in: type of a face */
	int Base_type,	   /** in: type of basis functions: */
			   /**   (APC_TENSOR) - tensor product */
			   /** 	(APC_COMPLETE) - complete polynomials */
	int *Pdeg_vec,	/** in: element degree of polynomial */
	int *Ngauss,	/** out: number of gaussian points */
	double *Xg,	/** out: coordinates of gaussian points */
	double *Wg	/** out: weights associated with points */
	);

/**--------------------------------------------------------
apr_L2_proj - to L2 project a function onto an element
---------------------------------------------------------*/
extern int apr_L2_proj(	/** returns: >0 - success, <=0 - failure */
		int Field_id,   /** in: field ID */
	int Mode,	/** in: mode of operation */
			/**    <-1 - projection from ancestor to father */
			/**          the value is the number of ancestors */
			/**     -1 - projection from father to son */
			/**     >0 - projection of function, the value */
			/**          is a flag for routine returning */
			/**          function value at point */
	int El,		/** in: element number */
	int *Pdeg_vec,	/** in: element degree of approximation */
	double* Dofs,	/** out: workspace for degress of freedom of El */
			/** 	NULL - write to  data structure */
	int *El_from,	/** in: list of elements to provide function */
	int *Pdeg_vec_from,	/** in: degree of polynomial for each El_from */
	/** if element Pdeg is a vector - its components must be suitably */
	/** placed in Pdeg_vec_from */
	double* Dofs_from, /** in: Dofs of El_from or...*/
	double (*Fun_p)(double*,double*,double*,double*)   /** in: pointer to */
			   /** function with field values and its derivatives */
	);

/// Flags controling apr_sol_xglob
enum ape_sol_xglob{
    /// the defaults are:
    /// not APE_SOL_XGLOB_MATCH_ALL_ELEMENTS
    /// = match first element only
    /// not APE_SOL_XGLOB_MATCH_WITH_ADAPTATION
    /// = match only at initial level
    /// not APE_SOL_XGLOB_CHECK_ONLY_GIVEN_ELEMENT
    /// = search for matching element(s)
    APE_SOL_XGLOB_DEFAULT = 0x00,
    APE_SOL_XGLOB_MATCH_ALL_ELEMENTS = 0x01,
    APE_SOL_XGLOB_MATCH_WITH_ADAPTATION = 0x02,
    /// if the element containing given point is already known, one
    /// can put its' id into El[0]
    /// and only this element will be considered
    APE_SOL_XGLOB_CHECK_ONLY_GIVEN_ELEMENT = 0x04

};

/**--------------------------------------------------------
apr_sol_xglob - to return the solution at a point with global
	coordinates specified. The procedure finds,
	for a given global point, the respective local coordinates
	within the proper initial mesh element, than computes
	corresponding local coordinates within an active ancestor of
	the initial element and finally finds the value of solution.
	There may be several initial mesh elements, several ancestors
	and several values due to the discontinuity of approximate solution.
---------------------------------------------------------*/
extern int apr_sol_xglob(/** returns: >0 - success, <=0 - failure */
  int Field_id,   /** in: field ID */
  double *Xglob,	/** in: global coordinates of a point */
  int Nb_sol,    /** in: which solution to take: 1 - sol_1, 2 - sol_2 */
  int* El,	/** out: list of element numbers,  */
                /**      El[0] - number of elements on the list */
  double* Xloc,	/** out: list of local coordinates within elements */
  double *Sol,	/** out: list of solutions at the point */
  double *Dxsol,  /** out: list of derivatives wrt x of solution */
  double *Dysol,  /** out: list of derivatives wrt y of solution */
  double *Dzsol,   /** out: list of derivatives wrt z of solution */
  double Close_proximity_precision, //< default to APC_CLOSE; how far from element point is considered to be within an element
  int Sol_xglob_flags //< see ape_sol_xlgob
);


/**-----------------------------------------------------------
  apr_create_constr_data - to create and fill constraint arrays
				   (assuming approximation and mesh data agree)
------------------------------------------------------------*/
extern int apr_create_constr_data(
		  /** returns: >=0 - success code, <0 - error code */
  int Field_id    /** in: approximation field ID  */
  );



/**-----------------------------------------------------------
  apr_get_el_constr_data - to return the number and the list of element's real
			   nodes, with the corresponding constraint coefficients
------------------------------------------------------------*/
extern int apr_get_el_constr_data(
		   /** returns: >=0 - success code, <0 - error code */
  int Field_id,   /** in: approximation field ID  */
  int El_id,        /** in: element ID */
  int* Nodes,  	  /** out: list of vertex node IDs */
			  /*	(Nodes[0] - number of nodes) */
  int* Nr_constr, /** out: list with the numbers of constraints for each vertex*/
  int* Constr_id, /** out: ID's of parent (constraining) nodes */
  int* Constr_type,   /** out: Type of parent (used in quadratic approximation and above) */
  double* Constr_val  /** out: the corresponding constraint coefficients */
  );

/**-----------------------------------------------------------
  apr_get_constr_data - to return constraints data for a node (dof entity)
------------------------------------------------------------*/
extern int apr_get_constr_data(/** returns: >=0 - success code, <0 - error code */
  int Field_id,     /** in: approximation field ID  */
  int Node_id,      /** in: id of node (can be vertex or edge) */
  int Node_type,    /** in: type of input node APC_VERTEX or APC_EDGE */

  int *Constr,      /** out: table of constraints data, 
		     *      Constr[0] - number of constraints 
		     *                  (in linear approximation
		     *                  2 - mid-edge node, 4 - mid-side node)
		     */
  int *Constr_type  /** out: table of constraints element type,
		     *      Constr_type[0] - id of coefficient vector
		     */
  );


/// apr_get_profile - to return values of solution at given points
//						to generate plots and diagrams

extern int apr_get_profile( // returns: >0 - success, <=0 - failure
	FILE* filePtr,	//in: pointer to file, where profile will be printed
	int fieldId,	//in: field ID
	int solNr,	//in: solution number
	int nSol,	//in: solution length (no. of components)
	double * pt1,	//in: begin point coords
	double * pt2,	//in: end point coords
	int nPoints	//in: number of points (including end points)
	);
    
/** @} */ // end of group

#ifdef __cplusplus
}
#endif

#endif
