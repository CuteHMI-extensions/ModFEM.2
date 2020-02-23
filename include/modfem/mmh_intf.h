/************************************************************************
File mmh_intf.h - interface with mesh manipulation module of the code.

Apart from some refinement functions (and possibly some special functions, 
e.g. for remesh module) mesh module functions are local (O(1)) and do not
have to be parallelized (they are even not assumed to be thread safe!!!) 

Contains declarations of constants and interface routines:
  mmr_module_introduce - to return the mesh name
  mmr_init_mesh - to initialize the mesh data structure
  mmr_export_mesh - to export (write to a file) mesh data in a specified format
  mmr_test_mesh - to test the integrity of mesh data
  OBSOLETE !!! mmr_get_mesh_i_params - to return mesh parameters
  mmr_get_nr_elem - to return the number of active elements
  mmr_get_max_elem_id - to return the maximal element id
  mmr_get_next_act_elem - to return the next active element id
  mmr_get_next_elem_all - to return the next element id (active or inactive)
  mmr_get_nr_face - to return the number of active faces
  mmr_get_max_face_id - to return the maximal face id
  mmr_get_next_act_face - to return the next active face id
  mmr_get_next_face_all - to return the next face id (active or inactive)
  mmr_get_nr_edge - to return the number of active edges
  mmr_get_max_edge_id - to return the maximal edge id
  mmr_get_next_edge_all - to return the next edgeent id (active or inactive)
  mmr_get_nr_node - to return the number of active nodes
  mmr_get_max_node_id - to return the maximal node id
  mmr_get_next_node_all - to return the next node id (active or inactive)
  mmr_get_max_gen - to get maximal allowed generation level for elements
  mmr_set_max_gen - to set maximal allowed generation level for elements
  mmr_get_max_gen_diff - to get maximal allowed generation difference between
                         neighboring elements
  mmr_set_max_gen_diff - to set maximal allowed generation difference between
                         neighboring elements

  mmr_get_max_elem_max - to return the maximal possible element id
  mmr_get_max_face_max - to return the maximal possible face id
  mmr_get_max_edge_max - to return the maximal possible edge id
  mmr_get_max_node_max - to return the maximal possible vertex id

  mmr_init_ref - to initialize the process of refinement
  mmr_refine_el - to refine an element  WITHOUT IRREGULARITY CHECK
  mmr_r_refine - to r-refine an element with given boundary conditions
  mmr_derefine_el - to derefine an element  WITHOUT IRREGULARITY CHECK
  mmr_refine_mesh - to refine the WHOLE mesh WITHOUT IRREGULARITY CHECK
  mmr_derefine_mesh - to derefine the WHOLE mesh WITHOUT IRREGULARITY CHECK
  mmr_final_ref - to finalize the process of refinement after rewriting DOFs
  mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready for dofs projection
  mmr_free_mesh - to free space allocated for mesh data structure

  mmr_gen_boundary_layer - to generate boundary layer (for tetrahedral meshes only)
  mmr_elem_structure - to return elem structure (e.g. for sending)
  mmr_el_status - to return element status (active, inactive, free space)
  mmr_el_type - to return element type
  mmr_el_groupID - to return material ID for element
  mmr_el_set_mate - to set material number for element
  mmr_el_type_ref - to return element's type of refinement
  mmr_el_faces - to get faces and big neighbors of an element
  mmr_el_edges - to get element's edges
  mmr_el_node_coor - to get the coordinates of element's nodes
  mmr_el_fam - to return family information for an element
  mmr_el_fam_all - to return all recursive family information for an element
  mmr_el_gen - to return generation level for an element
  mmr_el_ancestor - to find the ancestor of an element
  mmr_el_eq_neig - to get equal size (or larger) neighbors of an element
  mmr_el_hsize - to compute a characteristic linear size for an element
  mmr_set_elem_fath - to set family data for an elem
  mmr_set_elem_fam - to set family data for an elem
  mmr_el_fa_nodes - to get list of local indexes for face nodes in elem

  mmr_face_structure - to return face structure (e.g. for sending)
  mmr_fa_status - to return face status (active, inactive, free space)
  mmr_fa_type - to return face type (triangle, quad, free space)
  mmr_fa_bc - to get the boundary condition flag for a face
  mmr_fa_sub_bnd - to indicate face is on the subdomain boundary
  mmr_fa_set_sub_bnd - to indicate face is on the subdomain boundary
  mmr_fa_edges - to return a list of face's edges
  mmr_fa_eq_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (only equal size
		neighbors considered)
  mmr_fa_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (for active
		faces all neighbors are active)
  mmr_fa_node_coor - to get the list and coordinates of faces's nodes
  mmr_fa_elem_coor - to find coordinates within neighboring
		elements for a point on face
  mmr_fa_area - to compute the area of face and vector normal
  mmr_fa_fam - to return face's family information
  mmr_set_face_fam - to set family data for an face
  mmr_set_face_neig - to set neighbors data for a face

  mmr_edge_status - to return edge status (active, inactive, free space)
  mmr_edge_nodes - to return edge node's IDs
  mmr_edge_sons - to return edge son's numbers
  mmr_edge_elems - to return a list of elements containing the edge
  mmr_edge_structure - to return edge structure (e.g. for sending)
  mmr_set_edge_type - to set type for an edge
  mmr_set_edge_fam - to set family data for an edge

  mmr_node_status - to return node status (active, free space)
  mmr_node_coor - to return node coordinates
  mmr_set_node_coor - set node coordinates
  mmr_move_node - moves node by vector

  mmr_loc_loc - to compute local coordinates within an element,
	given local coordinates within an element of the same family

  mmr_create_element - to create an element and fill its data structure
  mmr_create_face - to create a new face
  mmr_create_edge - to create a new edge structure
  mmr_create_node - to create a node
  mmr_clust_face - to cluster back a family of faces
  mmr_clust_edge  - to cluster two edges back into their father

  mmr_divide_face - to divide a face into sons
  mmr_divide_edge - to divide an edge into two sons

  mmr_del_elem    - to free an element structure
  mmr_del_face    - to free a face structure
  mmr_del_edge    - to free an edge structure
  mmr_del_node    - to free a node structure


// The five next functions are not supported by mmd_prism, mmd_prism_2D and mmd_remesh modules!!!

  mmr_reserve - to inform mesh module about new target number of elems,faces,edges and nodes
  mmr_add_elem - to add (append) new element
  mmr_add_face - to add (append) new face
  mmr_add_edge - to add (append) new edge
  mmr_add_node - to add (append) new node

// Groups and contact BCs handling
extern int mmr_groups_number(const int Mesh_id);
extern int mmr_groups_ids(const int Mesh_id, int *tab, int l_tab);
extern int mmr_split_into_blocks_add_contact(const char *workdir,const int Mesh_id, int *tabMat,int l_mat,int *tabBlock,int l_block,int *war,int l_war,double *tempBlock,int l_tempBlock,int *groups_fluid,int l_groups_fluid);
int mmr_get_fa_el_bc_connect(int face_id,int &el_id) - return 0 if normal bc else return connected face_id


// Below are functions used by modules intarfacing with REMESH module ONLY !!!!

  mmr_init_dist2bound - to initialize and pre-compute resources for mmr_get_el_dist2bound
  mmr_get_el_dist2boundary - to get distance from point inside element to nearest boundary

  mmr_get_coor_from_element_motion(int idEl,int idLP,double *coor,int flagaSiatki);

  mmr_copyMESH(int flaga);
  mmr_create_mesh_Cube(const char *nazwa,int node_x,int node_y,int node_z,double size_x,double size_y,double size_z,double divide,int *warunki);
  mmr_init_all_change(int a);
  mmr_test_mesh_motion(int ileWarstw,int obecny_krok, int od_krok,int ileKrok,double minPoprawy,double px0,double py0,double pz0,double px1,double py1,double pz1,double endX,double endY,double endZ);
  mmr_test_weldpool(double obecnyKrok, double krok_start,double minPoprawy,double doX,double doY,double dl,double zmiejszaPrzes,double limit,double szerokosc);


------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
	06.2010 - Kazimierz Michalik, hybrid mesh
    10.2012 - Kazimierz Michalik, modyfied I/O
*************************************************************************/

#ifndef _mmh_intf_
#define _mmh_intf_

#include <stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

/** @defgroup MM Mesh
 *
 *  @{
 */

/** Constants */
#define MMC_MAXELFAC        5   /** maximal number of faces of an element */
#define MMC_MAXFAVNO        4   /** maximal number of vertices of a face */
#define MMC_MAXELVNO        6   /** maximal number of vertices of an element */
#define MMC_MAXELSONS       8   /** maximal number of sons of an element */
#define MMC_MAXELEDGES      9   /** maximal number of edges in an element */

extern const int MMC_MAX_EDGE_ELEMS;/*max number of elements containg an edge*/

extern const int MMC_CUR_MESH_ID; /** indicator for the current active mesh */
extern const int MMC_INIT_GEN_LEVEL;   // value of 'initial mesh' generation level

/** Data file types for mesh I/O data */
typedef enum {
    MMC_MOD_FEM_MESH_DATA   = 0, /** generic mesh for mesh module */
    MMC_MOD_FEM_PRISM_DATA  = 1, /** prism mesh read from dump files */
    MMC_MOD_FEM_TETRA_DATA  = 2, /** tetra mesh read from dump files */
    MMC_MOD_FEM_HYBRID_DATA = 3, /** hybrid mesh read from dump files */
    MMC_GRADMESH_DATA       = 10, /** mesh produced by 2D GRADMESH generator */
    MMC_HP_FEM_MESH_DATA    = 1, /// mesh read from dump files
    MMC_NASTRAN_DATA        = (int)'n',
    MMC_NASTRAN_FREE_DATA   = MMC_NASTRAN_DATA,
    MMC_NASTRAN_SHORT_DATA  = (int)'s',
    MMC_NASTRAN_LONG_DATA   = (int)'l',
    MMC_BINARY_DATA         = (int)'b', /// for binary dump
    MMC_BOUND_VERTS_DATA    = (int)'m', /// 'm' like prof. Milenin :)
    MMC_IN_ANSYS_DATA       = (int)'i', /// .in files
    MMC_PARAVIEW_VTK_DATA   = (int)'v', /// Paraview .vtk files
    MMC_MSH_DATA            = (int)'f', /// ANSYS meshing, TGrid .msh files
    MMC_MAX_FILE_TYPE=126
} mmt_file_type;


extern const int MMC_AUTO_GENERATE_ID;

/** Position and neighbors options */
extern const int MMC_BOUNDARY;         /** boundary indicator */
extern const int MMC_BIG_NGB;          /** big neghbor indicator */
extern const int MMC_SUB_BND;         /** inter-subdomain boundary indicator */

/** Identifiers of mesh entities */
/* Identifiers of mesh entities */
typedef enum {
    MMC_NODE    = 0,
    MMC_EDGE    = 1,
    MMC_FACE    = 2,
    MMC_ELEMENT = 3,
    MMC_ALL_N_TYPES
} mmt_mesh_entity;

/** Types of mesh entities */
typedef enum {
    MMC_QUAD            = 4, /** quadrilateral element or face */
    MMC_TRIA            = 3, /** triangular element or face */
    MMC_TETRA           = 7, /** tetrahedral element */
    MMC_PRISM           = 5, /** prismatic element */
    MMC_BRICK           = 6, /** hexahedral element */
    MMC_ENTITY_UNKNOWN  = -1 /// for handling unimplemented cases
} mmt_mesh_entity_types;



/** Refinement options */
extern const int MMC_DO_UNI_REF;   /** to perform uniform refinement */
extern const int MMC_DO_UNI_DEREF; /** to perform uniform derefinement */

/** Status indicators */
extern const int MMC_ACTIVE;           /** active mesh entity */
extern const int MMC_INACTIVE;         /** inactive (refined) mesh entity */
extern const int MMC_FREE;             /** free space in data structure */

/** Monitoring options */
extern const int MMC_PRINT_NOT;     /** not to print anything */
extern const int MMC_PRINT_ERRORS;  /** to print error messages only */
extern const int MMC_PRINT_INFO;    /** to print most important information */
extern const int MMC_PRINT_ALLINFO; /** to print all available information */

/** Refinement types */
extern const int MMC_NOT_REF;          /** not refined */
extern const int MMC_REF_ISO;          /** isotropic refinement */
extern const int MMC_REF_ANI;          /** anisotropic refinement */

/** Other */
extern const int MMC_INIT_GEN_LEVEL;   // value of 'initial mesh' generation level
extern const int MMC_NO_FATH;       /** no father indicator */

extern const int MMC_FACE_NODES_FOR_TETRA[4][3];
extern const int MMC_FACE_NODES_FOR_PRISM[5][4];

extern const int MMC_IS_SUPPORTING_NEW_MPI;

/** Declarations of interface routines: */



/**-----------------------------------------------------------
  mmr_module_introduce - to return the mesh name
------------------------------------------------------------*/
extern int mmr_module_introduce(
                  /** returns: >=0 - success code, <0 - error code */
  char* Mesh_name /** out: the name of the mesh */
  );

/**-----------------------------------------------------------
  mmr_module_internals - to return the mesh internal structures
------------------------------------------------------------*/
extern void* mmr_module_internals(
                  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id/** out: the name of the mesh */
  );


/**-----------------------------------------------------------
  mmr_init_mesh - to initialize the mesh data structure (and read data)
------------------------------------------------------------*/
extern int mmr_init_mesh(  /** returns: >0 - Mesh ID, <0 - error code */
  int Control,	    /** in: control variable to choose data format */
  /** MMC_MOD_FEM_MESH_DATA = 1 - mesh read from dump files */
  /** MMC_GRADMESH_DATA    = 2 - mesh produced by 2D GRADMESH mesh generator */
  char *Filename,    /** in: name of the file to read mesh data */
  FILE* interactive_output // in: name of the output file to write out messages
  );

//------------------------------------------------------------
//  mmr_init_mesh - to initialize the mesh data structure
//------------------------------------------------------------
// NOTE: this is NOT a resize function. It does NOT ADD any elements.
// It only creates new mesh and reserves resources for it.
// To add mesh entities use mmr_add_* functions.
//------------------------------------------------------------
extern int mmr_init_mesh2(  /** returns: >0 - Mesh ID, <0 - error code */
	FILE* Interactive_output, // in: name of the output file to write out messages
	const int N_nodes, //IN: target count of nodes
	const int N_edges, //IN: target count of edges
	const int N_faces, //IN: target count of faces
	const int N_elems  //IN: target count of elements
  );

/**--------------------------------------------------------
  mmr_export_mesh - to export (write to a file) mesh data in a specified format
---------------------------------------------------------*/
extern int mmr_export_mesh( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Control,	  /** in: control variable to choose data format */
        /** MMC_MOD_FEM_MESH_DATA = 1 - mesh written to standard dump files */
  char *Filename  /** in: name of the file to write mesh data */
  );

/**--------------------------------------------------------
  mmr_test_mesh - to test the integrity of mesh data
---------------------------------------------------------*/
extern int mmr_test_mesh( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_mesh_i_params - to return mesh parameters
---------------------------------------------------------*/
int mmr_mesh_i_params(  /** returns: >=0 - integer mesh parameter,
                                     <0 - error code  */
  int Mesh_id,  /** in: pointer to the mesh data structure */
  int Num          /** in: parameter number in control structure */
  );

/**--------------------------------------------------------
  mmr_get_nr_elem - to return the number of active elements
---------------------------------------------------------*/
extern int mmr_get_nr_elem(/** returns: >=0 - number of active elements, */
			   /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_max_elem_id - to return the maximal element id
---------------------------------------------------------*/
extern int mmr_get_max_elem_id(  /** returns: >=0 - maximal element id */
			         /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_next_act_elem - to return the next active element id
---------------------------------------------------------*/
extern int mmr_get_next_act_elem(/** returns: >=0 - the next active element ID */
			   /**                 (0 - input is the last element) */
		           /**                 <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Nel       /** in: input element (0 - return first active element) */
  );

/**--------------------------------------------------------
  mmr_get_next_elem_all - to return the next element id (active or inactive)
---------------------------------------------------------*/
extern int mmr_get_next_elem_all( /** returns: >=0 - ID of the next element */
			   /**                (0 - input is the last element) */
		           /**           <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Nel       /** in: input element (0 - return first element) */
  );

/**--------------------------------------------------------
  mmr_get_nr_face - to return the number of active faces
---------------------------------------------------------*/
extern int mmr_get_nr_face(/** returns: >=0 - number of active faces, */
		    /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_max_face_id - to return the maximal face id
---------------------------------------------------------*/
extern int mmr_get_max_face_id(  /** returns: >=0 - maximal face id */
			  /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_next_act_face - to return the next active face id
---------------------------------------------------------*/
extern int mmr_get_next_act_face( /** returns: >=0 - ID of the next active face */
			   /**                (0 - input is the last face) */
		           /**           <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Nfa       /** in: input face (0 - return first active face) */
  );

/**--------------------------------------------------------
  mmr_get_next_face_all - to return the next face id (active or inactive)
---------------------------------------------------------*/
extern int mmr_get_next_face_all( /** returns: >=0 - ID of the next face */
			   /**                (0 - input is the last face) */
		           /**           <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Nfa       /** in: input face (0 - return first face) */
  );

/**--------------------------------------------------------
  mmr_get_nr_edge - to return the number of active edges
---------------------------------------------------------*/
extern int mmr_get_nr_edge(/** returns: >=0 - number of active edges, */
		    /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_max_edge_id - to return the maximal edge id
---------------------------------------------------------*/
extern int mmr_get_max_edge_id(  /** returns: >=0 - maximal edge id */
			  /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_next_edge_all - to return the next edgeent id (active or inactive)
---------------------------------------------------------*/
extern int mmr_get_next_edge_all( /** returns: >=0 - ID of the next edge */
			   /**                (0 - input is the last edge) */
		           /**           <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ned       /** in: input edge (0 - return first edge) */
  );

/**--------------------------------------------------------
  mmr_get_nr_node - to return the number of active nodes
---------------------------------------------------------*/
extern int mmr_get_nr_node(/** returns: >=0 - number of active nodes, */
		    /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_max_node_id - to return the maximal node id
---------------------------------------------------------*/
extern int mmr_get_max_node_id(  /** returns: >=0 - maximal node id */
			  /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_next_node_all - to return the next node id (active or inactive)
---------------------------------------------------------*/
extern int mmr_get_next_node_all( /** returns: >=0 - ID of the next node */
			   /**                (0 - input is the last node) */
		           /**           <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Nno       /** in: input node (0 - return first node) */
  );

/**--------------------------------------------------------
  mmr_get_max_gen - to get maximal allowed generation level for elements
------------------------------------------------------------*/
extern int mmr_get_max_gen(
/** returns: >=0 - maximal generation, <0-error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_set_max_gen - to set maximal allowed generation level for elements
------------------------------------------------------------*/
extern int mmr_set_max_gen(/** returns: >=0-success code, <0-error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Max_gen   /** in: maximal generation difference */
  );

/**--------------------------------------------------------
  mmr_get_max_gen_diff - to get maximal allowed generation difference between
                         neighboring elements
------------------------------------------------------------*/
extern int mmr_get_max_gen_diff(
/** returns: >=0 - maximal generation difference, <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_set_max_gen_diff - to set maximal allowed generation difference between
                         neighboring elements
------------------------------------------------------------*/
extern int mmr_set_max_gen_diff(/** returns: >=0-success code, <0-error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Max_gen_diff    /** in: maximal generation difference */
  );

/**-----------------------------------------------------------
  mmr_init_ref - to initialize the process of refinement
------------------------------------------------------------*/
extern int mmr_init_ref(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
);

/**-----------------------------------------------------------
  mmr_refine_el - to refine an element WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_refine_el( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El   /** in: element ID  */
);

/**-----------------------------------------------------------
mmr_derefine_el - to derefine an element WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_derefine_el( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El  /** in: element ID  */
);

/**-----------------------------------------------------------
  mmr_refine_mesh - to refine the WHOLE mesh WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_refine_mesh( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
);

/**-----------------------------------------------------------
mmr_derefine_mesh - to derefine the WHOLE mesh WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_derefine_mesh( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
);

/**-----------------------------------------------------------
  mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready
                                  for dofs projection
------------------------------------------------------------*/
extern int mmr_is_ready_for_proj_dof_ref(
  int Mesh_id
);

/**-----------------------------------------------------------
mmr_r_refine - to r-refine an elements with given boundary condition
------------------------------------------------------------*/
extern int mmr_r_refine( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Bc,
    /** in: boundary conditon ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
  void (*reallocation_func)(double * x,double * y, double * z)
    /** function moving to new position*/
);

/**-----------------------------------------------------------
mmr_gen_boundary_layer - to generate boundary layer (for tetrahedral meshes only)
------------------------------------------------------------*/
extern int mmr_gen_boundary_layer(
			   /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Bc,
    /** in: boundary conditon ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
  int thicknessProc,
    // in: total thickness of boundary layers expressed as % of mesh size
  int noLayers, //in: number of layers
  int distribuiton, //in: layer distributin (linear=1, exp=2)
  double * ignoreVect
);

/**-----------------------------------------------------------
  mmr_final_ref - to finalize the process of refinement after rewriting DOFs
------------------------------------------------------------*/
extern int mmr_final_ref(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id 	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
);

/**--------------------------------------------------------
  mmr_free_mesh - to free space allocated for mesh data structure
---------------------------------------------------------*/
extern int mmr_free_mesh(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id 	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );


/**--------------------------------------------------------
mmr_elem_structure - to return elem structure (e.g. for sending)
---------------------------------------------------------*/
extern int mmr_elem_structure(/** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,  	   /** in: elem ID */
  int* Elem_struct /** out: elem structure in the form of integer array */
  );

/**--------------------------------------------------------
  mmr_el_status - to return element status (active, inactive, free space)
---------------------------------------------------------*/
extern int mmr_el_status( /** returns element status: */
			/** +1 (MMC_ACTIVE)   - active element */
			/**  0 (MMC_FREE)     - free space */
			/** -1 (MMC_INACTIVE) - inactive (refined) element */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El  	/** in: element ID */
  );


/**--------------------------------------------------------
  mmr_el_type - to return element type
---------------------------------------------------------*/
extern int mmr_el_type( /** returns element type or <0 - error code */
			/*	 7 (MMC_TETRA) - tetrahedron */
			/*	 5 (MMC_PRISM) - prism */
			/*	 6 (MMC_BRICK) - brick */
			/*	 0 (MMC_FREE)  - free space */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El  	/** in: element ID */
  );

/**--------------------------------------------------------
  mmr_el_groupID - to return group ID for element
---------------------------------------------------------*/
/** \brief Replaces old mmr_el_mat() routine.
 *  Returns ID of group for given element number.
 *  In order to get the material number, pass obtained groupID into \link utr_mat_get_matID \endlink.
 *  In order to get the block number, pass obtained groupID into
 *  \link utr_mat_get_blockID \endlink .
 *
 * \param Mesh_id mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
 * \param El element ID
 * \return returns group flag for element, <0 - error code
 *
 */
extern int mmr_el_groupID( /**  */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El  	/** in: element ID */
  );

/**--------------------------------------------------------
mmr_el_set_groupID - to set group number for element
---------------------------------------------------------*/
/** \brief Replaces old mmr_el_set_mat() routine.
 *  Assigns given element for group with given ID.
 *  Note, that this can effectively change element material number ( \link utr_mat_get_matID \endlink ),
 *  element block assignment ( \link utr_mat_get_blockID \endlink ),
 *  and is some circumstances also affects boundary conditions definitions ( \link utr_mesh_insert_BC_contact \endlink ).
 *
 * \param Mesh_id mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
 * \param El element ID
 * \param Group_id
 * \return group ID
 *
 */
extern int mmr_el_set_groupID(/** sets group flag for element */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,  	/** in: element ID */
  int Group_id    /** in: group ID */
  );

/**--------------------------------------------------------
  mmr_el_type_ref - to return element's type of refinement
---------------------------------------------------------*/
extern int mmr_el_type_ref(
		/** returns element's type of refinement, <0 - error code */
                /**         MMC_NOT_REF - not refined */
                /**         MMC_REF_ISO - isotropic refinement */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El  	/** in: element ID */
  );

/**--------------------------------------------------------
mmr_el_faces - to get faces of an element
---------------------------------------------------------*/
extern int mmr_el_faces( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,     	/** in: element ID */
  int* Faces,  	/** out: list of faces */
                /*	(Faces[0] - number of faces) */
  int* Orient	/** out: orientation for each face */
	      	/*	+1 (MMC_SAME_ORIENT) - the same as element */
	       	/*	-1 (MMC_OPP_ORIENT) - opposite */
  );

/**--------------------------------------------------------
  mmr_el_node_coor - to get the coordinates of element's nodes
---------------------------------------------------------*/
extern int mmr_el_node_coor( 	/** returns number of nodes */
  int Mesh_id,	  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,  	  /** in: element ID */
  int* Nodes,  	  /** out: list of vertex node IDs */
       		  /*	(Nodes[0] - number of nodes) */
  double* Node_coor /** out: coordinates of element vertices */
  );

/**--------------------------------------------------------
mmr_el_edges - to get the list of element's edges
---------------------------------------------------------*/
int mmr_el_edges( 	/** returns the number of edges or error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,  	  /** in: element ID */
  int* Edges   	  /** out : list of edge IDs */
       		  /*	(Edges[0] - number of edges) */
  );

/**--------------------------------------------------------
  mmr_el_fam - to return family information for an element
---------------------------------------------------------*/
extern int mmr_el_fam( /** returns: element father ID or 0 (MMC_NO_FATH) for */
		       /**          initial mesh element or <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,       /** in: element ID (1...) */
  int *Elsons,  /** out: list of element's sons */
               	/** 	Elsons[0] - number of sons */
  int *Type	/** out: type of refinement */
                /**         MMC_NOT_REF - not refined */
                /**         MMC_REF_ISO - isotropic refinement */
  );

//---------------------------------------------------------
// mmr_el_fam_all - to return all recursive family information for an element
//---------------------------------------------------------
extern int mmr_el_fam_all( /** returns: element father ID or 0 (MMC_NO_FATH) for */
		       /**          initial mesh element or <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,       /** in: element ID (1...) */
  int *Elsonscd   /** out: list of all recursive element's sons */
               	/** 	Elsons[0] - number of sons */
  );


/**--------------------------------------------------------
  mmr_el_gen - to return generation level for an element
---------------------------------------------------------*/
extern int mmr_el_gen(   /** returns: El's generation ID */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El        /** in: element ID (1...) */
  );

/**--------------------------------------------------------
  mmr_el_ancestor - to find the ancestor of an element
                    with generation level Ilev
---------------------------------------------------------*/
extern int mmr_el_ancestor( /** returns: >0 - ancestor ID, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,       /** in: element ID (1...) */
  int Ilev      /** in: level ID */
  );

/**--------------------------------------------------------
  mmr_el_hsize - to compute a characteristic linear size for an element
	         (for linear and multi-linear 3D elements)
---------------------------------------------------------*/
extern double mmr_el_hsize(   /** returns: element size */
  int Mesh_id,	  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,         /** in: element ID (1...) */
  double *Size_x, /** out: for anizotropic elements */
  double *Size_y,
  double *Size_z
  );

/**--------------------------------------------------------
mmr_el_eq_neig - to get equal size neighbors of an element
---------------------------------------------------------*/
extern int mmr_el_eq_neig( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,     	/** in: element ID */
  int* Neig,   	/** out: list of equal size neighbors */
       		/**  >0 - equal size neighbor */
       		/**  -1 (MMC_BIG_NGB) - big neighbor */
       		/**   0 (MMC_BOUNDARY) - boundary (always second neighbor)*/
  int* Neig_sides /** out: list of sides of neighbors */
  );

/**--------------------------------------------------------
mmr_face_structure - to return face structure (e.g. for sending)
---------------------------------------------------------*/
extern int mmr_face_structure(/** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,  	   /** in: face ID */
  int* Face_struct /** out: face structure in the form of integer array */
  );

/**--------------------------------------------------------
  mmr_fa_status - to return face status (active, inactive, free space)
---------------------------------------------------------*/
extern int mmr_fa_status( /** returns face status: */
			  /*	 1 (MMC_ACTIVE)   - active face */
			  /*	 0 (MMC_FREE)     - free space */
			  /*	-1 (MMC_INACTIVE) - inactive face */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa  	/** in: face ID */
  );

/**--------------------------------------------------------
  mmr_fa_type - to return face type (triangle, quad, free space)
---------------------------------------------------------*/
extern int mmr_fa_type( /** returns face type: */
			/*	 3 (MMC_TRIA) - triangle */
			/*	 4 (MMC_QUAD) - quadrilateral */
			/*	 0 (MMC_FREE) - free space */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa  	/** in: face ID */
  );

/**-------------------------------------------------------------------------
  mmr_fa_bc - to get the boundary condition flag for a face
---------------------------------------------------------------------------*/
extern int mmr_fa_bc( /** returns: bc flag for a face */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa	      /** in: global face ID */
  );

/**-------------------------------------------------------------------------
  mmr_fa_bc - to set the boundary condition flag for a face
---------------------------------------------------------------------------*/
extern int mmr_fa_set_bc( /** returns: bc flag for a face */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,	      /** in: global face ID */
  int Bc_num    //in: Boundary Condition number associated with face Fa
  );


/**--------------------------------------------------------
  mmr_fa_sub_bnd - to indicate face is on the boundary
---------------------------------------------------------*/
extern int mmr_fa_sub_bnd( /** returns: 1 - true, 0 - false */
  int Mesh_id,    /** in: mesh ID */
  int Face_id     /** in: face ID */
  );

/**--------------------------------------------------------
  mmr_fa_set_sub_bnd - to indicate face is on the boundary
---------------------------------------------------------*/
extern int mmr_fa_set_sub_bnd(/** returns: >=0 - success code, <0 - error code */
  int Mesh_id,    /** in: mesh ID */
  int Face_id,    /** in: face ID */
  int Side_id     /** in: side ID ??? */
  );

/**--------------------------------------------------------
mmr_fa_edges - to return a list of face's edges
---------------------------------------------------------*/
extern int mmr_fa_edges( /** returns: >=0 - number of edges, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,         /** in:  global face ID */
  int *Fa_edges,  /** out: face edges */
  int *Ed_orient  /** out: edges orientation: */
	          /*	+1 (MMC_SAME_ORIENT) - the same as face */
	          /*	-1 (MMC_OPP_ORIENT)  - opposite */
  );

/**--------------------------------------------------------
mmr_fa_eq_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (only equal size
		neighbors considered)
---------------------------------------------------------*/
extern void mmr_fa_eq_neig(
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,       /** in:  global face ID */
  int *Fa_neig, /** out: face neighbors */
       		/** (first neighbor has the same orientation */
	       	/**  as the face and ordering of its nodes */
	       	/**  defines ordering of nodes on the face) */
	       	/**  >0 - equal size neighbor ID */
       		/**  -1 (MMC_BIG_NGB) - big neighbor */
       		/**   0 (MMC_BOUNDARY) - boundary (always second neighbor)*/
       		/**  -2 (MMC_SUB_BND) - inter-subdomain boundary */
                /**  <-10 - ghost face ID for contact boundary */
  int *Neig_sides,      /** out: side local IDs for equal size neighbors */
  int *Node_shift	/** out: the difference in positions between */
			/*	first face's node for first and second */
			/** 	neighbor (usage in mmr_fa_node_coor) */
  );

/**--------------------------------------------------------
  mmr_fa_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (for active
		faces all neighbors are active)
---------------------------------------------------------*/
extern void mmr_fa_neig(
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,               /** in:  global face ID */
  int *Fa_neig,         /** out: face neighbors */
	       	        /** (first neighbor has the same orientation */
	       	        /**  as the face and ordering of its nodes */
	       	        /**  defines ordering of nodes on the face) */
	       	        /**  >0 - equal size neighbor ID */
	       	        /**  <0 - big neighbor ID */
       		        /**   0 (MMC_BOUNDARY) - boundary (always 2nd neighbor)*/
  int *Neig_sides,      /** out: side IDs for neighbors */
  int *Node_shift,      /** out: the difference in positions between */
       		        /*	first face's node for first and second */
       		        /** 	neighbor (usage in mmr_fa_node_coor) */
  int *Diff_gen,	/** out: generation difference between neighbors */
  double *Acoeff,       /** out: coefficients of linear transformation... */
  double *Bcoeff	/** 	between face coordinates and */
          		/** 	coordinates on an ancestor face  */
			/** 	being a side of big neighbor  */
  );

/**-------------------------------------------------------------------------
mmr_fa_fam - to return face's family information
---------------------------------------------------------------------------*/
extern int mmr_fa_fam( /** returns: face's father */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,      	/** in: global face ID */
  int* Fasons,	/** out: sons */
               	/** 	Fasons[0] - number of sons */
  int* Node_mid	/** out: node in the middle (if any) */
  );

/**--------------------------------------------------------
mmr_fa_node_coor - to get the list and coordinates of faces's nodes
---------------------------------------------------------*/
extern int mmr_fa_node_coor( /** returns: number of nodes for a face */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,  	     /** in: face ID */
  int* Nodes,  	     /** out: list of vertex node IDs */
  double* Node_coor  /** out: coordinates of face vertices */
  );


/**--------------------------------------------------------
  mmr_el_fa_nodes - to get list local face nodes indexes in elem
---------------------------------------------------------*/
extern int mmr_el_fa_nodes( // returns face type flag
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,		// in: global elem ID
  int Fa,		/** in: local face number in elem El */
  int *Fa_nodes		/** out: list of local indexes of face nodes */
  );

/**--------------------------------------------------------
  mmr_fa_elem_coor - to find coordinates within neighboring
	             elements for a point on face
---------------------------------------------------------*/
extern void mmr_fa_elem_coor(
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  double *Xloc,	        /** in: local coordinates on a face */
  int *Fa_neig,	        /** in: face neighbors (<0 - big) */
                        /*	first - same orientation */
	       	        /*	second - opposite orientation */
  int *Neig_sides,      /** in: which side face is for neighbors */
  int Node_shift,       /** in: the difference in positions between */
			/*	first face's node for first and second */
			/** 	neighbor (usage in mmr_fa_node_coor) */
  double *Acoeff,       /** in: coefficients of transformation between... */
  double *Bcoeff,	/** in: ...face coord and big neighb face coord */
  double *Xneig         /** out: local coordinates for neighbors */
  );

/**--------------------------------------------------------
  mmr_fa_area - to compute the area of face and vector normal
---------------------------------------------------------*/
extern void mmr_fa_area(
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,		/** in: global face ID */
  double *Area,	        /** out: face area */
  double *Vec_norm      /** out: normal vector */
  );

/**--------------------------------------------------------
mmr_edge_status - to return edge status (active, inactive, free space)
---------------------------------------------------------*/
extern int mmr_edge_status( /** returns edge status: */
			/** +1 (MMC_ACTIVE)   - active edge */
			/**  0 (MMC_FREE)     - free space */
			/** -1 (MMC_INACTIVE) - inactive (refined) edge */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed  	/** in: edge ID */
  );

/**--------------------------------------------------------
mmr_edge_nodes - to return edge node's IDs
---------------------------------------------------------*/
extern int mmr_edge_nodes( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed,		/** in: edge ID */
  int *Edge_nodes	/** out: IDs of edge nodes */
  );

/**--------------------------------------------------------
mmr_edge_elems - to return IDs of elements containing the edge
---------------------------------------------------------*/
extern int mmr_edge_elems( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed,		/** in: edge ID */
  int *Edge_elems	/** out: IDs of elements containing the edge */
                        /**      Edge_elems[0] - the number of elements */
  );

/**-----------------------------------------------------------
  mmr_create_edge_elems - to create (or recreate if Max_edge_id > 0) for each
                  edge a list of elements to which it belongs
------------------------------------------------------------*/
int mmr_create_edge_elems( /** returns: 1-success, <=0-failure */
  int  Mesh_id,	 /** in: ID of the mesh to be used or 0 for the current mesh */
  int Max_edge_id /** in: the range of edge IDs to consider or 0 for default */
  /** if Max_edge_id==0 it is assumed that there are no structures to free !!! */
  );

/**--------------------------------------------------------
mmr_edge_sons - to return edge son's numbers
---------------------------------------------------------*/
extern int mmr_edge_sons( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed,		/** in: edge ID */
  int *Edge_sons,	/** out: IDs of edge sons */
  int *Edge_mid_node	/** out: IDs of edge mid-node son */
  );

/**--------------------------------------------------------
mmr_edge_structure - to return edge structure (e.g. for sending)
---------------------------------------------------------*/
extern int mmr_edge_structure(/** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed,  	   /** in: edge ID */
  int* Edge_struct /** out: edge structure in the form of integer array */
  );

/**--------------------------------------------------------
mmr_set_edge_type - to set type for an edge
---------------------------------------------------------*/
extern int mmr_set_edge_type( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Edge_id,     /** in: edge ID */
  int Type         /** in: edge type (the number of attempted subdivisions !!!)*/
  );

/**--------------------------------------------------------
mmr_set_edge_fam - to set family data for an edge
---------------------------------------------------------*/
extern int mmr_set_edge_fam( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Edge_id,     /** in: edge ID */
  int Son1,        /** in: first son ID */
  int Son2         /** in: second son ID */
  );

/**--------------------------------------------------------
mmr_set_face_fam - to set family data for an face
---------------------------------------------------------*/
extern int mmr_set_face_fam( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Face_id,     /** in: face ID */
  int *Sons       /** in: sons IDs */
  );

/**--------------------------------------------------------
mmr_set_face_neig - to set neighbors data for a face
---------------------------------------------------------*/
extern int mmr_set_face_neig( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Face_id,     /** in: face ID */
  int Neig1,       /** in: first neig ID */
  int Neig2,        /** in: second neig ID */
  int Neig1Type,    /** in: first neig type */
  int Neig2Type     /** in: second neig type */
  );

/**--------------------------------------------------------
mmr_set_elem_fath - to set family data for an elem
---------------------------------------------------------*/
int mmr_set_elem_fath( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Elem_id,     /** in: elem ID */
  int Fath        /** in: father ID */
  );

/**--------------------------------------------------------
mmr_set_elem_fam - to set family data for an elem
---------------------------------------------------------*/
extern int mmr_set_elem_fam( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Elem_id,     /** in: elem ID */
  int Fath,        /** in: father ID */
  int *Sons        /** in: sons IDs */
  );

/**--------------------------------------------------------
mmr_node_status - to return node status (active, inactive, free space)
---------------------------------------------------------*/
extern int mmr_node_status( /** returns node status: */
			/** +1 (MMC_ACTIVE)   - active node */
			/**  0 (MMC_FREE)     - free space */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Node  	/** in: node ID */
  );

/**--------------------------------------------------------
mmr_node_coor - to return node coordinates
---------------------------------------------------------*/
extern int mmr_node_coor( /** returns success (>=0) or error (<0) code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Node,  	/** in: node ID */
  double *Coor  /** out: coordinates */
  );

/*---------------------------------------------------------
mmr_set_node_coor - set node coordinates
---------------------------------------------------------*/
extern int mmr_set_node_coor( /* returns success (>=0) or error (<0) code */
  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Node,  	/* in: node ID */
  double *Xcoor  	/* in: node coordinates */
  );

/*---------------------------------------------------------
mmr_move_node - moves node by vector
---------------------------------------------------------*/
extern int mmr_move_node( /* returns success (>=0) or error (<0) code */
  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Node,  	/* in: node ID */
  double *Dx  	/* in: displacement (vector) */
  );

/**--------------------------------------------------------
mmr_loc_loc - to compute local coordinates within an element,
	given local coordinates within an element of the same family
---------------------------------------------------------*/
extern int mmr_loc_loc(/** returns: 1 - success, 0 - failure */
  int Mesh_id,   /** in: field ID */
  int El_from, 	/** in: element number */
  double* X_from, /** in: local element coordinates */
  int El_to, 	/** in: another element number */
  double* X_to	/** out: local another element coordinates */
  );

extern int mmr_create_element(
		 /** returns: ID of the created element (<=0 - failure) */
  int  Mesh_id,	 /** in: ID of the mesh to be used or 0 for the current mesh */
  int  Type,     /** in: type for the face */
  int  Mate,     /** in: material indicator */
  int  Fath,     /** in: father element ID */
  int  Refi,     /** in: refinement type indicator */
  int* Faces,    /** in: list of faces' IDs */
  int* Sons      /** in: list of sons (only for inactive elements, Type<0) */
  );

extern int mmr_create_face( /** returns: ID of created face (<0 - error code) */
  int  Mesh_id,	 /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int  Type,     /** in: type for the face */
  int  Flag_bc,  /** in: bc flag for the face */
  int* Edges,    /** in: edges for the new face */
  int* Neig,     /** in (optional): neighbors (or NULL) */
  int* Sons      /** in (optional): sons (if Type<0) or NULL */
);

extern int mmr_create_edge( /** returns: ID of a new edge (<0 - error code)*/
  int  Mesh_id,	 /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Type,      /** in: type indicator (MMC_EDGE or */
                 /**     number of attempted divisions for inactive edges0 */
  int  Node1,    /** in: nodes for a new edge */
  int  Node2
);

extern int mmr_create_node( /** returns: node ID of created node (<0 - error) */
  int  Mesh_id,	 /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  double Xcoor,  /** in: coordinates of new node */
  double Ycoor,
  double Zcoor
);

extern int mmr_clust_face( /** returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Face,      /** in: face ID */
  int* Face_sons /** in (optional): face sons' IDs */
);

extern int mmr_clust_edge( /** returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Edge_id   /** in: clustered edge ID */
);


extern int mmr_divide_face(/** returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	   /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int  Fa,         /** in: face to be divided */
  int* Face_sons,  /** out: sons (sons[0] - the number of sons) */
  int* Sons_edges, /** out: created new edges (edges[0] - the number of edges)*/
  int* New_nodes   /** out: created new nodes (nodes[0] - the number of nodes)*/
);

extern int mmr_divide_edge( /** returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int  Edge,      /** in: edge to be divided */
  int* Edge_sons, /** out: two sons */
  int* Node_mid   /** out: node in the middle */
);

/**-----------------------------------------------------------
mmr_del_elem    - to free an element structure
------------------------------------------------------------*/
extern int mmr_del_elem( /** returns: 1-success, <=0-failure */
  int  Mesh_id,	/** in: ID of the mesh to be used or 0 for the current mesh */
  int Elem      /** in: element ID */
  );

/**-----------------------------------------------------------
  mmr_del_face    - to free a face structure
------------------------------------------------------------*/
extern int mmr_del_face( /** returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Face      /** in: face ID */
);

/**-----------------------------------------------------------
  mmr_del_edge    - to free an edge structure
------------------------------------------------------------*/
extern int mmr_del_edge( /** returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Edge      /** in: edge ID */
);

/**-----------------------------------------------------------
  mmr_del_node    - to free a node structure
------------------------------------------------------------*/
extern int mmr_del_node( /** returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Node      /** in: node ID */
);

/**--------------------------------------------------------
  mmr_get_max_elem_max - to return the maximal possible elem id
---------------------------------------------------------*/
extern int mmr_get_max_elem_max(  /** returns: >=0 - maximal possible elem id */
			  /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_max_face_max - to return the maximal possible face id
---------------------------------------------------------*/
extern int mmr_get_max_face_max(  /** returns: >=0 - maximal possible face id */
			  /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_max_edge_max - to return the maximal possible edge id
---------------------------------------------------------*/
extern int mmr_get_max_edge_max(  /** returns: >=0 - maximal possible edge id */
			  /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_max_node_max - to return the maximal possible vertex id
---------------------------------------------------------*/
extern int mmr_get_max_node_max(  /** returns: >=0 - maximal possible vertex id */
			  /**           <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );


//******************************************************************************//
// The five next functions are not supported by mmd_prism and mmd_prism_2D modules!!!!!!!
//******************************************************************************//

///---------------------------------------------------------
/// mmr_init_read - to inform mesh module about new target number of elems,faces,edges and nodes
///------------------------------------------------------------
/// NOTE: this is NOT a resize function. It does NOT ADD/REMOVE any elements.
/// It only reserves resources.
/// To add mesh entities use mmr_add_* functions.
///---------------------------------------------------------
extern int mmr_init_read( // return >=0 - success, <0 - error code
    const int Mesh_id, /// IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
    const int N_nodes, ///IN: target count of nodes
    const int N_edges, ///IN: target count of edges
    const int N_faces, ///IN: target count of faces
    const int N_elems  ///IN: target count of elements
  );


extern int mmr_finish_read( // return >=0 - success, <0 - error code
    const int Mesh_id // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
  );


//---------------------------------------------------------
//  mmr_add_elem - to add (append) new element
//---------------------------------------------------------
extern int mmr_add_elem(// returns: added element id
    const int Mesh_id, // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
    const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
    const int El_type, // IN: type of new element
    const int El_nodes[6], //IN: if known: nodes of new element, otherwise NULL
    const int El_faces[5]  //IN: if known: faces of new element, otherwise NULL
    , const int Material_idx); //IN: material id
// NOTE: if both El_nodes and El_faces are NULLs
// then an element is udefined and an error is returned!

//---------------------------------------------------------
//  mmr_add_face - to add (append) new face
//---------------------------------------------------------
extern int mmr_add_face(// returns: added face id
    const int Mesh_id,  // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
    const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
    const int Fa_type,  // IN: type of new face
    const int Fa_nodes[4], // IN: if known: nodes of new face, otherwise NULL
    const int Fa_edges[4], // IN: if known: edges of new face, otherwise NULL
    const int Fa_neigs[2]  // IN: if known: neigs of new face, otherwise NULL
    , const int B_cond_val);
// NOTE: if both Fa_nodes and Fa_edges are NULLs
// then a face is udefined and an error is returned!

//---------------------------------------------------------
//  mmr_add_edge - to add (append) new edge
//---------------------------------------------------------
extern int mmr_add_edge(// returns: added edge id
	const int Mesh_id,  //IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
	const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
 	const int nodes[2]  //IN: nodes for new edge
	);
// NOTE: if nodes==NULL edge is undefined and an error is returned!

//---------------------------------------------------------
//  mmr_add_node - to add (append) new node
//---------------------------------------------------------
extern int mmr_add_node( // returns: added node id
	const int Mesh_id,  //IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
	const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
	const double Coords[3] //IN: geometrical coordinates of node
	);
// NOTE: if Coords==NULL node is undefined and an error is returned!


///
/// \brief mmr_current_mesh_state_hash - returns an current unique value (state hash, state stamp)
/// that can be used to compare with value returned from this function at other time.
/// If values are different, then one can assume that mesh internal state has changed.
/// \param Mesh_id
/// \return
///
extern int mmr_current_mesh_state_hash( // returns a state-stamp for mesh
        const int Mesh_id); // Mesh_id for state-stamp


//******************************************************************************//
// Below are functions used by modules intarfacing with REMESH module ONLY !!!!
//******************************************************************************//

/**--------------------------------------------------------
  mmr_init_dist2boundary  - to return distance from given vertex to nearst boundary
---------------------------------------------------------*/
extern int mmr_init_dist2bound( // returns 0 if succesfull
  const int Mesh_id, //in: mesh ID or 0 (MMC_CUR_MESH_ID)
  const int* BCs, //in: c-array of accepted BC numbers
  const int nBCs); //in: length of BCs parameter


/**--------------------------------------------------------
  mmr_get_el_dist2boundary  - to return approx distance from point within element to nearst boundary
---------------------------------------------------------*/
extern double mmr_get_el_dist2boundary(  /** returns: dist>=0.0  */
			  /**          <0 - "far far away from boundary" code */
  const int Mesh_id, /*in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh*/
  const int El_id,	 // in: El_id wihch includes point Coord
  const double * Coord); // GLOBAL coords of point inside El_id



extern void mmr_get_coor_from_motion_element(
  const int Mesh_id,
  int idEl,
  int idLP,
  double *coor,
  int flagaSiatki
);

/**-----------------------------------------------------------
  mmr_copyMESH
------------------------------------------------------------*/

extern void mmr_copyMESH(const int Mesh_id,int flaga);


/**-----------------------------------------------------------
  mmr_create_mesh_Cube - test movement mesh
------------------------------------------------------------*/

extern void mmr_create_mesh_Cube(
  const int Mesh_id,
  const char *nazwa,
  int node_x,
  int node_y,
  int node_z,
  double size_x,
  double size_y,
  double size_z,
  double divide,
  int *warunki
);

/**-----------------------------------------------------------
  mmr_test_mesh_motion - test movement mesh
------------------------------------------------------------*/

extern void mmr_init_all_change(const int Mesh_id,int a);

extern void mmr_test_mesh_motion(const int Mesh_id,int ileWarstw,int obecny_krok, int od_krok,int ileKrok,double minPoprawy,double px0,double py0,double pz0,double px1,double py1,double pz1,double endX,double endY,double endZ);

extern double mmr_test_weldpool(const int Mesh_id,double obecnyKrok, double krok_start,double minPoprawy,double doX,double doY,double dl,double zmiejszaPrzes,double limit,double szerokosc);


extern int mmr_groups_number(const int Mesh_id);
extern int mmr_groups_ids(const int Mesh_id, int *tab, int l_tab);

///
/// \brief mmr_get_fa_el_bc_connect
/// \param face_id
/// \param el_id out: [0][1][2] - number of correspoding gauss points on the other (contact) 'connected face'
/// [3] - id of neighbouring element with face 'face_id'
/// [4] - id of the other element on the other side of 'connected face'
/// \return 0 if normal bc else return 'connected face_id'
///

extern void mmr_get_vec_norm_at_node(const int Mesh_id,int node_id,double *vec_norm);

extern int mmr_get_fa_el_bc_connect(
		const int Mesh_id,
        int face_id,
        int *el_id
        );

/** @} */ // end of group


#ifdef __cplusplus
}
#endif

#endif
