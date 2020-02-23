/******************************************************************************
File uth_intf.h - interface to utility procedures (and macros)
				  1. general purpose utilities (implementation in uts_util.c)
				  2. time measurments (implementation system dependent !!!)

Contains declarations of constants and interface routines:

General purpose utilities:

  utr_ctrl_c_signal_handler(int param) - to catch kill signals
  utr_set_interactive - to set up names for interactive input and output
		that control the run (possibly stdin and/or stdout for on-line control)
  utr_initialize_mesh - to initialize mesh of a specified type
  utr_initialize_field - to initialize field of a specified type
  utr_write_paraview_std_lin - to dump data in Paraview format (for std_lin only)
  utr_write_paraview_partmesh - to dump mesh in Paraview format with partition info
  !utr_export_mesh - to export mesh of a specified type in a specified format
  !utr_export_field - to export field of a specified type in a specified format

  utr_get_ent_dofs - to read geometry or solution DOFs for a given DOF entity
				 from suitable (approximation or mesh) data structures

in uts_adapt.c
  utr_create_patches - to create patches of elements containing a given node
					  (uses only mesh data - can be the same for all fields)
  utr_create_patches_small - the same as above but uses lists of elements nodes
							 instead of lists of elements dof entities
  utr_recover_derivatives - to recover deirvatives using patches of elements
					  (can be called separately for each field on the same mesh)
  utr_recover_derivatives_small-to recover derivatives using patches of elements
	(does not use apr_get_stiff_mat_data for getting lists of element nodes)
  utr_adapt - to enforce default adaptation strategy for SINGLE problem/field
  utr_test_refinements - to perform sequence of random refinements/unrefinements
  utr_manual_refinement - to perform manual refinement/unrefinement

in uts_ls_intf.c
  utr_get_list_ent - to return the list of:
					  1. integration entities - entities
						 for which stiffness matrices and load vectors are
						 provided by the FEM code to the solver module,
					  2. DOF entities - entities with which there are dofs
						 associated by the given approximation
  utr_get_list_dof_ent - to return the list  DOF entities
						 - entities with which there are dofs
			 associated by the given approximation
						 given the list of integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module
						 (not for DG approximation!!!)
  utr_get_list_int_ent - to return the list of integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module
						 (not for DG approximation!!!)
  utr_get_list_ent_coarse - to return the list of integration entities - entities
						 for which stiffness matrices and load vectors are
						 provided by the FEM code to the solver module,
						 and DOF entities - entities with which there are dofs
						 associated by the given approximation for COARSE level
						 given the corresponding lists from the fine level
			  (it calls separate implementations for dg and std?)
  utr_get_max_num_grid_levels - for limiting nr_levels in multigrid
							  based on mesh and field data
  utr_dof_ent_sons - to return info whether the entity is owned
					 and a list of dof entity sons for owned entity
  utr_comp_stiff_mat - to create a stiffness matrix
					  and a load vector corresponding to the specified
					  mesh entity, together with information on how to
					  assemble entries into the global stiffness matrix
					  and the global load vector
  utr_create_assemble_stiff_mat - to create element stiffness matrices
								 and assemble them to the global SM
  utr_renumber - to renumber (permute) graph vertices (used for bandwidth reduction)
  utr_color_int_ent_for_assembly_with_graph_creation - coloring elements for lock free
					  assembly (with creation of neighbourhood graphs)
  utr_color_int_ent_for_assembly - coloring elements for lock free assembly
  utr_rewr_sol - to rewrite solution from one vector to another

in uts_util:
  utr_average_sol_el - to compute the average of solution over element

  utr_dvector - to allocate a double vector: name[0..ncom-1]:
				  name=ut_dvector(ncom,error_text)
  utr_ivector - to allocate an integer vector: name[0..ncom-1]:
				  name=ut_ivector(ncom,error_text)
  utr_d_zero - to zero a double vector of length Num
  utr_chk_list - to check whether Num is on the list List
	with length Ll
  utr_put_list - to put Num on the list List with length Ll
	(filled with numbers and zeros at the end)
  utr_mat3_inv - to invert a 3x3 matrix (stored as a vector!)
  utr_vec3_add - to compute vector product of 3D vectors
  utr_vec3_subst - to compute sum of 3D vectors
  utr_vec3_prod - to compute difference of 3D vectors
  utr_vec3_mxpr - to compute mixed vector product of 3D vectors
  utr_vec3_length - to compute length of a 3D vector
  utr_mat3vec - to compute matrix vector product in 3D space
  utr_mat3mat - to compute matrix matrix product in 3D space
	(all matrices are stored by rows as vectors!)
  utr_mat_det - to compute determinant of a matrix
  utr_skip_rest_of_line - to allow for comments in input files
  utr_fprintf_double - to write double number to file with specified accuracy
					 (in the number of significant digits)
  utr_insert_sort - to sort an array by insertion (for small integer arrays)
  utr_heap_sort - to heap-sort an array (for large double precision arrays)

Input Output routines (IO):
 utr_io_read_mesh - to read mesh from given filename(s)


------------------------------
History:
	  02.2002 - Krzysztof Banas, kbanas@pk.edu.pl
	  10.2012 - Kazimierz.Michalik@agh.edu.pl
	  08.2012 - Kamil Wachala, kamil.wachala@gmail.com
	  2011 - Kazimierz Michalik, kamich@agh.edu.pl
	  2016    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
******************************************************************************/

#ifndef _uth_intf_
#define _uth_intf_

#include <stdio.h>

/** standard macro for max and min and abs */
#define utm_max(x,y) ((x)>(y)?(x):(y))
#define utm_min(x,y) ((x)<(y)?(x):(y))
#define utm_abs(x)   ((x)<0?-(x):(x))

#ifndef UTM_SAFE_FREE_PTR
# define UTM_SAFE_FREE_PTR(ptr) if((ptr)!=NULL) { free((ptr));(ptr)=NULL; };
#endif

#ifndef UTM_SAFE_FREE_PTR_ARRAY
# define UTM_SAFE_FREE_PTR_ARRAY(p_ar,n) if((p_ar)!=NULL){int i=0; for(;i<n;++i){UTM_SAFE_FREE_PTR(p_ar[i]);}};
#endif

#ifdef __cplusplus
extern "C"
{
#endif


/** @defgroup UTM Utilities
 *
 *  @{
 */

typedef enum
{
	UTE_POINT_DATA=0,
	UTE_CELL_DATA,
	UTE_SCALARS,
	UTE_VECTORS,
	UTE_TENSORS,
	UTE_DOUBLE,
	UTE_INT,
	UTE_VTK_LEGACY_FORMAT,
	UTE_VTK_XML,
	UTE_LAST
} ute_paraview_flags;

#define UTC_LIST_END_MARK -1

/** definition of type utt_patches - useful e.g. for error estimation */
#define UTC_MAXEL_PATCH 100 // maximal number of elements in the patch
#define UTC_MAXNO_PATCH 200 // maximal number of nodes in the patch

#define UTC_MAX_INT_ENT_NEIG 300 // maximal number of neighbours for int ent
  // i.e. maximal number of faces and elements that share some DOFS with
  // a given integration entity (face or element)

/*struct for patches - e.g. to recover derivatives by uts_recover_derivatives*/
typedef struct {
	int nr_elems;
	int elems[UTC_MAXEL_PATCH];
	int nr_nodes;
	int nodes[UTC_MAXNO_PATCH];
	double *deriv;
} utt_patches;

typedef struct {
  // practically every utr_ routine can have its own time counter
  double initialize_mesh;
  double initialize_field;
  double decompress_files;
  double problem_initialization;

  double create_int_ent_list;
  double create_dof_ent_list;
  double renumber_dofs;
  double total_coloring;
  double color_elements;
  double color_faces;
  double handle_solver_interface_data_structures;
  double handle_solver_interface_GPU_data_structures;
  double create_solver_data_structures;
  double total_sir_create;

  double read_initial_guess;
  double read_geo_dofs_vector;
  double read_single_sol_dofs_vector;
  double read_all_sol_dofs_vectors;
  double CPU_integration_elems;
  double CPU_assembly_elems;
  double CPU_integration_faces;
  double CPU_assembly_faces;
  double GPU_integration_and_assembly;
  double CPU_GPU_solver_ds_merging;
  double create_sm_uncon;
  double adding_constraints_info;
  double total_integration_and_assembly;
  double solve;
  double write_solution;
  double total_sir_solve;

  double export_mesh;
  double export_field;
  double compress_files;
  double create_patches;
  double recover_derivatives;
  double calculate_error;
  double refine_mesh;
  double project_field;
  double manual_refinemts;
  double total_outside_solver;

  double total;
} utt_time_measurement;

extern utt_time_measurement utv_time;

/**--------------------------------------------------------
  utr_ctrl_c_signal_handler(int param) - to catch kill signals
----------------------------------------------------------*/
  extern int utv_SIGINT_not_caught;
  extern void utr_ctrl_c_signal_handler(int param);

/**--------------------------------------------------------
  utr_reset_time_counters
----------------------------------------------------------*/
  extern void utr_reset_time_counters( int Control);

/**--------------------------------------------------------
  utr_display_time_counters
----------------------------------------------------------*/
  extern void utr_display_time_counters( int Control);


/**--------------------------------------------------------
  utr_set_interactive - to set up pointers for interactive input and output
		 that control the run (possibly stdin and/or stdout for on-line control)
----------------------------------------------------------*/
extern int utr_set_interactive(  /** returns: >0 - success code, <0 - error code */
  char* Workdir, // working directory
  int Argc, // passed from maine
  char **Argv,
  FILE **Interactive_input,
  FILE **Interactive_output
);

/**--------------------------------------------------------
utr_initialize_mesh - to initialize mesh of a specified type
---------------------------------------------------------*/
extern int utr_initialize_mesh( /** returns mesh_id */
  FILE *Interactive_output, /** file or stdout to write messages */
  const char* Work_dir, // path to working directory
  char Mesh_type, /** letter symbol denoting mesh_type (j, p, t or h) */
  const char* Mesh_file /** mesh file name - conforming to naming convention */
  );

/**--------------------------------------------------------
utr_initialize_field - to initialize field of a specified type
---------------------------------------------------------*/
extern int utr_initialize_field(  /** returns: >0 - field ID, <0 - error code */
  FILE *Interactive_output, /** file or stdout to write messages */
  char Field_type, /** options: s - standard, d - discontinuous */
  char Control,	 /** in: control variable: */
				 /**      z - to initialize the field to zero */
				 /**      r - to read field values from the file */
				 /**      i - to initialize the field using function */
				 /**                 provided by the problem dependent module */
  int Mesh_id,	 /** in: ID of the corresponding mesh */
  int Nreq,	 /** in: number of equations - solution vector components */
  int Nr_sol,	 /** in: number of solution vectors for each dof entity */
  int Pdeg,      /** in: degree of approximating polynomial */
  char *Filename,/** in: name of the file to read approximation data */
  double (*Fun_p)(int, double*, int) /** pointer to function that provides */
		 /** problem dependent initial condition data */
  );

/**--------------------------------------------------------
utr_write_paraview_std_lin - to dump data in Paraview format (for std_lin only)
---------------------------------------------------------*/
extern int utr_write_paraview_std_lin(
  int Mesh_id,       /** in: ID of the mesh */
  int Field_id,      /** in: ID of the field */
  const char* Workdir,
  const char *Filename,     /** in: name of the file to write data */
  ute_paraview_flags VTK_file_version
  );

/**--------------------------------------------------------
/// utr_write_paraview_mesh - to dump mesh in Paraview format
---------------------------------------------------------*/
extern int utr_write_paraview_mesh(int Mesh_id,
								   const char*   Work_dir,
								   const char *Filename,
								   ute_paraview_flags VTK_file_version);

/**--------------------------------------------------------
/// utr_write_paraview_partmesh - to dump mesh in Paraview format includiong partition info
---------------------------------------------------------*/
extern void utr_write_paraview_partmesh(
  int Mesh_id,
  const char*   Work_dir,
  const char *Filename,
  int *Parts, /** partition vector for the elements of the mesh */
  ute_paraview_flags VTK_file_version
);

typedef struct {
	int field_id;
	int dofs_write[20];		// dofs to write: dofs_write[0] - number of dofs
							//                dofs_write[i] - IDs of dofs to write
	ute_paraview_flags entity_type;
	ute_paraview_flags quantity_type;
	ute_paraview_flags value_type;
	char const *field_name;
} utt_paraview_field_descriptor;

/**--------------------------------------------------------
/// utr_write_paraview_field - to dump field in Paraview format
---------------------------------------------------------*/
extern int utr_write_paraview_fields(const char*   Work_dir,
  const char*   Filename,
  double Current_time,
  const int     N_descriptors,
  utt_paraview_field_descriptor* Field_descriptors,
  ute_paraview_flags    VTK_file_version
							 );
/**--------------------------------------------------------
/// utr_write_paraview_tensor - to dump tensor in Paraview format
---------------------------------------------------------*/
extern int utr_write_paraview_tensor(
  int Field_id,
  FILE *File,
  int* Dofs_write, // dofs to write: dofs_write[0] - number of dofs
				   //                dofs_write[i] - IDs of dofs to write
  double dt
							 );
/**--------------------------------------------------------
/// utr_write_paraview_bc - to dump boundary conditions 'field' in Paraview format
---------------------------------------------------------*/
extern int utr_write_paraview_bc(
  int Mesh_id,
  const char*   Work_dir,
  const char *Filename,
  ute_paraview_flags VTK_file_version
);

/*---------------------------------------------------------
  utr_get_ent_dofs - to read geometry or solution DOFs for a given DOF entity
				 from suitable (approximation or mesh) data structures
----------------------------------------------------------*/
int utr_get_ent_dofs(
  int Problem_id,     /** in: solver ID (used to identify the subproblem) */
  int Dof_ent_type,
  int Dof_ent_id,
  int Nrdofs,
  int Vect_id,       /** in: vector ID in case of multiple solution vectors */
					 /*      or geometrical DOFs */
  double* Vect_dofs  /** in: dofs to be written */
  );


/**--------------------------------------------------------
  utr_create_patches - to create patches of elements containing a given node
					  (uses only mesh data - can be the same for all fields)
---------------------------------------------------------*/
extern int utr_create_patches( /** returns: >0 - Nr_patches, <=0 - failure */
  int Problem_id,       /** in: problem data structure to be used */
  utt_patches **patches_p /** in/out - array of patches for real nodes */
);

/**--------------------------------------------------------
  utr_create_patches - to create patches of elements containing a given node
					  (uses only mesh data - can be the same for all fields)
---------------------------------------------------------*/
extern int utr_create_patches_small( /** returns: >0 - Nr_patches, <=0 - failure */
  int Problem_id,       /** in: problem data structure to be used */
  utt_patches **patches_p /** in/out - array of patches for real nodes */
);

/**--------------------------------------------------------
  utr_recover_derivatives - to recover derivatives using patches of elements
					  (can be called separately for each field on the same mesh)
---------------------------------------------------------*/
extern int utr_recover_derivatives( /** returns: >0 - success, <=0 - failure */
  int Field_id,         /** in: field data structure to be used */
  int Sol_vec_id,       /*in: which solution vector to take into account */
  int Nr_patches,
  utt_patches *patches /** in - array of patches for real nodes */
		  /** out: array of patches with computed derivatives for ALL nodes */
  );

/**--------------------------------------------------------
  utr_recover_derivatives_small-to recover derivatives using patches of elements
	(does not use apr_get_stiff_mat_data for getting lists of element nodes)
					  (can be called separately for each field on the same mesh)
---------------------------------------------------------*/
extern int utr_recover_derivatives_small( /** returns: >0 - success, <=0 - failure */
  int Field_id,         /** in: field data structure to be used */
  int Sol_vec_id,       /*in: which solution vector to take into account */
  int Nr_patches,
  utt_patches *patches /** in - array of patches for real nodes */
		  /** out: array of patches with computed derivatives for ALL nodes */
  );

/**--------------------------------------------------------
  utr_adapt - to enforce default adaptation strategy for SINGLE problem/field
---------------------------------------------------------*/
extern int utr_adapt( /** returns: >0 - success, <=0 - failure */
  int Problem_id,       /** in: problem data structure to be used */
  char* Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output
  );

/**--------------------------------------------------------
utr_test_refinements - to perform sequence of random refinements/unrefinements
---------------------------------------------------------*/
extern int utr_test_refinements(  /** returns: >0 - success, <=0 - failure */
  int Problem_id,       /** in: problem data structure to be used */
  char* Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output
  );

/**--------------------------------------------------------
utr_manual_refinement - to perform manual refinement/unrefinement
---------------------------------------------------------*/
extern int utr_manual_refinement(  /** returns: >0 - success, <=0 - failure */
  int Problem_id,       /** in: problem data structure to be used */
  char* Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output
  );

/**-----------------------------------------------------------
  utr_get_list_ent - to return the list of integration entities - entities
						 for which stiffness matrices and load vectors are
						 provided by the FEM code to the solver module,
						 and DOF entities - entities with which there are dofs
						 associated by the given approximation
------------------------------------------------------------*/
extern int utr_get_list_ent( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in:  problem (and solver) identification */
  int* Nr_int_ent,    /** out: number of integration entitites */
	/** GHOST ENTITIES HAVE NEGATIVE TYPE !!! */
  int** List_int_ent_type,/** out: list of types of integration entitites */
  int** List_int_ent_id,  /** out: list of IDs of integration entitites */
	/** GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!! */
  int* Nr_dof_ent,    /** out: number of dof entities (entities with which there
					  are dofs associated by the given approximation) */
  int** List_dof_ent_type,/** out: list of types of integration entitites */
  int** List_dof_ent_id,  /** out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/** out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /** out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent/** out: maximal number of dofs per dof entity */
  );



/*------------------------------------------------------------
  utr_get_list_dof_ent - to return the list  DOF entities
						 - entities with which there are dofs
			 associated by the given approximation
						 given the list of integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module
						 (not for DG approximation!!!)
------------------------------------------------------------*/
extern int utr_get_list_dof_ent( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in:  problem (and solver) identification */
  int* Nr_int_ent,    /* in: number of integration entitites */
  int** List_int_ent_type,/* in: list of types of integration entitites */
  int** List_int_ent_id,  /* in: list of IDs of integration entitites */
	/* GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!! */
  int* Nr_dof_ent,    /* out: number of dof entities (entities with which there
			  are dofs associated by the given approximation) */
  int** List_dof_ent_type,/* out: list of types of integration entitites */
  int** List_dof_ent_id,  /* out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /* out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent/* out: maximal number of dofs per dof entity */
			  );

/*------------------------------------------------------------
  utr_get_list_int_ent - to return the list of integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module
						 (not for DG approximation!!!)
------------------------------------------------------------*/
extern int utr_get_list_int_ent( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in:  problem (and solver) identification */
  int* Nr_int_ent,    /* out: number of integration entitites */
  int** List_int_ent_type,/* out: list of types of integration entitites */
  int** List_int_ent_id  /* out: list of IDs of integration entitites */
			  );

/**-----------------------------------------------------------
  utr_get_list_ent_coarse - to return the list of integration entities - entities
						 for which stiffness matrices and load vectors are
						 provided by the FEM code to the solver module,
						 and DOF entities - entities with which there are dofs
						 associated by the given approximation for COARSE level
						 given the corresponding lists from the fine level
			  (it calls separate implementations for dg and std?)
------------------------------------------------------------*/
extern int utr_get_list_ent_coarse( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in:  problem (and solver) identification */
  int Nr_int_ent_fine,	/* in: number of integration entitites */
  int *List_int_ent_type_fine,	/* in: list of types of integration entitites */
  int *List_int_ent_id_fine,	/* in: list of IDs of integration entitites */
  int Nr_dof_ent_fine,	/* in: number of dof entities (entities with which there
			   are dofs associated by the given approximation) */
  int *List_dof_ent_type_fine,	/* in: list of types of integration entitites */
  int *List_dof_ent_id_fine,	/* in: list of IDs of integration entitites */
  int *List_dof_ent_nrdof_fine,	/* in: list of no of dofs for 'dof' entity */
  int Nrdof_glob_fine,	/* in: global number of degrees of freedom (unknowns) */
  int Max_dof_per_ent_fine,	/* in: maximal number of dofs per dof entity */
  int *Pdeg_coarse_p,	/* in: degree of approximation for coarse space */
  int *Nr_int_ent_p,	/* out: number of integration entitites */
  int **List_int_ent_type,	/* out: list of types of integration entitites */
  int **List_int_ent_id,	/* out: list of IDs of integration entitites */
  int *Nr_dof_ent_p,	/* out: number of dof entities (entities with which there
			   are dofs associated by the given approximation) */
  int **List_dof_ent_type,	/* out: list of types of integration entitites */
  int **List_dof_ent_id,	/* out: list of IDs of integration entitites */
  int **List_dof_ent_nrdofs,	/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob_p,    /* out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent_p /* out: maximal number of dofs per dof entity */
  );


/*------------------------------------------------------------
utr_get_max_num_grid_levels - for limiting nr_levels in multigrid
							  based on mesh and field data
------------------------------------------------------------*/
extern int utr_get_max_num_grid_levels(
  int Problem_id
);


/*---------------------------------------------------------
  utr_dof_ent_sons - to return info whether the entity is owned
					 and a list of dof entity sons for owned entity
---------------------------------------------------------*/
extern int utr_dof_ent_sons( /* returns: info whether the entity is owned */
  int Problem_id,     /* in: problem ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,        /* in: mesh entity ID */
  int *Ent_sons     /* out: list of dof entity sons (for owned entities) */
					 /* 	Ent_sons[0] - number of sons */
  );

/*------------------------------------------------------------
  utr_comp_stiff_mat - to create a stiffness matrix
					  and a load vector corresponding to the specified
					  mesh entity, together with information on how to
					  assemble entries into the global stiffness matrix
					  and the global load vector
------------------------------------------------------------*/
extern int utr_comp_stiff_mat(	 /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in: approximation field ID  */
  int Int_ent_type,	/* in: unique identifier of the integration entity */
  int Int_ent_id,	/* in: unique identifier of the integration entity */
  int Comp_sm,	/* in: indicator for the scope of computations: */
  /*   PDC_NO_COMP  - do not compute anything */
  /*   PDC_COMP_SM - compute entries to stiff matrix only */
  /*   PDC_COMP_RHS - compute entries to rhs vector only */
  /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,	/* in: enforced degree of polynomial (if > 0 ) */
  int *Nr_dof_ent,	/* in: size of arrays, */
			/* out: number of mesh entities with which dofs and */
			/*      stiffness matrix blocks are associated */
  int *List_dof_ent_type,	/* out: list of types for 'dof' entities */
  int *List_dof_ent_id,	/* out: list of ids for 'dof' entities */
  int *List_dof_ent_nrdofs,	/* out: list of no of dofs for 'dof' entity */
  int *Nrdofs_loc,	/* in(optional): size of Stiff_mat and Rhs_vect */
		   /* out(optional): actual number of dofs per integration entity */
  double *Stiff_mat,	/* out(optional): stiffness matrix stored columnwise */
  double *Rhs_vect,	/* outpds_elast_ls_std_intf.c(optional): rhs vector */
  char *Rewr_dofs	/* out(optional): flag to rewrite or sum up entries */
			/*   'T' - true, rewrite entries when assembling */
			/*   'F' - false, sum up entries when assembling */
 );

/**-----------------------------------------------------------
 utr_create_assemble_stiff_mat - to create element stiffness matrices
								 and assemble them to the global SM
(implementation provided separately for different multithreaded environments)
------------------------------------------------------------*/
extern int utr_create_assemble_stiff_mat(
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


/**-----------------------------------------------------------
  utr_rewr_sol - to rewrite solution from one vector to another
------------------------------------------------------------*/
extern int utr_rewr_sol(
  int Field_id,      /** in: data structure to be used  */
  int Sol_from,      /** in: ID of vector to read solution from */
  int Sol_to         /** in: ID of vector to write solution to */
  );

/*------------------------------------------------------------
  utr_renumber - to renumber (permute) graph vertices
------------------------------------------------------------*/
extern void utr_renumber(
  int Nr_vertices, /* in: number of graph vertices */
  int * Nrneig,    /* in: list of numbers of neighbors */
  int ** L_neig,   /* in: adjacency list */
  int *Permutation_array  /* out: permutation vector */
		  );


/**--------------------------------------------------------
utr_average_sol_el - to compute the average of solution over element
----------------------------------------------------------*/
extern double utr_average_sol_el( /** returns: the average of solution over element */
		int Field_id,      /** in: data structure to be used  */
		int El          /** in: element number */
			   );

/**--------------------------------------------------------
utr_dvector - to allocate a double vector: name[0..ncom-1]:
				  name=ut_dvector(ncom,error_text)
---------------------------------------------------------*/
extern double *utr_dvector( /** return: pointer to allocated vector */
	int ncom,  	/** in: number of components */
	char error_text[]/** in: error text to be printed */
	);

/**--------------------------------------------------------
utr_ivector - to allocate an integer vector: name[0..ncom-1]:
				  name=ut_ivector(ncom,error_text)
---------------------------------------------------------*/
extern int *utr_ivector(    /** return: pointer to allocated vector */
	int ncom, 	/** in: number of components */
	char error_text[]/** in: error text to be printed */
	);


/**--------------------------------------------------------
utr_d_zero - to zero a double vector of length Num
---------------------------------------------------------*/
extern void utr_d_zero(
	double* Vec, 	/** in, out: vector to initialize */
	int Num		/** in: vector length */
	);

/**--------------------------------------------------------
utr_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
extern int utr_chk_list(	/** returns: */
			/** >0 - position on the list */
					/** 0 - not found on the list */
	int Num, 	/** number to be checked */
	int* List, 	/** list of numbers */
	int Ll		/** length of the list */
	);

/**--------------------------------------------------------
utr_put_list - to put Num on the list List with length Ll
	(filled with numbers and zeros at the end)
---------------------------------------------------------*/
extern int utr_put_list( /** returns*/
		/**  >0 - position already occupied on the list */
				/**  <0 - position at which put on the list */
				/**   0 - list full, not found on the list */
	int Num, 	/** in: number to put on the list */
	int* List, 	/** in: list */
	int Ll		/** in: total list's lengths */
	);

/**--------------------------------------------------------
utr_mat3_inv - to invert a 3x3 matrix (stored as a vector!)
---------------------------------------------------------*/
extern double utr_mat3_inv(
	/** returns: determinant of matrix to invert */
	const double *mat,	/** matrix to invert */
	double *mat_inv	/** inverted matrix */
	);

//----------------------------------------------------------
//utr_vec3_add - to compute sum of 3D vectors
//----------------------------------------------------------

extern void utr_vec3_add(
	const double vec_a[3], //IN
	const double vec_b[3], //IN
	double vec_c[3]);	//OUT

//----------------------------------------------------------
//utr_vec3_subst - to compute the difference of 3D vectors
//----------------------------------------------------------

extern void utr_vec3_subst(
	const double vec_a[3],	//IN
	const double vec_b[3],	// IN
	double vec_c[3]);	// OUT

/**--------------------------------------------------------
utr_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
extern void utr_vec3_prod(
	const double* vec_a, 	/** in: vector a */
	const double* vec_b, 	/** in: vector b */
	double* vec_c	/** out: vector product axb */
	);

/**--------------------------------------------------------
utr_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
extern double utr_vec3_mxpr( /** returns: mixed product [a,b,c] */
	const double* vec_a, 	/** in: vector a */
	const double* vec_b, 	/** in: vector b */
	const double* vec_c	/** in: vector c */
	);

/**--------------------------------------------------------
utr_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
extern double utr_vec3_length(	/** returns: vector length */
	const double* vec	/** in: vector */
	);

/**--------------------------------------------------------
utr_mat3vec - to compute matrix vector product in 3D space
---------------------------------------------------------*/
extern void utr_mat3vec(
	const double* m1, 	/** in: matrix (stored by rows as a vector!) */
	const double* v1, 	/** in: vector */
	double* v2	/** out: resulting vector */
	);

/**--------------------------------------------------------
utr_mat3mat - to compute matrix matrix product in 3D space
	(all matrices are stored by rows as vectors!)
---------------------------------------------------------*/
extern void utr_mat3mat(
	const double* m1,	/** in: matrix */
	const double* m2,	/** in: matrix */
	double* m3	/** out: matrix m1*m2 */
	);


//----------------------------------------------------------
//utr_vec3_dot - to compute vector dot product of 3D vectors
//----------------------------------------------------------
extern double utr_vec3_dot(
					const double vec_a[3],
					const double vec_b[3]);

//----------------------------------------------------------
//utr_vec3_angle - to compute vector dot product of 3D vectors
//----------------------------------------------------------
extern double utr_vec3_angle(
					const double vec_a[3],
					const double vec_b[3]);


/**--------------------------------------------------------
  utr_mat_det - to compute determinant of a matrix
---------------------------------------------------------*/
extern double utr_mat_det(const double *m, int n,char store, double * det);

/**--------------------------------------------------------
utr_skip_rest_of_line - to allow for comments in input files
---------------------------------------------------------*/
extern void utr_skip_rest_of_line(
			  FILE *Fp  /** in: input file */
			   );


/**------------------------------------------------------------
utr_fprintf_double -	prints double to stream in "%f" fashion
			with specified significant digits number
-------------------------------------------------------------*/
extern int utr_fprintf_double(
  FILE* stream,	// output stream
  int sdigits,	// number of significant digits to write - put 0 for
		// full accuracy (.12) written in "%g" fashion
  double value	// double to write
  );

/**--------------------------------------------------------
utr_heap_sort - to heap-sort an array (code taken from fortran...)
---------------------------------------------------------*/
extern void utr_heap_sort(
   int    *Ind_array,    /** in/out: index array for sorting */
						 /** Ind_array[0] - the length of both arrays */
   double *Val_array     /** in: array of values used for sorting */
   );

/**--------------------------------------------------------
utr_insert_sort - to sort an array by insertion (for small integer arrays)
---------------------------------------------------------*/
extern int utr_insert_sort(
			int* A, // array
			int p,  //first index
			int k   // last index
);


//--------------------------------------------------------
//  utr_io_read_mesh - to read mesh with given filename(s)
//  Reads all mesh files from current working directory
//  that matches with passed regular expression Mesh_files_regex.
//  All files HAVE TO be the same type (ex. nas,dat,in etc.).
//  NOTE: single filename ex. "mesh.dat" is also a valid regular expression.
//--------------------------------------------------------
extern int utr_io_read_mesh(
  // returns: > 0 - number of files correctly read; <=0 - error
  const int Mesh_id, // IN: id of mesh to read to
  const char * Working_dir, // IN: directory with mesh file(s)
  const char * Mesh_files_regex, // IN: regular expression pattern
  // NOTE: using regex 'POSIX grep' standard
  const char Mesh_type // IN: type of mesh
			  );


/**--------------------------------------------------------
utr_io_export_field - to export mesh of a specified type
---------------------------------------------------------*/
extern int utr_export_field(/** returns number of files exported */
	int Field_id,    /* in: field ID */
	int Nreq,        /* in: number of equations (scalar dofs) */
	int Select,      /* in: parameter to select written vectors */
	int Accuracy,    /* in: parameter specyfying accuracy - significant digits */
			 /* (put 0 for full accuracy in "%g" format) */
	char *Filename   /* in: name of the file to write field data */
  );

/**--------------------------------------------------------
utr_io_export_mesh - to export mesh of a specified type
---------------------------------------------------------*/
extern int utr_io_export_mesh(/** returns number of files exported */
  FILE *Interactive_output, /** file or stdout to write messages */
  const int Mesh_id, // id of mesh to export
  char Mesh_type, /** symbol denoting \link  mmt_file_type  mesh_type \endlink  */
  char* Mesh_file /** mesh file name - conforming to naming convention */
  );


/**--------------------------------------------------------
utr_io_export_mesh - to display menu for exporting mesh
---------------------------------------------------------*/
extern int utr_menu_export_mesh(/** returns number of files exported */
  FILE *Interactive_output, /** file or stdout to write messages */
  FILE *Interactive_input,
  const int Mesh_id, /// id of mesh to export
  const char *Work_dir, /// Working directory
  char* Mesh_file /** mesh file name - conforming to naming convention */
  );



extern int utr_io_write_img_to_pbm(// returns 0 if all ok.
		const char * Work_dir,  // Directory to write in.
		const char * Filename,  // without extension
		const char * Comment, // written into img file header, can be NULL if no comment.
		const int Width,  // >0
		const int Height, // >0
		const int Max_color_component_value,
		const int Magic_number, // {4,5,6}
		const unsigned char* Img_data, // pointer to the array with image data
		FILE* Interactive_output
		);






int utr_io_write_mesh_info_to_PAM(int Mesh_id
		, const char *Work_dir, const char *Filename, const char *Comment, FILE *Interactive_output);

// utr_io_write_img_to_pam - to write data into PAM image file.
// http://en.wikipedia.org/wiki/Netpbm
// Allowed combinations of parameters:
// -----------------------------------------------------------
// TUPLTYPE           |MAXVAL |DEPTH  |comment
// -----------------------------------------------------------
// BLACKANDWHITE        1       1       special case of GRAYSCALE
// GRAYSCALE            2…65535	1       2 bytes per pixel for MAXVAL > 255
// RGB                  1…65535	3       6 bytes per pixel for MAXVAL > 255
// BLACKANDWHITE_ALPHA	1       2       2 bytes per pixel
// GRAYSCALE_ALPHA      2…65535	2       4 bytes per pixel for MAXVAL > 255
// RGB_ALPHA            1…65535	4       8 bytes per pixel for MAXVAL > 255s of parameters:
// -----------------------------------------------------------
int utr_io_write_img_to_pam( // returns 0 if all ok.
		const char * Work_dir,  // Directory to write in.
		const char * Filename,  // without extension
		const int Width,  // >0
		const int Height, // >0
		const int Depth,  // <1:4>
		const int Maxval, // <1:65535>,
		const char* TUPLTYPE,
		const char* Img_data, // pointer to the array with image data
		FILE* Interactive_output
		);

// utr_io_write_img_to_pnm - to write data into PNM image file.
// The Portable Bit/Grey/PixMap formats PBM, PGM, PPM.
// They are collectively referred to as PNM (Portable aNy Map).
// http://en.wikipedia.org/wiki/Netpbm_format
// Allowed combinations:
// -----------------------------------------------------------
// Type                 Magic number	Extension	Colors
// -----------------------------------------------------------
// Portable BitMap		P4	binary      .pbm        0–1 (black & white)
// Portable GrayMap		P5	binary      .pgm        0–255 (gray scale)
// Portable PixMap		P6	binary      .ppm        0–255 (RGB)
// -----------------------------------------------------------
int utr_io_write_img_to_pnm(// returns 0 if all ok.
		const char * Work_dir,  // Directory to write in.
		const char * Filename,  // without extension
		const char * Comment, // written into img file header, can be NULL if no comment.
		const int Width,  // >0
		const int Height, // >0
		const int Max_color_component_value, // see Colors above
		const int Magic_number, // {4,5,6}
		const char* Img_data, // pointer to the array with image data
		FILE* Interactive_output
		);


int utr_color_int_ent_for_assembly(
  int Problem_id,           /* in: Problem id */
  int Level_id,             /* in: Level id */
  int nr_elems,             /* in: number of elements */
  int nr_faces,             /* in: number of faces */
  int* L_int_ent_type,      /* in: integration entities type,
							   out: coloured  */
  int* L_int_ent_id,        /* in: integration entities id,
							   out: coloured  */
  int nr_dof_ent,           /* number of dof entities dof_ents */
  int* nr_int_ent_loc,      /* number of int_ents for each dof_ent */
  int** l_int_ent_index,    /* list of indices in int_ent table for each dof_ent */
  int* Nr_colors_elems,     /* out: number of colours for elements  */
  int** L_color_index_elems,/* out: crs table with index of L_int_ent.. where colour starts and ends  */
  int* Nr_colors_faces,     /* out: number of colours for faces  */
  int** L_color_index_faces /* out: crs table with index of L_int_ent.. where colour starts and ends  */
				   );


// UTR_IO_RESULT (file with infomrations) should not be here!
// It is in uth_io_files.h

///
/// \brief utr_ctrl_pts_init - to initialize 'control points' functionality.
/// It is coupled with utr_io_results_... functions, but to simplify useage this function wraps up
/// several setup steps to make it work properly.
/// \param Work_dir - working directory in which Filename will be searched and Target_ctrl_pts_filename will be saved.
/// \param Filename - name of file from which 'control_points' block will be loaded (e.g. problem.dat).
/// \param N_solutions - number of values from all fields wich will be stored at each control point in each time step.
/// Usually this could be the sum of \link apr_get_nreq \endlink for all desired fields.
/// \param Target_ctrl_pts_filename - name of file in which control points values will be stored (e.g. control_points.csv).
/// \return
///
int utr_ctrl_pts_init(const char* Work_dir, const char* Filename, const int N_solutions,
					  const char* Target_ctrl_pts_filename);
///
/// \brief utr_ctrl_pts_add_values - this function computes values at previously defined
/// control points and stores them until they will we written to the file when
/// \link utr_io_result_write_values_and_proceed \endlink will be called.
/// NOTE that before this function can be called, the \link utr_ctrl_pts_init \endlink must be called!
///  This function can be called several times, to store more values from different fields,
/// but no more then first N_solutions given in \link utr_ctrl_pts_init \endlink will be saved to file.
/// \param Field_id - ID of field to compute values from at control points.
/// \return
///
int utr_ctrl_pts_add_values(const int Field_id);

int utr_generate_int_ent_coloring( // return number of colors
		const int nr_int_ent, // IN: mesh id to color
		int ** elem2elem, // IN: array entity neighborhood
		int *Elem_colors);  // OUT: array of size number of active elems containing colors (numbers)


///////////////////////////////////////////////////////


/** @} */ // end of group


#ifdef __cplusplus
}
#endif

#endif
