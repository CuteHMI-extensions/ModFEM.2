/************************************************************************
File mms_prism_intf.c - implementation of interface routines for
						meshes of prismatic elements

Contains definitions of interface routines:
  mmr_init_mesh - to initialize the mesh data structure (and read data)
  mmr_export_mesh - to export (write to a file) mesh data in a specified format
  mmr_test_mesh - to test the integrity of mesh data
  mmr_get_nr_elem - to return the number of active elements
  mmr_get_max_elem_id - to return the maximal element id
  mmr_get_max_elem_max - to return the maximal possible element id
  mmr_get_next_act_elem - to return the next active element id
  mmr_get_next_elem_all - to return the next element id (active or inactive)

  mmr_get_nr_face - to return the number of active faces
  mmr_get_max_face_id - to return the maximal face id
  mmr_get_next_act_face - to return the next active face id
  mmr_get_next_face_all - to return the next face id (active or inactive)
  mmr_get_nr_edge - to return the number of active edges
  mmr_get_max_edge_id - to return the maximal edge id
  mmr_get_nr_node - to return the number of active nodes
  mmr_get_max_node_id - to return the maximal node id
  mmr_get_max_node_max - to return the maximal possible vertex id

  mmr_get_max_gen - to get maximal allowed generation level for elements
  mmr_set_max_gen - to set maximal allowed generation level for elements
  mmr_get_max_gen_diff - to get maximal allowed generation difference between
						 neighboring elements
  mmr_set_max_gen_diff - to set maximal allowed generation difference between
						 neighboring elements

  mmr_write_mesh - to dump-out mesh data
  mmr_init_ref - to initialize the process of refinement
  mmr_refine_el - to refine an element
  mmr_derefine_el - to derefine an element
  mmr_final_ref - to finalize the process of refinement after rewriting DOFs
  mmr_free_mesh - to free space allocated for mesh data structure
  mmr_get_fa_el_bc_connect(int fa,int &el) - return 0 if normal bc else return connected face (el return el_id)
  mmr_init_mesh2


------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* interface for the mesh manipulation module */
#include <modfem/mmh_intf.h>

/* mesh manipulation data structure and headers for internal routines */
#include <modfem/mm_prism/mmh_prism.h>

#include <modfem/uth_log.h>
#include <modfem/uth_hash.h>

/* Initialization of constants */
const int MMC_MAX_EDGE_ELEMS = 20; /* max number of elements containg an edge */
const int MMC_CUR_MESH_ID = 0; /* indicator for the current active mesh */

const int MMC_AUTO_GENERATE_ID = 0;

/* Refinement options */
const int MMC_DO_UNI_REF   = -1; /* indicator to perform uniform refinement */
const int MMC_DO_UNI_DEREF = -2; /* indicator to perform uniform derefinement */

/* Status indicators */
const int MMC_ACTIVE        = 1;   /* active mesh entity */
const int MMC_INACTIVE      = -1;   /* inactive (refined) mesh entity */
const int MMC_FREE          = 0;   /* free space in data structure */

/* Monitoring options */
const int MMC_PRINT_NOT     = 0;  /* indicator not to print anything */
const int MMC_PRINT_ERRORS  = 1;  /* to print error messages only */
const int MMC_PRINT_INFO    = 2;  /* to print most important information */
const int MMC_PRINT_ALLINFO = 3;  /* to print all available information */

/* Position and neighbors options */
const int MMC_BOUNDARY      = 0;   /* boundary indicator */
const int MMC_BIG_NGB       = -1;  /* big neghbor indicator */
const int MMC_SUB_BND       = -2;  /* inter-subdomain boundary indicator */

/* Refinement types */
const int MMC_NOT_REF       = 0;   /* not refined */
const int MMC_REF_ISO       = 1;   /* isotropic refinement */
const int MMC_REF_ANI       = 2;   /* anisotropic refinement */

/* Other */
const int MMC_INIT_GEN_LEVEL = 0;  // initial elemenst generation number
const int MMC_NO_FATH       = 0;   /* no father indicator */
const int MMC_SAME_ORIENT   = 1;   /* indicator for the same orientation */
const int MMC_OPP_ORIENT   = -1;   /* indicator for the opposite orientation */

const int MMC_FACE_NODES_FOR_TETRA[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
const int MMC_FACE_NODES_FOR_PRISM[5][4] = {{0, 1, 2, -1}, {3, 4, 5, -1}, {0, 1, 3, 4}, {1, 2, 4, 5}, {2, 0, 5, 3}};

const int MMC_IS_SUPPORTING_NEW_MPI = 0;

/*** GLOBAL VARIABLES for the whole module ***/
int       mmv_nr_meshes = 0;   /* the number of meshes in the problem */
int       mmv_cur_mesh_id = 0;              /* ID of the current mesh */
mmt_mesh  mmv_meshes[MMC_MAX_NUM_MESH];            /* array of meshes */

/*------------------------------------------------------------
  mmr_module_introduce - to return the mesh name
------------------------------------------------------------*/
int mmr_module_introduce(
		/* returns: >=0 - success code, <0 - error code */
		char * Mesh_name /* out: the name of the mesh */
)
{
	char * string = "3D_PRISM";

	strcpy(Mesh_name, string);

	return (1);
}

/**-----------------------------------------------------------
  mmr_module_internals - to return the mesh internal structures
------------------------------------------------------------*/
void * mmr_module_internals(
		/** returns: >=0 - success code, <0 - error code */
		int Mesh_id/** out: the name of the mesh */
)
{
	return NULL;
}

/*------------------------------------------------------------
  mmr_init_mesh - to initialize the mesh data structure (and read data)
------------------------------------------------------------*/
int mmr_init_mesh(  /* returns: >0 - Mesh ID, <0 - error code */
		int Control,	    /* in: control variable to choose data format */
		/* MMC_MOD_FEM_MESH_DATA - generic mesh for this module */
		/* MMC_MOD_FEM_PRISM_DATA - prismatic mesh */
		/* MMC_GRADMESH_DATA - mesh produced by 2D GRADMESH mesh generator */
		char * Filename,   /* in: name of the file to read mesh data */
		FILE * interactive_output
)
{
	/* auxiliary variables */
	int iaux, ned;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* increase the counter for meshes */
	mmv_nr_meshes++;

	/* set the current mesh number */
	mmv_cur_mesh_id = mmv_nr_meshes;

	if (Control == MMC_MOD_FEM_PRISM_DATA || Control == MMC_MOD_FEM_MESH_DATA) {

		/* dump-in data stored by previous runs */
		iaux = mmr_read_mesh(mmv_cur_mesh_id, Filename);

	}
	else if (Control == MMC_GRADMESH_DATA) {

		/* read 2D unstructured mesh data and generate half-structured 3D mesh */
		iaux = mmr_import_mesh_grad(mmv_cur_mesh_id, Filename);
	}
	else {
		printf("Unknown mesh type in mmr_init_mesh.... ! Exiting\n");
		exit(-1);
	}

	/* allocate the space for elems structure for active and inactive edges */
	/* parameter Max_edge_id==0 - there are no structures to free !!! */
	mmr_create_edge_elems(mmv_cur_mesh_id, 0);

	if (iaux >= 0) return (mmv_cur_mesh_id);
	else return (iaux);
}

/*---------------------------------------------------------
  mmr_export_mesh - to export (write to a file) mesh data in a specified format
---------------------------------------------------------*/
extern int mmr_export_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Control,	  /* in: control variable to choose data format */
		/* MMC_MOD_FEM_MESH_DATA - dump files in standard format for this module */
		/* MMC_MOD_FEM_PRISM_DATA - dump files in standard format for prisms */
		char * Filename /* in: name of the file to dump mesh data */
)
{

	/* auxiliary variables */
	int iaux;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (Control == MMC_MOD_FEM_PRISM_DATA || Control == MMC_MOD_FEM_MESH_DATA) {

		/* dump-out data in the standard MOD_FEM format */
		iaux = mmr_write_mesh(Mesh_id, Filename);

	}

	return (iaux);
}

/*---------------------------------------------------------
  mmr_test_mesh - to test the integrity of mesh data
---------------------------------------------------------*/
int mmr_test_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* for each active element */

	/* get the list of faces */

	/* get the list of edges */

	/* get the list of vertices */

	/* for each element's face */

	/* check whether element is a face's neighbor */

	/* for each face's edge */

	/* check whether it is element's edge */

	/* for each face's vertex */

	/* check it is element's vertex */

	/* for each element's edge */

	/* for each edge's endpoint */

	/* check it is element's vertex */

	/* etc. */

	return (0);
}




/*---------------------------------------------------------
  mmr_mesh_i_params - to return mesh parameters
---------------------------------------------------------*/
int mmr_mesh_i_params(  /* returns: >=0 - integer mesh parameter,
									 <0 - error code  */
		int Mesh_id,  /* in: pointer to the mesh data structure */
		int Num          /* in: parameter number in control structure */
)
{
	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mmr_get_mesh_i_params(mesh, Num));
}

/*------------------------------------------------------------
mmr_split_into_blocks_add_contact
------------------------------------------------------------*/
int mmr_split_into_blocks_add_contact(
		const int Mesh_id,
		const char * workdir,
		int * tabMat,
		int l_mat,
		int * tabBlock,
		int l_block,
		int * war,
		int l_war,
		double * tempBlock,
		int l_tempBlock,
		int * groups_fluid,
		int l_groups_fluid) {

	printf("Method mmr_split_into_blocks_add_contact not implemented in mmd_prism!\n");
	exit(-1);
}

/*---------------------------------------------------------
  mmr_get_nr_elem - to return the number of active elements
---------------------------------------------------------*/
int mmr_get_nr_elem(/* returns: >=0 - number of active elements, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/*kbw
	printf("In get_nr_elem: Mesh_id %d, mesh %u, Nrel %d\n",
	 Mesh_id, mesh, mesh->parm.nrel);
	kew*/

	return (mesh->parm.nrel);
}


/*---------------------------------------------------------
  mmr_get_max_elem_id - to return the maximal element id
---------------------------------------------------------*/
int mmr_get_max_elem_id(  /* returns: >=0 - maximal element id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.nmel);
}

/*---------------------------------------------------------
  mmr_get_next_act_elem - to return the next active element id
---------------------------------------------------------*/
int mmr_get_next_act_elem( /* returns: >=0 - ID of the next active element */
		/*                (0 - input is the last element) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nel       /* in: input element (0 - return first active element) */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int iel, nmel;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);
	nmel = mesh->parm.nmel;

	/*kbw
	printf("In get_next_act_elem: Mesh_id %d, mesh %u, nmel %d, Nel %d\n",
	 Mesh_id, mesh, nmel, Nel);
	kew*/

	if (Nel < 0) return (-1);

	iel = Nel + 1;
	while (iel <= nmel && mmr_el_status(Mesh_id, iel) != MMC_ACTIVE) {
		iel++;
	}

	if (iel <= nmel) {
		/*kbw
		printf("returning nel=%d\n",iel);
		kew*/
		return (iel);
	}
	else if (iel == nmel + 1) return (0);
	else return (-1);
}

/*---------------------------------------------------------
  mmr_get_next_elem_all - to return the next element id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_elem_all( /* returns: >=0 - ID of the next element */
		/*                (0 - input is the last element) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nel       /* in: input element (0 - return first element) */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int iel, nmel;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);
	nmel = mesh->parm.nmel;

	if (Nel < 0) return (-1);

	iel = Nel + 1;
	while (iel <= nmel && mmr_el_status(Mesh_id, iel) == MMC_FREE) {
		iel++;
	}

	if (iel <= nmel) return (iel);
	else if (iel == nmel + 1) return (0);
	else return (-1);
}

void mmr_get_vec_norm_at_node(int Mesh_id, int node_id, double * vec_norm) {

	vec_norm[0] = -10000.0;
	//mmv_out_stream << "Method not implemented!: mmr_get_vec_norm_at_node";
	//exit(-1);

}

/*---------------------------------------------------------
  mmr_get_nr_face - to return the number of active faces
---------------------------------------------------------*/
int mmr_get_nr_face(/* returns: >=0 - number of active faces, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.nrfa);
}


/*---------------------------------------------------------
  mmr_get_max_face_id - to return the maximal face id
---------------------------------------------------------*/
int mmr_get_max_face_id(  /* returns: >=0 - maximal face id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.nmfa);
}

/*---------------------------------------------------------
  mmr_get_next_act_face - to return the next active face id
---------------------------------------------------------*/
int mmr_get_next_act_face( /* returns: >=0 - ID of the next active face */
		/*                (0 - input is the last face) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nfa       /* in: input face (0 - return first active face) */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int ifa, nmfa;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);
	nmfa = mesh->parm.nmfa;

	if (Nfa < 0) return (-1);

	ifa = Nfa + 1;
	while (ifa <= nmfa && mmr_fa_status(Mesh_id, ifa) != MMC_ACTIVE) {
		ifa++;
	}

	if (ifa <= nmfa) return (ifa);
	else if (ifa == nmfa + 1) return (0);
	else return (-1);
}

/*---------------------------------------------------------
  mmr_get_next_face_all - to return the next face id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_face_all( /* returns: >=0 - ID of the next face */
		/*                (0 - input is the last face) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nfa       /* in: input face (0 - return first face) */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int ifa, nmfa;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);
	nmfa = mesh->parm.nmfa;

	if (Nfa < 0) return (-1);

	ifa = Nfa + 1;
	while (ifa <= nmfa && mmr_fa_status(Mesh_id, ifa) == MMC_FREE) {
		ifa++;
	}

	if (ifa <= nmfa) return (ifa);
	else if (ifa == nmfa + 1) return (0);
	else return (-1);
}


/*---------------------------------------------------------
  mmr_get_nr_edge - to return the number of active edges
---------------------------------------------------------*/
int mmr_get_nr_edge(/* returns: >=0 - number of active edges, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.nred);
}


/*---------------------------------------------------------
  mmr_get_max_edge_id - to return the maximal edge id
---------------------------------------------------------*/
int mmr_get_max_edge_id(  /* returns: >=0 - maximal edge id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.nmed);
}

/*---------------------------------------------------------
  mmr_get_next_edge_all - to return the next edgeent id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_edge_all( /* returns: >=0 - ID of the next edge */
		/*                (0 - input is the last edge) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ned       /* in: input edge (0 - return first edge) */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int ied, nmed;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);
	nmed = mesh->parm.nmed;

	if (Ned < 0) return (-1);

	ied = Ned + 1;
	while (ied <= nmed && mmr_edge_status(Mesh_id, ied) == MMC_FREE) {
		ied++;
	}

	if (ied <= nmed) return (ied);
	else if (ied == nmed + 1) return (0);
	else return (-1);
}


/*---------------------------------------------------------
  mmr_get_nr_node - to return the number of active nodes
---------------------------------------------------------*/
int mmr_get_nr_node(/* returns: >=0 - number of active nodes, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.nrno);
}


/*---------------------------------------------------------
  mmr_get_max_node_id - to return the maximal node id
---------------------------------------------------------*/
int mmr_get_max_node_id(  /* returns: >=0 - maximal node id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.nmno);
}


/*---------------------------------------------------------
  mmr_get_next_node_all - to return the next node id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_node_all( /* returns: >=0 - ID of the next node */
		/*                (0 - input is the last node) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nno       /* in: input node (0 - return first node) */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int ino, nmno;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);
	nmno = mesh->parm.nmno;

	if (Nno < 0) return (-1);

	ino = Nno + 1;
	while (ino <= nmno && mmr_node_status(Mesh_id, ino) == MMC_FREE) {
		ino++;
	}

	if (ino <= nmno) return (ino);
	else if (ino == nmno + 1) return (0);
	else return (-1);
}


/*---------------------------------------------------------
  mmr_set_max_gen - to set maximal generation level for element refinements
------------------------------------------------------------*/
int mmr_set_max_gen(/* returns: >=0-success code, <0-error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Max_gen   /* in: maximal generation difference */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	mesh->parm.maxgen = Max_gen;

	return (0);
}

/*---------------------------------------------------------
  mmr_get_max_gen - to get maximal allowed generation level for elements
------------------------------------------------------------*/
int mmr_get_max_gen(/* returns: >=0-maximal generation
								 <0-error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.maxgen);
}


/*---------------------------------------------------------
  mmr_get_max_gen_diff - to get maximal allowed generation difference between
	neighboring elements
------------------------------------------------------------*/
int mmr_get_max_gen_diff(/* returns: >=0-maximal generation difference,
									  <0-error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.maxgendiff);
}

/*---------------------------------------------------------
  mmr_set_max_gen_diff - to set maximal generation difference between
	neighboring elements
------------------------------------------------------------*/
int mmr_set_max_gen_diff(/* returns: >=0-success code, <0-error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Max_gen_diff    /* in: maximal generation difference */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	mesh->parm.maxgendiff = Max_gen_diff;

	return (0);
}


/*------------------------------------------------------------
  mmr_init_ref - to initialize the process of refinement
------------------------------------------------------------*/
int mmr_init_ref(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	int iaux;

	/*++++++++++++++++ executable statements ++++++++++++++++*/



	return (0);
}


/*------------------------------------------------------------
  mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready
  for dofs projection // Added 06.2011 for compatibility reasons
------------------------------------------------------------*/
int mmr_is_ready_for_proj_dof_ref(
		int Mesh_id
)
{
	return 1;
}


/*------------------------------------------------------------
  mmr_refine_el - to refine an element WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_refine_el( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El   /* in: element ID  */
)
{
	return (mmr_divide_el8_p(Mesh_id, El));
}

/*------------------------------------------------------------
mmr_derefine_el - to derefine an element WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_derefine_el( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID  */
)
{
	return (mmr_clust_el8_p(Mesh_id, El));
}

/*------------------------------------------------------------
  mmr_refine_mesh - to refine the WHOLE mesh WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_refine_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	/* local variables */
	mmt_mesh * mesh;
	int i, iaux, iel, nmel_old, nrel_old, nrel_div, nelref, nrwait, listwait[100];
	int el_type, gen_el, max_gen;
	int ifa, num_face, face, face_sons[5], ifaneig;
	int iprint = MMC_PRINT_INFO + 1;
	int * list_temp;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	max_gen = mmr_get_max_gen(Mesh_id);

	nrel_div = 0;
	nmel_old = mesh->parm.nmel;
	nrel_old = mesh->parm.nrel;
	list_temp = (int *)malloc(nrel_old * sizeof(int));
	if (!list_temp) {
		printf("Not enough space for allocating list_temp ! Exiting\n");
		exit(-1);
	}

#ifdef DEBUG_MMM
	if (iprint > MMC_PRINT_ERRORS) {
		printf("\nStarting uniform refinement of the mesh.\n");
		printf("Parameters (number of active, maximal index, pointer to first free, dimension):\n");
		printf("Elements: nrel %d, nmel %d, pfel %d, mxel %d\n",
				mesh->parm.nrel, mesh->parm.nmel,
				mesh->parm.pfel, mesh->parm.mxel);
		printf("Faces:    nrfa %d, nmfa %d, pffa %d, mxfa %d\n",
				mesh->parm.nrfa, mesh->parm.nmfa,
				mesh->parm.pffa, mesh->parm.mxfa);
		printf("Edges:    nred %d, nmed %d, pfed %d, mxed %d\n",
				mesh->parm.nred, mesh->parm.nmed,
				mesh->parm.pfed, mesh->parm.mxed);
		printf("Nodes:    nrno %d, nmno %d, pfno %d, mxno %d\n",
				mesh->parm.nrno, mesh->parm.nmno,
				mesh->parm.pfno, mesh->parm.mxno);
		\
	}
#endif


	for (iel = 1; iel <= nmel_old; iel++) {

		if (mmr_el_status(Mesh_id, iel) == MMC_ACTIVE) {

			/* barrier */
			gen_el = mmr_el_gen(Mesh_id, iel);

			if (gen_el < max_gen) {

				list_temp[nrel_div] = iel;
				nrel_div++;

			}
		}
	}

	/*kbw
	  printf("Elements to refine: ");
	  for(iel=0;iel<nrel_div;iel++){
	  printf("%d, ",list_temp[iel]);
	  }
	  printf("\n");
	  kew*/


	for (iel = 0; iel < nrel_div; iel++) {
		iaux = mmr_divide_el8_p(Mesh_id, list_temp[iel]);
		if (iaux < 0) break;
	}

	free(list_temp);

#ifdef DEBUG_MMM
	if (iprint > MMC_PRINT_ERRORS && iaux > -1) {
		printf("\nAfter uniform refinement of the mesh.\n");
		printf("Parameters (number of active, maximal index, pointer to first free, dimension\n");
		printf("Elements: nrel %d, nmel %d, pfel %d, mxel %d\n",
				mesh->parm.nrel, mesh->parm.nmel,
				mesh->parm.pfel, mesh->parm.mxel);
		printf("Faces:    nrfa %d, nmfa %d, pffa %d, mxfa %d\n",
				mesh->parm.nrfa, mesh->parm.nmfa,
				mesh->parm.pffa, mesh->parm.mxfa);
		printf("Edges:    nred %d, nmed %d, pfed %d, mxed %d\n",
				mesh->parm.nred, mesh->parm.nmed,
				mesh->parm.pfed, mesh->parm.mxed);
		printf("Nodes:    nrno %d, nmno %d, pfno %d, mxno %d\n",
				mesh->parm.nrno, mesh->parm.nmno,
				mesh->parm.pfno, mesh->parm.mxno);
	}
#endif

	return (iaux);

}

/*------------------------------------------------------------
mmr_derefine_mesh - to derefine the WHOLE mesh WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_derefine_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	/* local variables */
	int i, iaux, iel, nmel_old, nrel_old, nrel_div;
	int iprint = 5;
	mmt_mesh * mesh;
	int * list_temp;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	nrel_div = 0;
	nmel_old = mesh->parm.nmel;
	nrel_old = mesh->parm.nrel;
	list_temp = (int *)malloc(nrel_old * sizeof(int));
	if (!list_temp) {
		printf("Not enough space for allocating list_temp ! Exiting\n");
		exit(-1);
	}

#ifdef DEBUG_MMM
	if (iprint > MMC_PRINT_ERRORS) {
		printf("\nStarting uniform derefinement of the mesh.\n");
		printf("Parameters (number of active, maximal index, pointer to first free, dimension):\n");
		printf("Elements: nrel %d, nmel %d, pfel %d, mxel %d\n",
				mesh->parm.nrel, mesh->parm.nmel,
				mesh->parm.pfel, mesh->parm.mxel);
		printf("Faces:    nrfa %d, nmfa %d, pffa %d, mxfa %d\n",
				mesh->parm.nrfa, mesh->parm.nmfa,
				mesh->parm.pffa, mesh->parm.mxfa);
		printf("Edges:    nred %d, nmed %d, pfed %d, mxed %d\n",
				mesh->parm.nred, mesh->parm.nmed,
				mesh->parm.pfed, mesh->parm.mxed);
		printf("Nodes:    nrno %d, nmno %d, pfno %d, mxno %d\n",
				mesh->parm.nrno, mesh->parm.nmno,
				mesh->parm.pfno, mesh->parm.mxno);
		\
	}
#endif

	for (iel = 1; iel <= nmel_old; iel++) {
		if (mmr_el_status(Mesh_id, iel) == MMC_ACTIVE) {

			list_temp[nrel_div] = iel;
			nrel_div++;

		}
	}

	for (iel = 0; iel < nrel_div; iel++) {

		if (mmr_el_status(Mesh_id, list_temp[iel]) == MMC_ACTIVE) {
			/* if element still active */
			iaux = mmr_clust_el8_p(Mesh_id, list_temp[iel]);
			if (iaux < 0) break;
		}
	}

	free(list_temp);

#ifdef DEBUG_MMM
	if (iprint > MMC_PRINT_ERRORS && iaux > -1) {
		printf("\nAfter uniform derefinement of the mesh.\n");
		printf("Parameters (number of active, maximal index, pointer to first free, dimension\n");
		printf("Elements: nrel %d, nmel %d, pfel %d, mxel %d\n",
				mesh->parm.nrel, mesh->parm.nmel,
				mesh->parm.pfel, mesh->parm.mxel);
		printf("Faces:    nrfa %d, nmfa %d, pffa %d, mxfa %d\n",
				mesh->parm.nrfa, mesh->parm.nmfa,
				mesh->parm.pffa, mesh->parm.mxfa);
		printf("Edges:    nred %d, nmed %d, pfed %d, mxed %d\n",
				mesh->parm.nred, mesh->parm.nmed,
				mesh->parm.pfed, mesh->parm.mxed);
		printf("Nodes:    nrno %d, nmno %d, pfno %d, mxno %d\n",
				mesh->parm.nrno, mesh->parm.nmno,
				mesh->parm.pfno, mesh->parm.mxno);
	}
#endif


	return (iaux);
}


/*------------------------------------------------------------
  mmr_final_ref - to finalize the process of refinement after rewriting DOFs
------------------------------------------------------------*/
int mmr_final_ref(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	/* local variables */
	mmt_mesh * mesh;
	int i, iel, nmel, num_sons, ison, son;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/* get number of elements */
	//i=14; nmel=mmr_get_mesh_i_params(mesh,i);
	nmel = mmr_get_max_elem_id(Mesh_id);

	/* in a loop over elements */
	for (iel = 1; iel <= nmel; iel++) {

		/* if element is active and still has sons */
		if (mmr_el_status(Mesh_id, iel) == MMC_ACTIVE) {
			if (mesh->elem[iel].refi == MMC_REF_ISO) {
				num_sons = 8;
				/* delete sons */
				for (ison = 0; ison < num_sons; ison++) {
					son = mesh->elem[iel].sons[ison];
					mesh->elem[son].fath = MMC_NO_FATH;
				}
				free(mesh->elem[iel].sons);
				mesh->elem[iel].refi = MMC_NOT_REF;
				mesh->elem[iel].sons = NULL;
			}
		}

	}

	return (0);
}


/*---------------------------------------------------------
  mmr_free_mesh - to free space allocated for mesh data structure
---------------------------------------------------------*/
int mmr_free_mesh(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int istr, i, nr_sol;
	int max_edge_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/* deallocate the space for elems structure for active and inactive edges */
	max_edge_id = mmr_get_max_edge_id(Mesh_id);
	mmr_delete_edge_elems(mmv_cur_mesh_id, max_edge_id); //Fix big memory leak

	free(mesh->node);
	free(mesh->edge);

	for (istr = 1; istr <= mesh->parm.nmfa; istr++) {
		if (mesh->face[istr].type < 0) {
			free(mesh->face[istr].sons);
		}
	}
	free(mesh->face);

	for (istr = 1; istr <= mesh->parm.nmel; istr++) {
		if (mesh->elem[istr].type != MMC_FREE) {
			free(mesh->elem[istr].face);
		}
		if (mesh->elem[istr].type < 0) {
			free(mesh->elem[istr].sons);
		}
	}
	free(mesh->elem);

	return (0);
}


/*---------------------------------------------------------
  mmr_get_max_elem_max - to return the maximal possible element id
---------------------------------------------------------*/
int mmr_get_max_elem_max(  /* returns: >=0 - maximal possible element id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.mxel);
}


/*---------------------------------------------------------
  mmr_get_max_face_max - to return the maximal possible vertex id
---------------------------------------------------------*/
int mmr_get_max_face_max(  /* returns: >=0 - maximal possible vertex id */
		/*           <0 - error code */
		int Mesh_id /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.mxfa);

}
/*---------------------------------------------------------
  mmr_get_max_edge_max - to return the maximal possible vertex id
---------------------------------------------------------*/
int mmr_get_max_edge_max(  /* returns: >=0 - maximal possible vertex id */
		/*           <0 - error code */
		int Mesh_id /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.mxed);

}
/*---------------------------------------------------------
  mmr_get_max_node_max - to return the maximal possible vertex id
---------------------------------------------------------*/
int mmr_get_max_node_max(  /* returns: >=0 - maximal possible vertex id */
		/*           <0 - error code */
		int Mesh_id /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->parm.mxno);

}

/*------------------------------------------------------------
  compatibility additions for hybris module
------------------------------------------------------------*/

//------------------------------------------------------------
//  mmr_init_mesh2 - to initialize the mesh data structure
//------------------------------------------------------------
// NOTE: this is NOT a resize function. It does NOT ADD any elements.
// It only creates new mesh and reserves resources for it.
// To add mesh entities use mmr_add_* functions.
//------------------------------------------------------------
extern int mmr_init_mesh2(  /* returns: >0 - Mesh ID, <0 - error code */
		FILE * Interactive_output, // in: name of the output file to write out messages
		const int N_nodes, //IN: target count of nodes
		const int N_edges, //IN: target count of edges
		const int N_faces, //IN: target count of faces
		const int N_elems  //IN: target count of elements
)
{
	fprintf(Interactive_output, "%s",
			"mmr_init_mesh2 function is not implemented in mmd_prism (check whetherMODFEM_NEW_MPI is set to FALSE, as it should be for prismatic meshes). Aborting.\n");
	exit(-1);
}


int mmr_finish_read( // return >=0 - success, <0 - error code
		const int Mesh_id // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
)
{
	mf_fatal_err("mmr_finish_read function is not implemented in mmd_prism.");
}

//---------------------------------------------------------
// mmr_reserve - to inform mesh module about new target number of elems,faces,edges and nodes
//------------------------------------------------------------
// NOTE: this is NOT a resize function. It does NOT ADD/REMOVE any elements.
// It only reserves resources.
// To add mesh entities use mmr_add_* functions.
//---------------------------------------------------------
extern int mmr_init_read( // return >=0 - success, <0 - error code
		const int Mesh_id, // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int N_nodes, //IN: target count of nodes
		const int N_edges, //IN: target count of edges
		const int N_faces, //IN: target count of faces
		const int N_elems  //IN: target count of elements
)
{
	mf_fatal_err("This function is not yet implemented. Aborting.");

	return -1;

}

//---------------------------------------------------------
//  mmr_add_elem - to add (append) new element
//---------------------------------------------------------
extern int mmr_add_elem(// returns: added element id
		const int Mesh_id, // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
		const int El_type, // IN: type of new element
		const int El_nodes[6], //IN: if known: nodes of new element, otherwise NULL
		const int El_faces[5],  //IN: if known: faces of new element, otherwise NULL
		const int Material_idx
)
// NOTE: if both El_nodes and El_faces are NULLs
// then an element is udefined and an error is returned!
{
	printf("%s", "This function is not yet implemented. Aborting.");
	exit(-1);
	return -1;

}

//---------------------------------------------------------
//  mmr_add_face - to add (append) new face
//---------------------------------------------------------
extern int mmr_add_face(// returns: added face id
		const int Mesh_id,  // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
		const int Fa_type,  // IN: type of new face
		const int Fa_nodes[4], // IN: if known: nodes of new face, otherwise NULL
		const int Fa_edges[4], // IN: if known: edges of new face, otherwise NULL
		const int Fa_neigs[2],  // IN: if known: neigs of new face, otherwise NULL
		const int B_cond_val
)
// NOTE: if both Fa_nodes and Fa_edges are NULLs
// then a face is udefined and an error is returned!
{
	printf("%s", "This function is not yet implemented. Aborting.");
	exit(-1);
	return -1;

}
//---------------------------------------------------------
//  mmr_add_edge - to add (append) new edge
//---------------------------------------------------------
extern int mmr_add_edge(// returns: added edge id
		const int Mesh_id,  //IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
		const int nodes[2]  //IN: nodes for new edge
)
// NOTE: if nodes==NULL edge is undefined and an error is returned!
{
	printf("%s", "This function is not yet implemented. Aborting.");
	exit(-1);
	return -1;

}

//---------------------------------------------------------
//  mmr_add_node - to add (append) new node
//---------------------------------------------------------
extern int mmr_add_node( // returns: added node id
		const int Mesh_id,  //IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
		const double Coords[3] //IN: geometrical coordinates of node
)
// NOTE: if Coords==NULL node is undefined and an error is returned!
{
	printf("%s", "This function is not yet implemented. Aborting.");
	exit(-1);
	return -1;
}

///
/// \brief mmr_current_mesh_state_hash - returns an current unique value (state hash, state stamp)
/// that can be used to compare with value returned from this function at other time.
/// If values are different, then one can assume that mesh internal state has changed.
/// \param Mesh_id
/// \return
///
int mmr_current_mesh_state_hash( // returns a state-stamp for mesh
		const int Mesh_id) // Mesh_id for state-stamp
{
	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return hashlittle(mesh->node, mesh->parm.nrno, Mesh_id);

}

/*------------------------------------------------------------
  compatibility additions for remesh module
------------------------------------------------------------*/


void mmr_get_coor_from_motion_element(
		const int Mesh_id,
		int idEl,
		int idLP,
		double * coor,
		int flagaSiatki
) {
	printf("Method mmr_get_coor_from_motion_element not implemented in mmd_prism!\n");
	exit(-1);
}

void mmr_copyMESH(int Mesh_id, int flaga) {
	printf("Method mmr_copyMESH not implemented in mmd_prism!\n");
	exit(-1);
}

void mmr_create_mesh_Cube(
		int Mesh_id,
		const char * nazwa,
		int node_x,
		int node_y,
		int node_z,
		double size_x,
		double size_y,
		double size_z,
		double divide,
		int * warunki
) {
	printf("Method mmr_create_mesh_Cube not implemented in mmd_prism!\n");
	exit(-1);
}

void mmr_init_all_change(
		int Mesh_id,
		int a
) {
	printf("Method mmr_init_all_change not implemented in mmd_prism!\n");
	exit(-1);
}

void mmr_test_mesh_motion(
		int Mesh_id,
		int ileWarstw,
		int obecny_krok,
		int od_krok,
		int ileKrok,
		double minPoprawy,
		double px0,
		double py0,
		double pz0,
		double px1,
		double py1,
		double pz1,
		double endX,
		double endY,
		double endZ
) {
	printf("Method mmr_test_mesh_motion not implemented in mmd_prism!\n");
	exit(-1);
}

double mmr_test_weldpool(
		int Mesh_id,
		double obecnyKrok,
		double krok_start,
		double minPoprawy,
		double doX,
		double doY,
		double dl,
		double zmiejszaPrzes,
		double limit,
		double szerokosc
) {
	printf("Method mmr_test_weldpool not implemented in mmd_prism!\n");
	exit(-1);
}


int mmr_groups_number(const int Mesh_id) {

	mf_fatal_err("Method mmr_groups_number not implemented in mmd_prism!");
}

int mmr_groups_ids(const int Mesh_id, int * tab, int l_tab) {

	mf_fatal_err("Method mmr_groups_ids not implemented in mmd_prism!");
}


/*------------------------------------------------------------
  mmr_get_fa_el_bc_connect(int fa,int &el) - return 0 if normal bc
					  else return connected face (el return el_id)
------------------------------------------------------------*/

int mmr_get_fa_el_bc_connect(const int Mesh_id, int face_id, int * el_id) {
	return 0;
}
