/************************************************************************
File mmh_intf.h - interface with mesh manipulation module of the code.

Contains declarations of constants and interface routines:
mmr_module_introduce - to return the mesh name
mmr_init_mesh - to initialize the mesh data structure (and read data)
mmr_export_mesh - to export (write to a file) mesh data in a specified format
mmr_test_mesh - to test the integrity of mesh data
mmr_get_mesh_i_params - to return mesh parameters
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
mmr_get_next_edge_all - to return the next edgeent id (active or inactive)
mmr_get_nr_node - to return the number of active nodes
mmr_get_max_node_id - to return the maximal node id
mmr_get_max_node_max - to return the maximal possible vertex id
mmr_get_next_node_all - to return the next node id (active or inactive)
mmr_set_max_gen - to set maximal generation level for element refinements
mmr_set_max_gen_diff - to set maximal generation difference between
neighboring elements

mmr_init_ref - to initialize the process of refinement
mmr_refine - to refine an element or the whole mesh
mmr_refine_el -
mmr_refine_mesh -

mmr_derefine - to derefine an element or the whole mesh
mmr_derefine_el -
mmr_derefine_mesh -h

mmr_final_ref - to finalize the process of refinement after rewriting DOFs
mmr_free_mesh - to free space allocated for mesh data structure

mmr_elem_structure - to return elem structure (e.g. for sending)
mmr_el_status - to return element status (active, inactive, free space)
mmr_el_type - to return element type
mmr_el_groupID - to return material ID for element
mmr_el_type_ref - to return element's type of refinement
mmr_el_faces - to get faces and big neighbors of an element
mmr_el_node_coor - to get the coordinates of element's nodes
mmr_el_fam - to return family information for an element
mmr_el_fam_all - to return all recursive family information for an element
mmr_el_gen - to return generation level for an element
mmr_el_ancestor - to find the ancestor of an element
mmr_el_eq_neig - to get equal size (or larger) neighbors of an element
mmr_el_hsize - to compute a characteristic linear size for an element
mmr_set_elem_fath - to set family data for an elem
mmr_set_elem_fam - to set family data for an elem
mmr_face_structure - to return face structure (e.g. for sending)
mmr_fa_status - to return face status (active, inactive, free space)
mmr_fa_type - to return face type (triangle, quad, free space)
mmr_fa_bc - to get the boundary condition flag for a face
mmr_fa_edges - to return a list of face's edges
mmr_fa_neig - to return a list of face's neighbors and
corresponding neighbors' sides (for active
faces all neighbors are active)
mmr_fa_node_coor - to get the list and coordinates of faces's nodes
mmr_fa_elem_coor - to find coordinates within neighboring
elements for a point on face
mmr_fa_area - to compute the area of face and vector normal
mmr_set_face_fam - to set family data for an face
mmr_set_face_neig - to set neighbors data for a face
mmr_edge_status - to return edge status (active, inactive, free space)
mmr_edge_nodes - to return edge node's IDs
mmr_edge_sons - to return edge son's numbers
mmr_edge_structure - to return edge structure (e.g. for sending)
mmr_set_edge_type - to set type for an edge
mmr_set_edge_fam - to set family data for an edge
mmr_node_status - to return node status (active, free space)
mmr_node_coor - to return node coordinates
mmr_el_fa_nodes - to get list local face nodes indexes in elem

mmr_loc_loc - to compute local coordinates within an element,
given local coordinates within an element of the same family
mmr_copyMESH(int flaga);

mmr_divide_face4_t - to break a triangular face into 4 sons
mmr_divide_face4_q - to divide a quadrilateral face into 4 sons
mmr_divide_edge - to divide an edge into two sons

mmr_del_elem    - to free an element structure
mmr_del_face    - to free a face structure
mmr_del_edge    - to free an edge structure
mmr_del_node    - to free a node structure


  mmr_reserve - to inform mesh module about new target number of elems,faces,edges and nodes
  mmr_add_elem - to add (append) new element
  mmr_add_face - to add (append) new face
  mmr_add_edge - to add (append) new edge
  mmr_add_node - to add (append) new node


mmr_test_mesh_motion - mesh motion test
mmr_set_new_and_get_old_uk_val(int idEl,int idGP,double *new_uk_val,double *old_uk_val);
mmr_get_coor_from_motion_element(int idEl,int idLP,double *coor,int flagaSiatki);

mmr_create_mesh_Cube(const char *nazwa,int node_x,int node_y,int node_z,double size_x,double size_y,double size_z,double divide,int *warunki);
mmr_init_all_change(int);

mmr_get_fa_el_bc_connect(int fa,int &el) - return 0 if normal bc else return connected face (el return el_id)
mmr_split_into_blocks_add_contact(int *tabMat,int l_mat,int *tabBlock,int l_block,int *war,int l_war,groups_fluid.data(),groups_fluid.size())

mmr_grups_number()
mmr_grups_ids(int *tab,int l_tab)
mmr_get_vec_norm_at_node(int node_id,int *vec_norm);

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
	07.2009 - Kazimierz Michalik, implementation
*************************************************************************/

//////////////////////////////////////////////////////////////////////////
/// NOTE: inside arrays are numbered 0..n-1, outside arrays are numered 1..n !
/// ASSUMPTION: interface ID == internal POS+1
///				interface ID != internal ID
//////////////////////////////////////////////////////////////////////////



#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cmath>
//#include <GL/gl.h>
#include <assert.h>

#include <modfem/mmh_intf.h>
#include <modfem/uth_log.h>

using namespace std;

#define _USE_MATH_DEFINES

#include <modfem/mm_remesh/ticooMesh3D.h>


TicooMesh3D mes3d(50);

TicooMesh3D * mmv_current_mesh = &mes3d;


#ifdef __cplusplus
extern "C"
{
#endif

const int FIRST = 1;

/* Initialization of constants */

const int MMC_CUR_MESH_ID = 0; /* indicator for the current active mesh */

const int MMC_AUTO_GENERATE_ID = 0;


/* Refinement options */
const int MMC_DO_UNI_REF   = -1; /* indicator to perform uniform refinement */
const int MMC_DO_UNI_DEREF = -2; /* indicator to perform uniform derefinement */
const int MMC_REF_ANI = 2;

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

/* Other */
const int MMC_NO_FATH       = 0;   /* no father indicator */
const int MMC_SAME_ORIENT   = 1;   /* indicator for the same orientation */
const int MMC_OPP_ORIENT   = -1;   /* indicator for the opposite orientation */
const int MMC_MAX_EDGE_ELEMS = 20;
const double SMALL = 1e-10;

//const int MMC_FACE_NODES_FOR_TETRA[4][3]={{3,1,0},{3,2,1},{3,0,2},{0,1,2}};
const int MMC_FACE_NODES_FOR_TETRA[4][3] = {{3, 0, 1}, {3, 1, 2}, {3, 2, 0}, {0, 2, 1}};

const int MMC_FACE_NODES_FOR_PRISM[5][4] = {{0, 1, 2, -1}, {3, 4, 5, -1}, {0, 1, 3, 4}, {1, 2, 4, 5}, {2, 0, 5, 3}};

const int MMC_IS_SUPPORTING_NEW_MPI = 0;
// Definition of global mesh variables.


int	mmv_current_mesh_id = 0;

void atExitClearMeshes()
{
	std::cout << "Method not implemented!: atExitClearMeshes() " << endl;
}



/*---------------------------------------------------------
mmr_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int mmr_chk_list(	/* returns: */
		/* >0 - position on the list */
		/* 0 - not found on the list */
		int Num, 	/* number to be checked */
		int * List, 	/* list of numbers */
		int Ll		/* length of the list */
) {

	int i, il;

	for (i = 0; i < Ll; i++) {
		if ((il = List[i]) == 0) break;
		/* found on the list on (i+1) position */
		if (Num == il) return (i + 1);
	}
	/* not found on the list */
	return (0);
}


/* Definitions of interface routines: */

inline void    select_mesh(const int mesh_id) {

	mmv_current_mesh = &mes3d;
	/*
	  //assert(mesh_id >= 0);
	  if (mesh_id != MMC_CUR_MESH_ID) {
		  //current_mesh = meshes.at(mesh_id);
		  //assert(current_mesh->meshId_ == mesh_id);

		  //un
		  //hObj::myMesh = current_mesh;
		}
		std::cout << "Method not implemented!: select_mesh ";
	*/
}

/*------------------------------------------------------------
  mmr_module_introduce - to return the mesh name
------------------------------------------------------------*/
int mmr_module_introduce(
		/* returns: >=0 - success code, <0 - error code */
		char * Mesh_name /* out: the name of the mesh */
) {
	//strcpy(Mesh_name,const_cast<char*>(current_mesh->name_.c_str()));    // must be, because of compatibility with C
	strcpy(Mesh_name, "3D_remesh");
	return 1;
}



void mmr_init_all_change(const int Mesh_id, int a) {

	mmv_current_mesh->init_all_change(a);

}


/**-----------------------------------------------------------
  mmr_module_internals - to return the mesh internal structures
------------------------------------------------------------*/
void * mmr_module_internals(
		/** returns: >=0 - success code, <0 - error code */
		int Mesh_id/** out: the name of the mesh */
)
{
	return mmv_current_mesh;
}

//------------------------------------------------------------
//  mmr_init_mesh - to initialize the mesh data structure
//------------------------------------------------------------
// NOTE: this is NOT a resize function. It does NOT ADD any elements.
// It only creates new mesh and reserves resources for it.
// To add mesh entities use mmr_add_* functions.
//------------------------------------------------------------
int mmr_init_mesh2(  /* returns: >0 - Mesh ID, <0 - error code */
		FILE * Interactive_output, // in: name of the output file to write out messages
		const int N_nodes, //IN: target count of nodes
		const int N_edges, //IN: target count of edges
		const int N_faces, //IN: target count of faces
		const int N_elems  //IN: target count of elements
)
{
	return 1;
}



int mmr_finish_read(const int Mesh_id)
{
	mmv_current_mesh->setPlansza(1.0, 1.0, 1.0, false, false);
	return 0;
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
//    mmv_current_mesh->points.czysc(100000,100000);
//    mmv_current_mesh->elements.czysc(600000,600000);
	return 0;
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
	mf_fatal_err("This function is not yet implemented.");
	return 0;
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
	mf_fatal_err("This function is not yet implemented. Aborting.");

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
	mf_fatal_err("This function is not yet implemented. Aborting.");

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
	mf_fatal_err("This function is not yet implemented. Aborting.");
}

/*------------------------------------------------------------
  mmr_set_new_and_get_old_uk_val(int idEl,double *new_uk_val,double *old_uk_val);
------------------------------------------------------------*/

void mmr_get_vec_norm_at_node(const int Mesh_id, int node_id, double * vec_norm) {


	mmv_current_mesh->mmr_get_vec_norm_at_node(node_id - 1, vec_norm);



}

void mmr_copyMESH(const int Mesh_id, int flaga) {

	mmv_current_mesh->copyMesh(flaga);


}

/*------------------------------------------------------------
  mmr_set_new_and_get_old_uk_val(int idEl,double *new_uk_val,double *old_uk_val);
------------------------------------------------------------
void mmr_set_new_and_get_old_uk_val(int idEl,int idGP,double *new_uk_val,double *old_uk_val){

	mmv_current_mesh->mmr_set_new_and_get_old_uk_val(idEl-1,idGP,new_uk_val,old_uk_val);

}
*/

/*------------------------------------------------------------
  mmr_test_mesh_motion - test movement mesh
------------------------------------------------------------*/
void mmr_test_mesh_motion(const int Mesh_id, int ileWarstw, int obecny_krok, int od_krok, int ileKrok, double minPoprawy, double px0, double py0, double pz0, double px1, double py1, double pz1, double endX, double endY, double endZ) {


	mmv_current_mesh->ruchPrzestrzen(ileWarstw, obecny_krok, od_krok, ileKrok, minPoprawy, px0, py0, pz0, px1, py1, pz1, endX, endY, endZ);
	//mmv_current_mesh->copyMesh(0);

}

/*------------------------------------------------------------
  mmr_test_weldpool - test movement mesh
------------------------------------------------------------*/
double mmr_test_weldpool(const int Mesh_id, double obecnyKrok, double krok_start, double minPoprawy, double doX, double doY, double dl, double zmiejszaPrzes, double limit, double szerokosc) {


	return mmv_current_mesh->ruchSpawanie(obecnyKrok, krok_start, minPoprawy, doX, doY, dl, zmiejszaPrzes, limit, szerokosc);
	//mmv_current_mesh->copyMesh(0);

}

/*------------------------------------------------------------
  mmr_get_bc_connect - return 0 if normal bc else return connected face_id
------------------------------------------------------------*/
int mmr_get_fa_el_bc_connect(const int Mesh_id, int face_id, int * el_id) {

	return mmv_current_mesh->mmr_get_fa_el_bc_connect(face_id - 1, el_id);


}


/*------------------------------------------------------------
  mmr_get_coor_from_motion_element - get coor from element
------------------------------------------------------------*/

void mmr_get_coor_from_motion_element(const int Mesh_id, int idEl, int idLP, double * coor, int flagaSiatki) {

	mmv_current_mesh->mmr_get_coor_from_motion_element(idEl - 1, idLP, coor, flagaSiatki);
}


/*------------------------------------------------------------
  mmr_create_mesh_Cube - warunki[6]={1,1,1,1,1,1}; //0...5 z0,z1,x0,x1,y0,y1
------------------------------------------------------------*/

void mmr_create_mesh_Cube(const int Mesh_id, const char * nazwa, int node_x, int node_y, int node_z, double size_x, double size_y, double size_z, double divide, int * warunki) {

	mmv_current_mesh->createCube(nazwa, node_x, node_y, node_z, size_x, size_y, size_z, divide, warunki);

}

///*------------------------------------------------------------
//mmr_split_into_blocks_add_contact
//------------------------------------------------------------*/

//int mmr_split_into_blocks_add_contact(const int Mesh_id,const char *workdir, int *tabMat,int l_mat,int *tabBlock,int l_block,int *war,int l_war,double *tempBlock,int l_tempBlock,int *groups_fluid,int l_groups_fluid){


//	mmv_current_mesh->mmr_split_into_blocks_add_contact(workdir,tabMat,l_mat,tabBlock,l_block,war,l_war,tempBlock,l_tempBlock,groups_fluid,l_groups_fluid);

//    return (l_war/4);
//}

int mmr_groups_number(const int Mesh_id) {

	return mmv_current_mesh->mmr_grups_number();

}
int mmr_groups_ids(const int Mesh_id, int * tab, int l_tab) {

	return mmv_current_mesh->mmr_grups_ids(tab, l_tab);

}

/*------------------------------------------------------------
  mmr_init_mesh - to initialize the mesh data structure (and read data)
------------------------------------------------------------*/
int mmr_init_mesh(  /* returns: >0 - Mesh ID, <0 - error code */
		int Control,	    /* in: control variable to choose data format */
		/* MMC_HP_FEM_MESH_DATA = 1 - mesh read from dump files */
		/* MMC_GRADMESH_DATA    = 2 - mesh produced by 2D GRADMESH mesh generator */
		char * Filename,   /* in: name of the file to read mesh data */
		// KM 06.2011
		// We assume that there is also a detailed boundary description file
		// which name is similar to the input file, but it ends with ".detail".
		FILE * interactive_output
) {

	switch (Control) {

		case MMC_MOD_FEM_PRISM_DATA:  /* prism mesh read from dump files */
		case MMC_MOD_FEM_MESH_DATA:  /* generic mesh for this module */
			//   importer = new MeshRead::DmpFileImporter(Filename);
			break;

		case MMC_MOD_FEM_TETRA_DATA: /* tetra mesh read from dump files */
		case MMC_MOD_FEM_HYBRID_DATA:  /* hybrid mesh read from dump files */
		{
			std::cout << "MMC_MOD_FEM_TETRA_DATA" << endl;
			mmv_current_mesh->wczytajNAS(Filename, 1, 0, 0, 0);
			mmv_current_mesh->setPlansza(1.0, 1.0, 1.0, false, false);
			//mmv_current_mesh->setPlansza(2,10,0.1,false,false);
			//mmv_current_mesh->setPlansza(0.1,0.005,0.1,false,false);
			//mmv_current_mesh->setPlansza(0.003,0.003,0.001,false,false);//spawanie

			//wczytywanie siatki
		}
		break;

		case MMC_NASTRAN_DATA:
			//std::cout << "MMC_NASTRAN_DATA"<<endl;

			mmv_current_mesh->wczytajNAS(Filename, 1, 0, 0, 0);
			mmv_current_mesh->setPlansza(1.0, 1.0, 1.0, false, false);
			//mmv_current_mesh->setPlansza(2,10,0.1,false,false);
			//mmv_current_mesh->setPlansza(0.1,0.005,0.1,false,false);
			//mmv_current_mesh->setPlansza(0.003,0.003,0.001,false,false);//spawanie


			//importer = new MeshRead::NasFileImporter(Filename);
			break;
		case MMC_GRADMESH_DATA:  /* mesh produced by 2D GRADMESH generator */
			std::cout << "Method not implemented!";
			exit(-1);
			break;
		case MMC_BINARY_DATA:
			//importer = new MeshRead::BinaryFileReader(Filename);
			break;
		case MMC_IN_ANSYS_DATA:

			mmv_current_mesh->wczytajNAS(Filename, 1, 0, 0, 0);





			/*
			mmv_current_mesh->wczytajAbaqus(Filename,1,0,0,0);
			string a = Filename;
			a+="_nas";
			mmv_current_mesh->ZapiszNAS(a.c_str(),1,0,0,0);
			std::cout << "you have new mesh file *_nas!";
			exit(-1);
			*/

//	  	importer = new MeshRead::InFileImporter(Filename);
			break;
	}//!switch(Control)

	//zwroc numer nowej siatki
	//std::cout << "Method not implemented!: mmr_init_mesh "<<endl;
	return 1;
}

/*---------------------------------------------------------
  mmr_export_mesh - to export (write to a file) mesh data in a specified format
---------------------------------------------------------*/
int mmr_export_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Control,	  /* in: control variable to choose data format */
		/* MMC_HP_FEM_MESH_DATA = 1 - mesh written to standard dump files */
		char * Filename /* in: name of the file to write mesh data */
) {
	if (Filename != NULL) {
		//select_mesh(Mesh_id);

		//MeshWrite::IMeshWriter * writer(NULL);
		switch (Control) {
			case MMC_BINARY_DATA:
				//writer = new MeshWrite::BinaryFileWriter(Filename);
				break;
			//case MMC_NASTRAN_DATA:
			//writer = new MeshWrite::Nas
			// break;
			case MMC_BOUND_VERTS_DATA:
				//writer = new MeshWrite::MileninExporter(Filename);
				break;
			case MMC_HP_FEM_MESH_DATA:
			default: {
				//writer = new MeshWrite::KazFileExporter(Filename);
				mmv_current_mesh->zapisParaView(Filename);
			}
			break;
		}

		//current_mesh->write(*writer);
		mf_log_info("mesh exported to %s", Filename);

		//delete writer;
		return 0;
	}

	std::cerr << "\n Error: export_mesh failed, no filename!" << endl;


	return -1;
}

/*---------------------------------------------------------
  mmr_test_mesh - to test the integrity of mesh data
---------------------------------------------------------*/
int mmr_test_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//select_mesh(Mesh_id);

	return mmv_current_mesh->getWielkoscTabX() ? 0 : -1;
}


/*---------------------------------------------------------
  mmr_get_nr_elem - to return the number of active elements
---------------------------------------------------------*/
int mmr_get_nr_elem(/* returns: >=0 - number of active elements, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//select_mesh(Mesh_id);

	return mmv_current_mesh->getNumberElements();
}

/*---------------------------------------------------------
  mmr_get_max_elem_id - to return the maximal element id
---------------------------------------------------------*/
int mmr_get_max_elem_id(  /* returns: >=0 - maximal element id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//select_mesh(Mesh_id);
	return mmv_current_mesh->getNumberElements();
}

/*---------------------------------------------------------
  mmr_get_max_elem_max - to return the maximal possible element id
---------------------------------------------------------*/
int mmr_get_max_elem_max(  /* returns: >=0 - maximal element id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//un
	//return hObj::maxPos();
	return mmv_current_mesh->getNumberElements();
}

/*---------------------------------------------------------
  mmr_get_next_act_elem - to return the next active element id
---------------------------------------------------------*/
int mmr_get_next_act_elem(/* returns: >=0 - the next active element ID */
		/*                 (0 - input is the last element) */
		/*                 <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nel       /* in: input element (0 - return first active element) */
) {

	// Because we have only 'active' elements, there is only need of checking is 'Nel' last element;
	// if not return next element id;
	// if Nel is last, next hObj is filled with zeros, so id will be zero

	if (mmv_current_mesh->getNumberElements() <= Nel) {
		return 0;
	}

	return ++Nel;

}

/*---------------------------------------------------------
  mmr_get_next_elem_all - to return the next element id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_elem_all( /* returns: >=0 - ID of the next element */
		/*                (0 - input is the last element) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nel       /* in: input element (0 - return first element) */
) {
	return mmr_get_next_act_elem(Mesh_id, Nel);
}

/*---------------------------------------------------------
  mmr_get_nr_face - to return the number of active faces
---------------------------------------------------------*/
int mmr_get_nr_face(/* returns: >=0 - number of active faces, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//un

	return mmv_current_mesh->getNumberFaces();
}

/*---------------------------------------------------------
  mmr_get_max_face_id - to return the maximal face id
---------------------------------------------------------*/
int mmr_get_max_face_id(  /* returns: >=0 - maximal face id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	return mmv_current_mesh->getNumberFaces();
}

/*---------------------------------------------------------
  mmr_get_next_act_face - to return the next active face id
---------------------------------------------------------*/
int mmr_get_next_act_face( /* returns: >=0 - ID of the next active face */
		/*                (0 - input is the last face) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nfa       /* in: input face (0 - return first active face) */
) {


	if (mmv_current_mesh->getNumberFaces() <= Nfa) {
		return 0;
	}
	return ++Nfa;
}


/*---------------------------------------------------------
  mmr_get_next_face_all - to return the next face id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_face_all( /* returns: >=0 - ID of the next face */
		/*                (0 - input is the last face) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nfa       /* in: input face (0 - return first face) */
) {
	return mmr_get_next_act_face(Mesh_id, Nfa);
}

/*---------------------------------------------------------
  mmr_get_nr_edge - to return the number of active edges
---------------------------------------------------------*/
int mmr_get_nr_edge(/* returns: >=0 - number of active edges, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//un
	return mmv_current_mesh->getNumberEdges();
}

/*---------------------------------------------------------
  mmr_get_max_edge_id - to return the maximal edge id
---------------------------------------------------------*/
int mmr_get_max_edge_id(  /* returns: >=0 - maximal edge id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	return mmv_current_mesh->getNumberEdges();
}

/*---------------------------------------------------------
  mmr_get_next_edge_all - to return the next edgeent id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_edge_all( /* returns: >=0 - ID of the next edge */
		/*                (0 - input is the last edge) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ned       /* in: input edge (0 - return first edge) */
) {

	if (mmv_current_mesh->getNumberEdges() <= Ned) {
		return 0;
	}
	return ++Ned;
}

/*---------------------------------------------------------
  mmr_get_nr_node - to return the number of active nodes
---------------------------------------------------------*/
int mmr_get_nr_node(/* returns: >=0 - number of active nodes, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {

	return mmv_current_mesh->getNumberPoints();
}

/*---------------------------------------------------------
  mmr_get_max_node_id - to return the maximal node id
---------------------------------------------------------*/
int mmr_get_max_node_id(  /* returns: >=0 - maximal node id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//select_mesh(Mesh_id);
	return mmv_current_mesh->getNumberPoints();
}

/*---------------------------------------------------------
mmr_get_max_node_max - to return the maximal possible vertex id
---------------------------------------------------------*/
extern int mmr_get_max_node_max(  /* returns: >=0 - maximal possible vertex id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {

	return mmv_current_mesh->getNumberPoints();
}

/*---------------------------------------------------------
  mmr_get_next_node_all - to return the next node id (active or inactive)
---------------------------------------------------------*/
int mmr_get_next_node_all( /* returns: >=0 - ID of the next node */
		/*                (0 - input is the last node) */
		/*           <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Nno       /* in: input node (0 - return first node) */
) {

	if (mmv_current_mesh->getNumberPoints() <= Nno) {
		return 0;
	}
	return ++Nno;
}

/*---------------------------------------------------------
  mmr_set_max_gen - to set maximal generation level for element refinements
------------------------------------------------------------*/
int mmr_set_max_gen(/* returns: >=0-success code, <0-error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Max_gen   /* in: maximal generation */
) {
	//select_mesh(Mesh_id);

	//current_mesh->maxGen_ = Max_gen;
	//std::cout << "Method not implemented! mmr_set_max_gen"<<endl;
	//exit(-1);
	return 0;
}

/*---------------------------------------------------------
  mmr_set_max_gen_diff - to set maximal generation difference between
						 neighboring elements
------------------------------------------------------------*/
int mmr_set_max_gen_diff(/* returns: >=0-success code, <0-error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Max_gen_diff    /* in: maximal generation difference */
) {
	//select_mesh(Mesh_id);
	// Mesh is 1-unregular

	//std::cout << "Method not implemented! mmr_set_max_gen_diff"<<endl;
	//exit(-1);
	return 0;
}

/*------------------------------------------------------------
  mmr_init_ref - to initialize the process of refinement
------------------------------------------------------------*/
int mmr_init_ref(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	//select_mesh(Mesh_id);
	std::cout << "Method not implemented! mmr_init_ref" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
  mmr_refine - to refine an element or the whole mesh
------------------------------------------------------------*/
int mmr_refine( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
) {
	/*
	select_mesh(Mesh_id);

	if(El==MMC_DO_UNI_REF) {
	  current_mesh->refine();
	}
	else {
	  current_mesh->refineElem(current_mesh->elements_[El]);
	}
	*/
	//std::cout << "\nDividing element " << El;
	std::cout << "Method not implemented! mmr_refine" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
  mmr_refine_el - to refine an element or the whole mesh
------------------------------------------------------------*/
int mmr_refine_el( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
) {
	/*
	select_mesh(Mesh_id);

	if(El==MMC_DO_UNI_REF) {
	  current_mesh->refine();
	}
	else {
	  current_mesh->refineElem(current_mesh->elements_[El]);
	}
	*/
	//std::cout << "\nDividing element " << El;
	std::cout << "Method not implemented! mmr_refine_el" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
  mmr_refine_mesh - to refine an element or the whole mesh
------------------------------------------------------------*/
int mmr_refine_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		/* in: element ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
) {
	/*
	select_mesh(Mesh_id);

	if(El==MMC_DO_UNI_REF) {
	  current_mesh->refine();
	}
	else {
	  current_mesh->refineElem(current_mesh->elements_[El]);
	}
	*/
	//std::cout << "\nDividing element " << El;
	std::cout << "Method not implemented! mmr_refine_mesh" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
mmr_gen_boundary_layer - generation boundary layer
------------------------------------------------------------*/
extern int mmr_gen_boundary_layer( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Bc, /* in: boundary conditon ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
		int thicknessProc,
		int noLayers,
		int distrib,
		double * ignoreVec
) {
	//select_mesh(Mesh_id);
	std::cout << "Method not implemented! mmr_gen_boundary_layer" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
mmr_r_refine - to r-refine an elements with given boundary condition
------------------------------------------------------------*/
extern int mmr_r_refine( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Bc,  /* in: boundary conditon ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
		void (*reallocation_func)(double * x, double * y, double * z)	/* function moving to new*/
) {
	//select_mesh(Mesh_id);
	std::cout << "Method not implemented! mmr_r_refine" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
  mmr_derefine - to derefine an element or the whole mesh
------------------------------------------------------------*/

int mmr_derefine( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID or -2 (MMC_DO_UNI_DEREF) for uniform derefinement */
) {
	std::cout << "Method not implemented! mmr_derefine" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
  mmr_derefine_el - to derefine an element or the whole mesh
------------------------------------------------------------*/

int mmr_derefine_el( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID or -2 (MMC_DO_UNI_DEREF) for uniform derefinement */
) {
	std::cout << "Method not implemented! mmr_derefine_el" << endl;
	exit(-1);
	return -1;
}

/*------------------------------------------------------------
  mmr_derefine_mesh - to derefine an element or the whole mesh
------------------------------------------------------------*/

int mmr_derefine_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		/* in: element ID or -2 (MMC_DO_UNI_DEREF) for uniform derefinement */
) {
	std::cout << "Method not implemented! mmr_derefine_mesh" << endl;
	exit(-1);
	return -1;
}


/*------------------------------------------------------------
  mmr_final_ref - to finalize the process of refinement after rewriting DOFs
------------------------------------------------------------*/
int mmr_final_ref(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	std::cout << "Method not implemented! mmr_final_ref" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
  mmr_free_mesh - to free space allocated for mesh data structure
---------------------------------------------------------*/
int mmr_free_mesh(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	std::cout << "Method not implemented! mmr_free_mesh" << endl;
	return 1;
}


/*---------------------------------------------------------
mmr_elem_structure - to return elem structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_elem_structure( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	   /* in: elem ID */
		int * Elem_struct /* out: elem structure in the form of integer array */
) {

	//const hObj & elem(current_mesh->elements_[El]);
	Elem_struct[0] = 7;
	//Elem_struct[1]=mesh->elem[El].ipid;
	Elem_struct[2] = 1;
	Elem_struct[3] = 0;
	Elem_struct[4] = 0;

	//un
	/*
	for (unsigned int i(0); i < elem.typeSpecyfic_.nComponents_; ++i) {
		Elem_struct[5+i]=hObj::posFromId(elem.components(i));
	}
	*/
	std::cout << "Method not implemented! mmr_elem_structure" << endl;
	exit(-1);
	return (-1);

}

/*---------------------------------------------------------
  mmr_el_status - to return element status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_el_status( /* returns element status: */
		/* +1 (MMC_ACTIVE)   - active element */
		/*  0 (MMC_FREE)     - free space */
		/* -1 (MMC_INACTIVE) - inactive (refined) element */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
) {
	if (mmv_current_mesh->getNumberElements() <= El - 1) {
		return 0;
	}
	return 1;
}


/*---------------------------------------------------------
  mmr_el_type - to return element type
---------------------------------------------------------*/

int mmr_el_type( /* returns element type or <0 - error code */
		/*	 7 (MMC_TETRA) - tetrahedron */
		/*	 5 (MMC_PRISM) - prism */
		/*	 6 (MMC_BRICK) - brick */
		/*	 0 (MMC_FREE)  - free space */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
) {


	return 7;
}

/*---------------------------------------------------------
  mmr_el_groupID - to return material ID for element
---------------------------------------------------------*/
int mmr_el_groupID( /* returns material flag for element, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
) {

	return mmv_current_mesh->mmr_el_groupID(El - 1);
	//return 1;
}
/*---------------------------------------------------------
mmr_el_set_mate - to set material number for element
---------------------------------------------------------*/
int mmr_el_set_groupID( /* sets material flag for element */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	/* in: element ID */
		int Mat_id    /* in: material ID */
) {
	std::cout << "Method not implemented! mmr_el_set_mate" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
  mmr_el_type_ref - to return element's type of refinement
---------------------------------------------------------*/
int mmr_el_type_ref(
		/* returns element's type of refinement, <0 - error code */
		/*         MMC_NOT_REF - not refined */
		/*         MMC_REF_ISO - isotropic refinement */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
) {

	return MMC_NOT_REF;
}

/*---------------------------------------------------------
mmr_el_faces - to get faces of an element
---------------------------------------------------------*/
int mmr_el_faces( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,     	/* in: element ID */
		int * Faces,  	/* out: list of faces */
		/*	(Faces[0] - number of faces) */
		int * Orient	/* out: orientation for each face */
		/*	+1 (MMC_SAME_ORIENT) - the same as element */
		/*	-1 (MMC_OPP_ORIENT) - opposite */
) {


	const int n_faces = mmv_current_mesh->mmr_el_faces(El - 1, Faces, Orient);
	if (Faces != NULL) {
		for (int i = 1, ileI = Faces[0] + 1; i < ileI; ++i) {
			Faces[i] += 1;
		}
	}
	return n_faces;
}

/*---------------------------------------------------------
  mmr_el_node_coor - to get the coordinates of element's nodes
---------------------------------------------------------*/
int mmr_el_node_coor( 	/* returns number of nodes */
		int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	  /* in: element ID */
		int * Nodes,  	 /* out: list of vertex node IDs */
		/*	(Nodes[0] - number of nodes) */
		double * Node_coor /* out: coordinates of element vertices */
) {


	if (El > 0) {
		mmv_current_mesh->mmr_el_node_coor(El - 1, Nodes, Node_coor);
		if (Nodes != NULL) {
			for (int i = 1, ileI = Nodes[0] + 1; i < ileI; ++i) {
				Nodes[i] += 1;
			}
		}
	}


	return 4;
}

/*---------------------------------------------------------
  mmr_el_fam - to return family information for an element
---------------------------------------------------------*/
int mmr_el_fam( /* returns: element father ID or 0 (MMC_NO_FATH) for */
		/*          initial mesh element or <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,       /* in: element ID (1...) */
		int * Elsons, /* out: list of element's sons */
		/* 	Elsons[0] - number of sons */
		int * Type	/* out: type of refinement */
		/*         MMC_NOT_REF - not refined */
		/*         MMC_REF_ISO - isotropic refinement */
) {

	Elsons[0] = 0;
	if (Type != NULL) {
		*Type = MMC_NOT_REF;
	}
	return 0;
}

//---------------------------------------------------------
// mmr_el_fam_all - to return all recursive family information for an element
//---------------------------------------------------------
int mmr_el_fam_all( /* returns: element father ID or 0 (MMC_NO_FATH) for */
		/*          initial mesh element or <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,       /* in: element ID (1...) */
		int * Elsons /* out: list of all recursive element's sons */
		/* 	Elsons[0] - number of sons */
)
{
	Elsons[0] = 0;

	return 0;
}


/*---------------------------------------------------------
  mmr_el_gen - to return generation level for an element
---------------------------------------------------------*/
int mmr_el_gen(   /* returns: El's generation ID */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El        /* in: element ID (1...) */
) {

	return 0;
}

/*---------------------------------------------------------
  mmr_el_ancestor - to find the ancestor of an element
					with generation level Ilev
---------------------------------------------------------*/
int mmr_el_ancestor( /* returns: >0 - ancestor ID, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,       /* in: element ID (1...) */
		int Ilev      /* in: level ID */
) {
	std::cout << "Method not implemented! mmr_el_ancestor " << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
  mmr_el_hsize - to compute a characteristic linear size for an element
			 (for linear and multi-linear 3D elements)
---------------------------------------------------------*/
double mmr_el_hsize(   /* returns: element size */
		int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,         /* in: element ID (1...) */
		double * Size_x, /* out: for anizotropic elements */
		double * Size_y,
		double * Size_z
) {

	if (El > 0) {
		return mmv_current_mesh->mmr_el_hsize(El - 1);
	}

	std::cout << "Method not implemented!: mmr_el_hsize El>0" << endl;
	exit(-1);
	return -1;

}

/*---------------------------------------------------------
mmr_el_eq_neig - to get equal size neighbors of an element
---------------------------------------------------------*/
int mmr_el_eq_neig( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,     	/* in: element ID */
		int * Neig,   	/* out: list of equal size neighbors */
		/*  >0 - equal size neighbor */
		/*  -1 (MMC_BIG_NGB) - big neighbor */
		/*   0 (MMC_BOUNDARY) - boundary (always second neighbor)*/
		int * Neig_sides /* out: list of sides of neighbors */
) {
	std::cout << "Method not implemented!: mmr_el_hsize mmr_el_eq_neig" << endl;
	return mmv_current_mesh->mmr_el_eq_neig(El - 1, Neig, Neig_sides);;
}

/*---------------------------------------------------------
mmr_face_structure - to return face structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_face_structure(/* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,  	   /* in: face ID */
		int * Face_struct /* out: face structure in the form of integer array */
) {

	/*
	select_mesh(Mesh_id);

	int noComponents = current_mesh->faces_[Fa].typeSpecyfic_.nComponents_;
	const hObj & face(current_mesh->faces_[Fa]);
	if (Face_struct != NULL) {
		Face_struct[0] = face.type_;
		//Face_struct[1]=mesh->face[Fa].ipid;
		Face_struct[2]=face.flags(B_COND);

		for (int i=0; i<noComponents; ++i) {
			Face_struct[3+i]=hObj::posFromId(current_mesh->faces_[Fa].components(i));
		}

		if(face.flags(F_TYPE)==F_FLIPPED) {
			std::swap(Face_struct[3+1],Face_struct[3+2]);
		}

		Face_struct[7]=face.neighs(0);
		Face_struct[8]=face.neighs(1);
		//  if(face.isBroken())
		//  {
		//	for(int i=0;i<4;i++)
		//	{
		//	  Face_struct[9+i]=mesh->face[Fa].sons[i];
		//	}
		//}
	}
	*/
	std::cout << "Method not implemented!: mmr_face_structure " << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
  mmr_fa_status - to return face status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_fa_status( /* returns face status: */
		/*	 1 (MMC_ACTIVE)   - active face */
		/*	 0 (MMC_FREE)     - free space */
		/*	-1 (MMC_INACTIVE) - inactive face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa  	/* in: face ID */
) {
	if (mmv_current_mesh->getNumberFaces() <= Fa - 1) {
		return 0;
	}
	return 1;
}

/*---------------------------------------------------------
  mmr_fa_type - to return face type (triangle, quad, free space)
---------------------------------------------------------*/
int mmr_fa_type( /* returns face type: */
		/*	 3 (MMC_TRIA) - triangle */
		/*	 4 (MMC_QUAD) - quadrilateral */
		/*	 0 (MMC_FREE) - free space */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa  	/* in: face ID */
) {

	return 3;
}

/*--------------------------------------------------------------------------
  mmr_fa_bc - to get the boundary condition flag for a face
---------------------------------------------------------------------------*/
int mmr_fa_bc( /* returns: bc flag for a face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa	      /* in: global face ID */
) {

	return mmv_current_mesh->getFace_bc(Fa - 1);
}

/*---------------------------------------------------------
mmr_fa_sub_bnd - to indicate face is on the boundary
---------------------------------------------------------*/
extern int mmr_fa_sub_bnd( /* returns: 1 - true, 0 - false */
		int Mesh_id,    /* in: mesh ID */
		int Face_id     /* in: face ID */
) {

	if (mmv_current_mesh->getFace_bc(Face_id - 1)) {
		return true;
	}
	return false;
}

/*---------------------------------------------------------
mmr_fa_set_sub_bnd - to indicate face is on the boundary
---------------------------------------------------------*/
extern int mmr_fa_set_sub_bnd(/* returns: >=0 - success code, <0 - error code */
		int Mesh_id,    /* in: mesh ID */
		int Face_id,    /* in: face ID */
		int Side_id     /* in: side ID ??? */
) {

	if (mmv_current_mesh->getNumberFaces() < Face_id - 1) {
		return 0;
	}
	mmv_current_mesh->setFace_bc(Face_id - 1, Side_id);
	return 1;
}
/*---------------------------------------------------------
mmr_fa_edges - to return a list of face's edges
---------------------------------------------------------*/
int mmr_fa_edges( /* returns: >=0 - number of edges, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,         /* in:  global face ID */
		int * Fa_edges, /* out: face edges */
		int * Ed_orient /* out: edges orientation: */
		/*	+1 (MMC_SAME_ORIENT) - the same as face */
		/*	-1 (MMC_OPP_ORIENT)  - opposite */
) {

	mmv_current_mesh->mmr_fa_edges(Fa - 1, Fa_edges, Ed_orient);
	for (int i = 1, ileI = Fa_edges[0] + 1; i < ileI; ++i) {
		Fa_edges[i] += 1;
	}
	return Fa_edges[0];
}

/*---------------------------------------------------------
mmr_fa_eq_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (only equal size
		neighbors considered)
---------------------------------------------------------*/
void mmr_fa_eq_neig(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,       /* in:  global face ID */
		int * Fa_neig, /* out: face neighbors */
		/* (first neighbor has the same orientation */
		/*  as the face and ordering of its nodes */
		/*  defines ordering of nodes on the face) */
		/*  >0 - equal size neighbor ID */
		/*  -1 (MMC_BIG_NGB) - big neighbor */
		/*   0 (MMC_BOUNDARY) - boundary (always second neighbor)*/
		int * Neig_sides,     /* out: side local IDs for equal size neighbors */
		int * Node_shift	/* out: the difference in positions between */
		/*	first face's node for first and second */
		/* 	neighbor (usage in mmr_fa_node_coor) */
) {

	if (Fa > 0) {
		mmv_current_mesh->getElwithFace(Fa - 1, Fa_neig, Neig_sides);
	}


	//++Fa_neig in fun() == getElwithFace
}

/*---------------------------------------------------------
  mmr_fa_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (for active
		faces all neighbors are active)
---------------------------------------------------------*/
void mmr_fa_neig(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,               /* in:  global face ID */
		int * Fa_neig,        /* out: face neighbors */
		/* (first neighbor has the same orientation */
		/*  as the face and ordering of its nodes */
		/*  defines ordering of nodes on the face) */
		/*  >0 - equal size neighbor ID */
		/*  <0 - big neighbor ID */
		/*   0 (MMC_BOUNDARY) - boundary (always 2nd neighbor)*/
		int * Neig_sides,     /* out: side IDs for neighbors */
		int * Node_shift,     /* out: the difference in positions between */
		/*	first face's node for first and second */
		/* 	neighbor (usage in mmr_fa_node_coor) */
		int * Diff_gen,	/* out: generation difference between neighbors */
		double * Acoeff,      /* out: coefficients of linear transformation... */
		double * Bcoeff	/* 	between face coordinates and */
		/* 	coordinates on an ancestor face  */
		/* 	being a side of big neighbor  */
) {


	//if(Acoeff!=NULL){std::cout << "Method not implemented!: Acoeff "<<endl;}
	//if(Bcoeff!=NULL){std::cout << "Method not implemented!: Bcoeff "<<endl;}
	if (Diff_gen != NULL) {
		std::cout << "Method not implemented!: Diff_gen " << endl;
	}
	/*
		assert(Fa > 0);
		select_mesh(Mesh_id);

		const hObj&	face(current_mesh->faces_.at(Fa));
		Fa_neig[0]=hObj::posFromId(face.neighs(0)) ;
		Fa_neig[1]=hObj::posFromId(face.neighs(1)) ;
		Neig_sides[0] = (Fa_neig[0] > B_COND) ? current_mesh->whichFaceAmI(face,current_mesh->elements_.getById(face.neighs(0) )) : 0;
		Neig_sides[1] = (Fa_neig[1] > B_COND) ? current_mesh->whichFaceAmI(face,current_mesh->elements_.getById(face.neighs(1) )) : 0;

		// HACK to conform with KB prism faces ordering
		if(Fa_neig[0] > B_COND) {
		if(current_mesh->elements_(face.neighs(0)).type_==ElemPrism::myType) {
			switch(Neig_sides[0]) {
			case 3: Neig_sides[0]=4; break;
			case 4: Neig_sides[0]=3; break;
			}
		}
		}
		if(Fa_neig[1] > B_COND) {
		if(current_mesh->elements_(face.neighs(1)).type_==ElemPrism::myType) {
			switch(Neig_sides[1]) {
			case 3: Neig_sides[1]=4; break;
			case 4: Neig_sides[1]=3; break;
			}
		}
		}
		// end of HACK

		if (Node_shift != NULL) *Node_shift = 0;
		if (Diff_gen != NULL) *Diff_gen = current_mesh->elements_[Fa_neig[0]].level_ - current_mesh->elements_[Fa_neig[1]].level_;
		assert(Fa_neig[0] != Fa_neig[1]);
		*/
	//std::cout << "Method not implemented!: mmr_fa_neig "<<endl;
	//exit(-1);


	if (Fa > 0) {
		mmv_current_mesh->getElwithFace(Fa - 1, Fa_neig, Neig_sides);
	}
	//++Fa_neig in fun() == getElwithFace

}

/*---------------------------------------------------------
mmr_fa_node_coor - to get the list and coordinates of faces's nodes
---------------------------------------------------------*/
extern int mmr_fa_node_coor( /* returns: number of nodes for a face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,  	     /* in: face ID */
		int * Nodes,  	    /* out: list of vertex node IDs */
		double * Node_coor /* out: coordinates of face vertices */
) {

	mmv_current_mesh->mmr_fa_node_coor(Fa - 1, Nodes, Node_coor);
	for (int i = 1, ileI = Nodes[0] + 1; i < ileI; ++i) {
		Nodes[i] += 1;
	}
	return Nodes[0];
}

/*--------------------------------------------------------------------------
mmr_fa_fam - to return face's family information
---------------------------------------------------------------------------*/
int mmr_fa_fam( /* returns: face's father */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,      	/* in: global face ID */
		int * Fasons,	/* out: sons */
		/* 	Fasons[0] - number of sons */
		int * Node_mid	/* out: node in the middle (if any) */
) {

	/*
		select_mesh(Mesh_id);

		assert(Fa >= 0);
		assert(Fa <  static_cast<int>(current_mesh->faces_.last()));
		assert(current_mesh->faces_[Fa].isBroken());

		Fasons[0] = current_mesh->faces_[Fa].nMyClassSons_;

		if (current_mesh->faces_[Fa].type_ == Face4::myType) {
		for (int i(1); i <= Fasons[0]; ++i) {
				Fasons[i]=current_mesh->faces_[Fa].getMyClassChild<Face4>(i-1)->pos_;
		}
			*Node_mid=current_mesh->faces_[Fa].sons(0) ;
		}
		else {
		for (int i(1); i <= Fasons[0]; ++i) {
				Fasons[i]=current_mesh->faces_[Fa].getMyClassChild<Face3>(i-1)->pos_;
		}
			*Node_mid=0;
		}
		*/

	std::cout << "Method not implemented!: mmr_fa_fam " << endl;
	exit(-1);
	return 0;
}


/*---------------------------------------------------------
  mmr_fa_elem_coor - to find coordinates within neighboring
				 elements for a point on face
---------------------------------------------------------*/
void mmr_fa_elem_coor(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		double * Xloc,	       /* in: local coordinates on a face */
		int * Fa_neig,	       /* in: face neighbors (<0 - big) */
		/*	first - same orientation */
		/*	second - opposite orientation */
		int * Neig_sides,     /* in: which side face is for neighbors */
		int, //Node_shift,       /* in: the difference in positions between */
		/*	first face's node for first and second */
		/* 	neighbor (usage in mmr_fa_node_coor) */
		double * Acoeff,  /* in: coefficients of transformation between... */
		double * Bcoeff,	 /* in: ...face coord and big neighb face coord */
		double * Xneig    /* out: local coordinates for neighbors */
) {
	/* auxiliary variables */
	int neig, ineig, el_type;
	double coor[2];

	/*++++++++++++++++ executable statements ++++++++++++++++*/



	//select_mesh(Mesh_id);

	/* loop over neighbors */
	for (ineig = 0; ineig < 2; ineig++) {
		if (Fa_neig[ineig] != 0) {
			neig = Fa_neig[ineig];
			neig *= neig < 0 ? -1 : 1;
			el_type = mmr_el_type(Mesh_id, neig);
			/* equal size neighbor - coordinate is just pased */
			if (Fa_neig[ineig] > 0) {
				coor[0] = Xloc[0];
				coor[1] = Xloc[1];
			} else {
				/* transform coeficients */
				coor[0] = Acoeff[0] * Xloc[0] + Acoeff[1] * Xloc[1] + Bcoeff[0];
				coor[1] = Acoeff[2] * Xloc[0] + Acoeff[3] * Xloc[1] + Bcoeff[1];
				// coor[0] = Xloc[1];
				// coor[1] = Xloc[0];

			}

			// if face is pointing inside elem, swap x and y
			//const int type=current_mesh->elements_[neig].flags(EL_TYPE);
			//int mask= 1<< (Neig_sides[ineig]);
			//if((type&mask) != 0) {
			//std::swap(coor[0],coor[1]);
			//}
			switch (el_type) {
				/*
							case MMC_PRISM: {
							assert(Xloc[0] >= -1.0-SMALL && Xloc[0] <= 1.0+SMALL);
						assert(Xloc[1] >= -1.0-SMALL && Xloc[1] <= 1.0+SMALL);
						assert(Xloc[0]+Xloc[1] <= 2.0+SMALL);
						assert(Xloc[0]+Xloc[1] >= -2.0-SMALL);

							switch(Neig_sides[ineig]) {
							case 0: // 0,1 - triangles (bottom, top)
							case 1: {
								Xneig[3*ineig  ]=coor[0];
								Xneig[3*ineig+1]=coor[1];
									Xneig[3*ineig+2]=(Neig_sides[ineig]==0? -1 : 1);
								} break;
								case 2: { // 2,4,3 - rectangles (sides)
									Xneig[3*ineig  ] = 0.5*(1+coor[0]);
									Xneig[3*ineig+1] = 0.0;
									Xneig[3*ineig+2] = coor[1];

							   } break;
							   case 3: {
									Xneig[3*ineig  ] = 1-0.5*(1+coor[0]);
									Xneig[3*ineig+1] = 0.5*(1+coor[0]);
									Xneig[3*ineig+2] = coor[1];
								} break;
								case 4: {
									Xneig[3*ineig  ] = 0.0;
									Xneig[3*ineig+1] = 0.5*(1+coor[0]);
									Xneig[3*ineig+2] = coor[1];
								} break;

								} //!switch(Neig_sides[ineig]
							}break;
				*/
				case MMC_TETRA: {
					assert(Xloc[0] >= -SMALL && Xloc[0] <= 1.0 + SMALL);
					assert(Xloc[1] >= -SMALL && Xloc[1] <= 1.0 + SMALL);
					assert(Xloc[0] + Xloc[1] <= 1.0);

					// all trans below are for faces pointing outside


					switch (Neig_sides[ineig]) {
						case 0: {

							Xneig[3 * ineig + 0] = coor[1];
							Xneig[3 * ineig + 1] = 0.0;
							Xneig[3 * ineig + 2] = coor[0];


							/*
								Xneig[3*ineig+0]=coor[1];
								Xneig[3*ineig+1]=coor[0];
								Xneig[3*ineig+2]=0.0;
							*/


						}
						break;
						case 1: {

							//coor[0]*=0.75;
							//coor[1]*=0.75;

							//coor[0]*=0.6445993;
							//coor[1]*=0.8165065;

							Xneig[3 * ineig + 0] = 1.0 - coor[0] - coor[1];
							Xneig[3 * ineig + 1] = coor[0];
							Xneig[3 * ineig + 2] = coor[1];

							/*
								Xneig[3*ineig+0]=coor[0];
								Xneig[3*ineig+1]=0.0;
								Xneig[3*ineig+2]=coor[1];
							*/

						}
						break;
						case 2: {

							Xneig[3 * ineig + 0] = 0.0;
							Xneig[3 * ineig + 1] = coor[1];
							Xneig[3 * ineig + 2] = coor[0];

							/*
								Xneig[3*ineig+0]=0.0;
								Xneig[3*ineig+1]=coor[1];
								Xneig[3*ineig+2]=coor[0];
							*/

						}
						break;
						case 3: {

							Xneig[3 * ineig + 0] = coor[1];
							Xneig[3 * ineig + 1] = coor[0];
							Xneig[3 * ineig + 2] = 0.0;

							/*
								coor[0]*=0.75;
								coor[1]*=0.75;
								Xneig[3*ineig+0]=1.0-coor[0]-coor[1];
								Xneig[3*ineig+1]=coor[0];
								Xneig[3*ineig+2]=coor[1];
							*/

						}
						break;
					} //!switch(Neig_sides[ineig]
				}
				break;
			}//!switch(el_type)
		}//!if
	}//!for

}

/*---------------------------------------------------------
  mmr_fa_area - to compute the area of face and vector normal
---------------------------------------------------------*/
void mmr_fa_area(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,		/* in: global face ID */
		double * Area,	       /* out: face area */
		double * Vec_norm     /* out: normal vector */
) {

	mmv_current_mesh->mmr_fa_area(Fa - 1, Area, Vec_norm);
}

/*---------------------------------------------------------
mmr_edge_status - to return edge status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_edge_status( /* returns edge status: */
		/* +1 (MMC_ACTIVE)   - active edge */
		/*  0 (MMC_FREE)     - free space */
		/* -1 (MMC_INACTIVE) - inactive (refined) edge */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed  	/* in: edge ID */
) {
	if (mmv_current_mesh->getNumberEdges() <= Ed - 1) {
		return 0;
	}
	return 1;
}

/*---------------------------------------------------------
mmr_edge_nodes - to return edge node's IDs
---------------------------------------------------------*/
int mmr_edge_nodes( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,		/* in: edge ID */
		int * Edge_nodes	/* out: IDs of edge nodes */
) {
	/*
		assert(Ed > 0);
		select_mesh(Mesh_id);
		Edge_nodes[0]=current_mesh->edges_[Ed].verts(0) ;
		Edge_nodes[1]=current_mesh->edges_[Ed].verts(1) ;
		assert(Edge_nodes[0] > 0);
		assert(Edge_nodes[1] > 0);
	*/
	//std::cout << "Method not implemented!: mmr_edge_nodes ";

	mmv_current_mesh->mmr_edge_nodes(Ed - 1, Edge_nodes);
	for (int i = 1, ileI = 3; i < ileI; ++i) {
		Edge_nodes[i] += 1;
	}
	return 1;
}

/*---------------------------------------------------------
mmr_edge_sons - to return edge son's numbers
---------------------------------------------------------*/
int mmr_edge_sons( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,		/* in: edge ID */
		int * Edge_sons,	/* out: IDs of edge sons */
		int * Edge_mid_node
) {

	std::cout << "Method not implemented!: mmr_edge_sons " << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_edge_structure - to return edge structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_edge_structure(/* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,  	   /* in: edge ID */
		int * Edge_struct /* out: edge structure in the form of integer array */
) {
	std::cout << "Method not implemented!: mmr_edge_structure " << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_set_edge_type - to set type for an edge
---------------------------------------------------------*/
int mmr_set_edge_type( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge_id,     /* in: edge ID */
		int Type         /* in: edge type (the number of attempted subdivisions !!!)*/
) {
	std::cout << "Method not implemented! mmr_set_edge_type" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_set_edge_fam - to set family data for an edge
---------------------------------------------------------*/
int mmr_set_edge_fam( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge_id,     /* in: edge ID */
		int Son1,        /* in: first son ID */
		int Son2         /* in: second son ID */
) {
	std::cout << "Method not implemented! mmr_set_edge_fam" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_set_face_fam - to set family data for an face
---------------------------------------------------------*/
int mmr_set_face_fam( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face_id,     /* in: face ID */
		int * Sons      /* in: sons IDs */
) {
	std::cout << "Method not implemented!" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_set_face_neig - to set neighbors data for a face
---------------------------------------------------------*/
int mmr_set_face_neig( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face_id,     /* in: face ID */
		int Neig1,       /* in: first neig ID */
		int Neig2,        /* in: second neig ID */
		int Neig1Type,    /** in: first neig type */
		int Neig2Type     /** in: second neig type */
) {
	std::cout << "Method not implemented!" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_set_elem_fath - to set family data for an elem
---------------------------------------------------------*/
int mmr_set_elem_fath( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Elem_id,     /* in: elem ID */
		int Fath        /* in: father ID */
) {
	std::cout << "Method not implemented!" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_set_elem_fam - to set family data for an elem
---------------------------------------------------------*/
int mmr_set_elem_fam( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Elem_id,     /* in: elem ID */
		int Fath,        /* in: father ID */
		int * Sons       /* in: sons IDs */
) {
	std::cout << "Method not implemented! mmr_set_elem_fam" << endl;
	exit(-1);
	return -1;
}

/*---------------------------------------------------------
mmr_node_status - to return node status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_node_status( /* returns node status: */
		/* +1 (MMC_ACTIVE)   - active node */
		/*  0 (MMC_FREE)     - free space */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node  	/* in: node ID */
) {

	return mmv_current_mesh->mmr_node_status(Node - 1);

}

/*---------------------------------------------------------
mmr_node_coor - to return node coordinates
---------------------------------------------------------*/
int mmr_node_coor( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Coor /* out: coordinates */
) {


	mmv_current_mesh->mmr_node_coor(Node - 1, Coor);

	return 3;
}

/* in file mms_prism_util.c:
mmr_loc_loc - to compute local coordinates within an element,
given local coordinates within an element of the same family
mmr_create_element - to create an element and fill its data structure
mmr_create_face - to create a new face
mmr_create_edge - to create a new edge structure
mmr_create_node - to create a node
mmr_clust_fa4_t - to cluster back a family of 4 triangular faces
mmr_clust_fa4_q - to cluster back a family of 4 quadrilateral faces
mmr_clust_edge  - to cluster two edges back into their father
mmr_divide_face4_t - to break a triangular face into 4 sons
mmr_divide_face4_q - to divide a quadrilateral face into 4 sons
mmr_divide_edge - to divide an edge into two sons

*/

/*---------------------------------------------------------
mmr_loc_loc - to compute local coordinates within an element,
given local coordinates within an element of the same family
---------------------------------------------------------*/
extern int mmr_loc_loc(/* returns: 1 - success, 0 - failure */
		int Mesh_id,   /* in: field ID */
		int El_from, 	/* in: element number */
		double * X_from, /* in: local element coordinates */
		int El_to, 	/* in: another element number */
		double * X_to	/* out: local another element coordinates */
) {
	std::cout << "Method not implemented! mmr_loc_loc" << endl;
	exit(-1);
	return -1;
}

extern int mmr_create_element(
		/* returns: ID of the created element (<=0 - failure) */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int  Type,     /* in: type for the face */
		int  Mate,     /* in: material indicator */
		int  Fath,     /* in: father element ID */
		int  Refi,     /* in: refinement type indicator */
		int * Faces,   /* in: list of faces' IDs */
		int * Sons     /* in: list of sons (only for inactive elements, Type<0) */
) {
	std::cout << "Method not implemented! mmr_create_element" << endl;
	exit(-1);
	return -1;
}


extern int mmr_create_face( /* returns: ID of created face (<0 - error code) */
		int  Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int  Type,     /* in: type for the face */
		int  Flag_bc,  /* in: bc flag for the face */
		int * Edges,   /* in: edges for the new face */
		int * Neig,    /* in (optional): neighbors (or NULL) */
		int * Sons     /* in (optional): sons (if Type<0) or NULL */
) {
	std::cout << "Method not implemented! mmr_create_face" << endl;
	exit(-1);
	return -1;
}


extern int mmr_create_edge( /* returns: ID of a new edge (<0 - error code)*/
		int  Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Type,      /* in: type indicator (MMC_EDGE or */
		/*     number of attempted divisions for inactive edges0 */
		int  Node1,    /* in: nodes for a new edge */
		int  Node2
) {
	std::cout << "Method not implemented! mmr_create_edge" << endl;
	exit(-1);
	return -1;
}


extern int mmr_create_node( /* returns: node ID of created node (<0 - error) */
		int  Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		double Xcoor,  /* in: coordinates of new node */
		double Ycoor,
		double Zcoor
) {
	std::cout << "Method not implemented! mmr_create_node" << endl;
	exit(-1);
	return -1;
}


extern int mmr_clust_fa4_t( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face,      /* in: face ID */
		int * Face_sons /* in (optional): face sons' IDs */
) {
	std::cout << "Method not implemented! mmr_clust_fa4_t" << endl;
	exit(-1);
	return -1;
}


extern int mmr_clust_fa4_q( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face,      /* in: face ID */
		int * Face_sons /* in (optional): face sons' IDs */
) {
	std::cout << "Method not implemented! mmr_clust_fa4_q" << endl;
	exit(-1);
	return -1;
}


extern int mmr_clust_edge( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge_id   /* in: clustered edge ID */
) {
	std::cout << "Method not implemented! mmr_clust_edge" << endl;
	exit(-1);
	return -1;
}



extern int mmr_divide_face4_t( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int  Fa,         /* in: face to be divided */
		int * Face_sons, /* out: face's sons */
		int * Sons_edges, /* out: created new edges */
		int * New_nodes  /* out: created new nodes */
) {
	std::cout << "Method not implemented! mmr_divide_face4_t" << endl;
	exit(-1);
	return -1;
}


extern int mmr_divide_face4_q(/* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int  Fa,         /* in: face to be divided */
		int * Face_sons, /* out: sons */
		int * Sons_edges, /* out: created new edges */
		int * New_nodes  /* out: created new nodes */
) {
	std::cout << "Method not implemented! mmr_divide_face4_q" << endl;
	exit(-1);
	return -1;
}


extern int mmr_divide_edge( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int  Edge,      /* in: edge to be divided */
		int * Edge_sons, /* out: two sons */
		int * Node_mid  /* out: node in the middle */
) {
	std::cout << "Method not implemented! mmr_divide_edge" << endl;
	exit(-1);
	return -1;
}



/*------------------------------------------------------------
mmr_del_elem    - to free an element structure
------------------------------------------------------------*/
int mmr_del_elem( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	/* in: ID of the mesh to be used or 0 for the current mesh */
		int Elem      /* in: element ID */
) {

	mmv_current_mesh->delElement(Elem - 1);
	return 1;
}

int mmr_del_face( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face      /* in: face ID */
) {

	//current_mesh->delFace(Face);
	std::cout << "Method not implemented! : mmr_del_face " << endl;
	exit(-1);
	return -1;
}

int mmr_del_edge( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge      /* in: edge ID */
) {
	//current_mesh->delEdge(Edge);
	std::cout << "Method not implemented! : mmr_del_face " << endl;
	exit(-1);
	return -1;

}

int mmr_del_node( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node      /* in: node ID */
) {
	mmv_current_mesh->delPoints(Node - 1);
	return 1;
}

/*---------------------------------------------------------
mmr_edge_elems - to return IDs of elements containing the edge
---------------------------------------------------------*/
extern int mmr_edge_elems( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,		/* in: edge ID */
		int * Edge_elems	/* out: IDs of elements containing the edge */
		/*      Edge_elems[0] - the number of elements */
) {

	int temp = mmv_current_mesh->mmr_edge_elems(Ed - 1, Edge_elems);
	for (int i = 1, ileI = Edge_elems[0] + 1; i < ileI; ++i) {
		Edge_elems[i] += 1;
	}

	return temp;

}

/*---------------------------------------------------------
mmr_el_edges - to get the list of element's edges
---------------------------------------------------------*/
int mmr_el_edges( 	/* returns the number of edges or error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	  /* in: element ID */
		int * Edges   	 /* out : list of edge IDs */
		/*	(Edges[0] - number of edges) */
) {
	//std::cout << "Method not implemented! : mmr_edge_elems ";

	int temp = mmv_current_mesh->mmr_el_edges(El - 1, Edges);
	for (int i = 1, ileI = Edges[0] + 1; i < ileI; ++i) {
		Edges[i] += 1;
	}
	return temp;
}



/*---------------------------------------------------------
  mmr_el_fa_nodes - to get list local face nodes indexes in elem
---------------------------------------------------------*/
extern int mmr_el_fa_nodes ( // returns face type flag
		int Mesh_id,	    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,		    // in: global elem ID
		int Fa,		    /* in: local face number in elem El */
		int * fa_nodes		/* out: list of local indexes of face nodes */
)
{

	int type = MMC_TRIA;

	switch (Fa) {
		case 0: {
			fa_nodes[0] = 3;
			fa_nodes[1] = 1;
			fa_nodes[2] = 0;
		}
		break;
		case 1: {
			fa_nodes[0] = 3;
			fa_nodes[1] = 2;
			fa_nodes[2] = 1;
		}
		break;
		case 2: {
			fa_nodes[0] = 3;
			fa_nodes[1] = 0;
			fa_nodes[2] = 2;
		}
		break;
		case 3: {
			fa_nodes[0] = 0;
			fa_nodes[1] = 1;
			fa_nodes[2] = 2;
		}
		break;
	}
	std::swap(fa_nodes[1], fa_nodes[2]);

	return type;
}

/*------------------------------------------------------------
  mmr_create_edge_elems - to create (or recreate if Max_edge_id > 0) for each
				  edge a list of elements to which it belongs
------------------------------------------------------------*/
int mmr_create_edge_elems( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int Max_edge_id /* in: the range of edge IDs to consider or 0 for default */
		/* if Max_edge_id==0 it is assumed that there are no structures to free !!! */
)
{

	std::cout << "Method not implemented! : mmr_create_edge_elems " << endl;
	exit(-1);
	return (-1);
}


/*---------------------------------------------------------
  mmr_get_max_gen - to set maximal allowed generation level for elements
------------------------------------------------------------*/
extern int mmr_get_max_gen(/* returns: >=0-success code, <0-error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	//std::cout << "Method not implemented! : mmr_get_max_gen ";
	//std::cout << "Method not implemented! : mmr_get_max_gen "<<endl;
	//exit(-1);
	return 0;
}

/*---------------------------------------------------------
  mmr_get_max_gen_diff - to get maximal allowed generation difference between
						 neighboring elements
------------------------------------------------------------*/
extern int mmr_get_max_gen_diff(
		/* returns: >=0 - maximal generation difference, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	//std::cout << "Method not implemented! : mmr_get_max_gen_diff "<<endl;
	//exit(-1);
	return 0;
}
/*
mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready for dofs projection
*/

extern int mmr_is_ready_for_proj_dof_ref(int Mesh_id)
{
	//std::cout << "Method not implemented! : mmr_is_ready_for_proj_dof_ref "<<endl;
	//exit(-1);
	return 0;
}

/*---------------------------------------------------------
  mmr_get_el_dist2boundary  - to return approx distance from point within element to nearst boundary
---------------------------------------------------------*/
double mmr_get_el_dist2boundary(  /* returns: dist>=0.0  */
		/*           <0 - "far far away from boundary" code */
		const int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		const int El_id,	 // in: El_id wihch includes point Coord
		const double * Coord) // GLOBAL coords of point belonging to El_id
{


	/*
	assert(Coord != NULL);
	assert(El_id >= FIRST);

	select_mesh(Mesh_id);
	const hHybridMesh & m(*current_mesh);
	const hObj & el(current_mesh->elements_[El_id]);

	double dist( -1.0 );

	// let us assume that near boundary are always prism elements
	if(el.type_ == ElemPrism::myType) {
	// if all vertices have the same "nearst boundary face"
	ID nearstFaces[ElemPrism::nVerts]={0};
	for(uTind i(0); i < el.typeSpecyfic_.nVerts_; ++i) {
	  nearstFaces[i]= m.vertices_[el.verts(i)].nearstFace_;
	}

	std::sort(nearstFaces,nearstFaces+ElemPrism::nVerts);
	const int nUniqueFaces(std::unique(nearstFaces,nearstFaces+ElemPrism::nVerts)-nearstFaces);

	assert(nUniqueFaces >= 1);

	dist=mmr_point_plane_dist(Coord,m.faces_(nearstFaces[0]).planeCoords_);
	for(int i(1); i < nUniqueFaces; ++i) {
	  double newDist = mmr_point_plane_dist(Coord,m.faces_(nearstFaces[i]).planeCoords_);
	  if(newDist < dist){
		dist = newDist;
	  }
	}
	}
	*/

	//un
	std::cout << "Method not implemented! : mmr_get_el_dist2boundary " << endl;
	exit(-1);
	return 0.01;
}


/*---------------------------------------------------------
  mmr_init_dist2boundary  - to return distance from given vertex to nearst boundary
---------------------------------------------------------*/
int mmr_init_dist2bound( // returns 0 if succesfull
		const int Mesh_id, //in: mesh ID or 0 (MMC_CUR_MESH_ID)
		const int * BCs, //in: c-array of accepted BC numbers
		const int nBCs) //in: length of BCs parameter
{
	//un
	//select_mesh(Mesh_id);
	//current_mesh->computeDist2Bound(BCs,nBCs);
	std::cout << "Method not implemented! : mmr_init_dist2bound " << endl;
	exit(-1);
	return 0;

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
//    /* auxiliary variables */
//      mmt_mesh* mesh;

//    /*++++++++++++++++ executable statements ++++++++++++++++*/

//    /* select the proper mesh data structure */
//      mesh = mmr_select_mesh(Mesh_id);

//    return hashlittle(mesh->node, mesh->parm.nrno, Mesh_id);

	mf_fatal_err("mmr_current_mesh_state_hash not implemented!");
	return 0;

}

int mmr_set_node_coor( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Xcoor  	/* in: node coordinates */
)
{
	mf_fatal_err("mmr_set_node_coor not implemented!");
	return 0;
}

#ifdef __cplusplus
}
#endif


//poniżej wylatuje
// -> i wylecialo.
