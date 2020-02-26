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
mmr_derefine - to derefine an element or the whole mesh
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
mmr_set_node_coor - set node coordinates
mmr_el_fa_nodes - to get list local face nodes indexes in elem

mmr_loc_loc - to compute local coordinates within an element,
given local coordinates within an element of the same family

mmr_divide_face4_t - to break a triangular face into 4 sons
mmr_divide_face4_q - to divide a quadrilateral face into 4 sons
mmr_divide_edge - to divide an edge into two sons

mmr_del_elem    - to free an element structure
mmr_del_face    - to free a face structure
mmr_del_edge    - to free an edge structure
mmr_del_node    - to free a node structure

/*------------------------------------------------------------
  mmr_test_mesh_motion - block of mmd_remesh interface
------------------------------------------------------------*/

/*------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
	07.2009 - Kazimierz Michalik, implementation
**	OC
***********************************************************************/

//////////////////////////////////////////////////////////////////////////
/// NOTE: inside arrays are numbered 0..n-1, outside arrays are numered 1..n !
/// ASSUMPTION: interface ID == internal POS+1
///				interface ID != internal ID
//////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>


/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */

#include "Common.h"
#include <modfem/mmh_intf.h>
#include "MeshModule/hHybridMesh.h"
#include "MeshModule/hHybridMeshWithContacts.h"
#include "MeshModule/ElemT4.h"
#include "MeshModule/ElemPrism.h"
//#include "../MeshWrite/IMeshWriter.h"
#include "MeshWrite/DmpFileExporter.h"
#include "MeshWrite/KazFileExporter.h"
#include "MeshWrite/BinaryFileWriter.h"
#include "MeshWrite/MileninExporter.h"
#include "MeshRead/DmpFileImporter.h"
#include "MeshRead/KazFileImporter.h"
#include "MeshRead/NasFileImporter.h"
#include "MeshRead/BinaryFileReader.h"
#include "MeshRead/InFileImporter.h"
#include "MeshModule/mmh_vec3.h"
#include <modfem/uth_log.h>
#include <modfem/uth_hash.h>


// Logging setup

std::ostringstream & mmv_out_stream = *( new std::ostringstream() );

void writeOutputStream()
{
	if (!mmv_out_stream.str().empty()) {
		mf_log_info("%s", mmv_out_stream.str().c_str());
	}
}

// Definition of global mesh variables.
std::vector<hHybridMeshWithContacts *>     mmv_meshes;

//hHybridMeshWithContacts * mmv_current_mesh = NULL;
int	mmv_current_mesh_id = 0;

void atExitClearMmv_Meshes()
{
	for (int i(0), end(mmv_meshes.size()); i < end; ++i) {
		safeDelete(mmv_meshes[i]);
	}
}


inline hHybridMeshWithContacts  & select_mesh(const int mesh_id) {
	//    assert(mesh_id >= 0);
	//  if (mesh_id != MMC_CUR_MESH_ID) {
	//      mmv_current_mesh = mmv_meshes.at(mesh_id);
	//      assert(m.meshId_ == mesh_id);
	//      //hObj::myMesh = mmv_current_mesh;
	//    }
	return *mmv_meshes[mesh_id];
}



#ifdef __cplusplus
extern "C" {
#endif

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
const int MMC_INIT_GEN_LEVEL = 0;  // value of 'initial mesh' generation level
const int MMC_NO_FATH       = 0;   /* no father indicator */
const int MMC_SAME_ORIENT   = 1;   /* indicator for the same orientation */
const int MMC_OPP_ORIENT   = -1;   /* indicator for the opposite orientation */
const int MMC_MAX_EDGE_ELEMS = 20;

const int MMC_FACE_NODES_FOR_TETRA[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
const int MMC_FACE_NODES_FOR_PRISM[5][4] = {{0, 1, 2, -1}, {3, 4, 5, -1}, {0, 1, 3, 4}, {1, 2, 4, 5}, {2, 0, 5, 3}};

const int MMC_IS_SUPPORTING_NEW_MPI = 1;

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

/*------------------------------------------------------------
  mmr_test_mesh_motion - test movement mesh
------------------------------------------------------------*/


int mmr_groups_number(const int Mesh_id) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.n_Groups();
}
int mmr_groups_ids(const int Mesh_id, int * tab, int l_tab) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.getGroups(tab);
}

int mmr_get_fa_el_bc_connect(const int Mesh_id, int face_id, int * el_id) {
	mf_fatal_err("Method not implemented!");
	return -1;
}

void mmr_get_coor_from_motion_element(const int Mesh_id, int idEl, int idLP, double * coor, int flagaSiatki) {
	mf_fatal_err("Method not implemented!");
}

void mmr_copyMESH(const int Mesh_id, int flaga) {
	mf_fatal_err("Method not implemented!");
}

void mmr_create_mesh_Cube(const int Mesh_id, const char * nazwa, int node_x, int node_y, int node_z, double size_x, double size_y, double size_z, double divide, int * warunki) {
	mf_fatal_err("Method not implemented!");
}

void mmr_init_all_change(const int Mesh_id, int a) {
	mf_fatal_err("Method not implemented!");
}

void mmr_test_mesh_motion(const int Mesh_id, int ileWarstw, int obecny_krok, int od_krok, int ileKrok, double minPoprawy, double px0, double py0, double pz0, double px1, double py1, double pz1, double endX, double endY, double endZ) {
	mf_fatal_err("Method not implemented!");
}

double mmr_test_weldpool(const int Mesh_id, double obecnyKrok, double krok_start, double minPoprawy, double doX, double doY, double dl, double zmiejszaPrzes, double limit, double szerokosc) {
	mf_fatal_err("Method not implemented!");
	return -1;
}


int mmr_finish_read(const int Mesh_id) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	m.normalizationProcessor();
	return 0;
}


/*------------------------------------------------------------
  mmr_test_mesh_motion - test movement mesh
------------------------------------------------------------*/

/*------------------------------------------------------------
  mmr_module_introduce - to return the mesh name
------------------------------------------------------------*/
int mmr_module_introduce(
		/* returns: >=0 - success code, <0 - error code */
		char * Mesh_name /* out: the name of the mesh */
) {
	//strcpy(Mesh_name,const_cast<char*>(current_mesh->name_.c_str()));    // must be, because of compatibility with C
	strcpy(Mesh_name, "3D_Hybrid");
	return 0;
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
		/* MMC_HP_FEM_MESH_DATA = 1 - mesh read from dump files */
		/* MMC_GRADMESH_DATA    = 2 - mesh produced by 2D GRADMESH mesh generator */
		char * Filename,   /* in: name of the file to read mesh data */
		/// KM 06.2011
		/// We assume that there is also a detailed boundary description file
		/// which name is similar to the input file, but it ends with ".detail".
		FILE * interactive_output
) {

	hHybridMeshWithContacts * m_ptr = NULL;

	// Create new mesh object.
	try {
		if (mmv_meshes.empty()) {
			atexit(writeOutputStream);      // flush info at exit
			atexit(atExitClearMmv_Meshes);  // assure freeing resources at program exit
		}
		m_ptr = new hHybridMeshWithContacts();
		if (m_ptr->meshId_ >= static_cast<int>(mmv_meshes.size())) {
			mmv_meshes.resize(m_ptr->meshId_ + 2, NULL);
		}
		assert(mmv_meshes[m_ptr->meshId_] == NULL);
		mmv_meshes[m_ptr->meshId_] = m_ptr;
		//hObj::myMesh = mmv_current_mesh;
	} catch (...) {
		std::cerr << "\n Error! Something gone wrong while creating hybrid mesh. Contact support.";
		return -1;
	}

	MeshRead::IMeshReader * importer(NULL);
	char * dot = strrchr(Filename, '.');
	++dot;

	switch (Control) {

		case MMC_MOD_FEM_PRISM_DATA:  /* prism mesh read from dump files */
		case MMC_MOD_FEM_MESH_DATA:  /* generic mesh for this module */
			importer = new MeshRead::DmpFileImporter(Filename);
			break;
		case MMC_MOD_FEM_TETRA_DATA: /* tetra mesh read from dump files */
		case MMC_MOD_FEM_HYBRID_DATA:  /* hybrid mesh read from dump files */
		{
			if (strcmp(dot, "nas") == 0) {
				importer = new MeshRead::NasFileImporter(Filename);
			}
			else {
				importer = new MeshRead::KazFileImporter(Filename);
			}
		}
		break;
		case MMC_NASTRAN_DATA:
			importer = new MeshRead::NasFileImporter(Filename);
			break;
		case MMC_GRADMESH_DATA:  /* mesh produced by 2D GRADMESH generator */
			mmv_out_stream << "Method not implemented!";
			exit(-1);
			break;
		case MMC_BINARY_DATA:
			importer = new MeshRead::BinaryFileReader(Filename);
			break;
		case MMC_IN_ANSYS_DATA:
			importer = new MeshRead::InFileImporter(Filename);
			break;
		default:
			throw "mmr_init_mesh: unknown mesh type!";
			break;
	}//!switch(Control)

	assert(importer != NULL);

	try {
		m_ptr->read(*importer);
	}
	catch (const std::string & e) {
		std::cerr << "\nError: "
				<< e
				<< " , while reading a mesh file: "
				<< importer->Name() << "\n"
				<< "Make sure that file format is correct and input files are also correct.\n"
				<< "Note that you are using hybrid mesh!" << std::endl;

	}
	catch (...) {
		std::cerr << "\nError! Something gone wrong while reading a mesh file: "
				<< importer->Name() << "\n"
				<< "Make sure that file format is correct and input files are also correct.\n"
				<< "Note that you are using hybrid mesh!" << std::endl;
		throw;
	}
	delete importer;

	assert(mmv_meshes[m_ptr->meshId_] == m_ptr);

	// to flush all of output (or errors)
	mmv_out_stream << std::endl;
	writeOutputStream();
	return m_ptr->meshId_;
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
	hHybridMeshWithContacts * m_ptr = NULL;
	// Create new mesh object.
	try {
		if (mmv_meshes.empty()) {
			atexit(atExitClearMmv_Meshes);  // assure freeing resources at program exit
		}
		m_ptr = new hHybridMeshWithContacts();
		if (m_ptr->meshId_ >= static_cast<int>(mmv_meshes.size())) {
			mmv_meshes.resize(m_ptr->meshId_ + 2, NULL);
		}
		assert(mmv_meshes[m_ptr->meshId_] == NULL);
		mmv_meshes[m_ptr->meshId_] = m_ptr;
		//hObj::myMesh = mmv_current_mesh;
	} catch (...) {
		std::cerr << "\n Error! Something gone wrong while creating hybrid mesh. Contact support.";
		return -1;
	}

	mmr_init_read(m_ptr->meshId_, N_nodes, N_edges, N_faces, N_elems);

	assert(mmv_meshes[m_ptr->meshId_] == m_ptr);

	// to flush all of output (or errors)
	mmv_out_stream << std::endl;
	writeOutputStream();
	return m_ptr->meshId_;
}

//---------------------------------------------------------
// mmr_reserve - to inform mesh module about new target number of elems,faces,edges and nodes
//------------------------------------------------------------
// NOTE: this is NOT a resize function. It does NOT ADD any elements. It only reserves resources.
// To add mesh entities use mmr_add_* functions.
//---------------------------------------------------------
int mmr_init_read( // return >=0 - success, <0 - error code
		const int Mesh_id, // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int N_nodes, //IN: target count of nodes
		const int N_edges, //IN: target count of edges
		const int N_faces, //IN: target count of faces
		const int N_elems  //IN: target count of elements
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	assert(N_nodes >= 0);
	assert(N_edges >= 0);
	assert(N_faces >= 0);
	assert(N_elems >= 0);

	if (static_cast<uTind>(N_nodes) > m.vertices_.capacity()) {
		m.vertices_.reserve(N_nodes);
	}

	if (static_cast<uTind>(N_edges) > m.edges_.capacity()) {
		m.edges_.reserve(N_edges, N_edges * sizeof(Edge));
	}

	if (static_cast<uTind>(N_faces) > m.faces_.capacity()) {
		m.faces_.reserve(N_faces, N_faces * sizeof(Face4));
	}

	if (static_cast<uTind>(N_elems) > m.elements_.capacity()) {
		m.elements_.reserve(N_elems, N_elems * sizeof(ElemPrism));
	}

	return 0;
}

//---------------------------------------------------------
//  mmr_add_elem - to add (append) new element
//---------------------------------------------------------
int mmr_add_elem(// returns: added element id
		const int Mesh_id, // IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
		const int El_type, // IN: type of new element
		const int El_nodes[6], //IN: if known: nodes of new element, otherwise NULL
		const int El_faces[5],  //IN: if known: faces of new element, otherwise NULL
		const int Material_idx)
// NOTE: if both El_nodes and El_faces are NULLs
// then an element is udefined and an error is returned!
{
	assert((El_nodes != NULL) || (El_faces != NULL));
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	Elem * el(NULL);

	if (El_nodes != NULL) {
		if (El_type == MMC_PRISM) {
			el = m.elements_.newObj<ElemPrism>(&m, reinterpret_cast<const uTind *>(El_nodes));
		}
		else if (El_type == MMC_TETRA) {
			el = m.elements_.newObj<ElemT4>(&m, reinterpret_cast<const uTind *>(El_nodes));
		}

		m.wasNormalized_ = false;
		mf_check(el != NULL, "Unknown element type at element creation (%d)", El_type);
	}
	else if (El_faces != NULL) {
		mf_log_err("Not implemented!");
	}
	else {
		mf_log_err("To few data to add element!");
	}



	el->flags(GROUP_ID) = Material_idx;

	return el->pos_;
}
//---------------------------------------------------------
//  mmr_add_face - to add (append) new face
//---------------------------------------------------------
int mmr_add_face(// returns: added face id
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
	assert( (Fa_nodes != NULL) || (Fa_edges != NULL) );

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	Face * face(NULL);
	if (Fa_nodes != NULL) {
		if ( Fa_type == MMC_QUAD ) {
			face = m.faces_.newObj<Face4>(&m, reinterpret_cast<const uTind *>(Fa_nodes));
		}
		else if (Fa_type == MMC_TRIA) {
			const uTind nodes[4] = {static_cast<uTind>(Fa_nodes[0]), static_cast<uTind>(Fa_nodes[1]), static_cast<uTind>(Fa_nodes[2]), UNKNOWN};
			face = m.faces_.newObj<Face3>(&m, nodes);
		}
		m.wasNormalized_ = false;
	}
	if (Fa_nodes == NULL && Fa_neigs != NULL) {
		mf_log_err("Not implented!");
	}

	if (Fa_neigs != NULL) {
		face->neighs(0) = Fa_neigs[0];
		face->neighs(1) = Fa_neigs[1];
	}
	else {
		face->neighs(0) = UNKNOWN;
		face->neighs(1) = UNKNOWN;
	}

	face->flags(B_COND) = B_cond_val;

	return face->pos_;
}

//---------------------------------------------------------
//  mmr_add_edge - to add (append) new edge
//---------------------------------------------------------
int mmr_add_edge(// returns: added edge id
		const int Mesh_id,  //IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
		const int nodes[2]  //IN: nodes for new edge
)
// NOTE: if nodes==NULL edge is undefined and an error is returned!
{
	assert(nodes != NULL);

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	Edge & edge = * m.edges_.newObj<Edge>(&m, reinterpret_cast<const uTind *>(nodes), NULL);
	m.wasNormalized_ = false;

	return edge.pos_;
}

//---------------------------------------------------------
//  mmr_add_node - to add (append) new node
//---------------------------------------------------------
int mmr_add_node( // returns: added node id
		const int Mesh_id,  //IN: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh
		const int Id, // IN: MMC_AUTO_GENERATE_ID or id numebr if known
		const double Coords[3] //IN: geometrical coordinates of node
)
// NOTE: if Coords==NULL node is undefined and an error is returned!
{
	assert(Coords != NULL);

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	int pos = 0;

	if ( (Id > 0) && (Id < m.vertices_.last()) ) {
		Vertex & node =  m.vertices_[Id];
		if ( node.coords_[0] != Coords[0]
				|| node.coords_[1] != Coords[1]
				|| node.coords_[2] != Coords[2]) {
			throw "Error in mmr_add_node: trying add node, but node with this id exist and it is diffrent!";
		}
	}
	else   if (Id == MMC_AUTO_GENERATE_ID) {
		Vertex & node = * m.vertices_.newObj<Vertex>(&m, NULL);
		m.wasNormalized_ = false;
		node.coords_[0] = Coords[0];
		node.coords_[1] = Coords[1];
		node.coords_[2] = Coords[2];
		pos = node.pos_;
	}
	mf_check(pos >= FIRST, "Elem number out of bounds!");
	return pos;
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
		hHybridMeshWithContacts & m = select_mesh(Mesh_id);

		MeshWrite::IMeshWriter * writer(NULL);
		switch (Control) {
			case MMC_BINARY_DATA:
				writer = new MeshWrite::BinaryFileWriter(Filename);
				break;
			//case MMC_NASTRAN_DATA:
			//writer = new MeshWrite::Nas
			// break;
			case MMC_BOUND_VERTS_DATA:
				writer = new MeshWrite::MileninExporter(Filename);
				break;
			case MMC_HP_FEM_MESH_DATA:
			default: {
				writer = new MeshWrite::KazFileExporter(Filename);
			}
			break;
		}

		m.write(*writer);
		mmv_out_stream << " \nMesh exported to " << Filename;

		delete writer;
		return 0;
	}

	std::cerr << "\n Error: export_mesh failed, no filename!";


	return -1;
}

/*---------------------------------------------------------
  mmr_test_mesh - to test the integrity of mesh data
---------------------------------------------------------*/
int mmr_test_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

//    mmv_out_stream <<"Creating hHybridMesh with sizeof(hObj)="<<sizeof(hObj)<<
//          ",\n" <<" sizeof(hBase)="<<sizeof(hBase)<<
//          ",\n" <<" sizeof(Vertex)="<<sizeof(Vertex)<<
//           ",\n" <<" sizeof(Edge)="<<sizeof(Edge)<<
//          ", sizeof(EdgeD)="<<sizeof(EdgeD)<<
//          ",\n"<<" sizeof(Face3)="<<sizeof(Face3)<<
//          ", sizeof(Face3D)="<<sizeof(Face3D)<<
//          ",\n"  << " sizeof(Face4)="<<sizeof(Face4)<<
//          ", sizeof(Face4D)="<<sizeof(Face4D)<<
//          ",\n" << " sizeof(ElemT4)="<<sizeof(ElemT4)<<
//          ", sizeof(ElemT4D)="<<sizeof(ElemT4D)<<
//          ",\n" <<" sizeof(ElemPrism)="<<sizeof(ElemPrism)<<
//          ", sizeof(ElemPrismD)="<<sizeof(ElemPrismD) <<
//                     "\n" << "sizeof(hHybridMesh)="<<m.totalSize();


	//m.printSetup();
	mf_log_info("Test mesh.");
	return m.checkAll() ? 0 : -1;
}


/*---------------------------------------------------------
  mmr_get_nr_elem - to return the number of active elements
---------------------------------------------------------*/
int mmr_get_nr_elem(/* returns: >=0 - number of active elements, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	return m.elements_.nonDividedObjs();
}

/*---------------------------------------------------------
  mmr_get_max_elem_id - to return the maximal element id
---------------------------------------------------------*/
int mmr_get_max_elem_id(  /* returns: >=0 - maximal element id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	assert(m.elements_.size() >= 0);
	return m.elements_.last();
}

/*---------------------------------------------------------
  mmr_get_max_elem_max - to return the maximal possible element id
---------------------------------------------------------*/
int mmr_get_max_elem_max(  /* returns: >=0 - maximal element id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	return mmr_get_max_elem_id(Mesh_id);
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
	assert(Nel >= 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	// Because we have only 'active' elements, there is only need of checking is 'Nel' last element;
	// if not return next element id;
	// if Nel is last, next hObj is filled with zeros, so id will be zero
	if (Nel > 0) {
		Nel = m.elements_[Nel].next()->pos_;
	}
	else { // Nel==0
		if (m.elements_.size() > 0) {
			Nel = m.elements_.first().pos_ ;
			//Nel =  FIRST;
		}
		// otherwise  Nel=0
	}

	while (Nel >= static_cast<int>(FIRST) && m.elements_[Nel].isBroken()) {
		Nel = m.elements_[Nel].next()->pos_;
	}
	assert(Nel >= 0 && Nel <= static_cast<int>(m.elements_.last()));
	return Nel;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.faces_.nonDividedObjs();
}

/*---------------------------------------------------------
  mmr_get_max_face_id - to return the maximal face id
---------------------------------------------------------*/
int mmr_get_max_face_id(  /* returns: >=0 - maximal face id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.faces_.last();
}

/*---------------------------------------------------------
  mmr_get_max_face_max - to return the maximal possible face id
---------------------------------------------------------*/
int mmr_get_max_face_max(  /* returns: >=0 - maximal possible face id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	return mmr_get_max_face_id(Mesh_id);
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
	assert(Nfa >= 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	// Because we have only 'active' elements, there is only need of checking is 'Nel' last element;
	// if not return next element id;
	// if Nfa is last, next hObj is filled with zeros, so id will be zero
	if (Nfa > 0) {
		Nfa = m.faces_[Nfa].next()->pos_;
	}
	else {// Nfa==0
		if (m.faces_.size() > 0) {
			Nfa = m.faces_.first().pos_ ;
			//Nfa = FIRST;
		}
		// otherwise  Nfa=0
	}

	while (Nfa >= static_cast<int>(FIRST) && m.faces_[Nfa].isBroken()) {
		Nfa = m.faces_[Nfa].next()->pos_;
	}
	assert(Nfa >= 0 && Nfa <= static_cast<int>(m.faces_.last()));
	return Nfa;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.edges_.nonDividedObjs();
}

/*---------------------------------------------------------
  mmr_get_max_edge_id - to return the maximal edge id
---------------------------------------------------------*/
int mmr_get_max_edge_id(  /* returns: >=0 - maximal edge id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.edges_.last();
}

/*---------------------------------------------------------
  mmr_get_max_edge_max - to return the maximal possible edge id
---------------------------------------------------------*/
int mmr_get_max_edge_max(  /* returns: >=0 - maximal possible edge id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	return mmr_get_max_edge_id(Mesh_id);
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
	assert(Ned >= 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	// Because we have only 'active' elements, there is only need of checking is 'Nel' last element;
	// if not return next element id;
	// if Nel is last, next hObj is filled with zeros, so id will be zero
	if (Ned > 0) {
		Ned = m.edges_[Ned].next()->pos_;
	}
	else {// Ned==0
		if (m.edges_.size() > 0) {
			//Ned = FIRST;
			Ned = m.edges_.first().pos_ ;
		}
		// otherwise  Nfa=0
	}

	while (Ned >= static_cast<int>(FIRST) && m.edges_[Ned].isBroken()) {
		Ned = m.edges_[Ned].next()->pos_;
	}
	assert(Ned >= 0 && Ned <= static_cast<int>(m.edges_.last()));
	return Ned;
}

/*---------------------------------------------------------
  mmr_get_nr_node - to return the number of active nodes
---------------------------------------------------------*/
int mmr_get_nr_node(/* returns: >=0 - number of active nodes, */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	return static_cast<int>(m.vertices_.size());
}

/*---------------------------------------------------------
  mmr_get_max_node_id - to return the maximal node id
---------------------------------------------------------*/
int mmr_get_max_node_id(  /* returns: >=0 - maximal node id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return static_cast<int>(m.vertices_.last());
}

/*---------------------------------------------------------
mmr_get_max_node_max - to return the maximal possible vertex id
---------------------------------------------------------*/
extern int mmr_get_max_node_max(  /* returns: >=0 - maximal possible vertex id */
		/*           <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.vertices_.last() * 1000;
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
	assert(Nno >= 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	if (Nno > 0) {
		Nno = m.vertices_[Nno].next()->pos_;
	}
	else {// Nno==0
		if (m.vertices_.size() > 0) {
			//Nno = FIRST;
			Nno = m.vertices_.first().pos_ ;
		}
		// otherwise  Nno=0
	}

	while (Nno >= static_cast<int>(FIRST) && m.vertices_[Nno].isBroken()) {
		Nno = m.vertices_[Nno].next()->pos_;
	}

	assert(Nno >= 0 && Nno <= static_cast<int>(m.vertices_.last()));
	return Nno;
}

/*---------------------------------------------------------
  mmr_set_max_gen - to set maximal generation level for element refinements
------------------------------------------------------------*/
int mmr_set_max_gen(/* returns: >=0-success code, <0-error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Max_gen   /* in: maximal generation */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	m.maxGen_ = Max_gen;

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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	// Mesh is 1-unregular
	return m.maxGenDiff_ == 1 ? 0 : -1;
}

/*------------------------------------------------------------
  mmr_init_ref - to initialize the process of refinement
------------------------------------------------------------*/
int mmr_init_ref(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.initRefine();
}


/*------------------------------------------------------------
  mmr_refine_el - to refine an element WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_refine_el( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El   /* in: element ID  */
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	m.refineElem(El);
	return 0;
}

/*------------------------------------------------------------
mmr_derefine_el - to derefine an element WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_derefine_el( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID  */
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	m.derefineElem(El);
	return 0;
}

/*------------------------------------------------------------
  mmr_refine_mesh - to refine the WHOLE mesh WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_refine_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	m.refine();
	return 0;
}

/*------------------------------------------------------------
mmr_derefine_mesh - to derefine the WHOLE mesh WITHOUT IRREGULARITY CHECK
------------------------------------------------------------*/
int mmr_derefine_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	m.derefine();
	return 0;
}


/*------------------------------------------------------------

  mmr_refine - to refine an element or the whole mesh
------------------------------------------------------------*/
// int mmr_refine( /* returns: >=0 - success code, <0 - error code */
//     int Mesh_id,/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
//     int El  /* in: element ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
// ) {


//     if(El==MMC_DO_UNI_REF) {
// 	  current_mesh->refine();
// 	}
// 	else {
// 	  current_mesh->refineElem(current_mesh->elements_[El]);
// 	}
// 	//out_stream << "\nDividing element " << El;
//     return 0;
// }

/*------------------------------------------------------------
mmr_r_refine - to r-refine an elements with given boundary condition
------------------------------------------------------------*/
extern int mmr_gen_boundary_layer( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Bc, /* in: boundary conditon ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
		int thicknessProc,
		int noLayers,
		int distrib,
		double * ignoreVec
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	return m.createBoundaryLayer(thicknessProc, noLayers, distrib == 2, ignoreVec) ? 0 : -1;
}

/*------------------------------------------------------------
mmr_r_refine - to r-refine an elements with given boundary condition
------------------------------------------------------------*/
extern int mmr_r_refine( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Bc,  /* in: boundary conditon ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
		void (*reallocation_func)(double * x, double * y, double * z)	/* function moving to new*/
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	return m.rRefine(Bc, reallocation_func);
}

/*------------------------------------------------------------
  mmr_derefine - to derefine an element or the whole mesh
------------------------------------------------------------*/

int mmr_derefine( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID or -2 (MMC_DO_UNI_DEREF) for uniform derefinement */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	bool result = (El == MMC_DO_UNI_DEREF ? m.derefine() : m.derefineElem(m.elements_.at(El)));
	return result ? 0 : -1;
}


/*------------------------------------------------------------
  mmr_final_ref - to finalize the process of refinement after rewriting DOFs
------------------------------------------------------------*/
int mmr_final_ref(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	// mark all inactive as active  -ALL ARE ALWAYS ACTIVE
	//const Vertex * const end(current_mesh->vertices_.end());
	//for(hHybridMesh::VertexPool::Iterator it(& current_mesh->vertices_.first()); it < end; ++it) {
	//    if(it->status() == MMC_INACTIVE) {
	//        it->status(MMC_ACTIVE);
	//    }
	//}

	for (hHybridMesh::EdgePool::Iterator<> it(& m.edges_); !it.done(); ++it) {
		assert(it->test(&m));
//        if(it->status() == MMC_INACTIVE && (!it->isBroken())) {
//            it->status(MMC_ACTIVE);
//        }
	}

	for (hHybridMesh::FacePool::Iterator<> it(& m.faces_); !it.done(); ++it) {
		assert(it->test(&m));
//        if(it->status() == MMC_INACTIVE && (!it->isBroken())) {
//            it->status(MMC_ACTIVE);
//        }
	}

	for (hHybridMesh::ElemPool::Iterator<> it(& m.elements_); !it.done(); ++it) {
		assert(it->test(&m));
//        if(it->status() == MMC_INACTIVE && (!it->isBroken())) {
//            it->status(MMC_ACTIVE);
//        }
	}

	return mmr_test_mesh(Mesh_id);
}

/*---------------------------------------------------------
  mmr_free_mesh - to free space allocated for mesh data structure
---------------------------------------------------------*/
int mmr_free_mesh(  /* returns: >=0 - success code, <0 - error code */
		int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.free() ? 0 : -1;
}


/*---------------------------------------------------------
mmr_elem_structure - to return elem structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_elem_structure( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	   /* in: elem ID */
		int * Elem_struct /* out: elem structure in the form of integer array */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	const hObj & elem(m.elements_[El]);
	Elem_struct[0] = elem.type_;
	//Elem_struct[1]=mesh->elem[El].ipid;
	Elem_struct[2] = elem.flags(GROUP_ID);
	Elem_struct[3] = elem.parent_;
	Elem_struct[4] = elem.flags(REFINEMENT);
	for (unsigned int i(0); i < elem.typeSpecyfic_.nComponents_; ++i) {
		Elem_struct[5 + i] = hObj::posFromId(elem.components(i));
	}

	// if(elem.isBroken())
	// {
	//for(i=0;i<8;i++)
	//{
	//  Elem_struct[11+i]=mesh->elem[El].sons[i];
	//}
	// }

	return (0);

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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	assert(El >= FIRST);
	const Elem & el(m.elements_.at(El));
	return (el.isBroken() || el.isMarkedBreak()) ? MMC_INACTIVE : MMC_ACTIVE;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	assert(El >= FIRST);
	return m.elements_.at(El).type_;
}

/*---------------------------------------------------------
  mmr_el_groupID - to return material ID for element
---------------------------------------------------------*/
int mmr_el_groupID( /* returns material flag for element, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
) {
	assert(El > 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return (El > 0 ? m.elements_.at(El).flags(GROUP_ID) : 0);
}
/*---------------------------------------------------------
mmr_el_set_mate - to set material number for element
---------------------------------------------------------*/
int mmr_el_set_groupID( /* sets group flag for element */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	/* in: element ID */
		int Group_id    /* in: group ID */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	m.elements_.at(El).flags(GROUP_ID) = Group_id;
	return 1;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.elements_.at(El).isBroken() ? MMC_REF_ISO : MMC_NOT_REF;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	const hObj & el(m.elements_[El]);
	if (Faces != NULL) {
		Faces[0] = el.typeSpecyfic_.nComponents_;
		for (int i(0); i < *Faces; ++i) {
			Faces[1 + i] = hObj::posFromId(el.components(i));
		}
	}



	if (Orient != NULL) {
		BYTE mask(1),
			 type(static_cast<BYTE>(el.flags(EL_TYPE)));
		for (uTind i(0); i < el.typeSpecyfic_.nComponents_; ++i) {
			// if bit for given face is 1 then face points INSIDE and should be flipped.
			Orient[i] = (type & mask) == 0 ? 1 : -1;
			mask <<= 1; // shift mask to next face
		}
	}

	if (el.type_ == ElemPrism::myType) {
		std::swap(Faces[4], Faces[5]); // normalize face order outside mmr
		if (Orient != NULL) {
			std::swap(Orient[3], Orient[4]);
		}
	}

	return 0;
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
	assert(El > 0);

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	const uTind & no_nodes  = m.elements_.at(El).typeSpecyfic_.nVerts_;
	if (Nodes != NULL) {
		Nodes[0] = no_nodes;
		memcpy(Nodes + 1, m.elements_[El].verts(), sizeof(ID)*no_nodes);
	}
	if (Node_coor != NULL) {
		const ID * vts_ptr = m.elements_[El].verts();
		for (int i = 0; i < no_nodes; ++i, ++vts_ptr) {
			memcpy(Node_coor + (3 * i), m.vertices_.at(*vts_ptr).coords_.getArray(), 3 * sizeof(double));
		}
	}


	return no_nodes;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	const hObj & el(m.elements_[El]);
	if (Elsons != NULL)
		if (el.isBroken()) {
			Elsons[0] = el.nMyClassSons_;
			if (el.type_ == ElemT4::myType) {
				for (int i(0); i < el.nMyClassSons_; ++i) {
					Elsons[i + 1] = el.getMyClassChild<ElemT4>(i)->pos_;
				}
			}
			else {
				for (int i(0); i < el.nMyClassSons_; ++i) {
					Elsons[i + 1] = el.getMyClassChild<ElemPrism>(i)->pos_;
				}
			}
		}

	if (Type != NULL) {
		*Type = el.isBroken() ? MMC_REF_ISO : MMC_NOT_REF;
	}
	const int parentPos = (el.level_ > 0 ? static_cast<int>(hObj::posFromId(el.parent_)) : MMC_NO_FATH);
	assert(parentPos >= 0);
	return parentPos;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	const Elem * el(& m.elements_[El]);
	const int parent(el->parent_);

	if ( Elsons != NULL ) {
		Elsons[0] = 0;

		const int level(el->level_);


		el = el->next();
		for (int i = 1; el->level_ > level; ++i)  {
			++(Elsons[0]);
			Elsons[i] = el->pos_;
			el = el->next();
		}
	}
	return parent;
}


/*---------------------------------------------------------
  mmr_el_gen - to return generation level for an element
---------------------------------------------------------*/
int mmr_el_gen(   /* returns: El's generation ID */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El        /* in: element ID (1...) */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.elements_[El].level_;
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
	mmv_out_stream << "Method not implemented!";
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

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	double hsize(0.0);
	if (El > 0) {
		hsize = m.elemHSize(m.elements_[El]);
	}
	return (hsize);
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
	assert(Neig != NULL);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	const hObj & elem( m.elements_[El]);
	Neig[0] = elem.typeSpecyfic_.nComponents_;
	for (register int n(1); n <= Neig[0]; ++n) {
		Neig[n] = hObj::posFromId(m.elemNeigh(elem, n - 1));
	}
	assert(Neig[0] > 0);
	return Neig[0];
}

/*---------------------------------------------------------
mmr_face_structure - to return face structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_face_structure(/* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,  	   /* in: face ID */
		int * Face_struct /* out: face structure in the form of integer array */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	int noComponents = m.faces_[Fa].typeSpecyfic_.nComponents_;
	const hObj & face(m.faces_[Fa]);
	if (Face_struct != NULL) {
		Face_struct[0] = face.type_;
		//Face_struct[1]=mesh->face[Fa].ipid;
		Face_struct[2] = face.flags(B_COND);

		for (int i = 0; i < noComponents; ++i) {
			Face_struct[3 + i] = hObj::posFromId(m.faces_[Fa].components(i));
		}

		if (face.flags(F_TYPE) == F_FLIPPED) {
			std::swap(Face_struct[3 + 1], Face_struct[3 + 2]);
		}

		Face_struct[7] = face.neighs(0);
		Face_struct[8] = face.neighs(1);
		//  if(face.isBroken())
		//  {
		//	for(int i=0;i<4;i++)
		//	{
		//	  Face_struct[9+i]=mesh->face[Fa].sons[i];
		//	}
		//}
	}

	return 0;
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
	assert(Fa > 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	int result(MMC_ACTIVE);
	if (Fa <= static_cast<int>(m.faces_.last()))
		if (m.faces_[Fa].isBroken())
			result = MMC_INACTIVE;

	return result;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.faces_[Fa].type_;
}

/*--------------------------------------------------------------------------
  mmr_fa_bc - to get the boundary condition flag for a face
---------------------------------------------------------------------------*/
int mmr_fa_bc( /* returns: bc flag for a face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa	      /* in: global face ID */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.faces_[Fa].flags(B_COND);
}

/*--------------------------------------------------------------------------
  mmr_fa_bc - to set the boundary condition flag for a face
---------------------------------------------------------------------------*/
int mmr_fa_set_bc( /* returns: bc flag for a face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,	      /* in: global face ID */
		int BC_num
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	return m.faces_[Fa].flags(B_COND) = BC_num;
}


/*---------------------------------------------------------
mmr_fa_sub_bnd - to indicate face is on the boundary
---------------------------------------------------------*/
extern int mmr_fa_sub_bnd( /* returns: 1 - true, 0 - false */
		int Mesh_id,    /* in: mesh ID */
		int Face_id     /* in: face ID */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.faces_[Face_id].neighs(1) == MMC_SUB_BND ? 1 : 0;
}

/*---------------------------------------------------------
mmr_fa_set_sub_bnd - to indicate face is on the boundary
---------------------------------------------------------*/
extern int mmr_fa_set_sub_bnd(/* returns: >=0 - success code, <0 - error code */
		int Mesh_id,    /* in: mesh ID */
		int Face_id,    /* in: face ID */
		int Side_id     /* in: side ID ??? */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	mfp_check_debug( ( (Side_id == 0) || (Side_id == 1) ), "Wrong side(%d) id in mmr_fa_set_sub_bnd", Side_id);
	//mfp_debug("Setting MMC_SUB_BND on %d (instead of %d)", Face_id, hObj::posFromId(m.faces_[Face_id].neighs(Side_id)));
	m.faces_[Face_id].neighs(Side_id) = MMC_SUB_BND;
	//m.faces_[Face_id].flags(B_COND)=Side_id;
	return 0;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	const hObj & face( m.faces_[Fa]);
	const int noEdges( face.typeSpecyfic_.nComponents_ );

	if (Fa_edges != NULL) {
		for (int i(0); i < noEdges; ++i) {
			Fa_edges[i] = hObj::posFromId(face.components(i));
		}
	}

	if (Ed_orient != NULL) {
		for (int i(0); i < noEdges; ++i) {
			Ed_orient[i] = 1;
		}
	}

	return noEdges;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	assert(Fa >= FIRST);
	assert(Fa_neig != NULL);

	const Face & f = m.faces_[Fa];
	Fa_neig[0] = hObj::posFromId(f.neighs(0));
	Fa_neig[1] = hObj::posFromId(f.neighs(1));

	if (Neig_sides != NULL) {
		Neig_sides[0] = m.whichFaceAmI(f, m.elements_[Fa_neig[0]]);
		assert(Neig_sides[0] >= 0);
		assert(Neig_sides[0] <= MMC_MAXELFAC);
		Neig_sides[1] = m.whichFaceAmI(f, m.elements_[Fa_neig[1]]);
		assert(Neig_sides[1] >= 0);
		assert(Neig_sides[1] <= MMC_MAXELFAC);
	}

	if (Node_shift != NULL) {
		mf_log_err("Not implemented!");
	}

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
	assert(Fa > 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	assert(Fa_neig != NULL);

	const hObj	& face(m.faces_.at(Fa));
	Fa_neig[0] = hObj::posFromId(face.neighs(0)) ;
	Fa_neig[1] = hObj::posFromId(face.neighs(1)) ;

	if (Neig_sides != NULL) {
		Neig_sides[0] = (Fa_neig[0] > B_COND) ? m.whichFaceAmI(face, m.elements_.getById(face.neighs(0) )) : 0;
		Neig_sides[1] = (Fa_neig[1] > B_COND) ? m.whichFaceAmI(face, m.elements_.getById(face.neighs(1) )) : 0;

		// HACK to conform with KB prism faces ordering
		if (Fa_neig[0] > B_COND) {
			if (m.elements_(face.neighs(0)).type_ == ElemPrism::myType) {
				switch (Neig_sides[0]) {
					case 3:
						Neig_sides[0] = 4;
						break;
					case 4:
						Neig_sides[0] = 3;
						break;
				}
			}
		}
		if (Fa_neig[1] > B_COND) {
			if (m.elements_(face.neighs(1)).type_ == ElemPrism::myType) {
				switch (Neig_sides[1]) {
					case 3:
						Neig_sides[1] = 4;
						break;
					case 4:
						Neig_sides[1] = 3;
						break;
				}
			}
		}
		// end of HACK
	}

	if (Node_shift != NULL) *Node_shift = 0;
	if (Diff_gen != NULL) *Diff_gen = m.elements_[Fa_neig[0]].level_ - m.elements_[Fa_neig[1]].level_;
	assert(Fa_neig[0] != Fa_neig[1]);
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
	assert(Nodes != NULL);
	// transition from internal to external order of nodes
	static const int face3Order[3] = {0, 1, 2},
			face4Order[4] = {0, 1, 3, 2};
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	hObj & face(m.faces_[Fa]);
	Nodes[0] = face.typeSpecyfic_.nVerts_;
	const int * order(face.type_ == Face3::myType ? face3Order : face4Order);
	for (register int i(0); i < Nodes[0]; ++i) {
		Nodes[i + 1] = face.verts(order[i]);
	}

	//if(face.flags(F_TYPE) == F_FLIPPED) {
	//	std::swap(Nodes[2],Nodes[3]);
	//}
	if (Node_coor != NULL) {
		for (register int c(0); c < Nodes[0]; ++c) {
			Node_coor[3 * c  ] = m.vertices_[Nodes[c]].coords_[0];
			Node_coor[3 * c + 1] = m.vertices_[Nodes[c]].coords_[1];
			Node_coor[3 * c + 2] = m.vertices_[Nodes[c]].coords_[2];
		}
	}

#ifdef _DEBUG
	for (int i(1), j(1); i < Nodes[0] - 1; ++j) {
// ?? what's that??
		assert( i == j ? true : Nodes[i] != Nodes[j] );

		if (j == Nodes[0] - 1) {
			j = 1;
			++i;
		}
	}
#endif
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	assert(Fa >= 0);
	assert(Fa <  static_cast<int>(m.faces_.last()));
	assert(m.faces_[Fa].isBroken());

	if (Fasons != NULL) {
		Fasons[0] = m.faces_[Fa].nMyClassSons_;
	}

	if (m.faces_[Fa].type_ == Face4::myType) {
		if (Fasons != NULL) {
			for (int i(1); i <= Fasons[0]; ++i) {
				Fasons[i] = m.faces_[Fa].getMyClassChild<Face4>(i - 1)->pos_;
			}
		}
		if (Node_mid != NULL) {
			*Node_mid = m.faces_[Fa].sons(0) ;
		}
	}
	else
	{
		if (Fasons != NULL) {
			for (int i(1); i <= Fasons[0]; ++i) {
				Fasons[i] = m.faces_[Fa].getMyClassChild<Face3>(i - 1)->pos_;
			}
		}
		if (Node_mid != NULL) {
			*Node_mid = 0;
		}
	}
	return m.faces_[Fa].parent_;
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
		double * Acoeff,      /* in: coefficients of transformation between... */
		double * Bcoeff,	/* in: ...face coord and big neighb face coord */
		double * Xneig        /* out: local coordinates for neighbors */
) {
	/* auxiliary variables */
	int neig, ineig, el_type;
	double coor[2];

	/*++++++++++++++++ executable statements ++++++++++++++++*/



	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

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
			}

			// if face is pointing inside elem, swap x and y
			const int type = m.elements_[neig].flags(EL_TYPE);
			int mask = 1 << (Neig_sides[ineig]);
			if ((type & mask) != 0) {
				std::swap(coor[0], coor[1]);
			}
			switch (el_type) {
				case MMC_PRISM: {
					assert(Xloc[0] >= -1.0 - SMALL && Xloc[0] <= 1.0 + SMALL);
					assert(Xloc[1] >= -1.0 - SMALL && Xloc[1] <= 1.0 + SMALL);
					assert(Xloc[0] + Xloc[1] <= 2.0 + SMALL);
					assert(Xloc[0] + Xloc[1] >= -2.0 - SMALL);

					switch (Neig_sides[ineig]) {
						case 0: // 0,1 - triangles (bottom, top)
						case 1: {
							Xneig[3 * ineig  ] = coor[0];
							Xneig[3 * ineig + 1] = coor[1];
							Xneig[3 * ineig + 2] = (Neig_sides[ineig] == 0 ? -1 : 1);
						}
						break;
						case 2: { // 2,4,3 - rectangles (sides)
							Xneig[3 * ineig  ] = 0.5 * (1 + coor[0]);
							Xneig[3 * ineig + 1] = 0.0;
							Xneig[3 * ineig + 2] = coor[1];

						}
						break;
						case 3: {
							Xneig[3 * ineig  ] = 1 - 0.5 * (1 + coor[0]);
							Xneig[3 * ineig + 1] = 0.5 * (1 + coor[0]);
							Xneig[3 * ineig + 2] = coor[1];
						}
						break;
						case 4: {
							Xneig[3 * ineig  ] = 0.0;
							Xneig[3 * ineig + 1] = 0.5 * (1 + coor[0]);
							Xneig[3 * ineig + 2] = coor[1];
						}
						break;

					} //!switch(Neig_sides[ineig]
				}
				break;

				case MMC_TETRA: {
					assert(Xloc[0] >= -SMALL && Xloc[0] <= 1.0 + SMALL);
					assert(Xloc[1] >= -SMALL && Xloc[1] <= 1.0 + SMALL);
					assert(Xloc[0] + Xloc[1] <= 1.0);

					// all trans below are for faces pointing outside
					switch (Neig_sides[ineig]) {
						case 0: {
							Xneig[3 * ineig  ] = coor[1];
							Xneig[3 * ineig + 1] = coor[0];
							Xneig[3 * ineig + 2] = 0.0;
						}
						break;
						case 1: {
							Xneig[3 * ineig  ] = coor[0];
							Xneig[3 * ineig + 1] = 0.0;
							Xneig[3 * ineig + 2] = coor[1];
						}
						break;
						case 2: {
							Xneig[3 * ineig  ] = 0.0;
							Xneig[3 * ineig + 1] = coor[1];
							Xneig[3 * ineig + 2] = coor[0];
						}
						break;
						case 3: {
							coor[0] *= 0.75;
							coor[1] *= 0.75;
							Xneig[3 * ineig + 0] = 1.0 - coor[0] - coor[1];
							Xneig[3 * ineig + 1] = coor[0];
							Xneig[3 * ineig + 2] = coor[1];

							//mfp_check_debug(Xneig[3*ineig+0] + Xneig[3*ineig+1] + Xneig[3*ineig+2] -1 < abs(SMALL), "Checking plane coords failed.");

						}
						break;
					} //!switch(Neig_sides[ineig]



				}
				break;
			}//!switch(el_type)
		}//!if
	}//!for

	/*kbw
	#ifdef DEBUG_MMM
	printf("In elem coord for a point %lf, %lf on a face\n",
	coor[0],coor[1]);
	printf("Neigh1 %d, side %d\n", Fa_neig[0], Neig_sides[0]);
	printf("Coor: %lf, %lf, %lf\n", Xneig[0], Xneig[1], Xneig[2]);
	if(Fa_neig[1]!=0){
	printf("Neigh2 %d, side %d, node_shift %d\n",
	Fa_neig[1], Neig_sides[1], Node_shift);
	printf("Coor: %lf, %lf, %lf\n", Xneig[3], Xneig[4], Xneig[5]);
	}
	#endif
	kbw*/

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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	double  vecNormArr[3];
	m.faceNormal(m.faces_[Fa], vecNormArr, Area);

	if (Vec_norm != NULL) {
		memcpy(Vec_norm, vecNormArr, sizeof(double) * 3);
	}
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	if (Ed > static_cast<int>(m.edges_.size())) {
		return MMC_FREE;
	} else if (Ed <= static_cast<int>(m.edges_.last())) {
		if (m.edges_.at(Ed).isBroken())
			return MMC_INACTIVE;
	}
	return MMC_ACTIVE;
}

/*---------------------------------------------------------
mmr_edge_nodes - to return edge node's IDs
---------------------------------------------------------*/
int mmr_edge_nodes( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,		/* in: edge ID */
		int * Edge_nodes	/* out: IDs of edge nodes */
) {
	assert(Ed > 0);
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	Edge_nodes[0] = m.edges_[Ed].verts(0) ;
	Edge_nodes[1] = m.edges_[Ed].verts(1) ;
	assert(Edge_nodes[0] > 0);
	assert(Edge_nodes[1] > 0);
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
	assert(Ed > 0);

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	if (Edge_sons != NULL) {
		Edge_sons[0] = m.edges_.at(Ed).getMyClassChild<Edge>(0)->pos_ ;
		Edge_sons[1] = m.edges_.at(Ed).getMyClassChild<Edge>(1)->pos_ ;

		assert(Edge_sons[0] > 0);
		assert(Edge_sons[1] > 0);
	}

	if (Edge_mid_node != NULL) {
		*Edge_mid_node = m.edges_.at(Ed).sons(0);
	}

	return 1;
}

/*---------------------------------------------------------
mmr_edge_structure - to return edge structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_edge_structure(/* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,  	   /* in: edge ID */
		int * Edge_struct /* out: edge structure in the form of integer array */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	const hObj & edge(m.edges_[Ed]);
	Edge_struct[0] = edge.type_;
	Edge_struct[1] = edge.verts(0);
	Edge_struct[2] = edge.verts(1);
	return 0;
}

/*---------------------------------------------------------
mmr_set_edge_type - to set type for an edge
---------------------------------------------------------*/
int mmr_set_edge_type( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge_id,     /* in: edge ID */
		int Type         /* in: edge type (the number of attempted subdivisions !!!)*/
) {
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

//    return static_cast<int>(m.vertices_[Node].status());
	int stat = MMC_ACTIVE;

	if ( (Node < FIRST) || (Node > m.vertices_.last()) ) {
		stat = MMC_INACTIVE;
	}
	else if ( & m.vertices_[Node] == NULL ) {
		stat = MMC_FREE;
	}

	return stat;
}

/*---------------------------------------------------------
mmr_node_coor - to return node coordinates
---------------------------------------------------------*/
int mmr_node_coor( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Coor /* out: coordinates */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	Coor[0] = m.vertices_[Node].coords_[0];
	Coor[1] = m.vertices_[Node].coords_[1];
	Coor[2] = m.vertices_[Node].coords_[2];

	return 3;
}

/*---------------------------------------------------------
mmr_move_node - moves node by vector
---------------------------------------------------------*/
int mmr_move_node( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Dx  	/* in: displacement (vector) */
)
{

	/* select the proper mesh data structure */
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	m.vertices_[Node].coords_[0] += Dx[0];
	m.vertices_[Node].coords_[1] += Dx[1];
	m.vertices_[Node].coords_[2] += Dx[2];

	return (1);

}

/*---------------------------------------------------------
mmr_set_node_coor - set node coordinates
---------------------------------------------------------*/
int mmr_set_node_coor( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Xcoor  	/* in: node coordinates */
)
{

	/* select the proper mesh data structure */
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	m.vertices_[Node].coords_[0] = Xcoor[0];
	m.vertices_[Node].coords_[1] = Xcoor[1];
	m.vertices_[Node].coords_[2] = Xcoor[2];

	return (1);

}

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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
	exit(-1);
	return -1;
}


extern int mmr_create_node( /* returns: node ID of created node (<0 - error) */
		int  Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		double Xcoor,  /* in: coordinates of new node */
		double Ycoor,
		double Zcoor
) {
	mmv_out_stream << "Method not implemented!";
	exit(-1);
	return -1;
}


extern int mmr_clust_fa4_t( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face,      /* in: face ID */
		int * Face_sons /* in (optional): face sons' IDs */
) {
	mmv_out_stream << "Method not implemented!";
	exit(-1);
	return -1;
}


extern int mmr_clust_fa4_q( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face,      /* in: face ID */
		int * Face_sons /* in (optional): face sons' IDs */
) {
	mmv_out_stream << "Method not implemented!";
	exit(-1);
	return -1;
}


extern int mmr_clust_edge( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge_id   /* in: clustered edge ID */
) {

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
	mmv_out_stream << "Method not implemented!";
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
	mmv_out_stream << "Method not implemented!";
	exit(-1);
	return -1;
}


extern int mmr_divide_edge( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int  Edge,      /* in: edge to be divided */
		int * Edge_sons, /* out: two sons */
		int * Node_mid  /* out: node in the middle */
) {
	mmv_out_stream << "Method not implemented!";
	exit(-1);
	return -1;
}



/*------------------------------------------------------------
mmr_del_elem    - to free an element structure
------------------------------------------------------------*/
int mmr_del_elem( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	/* in: ID of the mesh to be used or 0 for the current mesh */
		int Elem_pos      /* in: element ID */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	Elem & el = m.elements_.at(Elem_pos);
	el.mark2Del(&m);

	for (int i = 0; i < el.typeSpecyfic_.nFaces_; ++i) {
		Face & face = m.faces_.getById(el.components(i));
		if (face.neighs(1) == el.id_) {
			face.neighs(1) = B_COND;
		}
		else if (face.neighs(0) == el.id_) { //
			face.neighs(0) = face.neighs(1);
			face.neighs(1) = B_COND;
		}

		if ( ((face.neighs(0) == B_COND) || (face.neighs(0) == MMC_SUB_BND))
				&& ((face.neighs(1) == B_COND) || (face.neighs(1) == MMC_SUB_BND)) ) {
			face.mark2Del(&m);
		}
	}

	return 1;
}

int mmr_del_face( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face      /* in: face ID */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	m.faces_.at(Face).mark2Del(&m);
	return 1;
}

int mmr_del_edge( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge      /* in: edge ID */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	m.edges_.at(Edge).mark2Del(&m);
	return 1;

}

int mmr_del_node( /* returns: >=0 - success code, <0 - error code */
		int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node      /* in: node ID */
) {
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	m.vertices_.at(Node).mark2Del(&m);
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

	if (Edge_elems != NULL) {
		hHybridMeshWithContacts & m = select_mesh(Mesh_id);
		Edge_elems[0] = m.edge_elems[Ed].size();
		for (int i = 0; i < Edge_elems[0]; i++) {
			Edge_elems[i + 1] = m.edge_elems[Ed][i];
		}
	}

	return (0);
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	const Elem & el(m.elements_[El]);
	if (Edges != NULL) {
		Edges[0] = el.typeSpecyfic_.nEdges_;
		switch (el.type_) {
			case ElemT4::myType:
				Edges[1] = m.edge(el.verts(0), el.verts(1)).pos_;
				Edges[2] = m.edge(el.verts(0), el.verts(2)).pos_;
				Edges[3] = m.edge(el.verts(0), el.verts(3)).pos_;
				Edges[4] = m.edge(el.verts(1), el.verts(2)).pos_;
				Edges[5] = m.edge(el.verts(1), el.verts(3)).pos_;
				Edges[6] = m.edge(el.verts(2), el.verts(3)).pos_;
				break;
			case ElemPrism::myType:
				Edges[1] = m.edge(el.verts(0), el.verts(1)).pos_;
				Edges[2] = m.edge(el.verts(1), el.verts(2)).pos_;
				Edges[3] = m.edge(el.verts(2), el.verts(0)).pos_;
				Edges[4] = m.edge(el.verts(3), el.verts(4)).pos_;
				Edges[5] = m.edge(el.verts(4), el.verts(5)).pos_;
				Edges[6] = m.edge(el.verts(5), el.verts(3)).pos_;
				Edges[7] = m.edge(el.verts(0), el.verts(3)).pos_;
				Edges[8] = m.edge(el.verts(1), el.verts(4)).pos_;
				Edges[9] = m.edge(el.verts(2), el.verts(5)).pos_;
				break;
		}
	}
	return el.typeSpecyfic_.nEdges_;
}

/*---------------------------------------------------------
  mmr_el_fa_nodes - to get list local face nodes indexes in elem
---------------------------------------------------------*/
extern int mmr_el_fa_nodes( // returns face type flag
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,		// in: global elem ID
		int Fa,		/* in: local face number in elem El */
		int * fa_nodes		/* out: list of local indexes of face nodes */
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	int type = MMC_QUAD;
	if (m.elements_[El].type_ == MMC_PRISM) {
		switch (Fa) {
			case 0: {
				fa_nodes[0] = 0;
				fa_nodes[1] = 1;
				fa_nodes[2] = 2;
				type = MMC_TRIA;
			}
			break;
			case 1: {
				fa_nodes[0] = 3;
				fa_nodes[1] = 4;
				fa_nodes[2] = 5;
				type = MMC_TRIA;
			}
			break;
			case 2: {
				fa_nodes[0] = 0;
				fa_nodes[1] = 1;
				fa_nodes[2] = 4;
				fa_nodes[3] = 3;
			}
			break;
			case 3: {
				fa_nodes[0] = 1;
				fa_nodes[1] = 2;
				fa_nodes[2] = 5;
				fa_nodes[3] = 4;
			}
			break;
			case 4: {
				fa_nodes[0] = 0;
				fa_nodes[1] = 2;
				fa_nodes[2] = 5;
				fa_nodes[3] = 3;
			}
			break;
		} //!switch(neig..)
	}
	else { // MMC_TETRA
		type = MMC_TRIA;
		switch (Fa) {
			case 0: {
				fa_nodes[0] = 0;
				fa_nodes[1] = 1;
				fa_nodes[2] = 2;
			}
			break;
			case 1: {
				fa_nodes[0] = 0;
				fa_nodes[1] = 1;
				fa_nodes[2] = 3;
			}
			break;
			case 2: {
				fa_nodes[0] = 0;
				fa_nodes[1] = 2;
				fa_nodes[2] = 3;
			}
			break;
			case 3: {
				fa_nodes[0] = 1;
				fa_nodes[1] = 2;
				fa_nodes[2] = 3;
			}
			break;
		}//!switch(neig..)
		int flag = m.elements_[El].flags(EL_TYPE);
		if ((flag & 1 << Fa) != 0) { // if face is pointing interior of elem, make it pointing outside
			std::swap(fa_nodes[1], fa_nodes[2]);
		}
	}//!if(mmr_el_type)

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

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	m.findEdgeElems();

	return (1);
}


/*---------------------------------------------------------
  mmr_get_max_gen - to set maximal allowed generation level for elements
------------------------------------------------------------*/
extern int mmr_get_max_gen(/* returns: >=0-success code, <0-error code */
		int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);

	return m.maxGen_;
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.maxGenDiff_;
}
/*
mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready for dofs projection
*/

extern int mmr_is_ready_for_proj_dof_ref(int Mesh_id)
{
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return m.finalRef() ? 1 : 0;
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
	assert(Coord != NULL);
	assert(El_id >= static_cast<int>(FIRST));

	double dist( -1.0 );


#ifdef TURBULENTFLOW

	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	const hHybridMesh & m(*current_mesh);
	const hObj & el(current_mesh->elements_[El_id]);

	// let us assume that near boundary are always prism elements
	if (el.type_ == ElemPrism::myType) {
		// if all vertices have the same "nearst boundary face"
		ID nearstFaces[ElemPrism::nVerts] = {0};
		for (uTind i(0); i < el.typeSpecyfic_.nVerts_; ++i) {
			nearstFaces[i] = m.vertices_[el.verts(i)].nearstFace_;
		}

		std::sort(nearstFaces, nearstFaces + ElemPrism::nVerts);
		const int nUniqueFaces(std::unique(nearstFaces, nearstFaces + ElemPrism::nVerts) - nearstFaces);

		assert(nUniqueFaces >= 1);

		dist = mmr_point_plane_dist(Coord, m.faces_(nearstFaces[0]).planeCoords_);
		for (int i(1); i < nUniqueFaces; ++i) {
			double newDist = mmr_point_plane_dist(Coord, m.faces_(nearstFaces[i]).planeCoords_);
			if (newDist < dist) {
				dist = newDist;
			}
		}
	}
#endif
	return dist;
}


/*---------------------------------------------------------
  mmr_init_dist2boundary  - to return distance from given vertex to nearst boundary
---------------------------------------------------------*/
int mmr_init_dist2bound( // returns 0 if succesfull
		const int Mesh_id, //in: mesh ID or 0 (MMC_CUR_MESH_ID)
		const int * BCs, //in: c-array of accepted BC numbers
		const int nBCs) //in: length of BCs parameter
{
#ifdef TURBULENTFLOW
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	current_mesh->computeDist2Bound(BCs, nBCs);
#endif
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
	hHybridMeshWithContacts & m = select_mesh(Mesh_id);
	return hashlittle(&m.vertices_.first(), m.vertices_.size(), Mesh_id);

}

/** @}*/

#ifdef __cplusplus
}
#endif

