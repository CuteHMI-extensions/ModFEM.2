
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cmath>
//#include <GL/gl.h>
#include <assert.h>
#include <vector>
using namespace std;

#include <modfem/uth_mesh.h>
#include <modfem/uth_log.h>
#include <modfem/uth_mat.h>
#include <modfem/mmh_intf.h>
#include <modfem/uth_bc.h>

//namespace MESH_UTILS {

#include <modfem/mm_remesh/ticooMesh3D.h>
//}


TicooMesh3D & utr_mesh_read_ticoo_from_intf(const int Mesh_id)
{
	char mm_name[255] = {0};
	mmr_module_introduce(mm_name);
	if ( 0 == strcmp("3D_remesh", mm_name) ) {
		// so utl_mesh is compiled with mml_remesh,
		// no need of mesh conversion, can use mesh from mml_remesh
		return *reinterpret_cast<TicooMesh3D *>( mmr_module_internals(Mesh_id) );

	}

//	int inode = 0;
//	while (0 != (mmr_get_next_node_all(Mesh_id, inode)) {

//    }

//	int iface = 0;
//	while (0 != (mmr_get_next_face_all(Mesh_id, iface)) {

//    }
	mf_log_err("Mesh conversion to remesh(PC) not implemented!");
	static TicooMesh3D work_mesh_ptr(50);

	//mf_log_info("Ticoo elemes: %d", work_mesh_ptr->getNumberElements() );

	return work_mesh_ptr;
}

int utr_mesh_write_ticoo_to_intf(const int Mesh_id, const  TicooMesh3D & t_mesh)
{
	char mm_name[255] = {0};
	mmr_module_introduce(mm_name);
	if ( 0 == strcmp("3D_remesh", mm_name) ) {
		// so utl_mesh is compiled with mml_remesh,
		// no need of mesh conversion, can use mesh from mml_remesh
		return 0;
	}

	/// MAGIC HERE
	///
	mf_log_err("Mesh conversion from remesh(PC) not implemented!");
	return Mesh_id; // not working yet
}

#ifdef __cplusplus
extern "C" {
#endif

/**
/** \defgroup UTM_MESH Mesh Utilities
/** \ingroup UTM
/** @{
/** */


///
/// \brief utr_mesh_insert_BC_contact To insert contact faces into the given mesh.
///
/// \param Mesh_id mesh in: to update
/// \param N_contacts in: number of contacts
/// \param BC_nums[N_contacts] in: array of BC numbers
/// \param groupIDs[2*N_contacts] in: ID between which contact BC will be implemeted
/// \param Types[N_contacts] in: type specyfing IDs type
/// \return number of inserted BC contacts
///




int utr_mesh_insert_BC_contact(
		const char * Workdir,
		const int Mesh_id,
		const int N_contacts,
		const int * BC_nums,
		const int * IDs,
		const utt_mesh_bc_type * Types)
{
	int n_inserted_contacts = 0;

	if (N_contacts > 0) {
		mf_check_mem(BC_nums);
		mf_check_mem(IDs);
		mf_check_mem(Types);

		int n_groups = mmr_groups_number(Mesh_id);

		mf_log_info("Found %d groups in mesh.", n_groups);

		std::vector<int>    groups(n_groups), blocks, mats, bcs, groups_fluid;
		std::vector<double> tempB;

		mmr_groups_ids(Mesh_id, groups.data(), n_groups);
		//mf_log_info("Groups: %d, %d",groups[0],groups[1]);

		for (int i = 0, ile = utr_temp2block_size(); i < ile; ++i) {

			tempB.push_back(utr_temp2block_id(i));
		}

		for (int i = 0; i < n_groups; ++i) {
			blocks.push_back(groups[i]);
			blocks.push_back(utr_bc_get_blockID(groups[i]));

			mats.push_back(groups[i]);
			mats.push_back(utr_mat_get_matID(groups[i]));
		}


		for (int i = 0; i < N_contacts; ++i) {
			bcs.push_back(Types[i]);
			bcs.push_back(IDs[2 * i]);
			bcs.push_back(IDs[2 * i + 1]);
			bcs.push_back(BC_nums[i]);
		}

		int * fluidT, iterFluid = 0, iterNOTFluid = 0;

		for (int i = 0; i < mats.size(); i += 2) {

			if (utr_mat_get_material( mats[i] )->is_fluid == FLUID) {
				groups_fluid.push_back(mats[i]);
			}

			//printf("bleeeeeeeeeeeeee %lf",mat1->atT_dynamic_viscosity[0]);

		}

		mf_log_info("Inserting new BCs.");



		TicooMesh3D & t_mesh = utr_mesh_read_ticoo_from_intf(Mesh_id);

		if (t_mesh.getNumberElements() > 0) {

			t_mesh.mmr_split_into_blocks_add_contact(Workdir, mats.data(), mats.size(), blocks.data(), blocks.size(), bcs.data(), bcs.size(), tempB.data(), tempB.size(), groups_fluid.data(), groups_fluid.size());

			utr_mesh_write_ticoo_to_intf(Mesh_id, t_mesh);
		}

//mf_log_info("Groups: %f, %f, %f, %f",tempB[0],tempB[1],tempB[2],tempB[3]);



		//mf_log_info("Groups: %f, %f",tempB[2],tempB[3]);
		//mf_log_info("Groups: %d, %d",blocks[2],blocks[3]);
		//TicooMesh3D mes3d_(50);

//		// wczytanie?


		//double tB[4] ={1.0,1102.0,2.0,295.0};
		//int lB=4;

		/*
		mes3d_.mmr_split_into_blocks_add_contact(
			mats.data(), mats.size(),
			blocks.data(), blocks.size(),
			bcs.data(), bcs.size(),tempB.data(),tempB.size());
		*/

		// zapis

		n_inserted_contacts = N_contacts;
	}
	return n_inserted_contacts;
}



int utr_mesh_insert_new_bcnums(const char * Workdir, const int Mesh_id)
{
	const int count = utr_bc_to_insert_n_assigments();

	/// Initialziation of new bc data
	if (count > 0) {
		const utt_bc_assignment * BCs2insert = utr_bc_to_insert_get_assigments();

		std::vector<int> bc_nums(count, -1);
		std::vector<int> ids(2 * count, -1);
		std::vector<utt_mesh_bc_type> types(count, UTE_MESH_BC_UNDEFINED);

		for (int i = 0; i < count; ++i) {
			bc_nums[i] = BCs2insert[i].BC_num;
			types[i] = BCs2insert[i].Id_type;
			ids[2 * i] = BCs2insert[i].Ids[0];
			ids[2 * i + 1] = BCs2insert[i].Ids[1];
		}

		int n_inserted = utr_mesh_insert_BC_contact(Workdir,
						Mesh_id,
						count,
						bc_nums.data(),
						ids.data(),
						types.data());

		mf_check(n_inserted == count,
				"Not all contact BCs ware inserted correctly (%d from %d)",
				n_inserted, count);

		utr_mesh_compute_vec_norm_at_nodes(Mesh_id);
	}

	return count;	// Function must return value, otherwise UB occurs.
}

/// Normal vector at node computations.
///
/// \brief utv_mesh_vec_norm_at_nodes - holds normal vector at node
///
std::vector< std::vector<double>> utv_mesh_vec_norm_at_nodes;

void utr_mesh_compute_vec_norm_at_nodes(const int Mesh_id)
{
	if (Mesh_id >= utv_mesh_vec_norm_at_nodes.size()) {
		utv_mesh_vec_norm_at_nodes.resize(Mesh_id + 1);
	}

	std::vector<double> & vec_norm_at_nodes = utv_mesh_vec_norm_at_nodes[Mesh_id];
	const int n_nodes = mmr_get_max_node_id(Mesh_id);
	vec_norm_at_nodes.clear();
	vec_norm_at_nodes.resize(3 * (n_nodes + 1) );

	TicooMesh3D & t_mesh = utr_mesh_read_ticoo_from_intf(Mesh_id);

	if (! t_mesh.empty()) {

		//t_mesh.wyszukajNormalnePowierzchniDlaPunktowGranicznych();

		int node = 0;
		double * v_begin = vec_norm_at_nodes.data();
		while ( node = mmr_get_next_node_all(Mesh_id, node)) {
			t_mesh.mmr_get_vec_norm_at_node( node - 1, v_begin + (3 * node) );
			//        mf_log_info("Vec norm at node %d (%lf,%lf,%lf)", node,
			//                    v_begin[3*node],
			//                    v_begin[3*node+1],
			//                    v_begin[3*node+2]);


		}
		utr_mesh_write_ticoo_to_intf(Mesh_id, t_mesh);

	}
}

void utr_mesh_get_vec_norm_at_node(const int Mesh_id, int Node_id, double * Vec_norm)
{
	memcpy(Vec_norm, & utv_mesh_vec_norm_at_nodes[Mesh_id][3 * Node_id], 3 * sizeof(double) );
}

int utr_mesh_smoothing(const int Mesh_id, double l_step, int repetition)
{
	TicooMesh3D & t_mesh = utr_mesh_read_ticoo_from_intf(Mesh_id);

	if (! t_mesh.empty()) {

		mf_log_info("Smoothing mesh with parameters %lf, %d", l_step, repetition);
		t_mesh.utr_smoothing(l_step, repetition);

		utr_mesh_write_ticoo_to_intf(Mesh_id, t_mesh);

	}

	return 0;
}

int utr_mesh_smoothing_menu(const int Mesh_id, FILE * Interactive_input, FILE * Interactive_output)
{
	if (Interactive_input == stdin) {
		fprintf(Interactive_output, "\ts - mesh quality statistics\n");
		fprintf(Interactive_output, "\tc - coarse mesh smoothing\n");
		fprintf(Interactive_output, "\ta - average mesh smoothing (more time consuming)\n");
		fprintf(Interactive_output, "\tf - fine mesh smoothing (even more time consuming)\n");
		fprintf(Interactive_output, "\tq - go back to prev. menu\n");
		fprintf(Interactive_output, "\te - utr_count_bc_surface\n");
	}


	int rc = 0;
	char c = 0;
	do {
		fscanf(Interactive_input, "%c", &c);
		switch (c) {
			case 0:
			case '\n': //ignore Enter
				break;
			case 'c':
				utr_mesh_smoothing(Mesh_id, 0.1, 10);
				utr_mesh_quality_statistics(Mesh_id, Interactive_output);
				rc = 1;
				break;
			case 'a':
				utr_mesh_smoothing(Mesh_id, 0.05, 15);
				utr_mesh_quality_statistics(Mesh_id, Interactive_output);
				rc = 1;
				break;
			case 'f':
				utr_mesh_smoothing(Mesh_id, 0.1, 10);
				utr_mesh_smoothing(Mesh_id, 0.01, 10);
				utr_mesh_quality_statistics(Mesh_id, Interactive_output);
				rc = 1;
				break;
			case 's':
				utr_mesh_quality_statistics(Mesh_id, Interactive_output);
				rc = 1;
				break;
			case 'q':
				rc = 1;
				break;
			case 'e':
				utr_count_bc_surface(Mesh_id, Interactive_output);
				rc = 1;
				break;
			default:
				rc = -1;
				mf_log_err("Unknown command %c", c);
				break;
		}

	} while (rc == 0);

	return rc;
}


int utr_count_bc_surface(const int Mesh_id, FILE * Interactive_output)
{
	TicooMesh3D & t_mesh = utr_mesh_read_ticoo_from_intf(Mesh_id);

	if (! t_mesh.empty()) {

		mf_log_info("Mesh count_bc_surface(%d)", Mesh_id);
		t_mesh.count_bc_surface(Interactive_output);

		utr_mesh_write_ticoo_to_intf(Mesh_id, t_mesh);
	}
	return 0;
}

int utr_mesh_quality_statistics(const int Mesh_id, FILE * Interactive_output)
{
	TicooMesh3D & t_mesh = utr_mesh_read_ticoo_from_intf(Mesh_id);

	if (! t_mesh.empty()) {

		mf_log_info("Mesh quality statistics for mesh_id=%d", Mesh_id);
		t_mesh.mesh_quality_stats(Interactive_output);

		utr_mesh_write_ticoo_to_intf(Mesh_id, t_mesh);

	}

	return 0;
}

/** @} */ // end of group


#ifdef __cplusplus
}
#endif
