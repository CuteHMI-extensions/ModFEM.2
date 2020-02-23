/************************************************************************
File uts_write_paraview.c - paraview output generation

Contains definition of routines:
  utr_write_paraview_mesh - to dump mesh in Paraview format
  utr_write_paraview_partmesh - to dump mesh in Paraview format includiong partition info
  utr_write_paraview_field - to dump field in Paraview format

------------------------------
History:
	08.2012 - Kamil Wachala, kamil.wachala@gmail.com
	05.2011 - Kazimierz Michalik, kamich@agh.edu.pl, initial version
	06.2011 - KGB, pobanas@cyf-kr.edu.pl
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
	05.2016 - Kazimierz Michalik, kamich@agh.edu.pl, XML-VTK version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string>
#include<math.h>
#include<fstream>
#include<map>
#include<locale>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <zlib.h>

#include <modfem/pdh_intf.h>
#include <modfem/mmh_intf.h>
#include <modfem/aph_intf.h>
#include <modfem/uth_intf.h>
#include <modfem/uth_mat.h>
#include <modfem/uth_bc.h>
#include <modfem/uth_log.h>
#include "base64.h"

namespace pt = boost::property_tree;

typedef struct {
	double dofs[PDC_MAXEQ];
} PVSolInfo;


static const char * utv_paraview_flags_strings_legacy[UTE_LAST] = {
	"POINT_DATA",
	"CELL_DATA",
	"SCALARS",
	"VECTORS",
	"TENSORS",
	"double",
	"int",
	"# vtk DataFile Version 3.1\n",
	"<?xml version=\"1.0\"?>"
};

static const char * utv_paraview_flags_strings_XML[UTE_LAST] = {
	"PointData",
	"CellData",
	"Scalars",
	"Vectors",
	"Tensors",
	"Float64",
	"Int32",
	"# vtk DataFile Version 3.1\n",
	"<?xml version=\"1.0\"?>"
};


struct utt_paraview_vtk_xml {
	std::string save_dir;
	std::string main_xml_filename;
	std::string byte_order; // this two should be changed together
	std::string format;     // this two should be changed together
	std::string compressor;
	std::string partmesh_file;
	int mesh_id, n_elems, n_nodes, mesh_hash;
	std::map<double, std::string>    time2file;
	int  version;

	utt_paraview_vtk_xml()
		: mesh_id(-1), n_elems(-1),
		  byte_order("LittleEndian"),
		  save_dir("paraview"),
		  format("binary"), // ascii or binary // this two should be changed together
		  compressor("vtkZLibDataCompressor") //"" or vtkZLibDataCompressor // this two should be changed together
	{}

} ;

utt_paraview_vtk_xml utv_paraview_vtk_xml_data;

std::string  &  utr_compress_and_encode64(char const * data, size_t size, std::string & ret, std::string & header)
{
	uint32_t header_uints[4];
	header_uints[0] = 1; // no. of blocks
	header_uints[1] = size; // length of data/block
	header_uints[2] = size; // length of data

	uLongf compress_size = compressBound(size);

	std::vector<Bytef> comp_data;
	comp_data.resize(compress_size);
	compress(comp_data.data(), &compress_size,
			reinterpret_cast<const Bytef *>(data), size);
	comp_data.resize(compress_size);    //shrink to fit actual compressed size

	header_uints[3] = compress_size; // length of compressed data

//    mf_log_info("Compression header: [%d,%d,%d,%d] -> [%s]",
//                header_uints[0],header_uints[1],header_uints[2],header_uints[3],
//                reinterpret_cast<char*>(header_uints) );
//    mf_log_info("Compression ratio: %.2f\%",100. - 100.*(float)compress_size/(float)size);

	/**need to encode here **/
	base64_encode(reinterpret_cast<const unsigned char *>(header_uints), 4 * sizeof(uint32_t), header);

	/**need to encode here **/
	base64_encode(comp_data.data(),
			compress_size,
			ret);

	return ret;
}

template<typename T>
void    utr_paraview_do_write(std::ofstream & file, T const * data, const int size)
{
	if (utv_paraview_vtk_xml_data.format == "binary"
			&& utv_paraview_vtk_xml_data.compressor == "vtkZLibDataCompressor")
	{
		std::string compressed_n_encoded_data, header;
		utr_compress_and_encode64(reinterpret_cast<const char *>(data), size * sizeof(T), compressed_n_encoded_data, header);
		file << header << compressed_n_encoded_data;
	}
	else if (utv_paraview_vtk_xml_data.format == "ascii") {
		for (int i = 0; i < size; ++i) {
			file << data[i] << " ";
		}
	}
	else {
		mf_fatal_err("Ill-formed parameters of paraview writer!");
	}
}

template<>
void    utr_paraview_do_write<uint8_t>(std::ofstream & file, uint8_t const * data, const int size)
{
	if (utv_paraview_vtk_xml_data.format == "binary"
			&& utv_paraview_vtk_xml_data.compressor == "vtkZLibDataCompressor")
	{
		std::string compressed_n_encoded_data, header;
		utr_compress_and_encode64(reinterpret_cast<const char *>(data), size * sizeof(uint8_t), compressed_n_encoded_data, header);
		file << header << compressed_n_encoded_data;
	}
	else if (utv_paraview_vtk_xml_data.format == "ascii") {
		for (int i = 0; i < size; ++i) {
			file << (int)data[i] << " ";
		}
	}
	else {
		mf_fatal_err("Ill-formed parameters of paraview writer!");
	}
}

std::string & utr_vtk_xml_get_timestep_filename(utt_paraview_vtk_xml * Data,
		const double Time,
		const int Group,
		const int Part,
		const char * Workdir,
		const char * Filename
)
{


	mf_check( Filename != NULL, "Paraview VTK filename not set!");

	static bool was_init = false;
	// First-time initialization.
	if (was_init == false) {
		was_init = true;
		// reverse engeneer name
		if (Workdir != NULL
				&& strlen(Workdir) > 0) {
			Data->main_xml_filename = Workdir;
			Data->main_xml_filename.append("/").append(Filename);
		}
		else {
			Data->main_xml_filename = Filename;
		}

		Data->main_xml_filename.resize(Data->main_xml_filename.size() - 2);
		Data->main_xml_filename.append(".pvd");

		// write xml premable and usual stuff
		std::ofstream xml_file(Data->main_xml_filename.c_str());
		xml_file << "<?xml version=\"1.0\"?>\n"
				"<VTKFile type=\"Collection\" version=\"0.1\"\n"
				"byte_order=\"" << Data->byte_order << "\"\n"
				" compressor=\"" << Data->compressor << "\"\n"
				">"
				"<Collection>\n"
				"</Collection>\n"
				"</VTKFile>\n";

		xml_file.close();
	}

	mf_check( !Data->main_xml_filename.empty(), "Paraview VTK filename not set!");

	if (Data->time2file.find(Time) == Data->time2file.end()) {
		// Create appropriate XML node.
		std::string file_name(Filename);
		file_name.append(".vtu");
		std::string time_full_file_name;
		if (Workdir != NULL
				&& strlen(Workdir) > 0) {
			time_full_file_name = Workdir;
			time_full_file_name.append("/").append(file_name);
		}
		else {
			time_full_file_name = file_name;
		}


		Data->time2file[Time] = time_full_file_name;

		// Updateing file.
		// Create empty property tree object
		// Locale - to correct decimal point.
		std::locale::global( std::locale( "" ) );
		pt::ptree tree;
		// Parse the XML into the property tree.
		pt::read_xml(Data->main_xml_filename, tree);
		pt::ptree & new_data_set = tree.add("VTKFile.Collection.DataSet", "");
		new_data_set.add("<xmlattr>.timestep", Time);
		new_data_set.add("<xmlattr>.group", Group);
		new_data_set.add("<xmlattr>.part", Part);
		new_data_set.add("<xmlattr>.file", file_name);
		pt::write_xml(Data->main_xml_filename, tree);

		std::locale::global( std::locale::classic() );
	}

	return  Data->time2file[Time];
}

void utr_write_paraview_BC_XML_VTK(int Mesh_id, const char * Workdir, const char * Filename, int * Part, ute_paraview_flags VTK_file_version)
{
	std::string bcmesh_filename = Workdir;
	bcmesh_filename.append("/").append(Filename).append("_BC.vtu");

	const int max_node_id = mmr_get_nr_node(Mesh_id);
	const int n_faces = mmr_get_nr_face(Mesh_id);

	std::vector<int32_t> offsets;
	offsets.reserve(n_faces);
	std::vector<uint8_t> fa_types;
	fa_types.reserve(n_faces);
	std::vector<int32_t>	bc_nums;
	bc_nums.reserve(n_faces / 10);


	std::vector<int> nodes_data;
	nodes_data.reserve(MMC_MAXFAVNO * n_faces);
	int fa_nodes[MMC_MAXFAVNO + 1] = { 0 };
	int offset = 0, bc = 0;
	for (int fa_id = 0; (fa_id = mmr_get_next_act_face(Mesh_id, fa_id)) != 0;) {
		if ((bc = mmr_fa_bc(Mesh_id, fa_id)) != 0) {
			bc_nums.push_back(bc);
			mmr_fa_node_coor(Mesh_id, fa_id, fa_nodes, NULL);
			offsets.push_back(offset += fa_nodes[0]);
			for (int i = 1; i <= fa_nodes[0]; ++i) {
				nodes_data.push_back(fa_nodes[i] - 1);
			}
			switch (mmr_fa_type(Mesh_id, fa_id))
			{
				case MMC_TRIA:
					fa_types.push_back(5);
					break;
				case MMC_QUAD:
					fa_types.push_back(9);
					break;
			}
		}
	}

	std::ofstream bcmesh_file(bcmesh_filename.c_str());
	bcmesh_file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\""
			<< utv_paraview_vtk_xml_data.byte_order
			<< "\" compressor=\""
			<< utv_paraview_vtk_xml_data.compressor << "\" >\n"
			<< "<UnstructuredGrid>\n";

	//////////////////////////////////////// SURFACE MESH - BC ///////////////////////////////////////////////
	bcmesh_file << "<Piece NumberOfPoints=\"" << max_node_id << "\" NumberOfCells=\"" << bc_nums.size() << "\">\n"
			<< "<Points>\n"
			<< "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";
	// writing out points
	double * node_coor = new double[3 * max_node_id];
	for (int ino = 1; ino <= max_node_id; ino++) {
		if (mmr_node_status(Mesh_id, ino) == MMC_ACTIVE) {
			mmr_node_coor(Mesh_id, ino, node_coor + (ino - 1) * 3);
			//partmesh_file << node_coor[0] << " " << node_coor[1] << " " << node_coor[2] << " ";
		}
		else {
			node_coor[(ino - 1) * 3] = 0.0;
			node_coor[(ino - 1) * 3 + 1] = 0.0;
			node_coor[(ino - 1) * 3 + 2] = 0.0;
			//partmesh_file << "0.0 0.0 0.0";
		}
	}
	utr_paraview_do_write(bcmesh_file, node_coor, 3 * max_node_id);
	delete[] node_coor;
	bcmesh_file << "</DataArray>\n"
			<< "</Points>\n";

	//////////////////////////////////////// CELLS - SURFACE TRIANGLES ///////////////////////////////////////////////

	bcmesh_file << "<Cells>\n";
	/// Element nodes.
	bcmesh_file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << utv_paraview_vtk_xml_data.format << "\" >\n";

	utr_paraview_do_write(bcmesh_file, nodes_data.data(), nodes_data.size());
	nodes_data.clear();
	bcmesh_file << "</DataArray>\n";

	/// Element offsets of nodes.
	bcmesh_file << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";
	utr_paraview_do_write(bcmesh_file, offsets.data(), offsets.size());
	offsets.clear();
	bcmesh_file << "</DataArray>\n";

	/// Element types
	bcmesh_file << "<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";
	utr_paraview_do_write(bcmesh_file, fa_types.data(), fa_types.size());
	fa_types.clear();
	bcmesh_file << "</DataArray>\n"
			<< "</Cells>\n";


	bcmesh_file << "<CellData>\n"
			<< "<DataArray Name=\"BC_num\" type=\"Int32\" NumberOfComponents=\"1\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";

	utr_paraview_do_write(bcmesh_file, bc_nums.data(), bc_nums.size());

	bcmesh_file << "</DataArray>\n";
	bcmesh_file << "</CellData>\n";
	bcmesh_file << "</Piece>\n";
	bcmesh_file << "</UnstructuredGrid>\n"
			<< "</VTKFile>\n";
	bcmesh_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void utr_write_paraview_partmesh_XML_VTK(int Mesh_id, const char * Workdir, const char * Filename, int * Part, ute_paraview_flags VTK_file_version)
{

	const int max_node_id = mmr_get_max_node_id(Mesh_id);
	const int n_elems = mmr_get_nr_elem(Mesh_id);
	const int mesh_hash = mmr_current_mesh_state_hash(Mesh_id);

	if ( (utv_paraview_vtk_xml_data.mesh_id == Mesh_id)
			&&  (utv_paraview_vtk_xml_data.n_elems == n_elems)
			&&  (utv_paraview_vtk_xml_data.n_nodes == max_node_id)
			&&  (utv_paraview_vtk_xml_data.mesh_hash == mesh_hash) ) {
		mf_log_info("Reusing prev. partmesh file (%s), because it looks that mesh is unchanged.",
				utv_paraview_vtk_xml_data.partmesh_file.c_str());
		return;
	}

	utr_write_paraview_BC_XML_VTK(Mesh_id, Workdir, Filename, Part, VTK_file_version);

	std::string partmesh_filename = Workdir;
	partmesh_filename.append("/").append(Filename).append(".partmesh");
	utv_paraview_vtk_xml_data.mesh_id = Mesh_id;
	utv_paraview_vtk_xml_data.n_elems = n_elems;
	utv_paraview_vtk_xml_data.n_nodes = max_node_id;
	utv_paraview_vtk_xml_data.mesh_hash = mesh_hash;
	// Writing mesh as "UnstructuredGrid (. vtu) —Serial vtkUnstructuredGrid (unstructured)."


	std::ofstream partmesh_file(partmesh_filename.c_str());
	partmesh_file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\""
			<< utv_paraview_vtk_xml_data.byte_order
			<< "\" compressor=\""
			<< utv_paraview_vtk_xml_data.compressor << "\" >\n"
			<< "<UnstructuredGrid>\n";

	//////////////////////////////////////// VOLUME MESH  ///////////////////////////////////////////////
	partmesh_file << "<Piece NumberOfPoints=\"" << max_node_id << "\" NumberOfCells=\"" << n_elems << "\">\n"

			//////////////////////////////////////// POINTS ///////////////////////////////////////////////
			<< "<Points>\n"
			<< "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";
	// writing out points
	double * node_coor = new double[3 * max_node_id];
	for (int ino = 1; ino <= max_node_id; ino++) {
		if (mmr_node_status(Mesh_id, ino) == MMC_ACTIVE) {
			mmr_node_coor(Mesh_id, ino, node_coor + (ino - 1) * 3);
			//partmesh_file << node_coor[0] << " " << node_coor[1] << " " << node_coor[2] << " ";
		}
		else {
			node_coor[(ino - 1) * 3] = 0.0;
			node_coor[(ino - 1) * 3 + 1] = 0.0;
			node_coor[(ino - 1) * 3 + 2] = 0.0;
			//partmesh_file << "0.0 0.0 0.0";
		}
	}
	utr_paraview_do_write(partmesh_file, node_coor, 3 * max_node_id);
	delete [] node_coor;
	partmesh_file << "</DataArray>\n"
			<< "</Points>\n";

	//////////////////////////////////////// CELLS ///////////////////////////////////////////////
	partmesh_file << "<Cells>\n";
	std::vector<int32_t> offsets;
	offsets.reserve(n_elems);
	std::vector<uint8_t> el_types;
	el_types.reserve(n_elems);

	/// Element nodes.
	partmesh_file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << utv_paraview_vtk_xml_data.format << "\" >\n";
	std::vector<int> nodes_data;
	nodes_data.reserve(MMC_MAXELVNO * mmr_get_nr_elem(Mesh_id));
	int el_nodes[MMC_MAXELVNO + 1] = {0};
	int offset = 0;
	for (int el_id = 0; (el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0;) {
		mmr_el_node_coor(Mesh_id, el_id, el_nodes, NULL);
		offsets.push_back(offset += el_nodes[0]);
		for (int i = 1; i <= el_nodes[0]; ++i) {
			nodes_data.push_back(el_nodes[i] - 1);
		}
		switch (mmr_el_type(Mesh_id, el_id))
		{
			case MMC_PRISM:
				el_types.push_back(13);
				break;
			case MMC_TETRA:
				el_types.push_back(10);
				break;
		}
	}
	utr_paraview_do_write(partmesh_file, nodes_data.data(), nodes_data.size());
	nodes_data.clear();
	partmesh_file << "</DataArray>\n";

	/// Element offsets of nodes.
	partmesh_file << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";
	utr_paraview_do_write(partmesh_file, offsets.data(), offsets.size());
	offsets.clear();
	partmesh_file << "</DataArray>\n";

	/// Element types
	partmesh_file << "<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";
	utr_paraview_do_write(partmesh_file, el_types.data(), el_types.size());
	el_types.clear();
	partmesh_file << "</DataArray>\n"
			<< "</Cells>\n"

			//////////////////////////////////////// CELL DATA ///////////////////////////////////////////////
			<< "<CellData>\n"
			<< "<DataArray Name=\"GroupID\" type=\"Int32\" NumberOfComponents=\"1\" format=\""
			<< utv_paraview_vtk_xml_data.format << "\">\n";
	nodes_data.reserve(mmr_get_nr_elem(Mesh_id));
	for (int el_id = 0; (el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0;) {
		nodes_data.push_back(mmr_el_groupID(Mesh_id, el_id));
	}
	utr_paraview_do_write(partmesh_file, nodes_data.data(), nodes_data.size());
	partmesh_file << "</DataArray>\n";
	nodes_data.clear();

	if (utr_mat_get_n_materials() > 0) {

		partmesh_file << "<DataArray Name=\"MaterialID\" type=\"Int32\" NumberOfComponents=\"1\" format=\""
				<< utv_paraview_vtk_xml_data.format << "\">\n";
		for (int el_id = 0; (el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0;) {
			nodes_data.push_back(utr_mat_get_matID(mmr_el_groupID(Mesh_id, el_id)));
		}
		utr_paraview_do_write(partmesh_file, nodes_data.data(), nodes_data.size());
		partmesh_file << "</DataArray>\n";
		nodes_data.clear();
	}

	if (utr_bc_get_n_block_assignments() > 0) {
		partmesh_file << "<DataArray Name=\"BlockID\" type=\"Int32\" NumberOfComponents=\"1\" format=\""
				<< utv_paraview_vtk_xml_data.format << "\">\n";
		for (int el_id = 0; (el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0;) {
			nodes_data.push_back(utr_bc_get_blockID(mmr_el_groupID(Mesh_id, el_id)));
		}
		utr_paraview_do_write(partmesh_file, nodes_data.data(), nodes_data.size());
		partmesh_file << "</DataArray>\n";
		nodes_data.clear();
	}
	partmesh_file << "</CellData>\n";

	partmesh_file.close();

	utv_paraview_vtk_xml_data.partmesh_file = partmesh_filename;

	std::string mesh_filename(Workdir);
	mesh_filename.append("/").append(Filename).append(".vtu");

	//printf("REFERENCES TO BOOST_FILESYSTEM TEMPORARILY SWITCHED OFF DUE TO SOME ERRORS IN GCC\n");
	const boost::filesystem::path p_partmesh(partmesh_filename),
		  p_mesh(mesh_filename);
	boost::filesystem::copy_file(p_partmesh, p_mesh, boost::filesystem::copy_option::overwrite_if_exists);

	std::ofstream mesh_file(mesh_filename.c_str(), std::ios_base::app);
	mesh_file << "</Piece>\n"
			<< "</UnstructuredGrid>\n"
			<< "</VTKFile>\n";
	mesh_file.close();

}

/*---------------------------------------------------------
/// utr_write_paraview_partmesh - to dump mesh in Paraview format includiong partition info
---------------------------------------------------------*/
void utr_write_paraview_partmesh_legacy_VTK(int Mesh_id, const char * Workdir, const char * Filename, int * Part, ute_paraview_flags VTK_file_version)
{
	int i, el_id = 0, el_type, el_group;
	int max_node_id, ino;
	int numElems;
	double nodeCoor[3];
	int numConn = 0;
	int el_nodes[10];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	max_node_id = mmr_get_max_node_id(Mesh_id);
	numElems = mmr_get_nr_elem(Mesh_id);

	std::string filename(Workdir);
	filename.append("/").append(Filename);

	FILE * fp = fopen(filename.c_str(), "w");
	mf_check(fp != NULL, "Unable to open file to save(%s)", Filename);

	fprintf(fp, "# vtk DataFile Version 3.1\n");
	fprintf(fp, "crated by MODFEM\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d double\n", max_node_id);

	/* loop over vertices-nodes */
	for (ino = 1; ino <= max_node_id; ino++) {
		if (mmr_node_status(Mesh_id, ino) == MMC_ACTIVE) {
			mmr_node_coor(Mesh_id, ino, nodeCoor);
			fprintf(fp, "%.12lg %.12lg %.12lg\n", nodeCoor[0], nodeCoor[1], nodeCoor[2]);
		}
		else {
			fprintf(fp, "0.0 0.0 0.0\n");
		}
	}

	/* loop over elements */
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		numConn += mmr_el_node_coor(Mesh_id, el_id, NULL, NULL);
	}

	fprintf(fp, "CELLS %d %d\n", numElems, numConn + numElems);

	/* loop over elements */
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		mmr_el_node_coor(Mesh_id, el_id, el_nodes, NULL);
		if (el_nodes[0] == 6) { // PRISM
			fprintf(fp, "6 %d %d %d %d %d %d\n",
					el_nodes[1] - 1, el_nodes[2] - 1, el_nodes[3] - 1,
					el_nodes[4] - 1, el_nodes[5] - 1, el_nodes[6] - 1);
		}
		else if (el_nodes[0] == 4) { //TETRA
			fprintf(fp, "4 %d %d %d %d\n", el_nodes[1] - 1, el_nodes[2] - 1,
					el_nodes[3] - 1, el_nodes[4] - 1);
		}
		else {
			fprintf(fp, "ERROR: WRONG ELEMENT TYPE");
		}
	}

	fprintf(fp, "CELL_TYPES %d\n", numElems);
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		el_type = mmr_el_type(Mesh_id, el_id);
		if (el_type == MMC_PRISM) fprintf(fp, "13\n"); // PRISM
		else if (el_type == MMC_TETRA) fprintf(fp, "10\n"); // TETRA
	}

	// CELL DATA BEGIN
	fprintf(fp, "CELL_DATA %d\n", numElems);
	// GROUPS ID
	fprintf(fp, "SCALARS GroupID int\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		el_group = mmr_el_groupID(Mesh_id, el_id);
		fprintf(fp, "%d\n", el_group);
	}
	// MATERIAL IDS
	fprintf(fp, "SCALARS MateriaID int\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		el_group = utr_mat_get_matID(mmr_el_groupID(Mesh_id, el_id));
		fprintf(fp, "%d\n", el_group);
	}

	// BLOCK IDS
	if (utr_bc_get_n_block_assignments() > 0) {
		fprintf(fp, "SCALARS BlockID int\n");
		fprintf(fp, "LOOKUP_TABLE default\n");
		while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
			el_group = utr_bc_get_blockID(mmr_el_groupID(Mesh_id, el_id));
			fprintf(fp, "%d\n", el_group);
		}
	}

	if (Part != NULL)
	{
		fprintf(fp, "SCALARS ProcessorsDistribuation float\n");
		fprintf(fp, "LOOKUP_TABLE default\n");

		for (i = 0; i < numElems; i++)
		{
			fprintf(fp, "%d\n", Part[i]);
		}
	}

	fclose(fp);
}

void utr_get_solinfos(const int Field_id, PVSolInfo * solInfos)
{
	const int mesh_id = apr_get_mesh_id(Field_id);
	const int max_node_id = mmr_get_max_node_id(mesh_id);
	const int nreq = apr_get_nreq(Field_id);
	for (int ino = 0; ino <= max_node_id; ino++) {
		for (int idof = 0; idof < nreq; idof++) {
			solInfos[ino].dofs[idof] = 0.0;
		}
	}

	int el_nodes[MMC_MAXELVNO + 1] = {0};

	/* solution dofs in an element - linear aprox.*/
	//double el_dofs[APC_MAXELSD] = {0};
	double el_dofs[APC_DOUBLE_MAXELSD] = {0}; // for satisfy plast_flow module in some situation

	/* loop over elements */
	int el_id = 0;
	while ((el_id = mmr_get_next_act_elem(mesh_id, el_id)) != 0) {
		int dof_counter = 0;
		int nrnodes = mmr_el_node_coor(mesh_id, el_id, el_nodes, NULL);
		apr_get_el_dofs(Field_id, el_id, 1, el_dofs);
		for (int ino = 1; ino <= nrnodes; ino++) {
			int nrdofs = apr_get_ent_nrdofs(Field_id, APC_VERTEX, ino);
			for (int idof = 0; idof < nreq; idof++) {
				int node_id = el_nodes[ino];
				/*kbw
				printf("node %d, value %.12lg (dof_counter %d, idof %d)\n",node_id, solInfos[node_id].dofs[0], dof_counter, idof);
				/*kew*/
				mf_check(dof_counter < (sizeof el_dofs) / sizeof(double),
						"Dof counter (%d) exceeds el_dofs size (%d)",
						dof_counter, (sizeof el_dofs) / sizeof(double) );
				solInfos[node_id].dofs[idof] = el_dofs[dof_counter];
				//Fk-for testing material
				//solInfos[node_id].dofs[idof]=mmr_el_groupID(mesh_id, el_id);

				dof_counter++;
// 				if(dof_counter > APC_MAXELSD) {
// 				  printf("Number of dofs %d > maximum allowed %d in utr_get_solinfos\n",
// 					dof_counter,APC_MAXELSD);
// 				  printf(" Recompile utr_get_solinfos with greater APC_MAXELSD. Exiting!\n");
// 				  exit(0);
// 				}
			}
		}
	}
}

typedef std::map<int, PVSolInfo *> MAP_F2S;

void utr_write_paraview_fields_XML_VTU(  const char * Filename,
		const int N_desc,
		utt_paraview_field_descriptor * Desc,
		MAP_F2S & fieldID2solInfos)
{
	//printf("REFERENCES TO BOOST_FILESYSTEM TEMPORARILY SWITCHED OFF DUE TO SOME ERRORS IN GCC\n");
	boost::filesystem::path partmesh_path(utv_paraview_vtk_xml_data.partmesh_file),
		  fullfile_path(Filename);

	boost::filesystem::copy_file(partmesh_path, fullfile_path, boost::filesystem::copy_option::overwrite_if_exists);

	// Creating content of the time file.
	std::ofstream file(Filename, std::ios_base::app );
	file << "<PointData>\n";

	std::vector<double> node_dofs;
	for (int i = 0; i < N_desc; ++i) {
		mf_check(Desc[i].entity_type == UTE_POINT_DATA, "CellData - to implement!");

		// Write XML node header
		file << "<DataArray "
				"Name=\"" << Desc[i].field_name << "\" "
				"type=\"" << utv_paraview_flags_strings_XML[Desc[i].value_type] << "\" "
				"NumberOfComponents=\"" << Desc[i].dofs_write[0] << "\" "
				"format=\"" << utv_paraview_vtk_xml_data.format << "\">\n";

		const int mesh_id = apr_get_mesh_id(Desc[i].field_id);
		const int max_node_id = mmr_get_max_node_id(mesh_id);

		// Write content

		node_dofs.reserve(max_node_id * Desc[i].dofs_write[0]);
		const PVSolInfo * solInfos = fieldID2solInfos[Desc[i].field_id];
		for (int n = 1; n <= max_node_id; n++) {
			for (int idofs = 1; idofs <= Desc[i].dofs_write[0]; idofs++) {
				node_dofs.push_back(solInfos[n].dofs[Desc[i].dofs_write[idofs]]);
			}
		}
		utr_paraview_do_write(file, node_dofs.data(), node_dofs.size());
		node_dofs.clear();
		// Write end stuff.
		file << "</DataArray>\n";
	}

	// End whole file.
	file << "</PointData>"
			"</Piece>\n"
			"</UnstructuredGrid>\n"
			"</VTKFile>\n";
	file.close();
}

void utr_write_paraview_fields_legacy_VTK(  const char * Filename,
		const int N_desc,
		utt_paraview_field_descriptor * Desc,
		MAP_F2S & fieldID2solInfos)
{
	FILE * fp = fopen(Filename, "a");
	mf_check(fp != NULL, "Unable to open file to save(%s)", Filename);

	for (int i = 0; i < N_desc; ++i) {
		const int mesh_id = apr_get_mesh_id(Desc[i].field_id);
		const int max_node_id = mmr_get_max_node_id(mesh_id);

		/// Paraview - specyfic statements.
		fprintf(fp, "%s %d\n", utv_paraview_flags_strings_legacy[Desc[i].entity_type], max_node_id);
		fprintf(fp, "%s %s %s 1\n",
				utv_paraview_flags_strings_legacy[Desc[i].quantity_type],
				Desc[i].field_name,
				utv_paraview_flags_strings_legacy[Desc[i].value_type]);
		fprintf(fp, "LOOKUP_TABLE default\n");

		PVSolInfo * solInfos = fieldID2solInfos[Desc[i].field_id];

		for (int n = 1; n <= max_node_id; n++) {

			for (int idofs = 1; idofs <= Desc[i].dofs_write[0]; idofs++) {
				fprintf(fp, "%.12lg\n", solInfos[n].dofs[Desc[i].dofs_write[idofs]]);
			}
			fprintf(fp, "\n");
		}

	}
	fclose(fp);
}

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------
/// utr_write_paraview_field - to dump field in Paraview format
---------------------------------------------------------*/
int utr_write_paraview_fields(const char * Work_dir,
		const char * Filename,
		double Current_time,
		const int N_desc,
		utt_paraview_field_descriptor * Desc,
		ute_paraview_flags VTK_file_version)
{

	int i, nrdofs, idofs;

	int numNodes, max_node_id, numElems, nreq, ino, idof, el_id, nrnodes;
	int mesh_id, el_type;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	mf_check(Filename != NULL, "Filename not set!");

	int part = 0;

#ifdef PARALLEL
	part = pcv_my_proc_id;
#endif


	/// 1. GATHERING FIELD DATA
	/// Mapping from Field_id to PVsolInfo array.
	MAP_F2S   fieldID2solInfos;

	// Find sol_infos for all unique desired fields.
	for (int i = 0; i < N_desc; ++i) {
		if (fieldID2solInfos.find(Desc[i].field_id) == fieldID2solInfos.end()) {
			int mesh_id = apr_get_mesh_id(Desc[i].field_id);
			int max_node_id = mmr_get_max_node_id(mesh_id);
			//int numElems = mmr_get_nr_elem(mesh_id);

			int nreq = apr_get_nreq(Desc[i].field_id);
			if (nreq > PDC_MAXEQ) {
				printf("Number of equations %d > maximum allowed in utr_write_paraview_paraview %d (field_id %d, %s)\n",
						nreq, PDC_MAXEQ, Desc[i].field_id, Desc[i].field_name );
				printf(" Recompile uts_write_paraview with greater PDC_MAXEQ. Exiting!\n");
				exit(0);
			}

			// dof_entities (nodes) are numbered from 1 to max_node_id
			PVSolInfo * solInfos = (PVSolInfo *)malloc((max_node_id + 1) * sizeof(PVSolInfo));
			mf_check_mem(solInfos);
			utr_get_solinfos( Desc[i].field_id, solInfos);
			fieldID2solInfos[ Desc[i].field_id ] = solInfos;
		}
	}


	/// 2. Writing initial information etc., but not fields.
	std::string data_filename(Work_dir);
	data_filename.append("/").append(Filename);
	if (VTK_file_version == UTE_VTK_XML) {
		// writing XML description into time file
		data_filename = utr_vtk_xml_get_timestep_filename(&utv_paraview_vtk_xml_data,
						Current_time,
						0,
						part,
						Work_dir,
						Filename);
	}
	else { // legacy VTK
		data_filename.append(".vtk");
	}

	/// 3. Writing all fields.
	if (VTK_file_version == UTE_VTK_XML) {
		utr_write_paraview_fields_XML_VTU(data_filename.c_str(), N_desc, Desc, fieldID2solInfos);
	}
	else { // legacy VTK
		utr_write_paraview_fields_legacy_VTK(data_filename.c_str(), N_desc, Desc, fieldID2solInfos);
	}

	/// 4. Ending and cleaning.
	/// Freeing memory.
	for (MAP_F2S::iterator it = fieldID2solInfos.begin();
			it != fieldID2solInfos.end();
			++it) {
		free(it->second);
	}
	return 0;
}






/*---------------------------------------------------------
/// utr_write_paraview_mesh - to dump mesh in Paraview format
---------------------------------------------------------*/
int utr_write_paraview_mesh(int Mesh_id,
		const char * Work_dir,
		const char * Filename,
		ute_paraview_flags VTK_file_version)
{
	utr_write_paraview_partmesh(Mesh_id, Work_dir, Filename, NULL, VTK_file_version);
	return 0;
}

void utr_write_paraview_partmesh(int Mesh_id, const char * Work_dir, const char * Filename, int * Part, ute_paraview_flags VTK_file_version)
{
	mf_check(Filename != NULL, "Filename not set!");
	if (VTK_file_version == UTE_VTK_XML) {
		utr_write_paraview_partmesh_XML_VTK(Mesh_id, Work_dir, Filename, Part, VTK_file_version);
	}
	else {
		utr_write_paraview_partmesh_legacy_VTK(Mesh_id, Work_dir, Filename, Part, VTK_file_version);
	}
}




/*---------------------------------------------------------
/// utr_write_paraview_bc - to dump boundary conditions 'field' in Paraview format
---------------------------------------------------------*/
int utr_write_paraview_bc(int Mesh_id,
		const char * Work_dir,
		const char * Filename
		, ute_paraview_flags VTK_file_version)
{

	int i, faceId, bc = 0;
	int * nodeBCs;
	int numNodes, max_node_id, numFaces;
	int fNodes[5] = {0};

	double el_dofs[APC_MAXELSD];  /* solution dofs in an element */

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	max_node_id = mmr_get_max_node_id(Mesh_id);
	numFaces = mmr_get_nr_face(Mesh_id);

	std::string filename(Work_dir);
	filename.append("/").append(Filename);

	FILE * fp = fopen(filename.c_str(), "a");

	mf_check(fp != NULL, "Unable to open file to save(%s)", Filename);

	fprintf(fp, "POINT_DATA %d\n", max_node_id);
	fprintf(fp, "SCALARS bound_cond double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	// dof_entities (nodes) are numbered from 1 to max_node_id
	nodeBCs = (int *)malloc((max_node_id + 1) * sizeof(int));
	memset(nodeBCs, 0, (max_node_id + 1)* sizeof(int));
	// gather bc info for nodes
	faceId = 0;
	while (faceId = mmr_get_next_face_all(Mesh_id, faceId)) {
		bc = mmr_fa_bc(Mesh_id, faceId);
		if (bc != 0) {
			mmr_fa_node_coor(Mesh_id, faceId, fNodes, NULL);
			for (i = 1; i <= fNodes[0]; ++i) {
				if (nodeBCs[fNodes[i]] < bc) {
					nodeBCs[fNodes[i]] = bc;
				}
			}
		}
	}

	for (i = 1; i <= max_node_id; i++) {
		fprintf(fp, "%d\n\n", nodeBCs[i]);
	}

	fclose(fp);
	free(nodeBCs);

	return 0;
}


/******************* OLD INTERFACE *********************************/


/// \file uts_dump_paraview
/// \brief generic routines for dumping mesh and/or field into .vtu files (ParaView format).
///
/// int utr_write_paraview_mesh - to dump mesh in Paraview format
/// int utr_write_paraview_field - to dump field in Paraview format
///-----------------------------------------------
/// History:
/// 05.2011 - Kazimierz Michalik, kamich@agh.edu.pl, initial version
/// 06.2011 - KGB, pobanas@cyf-kr.edu.pl - utr_write_paraview_std_lin



int utr_write_paraview_std_lin_VTK_XML(int Mesh_id, int Field_id, const char * Workdir, const char * Filename, ute_paraview_flags VTK_file_version)
{

	utr_write_paraview_partmesh_XML_VTK(Mesh_id, Workdir, Filename, NULL, VTK_file_version);
	utt_paraview_field_descriptor field_desc;
	field_desc.dofs_write[0] = 1;
	field_desc.dofs_write[1] = 0;
	field_desc.entity_type = UTE_POINT_DATA;
	field_desc.field_id = Field_id;
	field_desc.field_name = "conf_diff_val";
	field_desc.quantity_type = UTE_SCALARS;
	field_desc.value_type = UTE_DOUBLE;
	return utr_write_paraview_fields(Workdir, Filename, 0.0, 1, & field_desc, VTK_file_version);
}

int utr_write_paraview_std_lin_legacy_VTK(int Mesh_id, int Field_id,
		const char * Workdir,
		const char * Filename,
		ute_paraview_flags VTK_file_version)
{
	FILE * fp;
	int i, ino, el_type, el_id, nrnodes, nreq, nrdofs, idof; //, node_id
	int numElems, max_node_id;
	double nodeCoor[3];
	int numConn = 0;
	int el_nodes[10];

	PVSolInfo * solInfos;
	double el_dofs[APC_MAXELSD];  /* solution dofs in an element */

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* open the output file */
	std::string filename(Workdir);
	filename.append("/").append(Filename);
	fp = fopen(filename.c_str(), "w");
	if (fp == NULL) {
		printf("Cannot open file '%s' for Paraview mesh data\n", Filename);
		return (-1);
	}

	max_node_id = mmr_get_max_node_id(Mesh_id);
	numElems = mmr_get_nr_elem(Mesh_id);

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "Navier-Stokes: pdd_navstokes\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d double\n", max_node_id);

	/* loop over vertices-nodes */
	for (ino = 1; ino <= max_node_id; ino++) {
		if (mmr_node_status(Mesh_id, ino) == MMC_ACTIVE) {
			mmr_node_coor(Mesh_id, ino, nodeCoor);
			fprintf(fp, "%.12lg %.12lg %.12lg\n", nodeCoor[0], nodeCoor[1], nodeCoor[2]);
		}
		else {
			fprintf(fp, "0.0 0.0 0.0\n");
		}
	}

	/* loop over elements */
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		numConn += mmr_el_node_coor(Mesh_id, el_id, NULL, NULL);
	}

	fprintf(fp, "CELLS %d %d\n", numElems, numConn + numElems);

	/* loop over elements */
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		mmr_el_node_coor(Mesh_id, el_id, el_nodes, NULL);
		if (el_nodes[0] == 6) { // PRISM
			fprintf(fp, "6 %d %d %d %d %d %d\n",
					el_nodes[1] - 1, el_nodes[2] - 1, el_nodes[3] - 1,
					el_nodes[4] - 1, el_nodes[5] - 1, el_nodes[6] - 1);
		}
		else if (el_nodes[0] == 4) { //TETRA
			fprintf(fp, "4 %d %d %d %d\n", el_nodes[1] - 1, el_nodes[2] - 1,
					el_nodes[3] - 1, el_nodes[4] - 1);
		}
		else {
			fprintf(fp, "ERROR: WRONG ELEMENT TYPE");
		}
	}

	fprintf(fp, "CELL_TYPES %d\n", numElems);
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		el_type = mmr_el_type(Mesh_id, el_id);
		if (el_type == MMC_PRISM) fprintf(fp, "13\n"); // PRISM
		if (el_type == MMC_TETRA) fprintf(fp, "10\n"); // TETRA
	}

	/* writing field data */

	nreq = apr_get_nreq(Field_id);
	if (nreq > PDC_MAXEQ) {
		printf("Number of equations %d > maximum allowed in uts_dump_paraview %d\n",
				nreq, PDC_MAXEQ);
		printf(" Recompile uts_dump_paraview with greater PDC_MAXEQ. Exiting!\n");
		exit(0);
	}

	solInfos = (PVSolInfo *)malloc((max_node_id + 1) * sizeof(PVSolInfo));
	for (ino = 0; ino <= max_node_id; ino++) {
		for (idof = 0; idof < nreq; idof++) {
			solInfos[ino].dofs[idof] = 0.0;
		}
	}

	/* loop over elements */
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(Mesh_id, el_id)) != 0) {
		int dof_counter = 0;
		nrnodes = mmr_el_node_coor(Mesh_id, el_id, el_nodes, NULL);
		apr_get_el_dofs(Field_id, el_id, 1, el_dofs);
		for (ino = 1; ino <= nrnodes; ino++) {
			nrdofs = apr_get_ent_nrdofs(Field_id, APC_VERTEX, ino);
			for (idof = 0; idof < nreq; idof++) {
				int node_id = el_nodes[ino];
				solInfos[node_id].dofs[idof] = el_dofs[dof_counter];
				//Fk-for testing material
				//solInfos[node_id].dofs[idof]=mmr_el_groupID(Mesh_id, el_id);
				/*kbw
					printf("node %d, value %.12lg\n",node_id, solInfos[node_id].dofs[0]);
				/*kew*/
				dof_counter++;
			}
		}
	}

	// NASTY HACKS - SHOULD BE MUCH MORE INPUT CONTROL PARAMETERS

	fprintf(fp, "POINT_DATA %d\n", max_node_id);

	if (nreq == 1) { // Laplace

		fprintf(fp, "SCALARS field_value double 1\n");
		fprintf(fp, "LOOKUP_TABLE default\n");

		for (i = 1; i <= max_node_id; i++) {
			fprintf(fp, "%.12lg\n", solInfos[i].dofs[0]);
			/*kbw
				  printf("node %d, value %.12lg\n",i, solInfos[i].dofs[0]);
			/*kew*/
		}

	}
	else if (nreq == 3) { // elasticity

		fprintf(fp, "VECTORS displacement double\n");
		for (i = 1; i <= max_node_id; i++) {
			fprintf(fp, "%.12lg %.12lg %.12lg\n",
					solInfos[i].dofs[0], solInfos[i].dofs[1], solInfos[i].dofs[2]);
		}

	}
	else if (nreq == 4) { // NS

		fprintf(fp, "SCALARS pressure double 1\n");
		fprintf(fp, "LOOKUP_TABLE default\n");
		for (i = 1; i <= max_node_id; i++) {
			fprintf(fp, "%.12lg\n", solInfos[i].dofs[nreq - 1]);
		}

		fprintf(fp, "VECTORS velocity double\n");
		for (i = 1; i <= max_node_id; i++) {
			fprintf(fp, "%.12lg %.12lg %.12lg\n",
					solInfos[i].dofs[0], solInfos[i].dofs[1], solInfos[i].dofs[2]);
		}

	}
	else if (nreq == 5) { // NS+THERMO

		fprintf(fp, "SCALARS pressure double 1\n");
		fprintf(fp, "LOOKUP_TABLE default\n");
		for (i = 1; i <= max_node_id; i++) {
			fprintf(fp, "%.12lg\n", solInfos[i].dofs[nreq - 2]);
		}

		fprintf(fp, "VECTORS velocity double\n");
		for (i = 1; i <= max_node_id; i++) {
			fprintf(fp, "%.12lg %.12lg %.12lg\n",
					solInfos[i].dofs[0], solInfos[i].dofs[1], solInfos[i].dofs[2]);
		}

		fprintf(fp, "SCALARS temperature double 1\n");
		fprintf(fp, "LOOKUP_TABLE default\n");
		for (i = 1; i <= max_node_id; i++) {
			fprintf(fp, "%.12lg\n", solInfos[i].dofs[nreq - 1]);
		}

	}
	else {

		printf("write_Paraview too stupid to handle %d solution vector components\n",
				nreq);

	}

	fclose(fp);

	free(solInfos);

	return 0;
}

/// utr_write_paraview_std_lin - to dump std_lin mesh and field in Paraview format
/** \param Mesh_id - id of the mesh associated with given field Field_id
	\param Field_id - id of the field to dump
	\param Filename - c-string with name of the file to write on disk
	\param Desc - c-string c-array with name of the field values
	(at least Desc[0]!= NULL should be passed)
*/
int utr_write_paraview_std_lin(int Mesh_id, int Field_id, const char * Workdir, const char * Filename, ute_paraview_flags VTK_file_version)
{
	int rc = 0;
	if (VTK_file_version == UTE_VTK_XML) {
		rc = utr_write_paraview_std_lin_VTK_XML(Mesh_id, Field_id, Workdir, Filename, VTK_file_version);
	}
	else {
		rc = utr_write_paraview_std_lin_legacy_VTK(Mesh_id, Field_id, Workdir, Filename, VTK_file_version);
	}
	return rc;
}



//typedef struct{
//  int elNodes[7];
//  int numNodes;
//  int elType;
//}PVElInfo;

///// utr_write_paraview_mesh_old - to dump mesh in Paraview format
///** \param Mesh_id - id of the mesh to dump
//	\param Filename - c-string with name of the file to write(dump) on disk
//*/
//int utr_write_paraview_mesh_old(int Mesh_id, char *Filename)
//{
//  FILE *fp;
//  int i,el_id=0;
//  PVElInfo * elInfos;
//  int numNodes;
//  int numElems;
//  double nodeCoor[3];
//  int numConn = 0;


///*++++++++++++++++ executable statements ++++++++++++++++*/

///* open the output file */
//  fp = fopen(Filename, "w");
//  if(fp==NULL) {
//	printf("Cannot open file '%s' for Paraview mesh data\n",Filename);
//	return(-1);
//  }


//  numNodes = mmr_get_nr_node(Mesh_id);
//  numElems = mmr_get_nr_elem(Mesh_id);
//  elInfos = (PVElInfo*)malloc((numElems+1) * sizeof(PVElInfo));



//  fprintf(fp,"# vtk DataFile Version 2.0\n");
//  fprintf(fp,"Navier-Stokes: pdd_navstokes\n");
//  fprintf(fp,"ASCII\n");
//  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
//  fprintf(fp,"POINTS %d double\n",numNodes);
//  for(i=1;i<=numNodes;i++){
//	mmr_node_coor(Mesh_id, i, nodeCoor);
//	fprintf(fp,"%.12lg %.12lg %.12lg\n",nodeCoor[0], nodeCoor[1], nodeCoor[2]);
//  }
//  for(i=1,el_id=0;(el_id=mmr_get_next_act_elem(Mesh_id,el_id)) > 0;i++){
//	mmr_el_node_coor(Mesh_id, el_id, elInfos[i].elNodes, NULL);
//	//apr_get_el_dofs(Field_id, i, 1, elInfos[i].solDofs);
//	elInfos[i].numNodes = elInfos[i].elNodes[0];
//	if(elInfos[i].numNodes == 6)
//	  elInfos[i].elType = 13;
//	else if(elInfos[i].numNodes == 4)
//	  elInfos[i].elType = 10;
//	else{
//	  elInfos[i].elType = -1;
//	}
//	numConn += elInfos[i].numNodes;
//  }

//  fprintf(fp,"CELLS %d %d\n",numElems, numConn+numElems);

//  for(i=1;i<=numElems;i++){
//	  if(elInfos[i].elType == 13) //PRISM
//	fprintf(fp,"6 %d %d %d %d %d %d\n", elInfos[i].elNodes[1]-1, elInfos[i].elNodes[2]-1, elInfos[i].elNodes[3]-1, elInfos[i].elNodes[4]-1, elInfos[i].elNodes[5]-1, elInfos[i].elNodes[6]-1);
//	  else if(elInfos[i].elType == 10) //TETRA
//	fprintf(fp,"4 %d %d %d %d\n", elInfos[i].elNodes[1]-1, elInfos[i].elNodes[2]-1, elInfos[i].elNodes[3]-1, elInfos[i].elNodes[4]-1);
//	  else
//	fprintf(fp, "ERROR: WRONG ELEMENT TYPE");
//  }

//  fprintf(fp,"CELL_TYPES %d\n",numElems);
//  for(i=1;i<=numElems;i++){
//	fprintf(fp,"%d\n", elInfos[i].elType);
//  }

//  fclose(fp);
//  free(elInfos);

//  return 0;
//}

///// utr_write_paraview_field_old - to dump field in Paraview format
///** \param Mesh_id - id of the mesh associated with given field Field_id
//	\param Field_id - id of the field to dump
//	\param Filename - c-string with name of the file to write on disk
//	\param Desc - c-string c-array with name of the field values
//	(at least Desc[0]!= NULL should be passed)
// */
//int utr_write_paraview_field_old(int Mesh_id, int Field_id, char *Filename, char ** Desc)
//{
//  FILE *fp;
//  int i,nrdofs;
//  PVSolInfo * solInfos;
//  int numNodes;


///*++++++++++++++++ executable statements ++++++++++++++++*/

///* open the output file */
//  fp = fopen(Filename, "w");
//  if(fp==NULL) {
//	printf("Cannot open file '%s' for Paraview field data\n",Filename);
//	return(-1);
//  }

//  numNodes = mmr_get_nr_node(Mesh_id);
//  solInfos = (PVSolInfo*)malloc((numNodes+1) * sizeof(PVSolInfo));

//  fprintf(fp,"POINT_DATA %d\n", numNodes);
//  fprintf(fp,"%s",Desc[0]);
//  fprintf(fp,"LOOKUP_TABLE default\n");



//  for(i=1;i<=numNodes;i++){
//	nrdofs=apr_get_ent_nrdofs(Field_id,APC_VERTEX,i);
//	apr_read_ent_dofs(Field_id,APC_VERTEX, i, nrdofs, 1, solInfos[i].dofs);
//  }

//  for(i=1;i<=numNodes;i++){
//	fprintf(fp,"%.12lg\n",solInfos[i].dofs[0]);
//  }

//  if(nrdofs > 1) {
//	fprintf(fp,"%s",Desc[1]);
//	for(i=1;i<=numNodes;i++){
//	fprintf(fp,"%.12lg %.12lg %.12lg\n",solInfos[i].dofs[0], solInfos[i].dofs[1], solInfos[i].dofs[2]);
//	}
//  }
//  fclose(fp);

//  free(solInfos);

//  return 0;
//}
///// utr_merge_paraview_mesh_field - to merge previously dumped Paraview mesh and Paraview field files
///// into one .vtu file.
///** \param meshFile - c-string with Paraview mesh filename
//	\param fieldFile - c-string with Paraview field filename
//	\param mergedFile - c-string with name of the file to write on disk
// */
//int utr_merge_paraview_mesh_field(char * meshFile, char * fieldFile, char * mergedFile)
//{
//	char call[255]={0};
//	sprintf(call, "cat %s %s > %s",meshFile,fieldFile,mergedFile);
//	system(call);
//	return 0;
//}


#ifdef __cplusplus
}
#endif
