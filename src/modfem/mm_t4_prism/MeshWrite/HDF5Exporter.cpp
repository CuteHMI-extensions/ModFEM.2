
#include "HDF5Exporter.h"
#include "../MeshModule/hHybridMeshTypes.h"

const hsize_t HDF5Exporter::dims[25] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};


HDF5Exporter::HDF5Exporter(const std::string & file_name)
:	coords_t(PredType::NATIVE_DOUBLE,1,  dims + 3),
	sons3_t(PredType::NATIVE_INT,1,  dims + 3),
	sonsFace_t(PredType::NATIVE_INT,1, dims + Face3::sons_no),
	hObj_t(sizeof(hObj)),
	vertex_t( sizeof(Vertex1)),
	edge_t(sizeof(Edge2)),
	side4_t(sizeof(Face3)),
	elemT4_t(sizeof(ElemT4)),
	elemPrism_t(sizeof(ElemPrism))
{
	// Creating compud datatypes for operating on HDF5 files.

//	hObj_t.insertMember( "id", HOFFSET(hObj, id), PredType::NATIVE_INT);
//	hObj_t.insertMember( "parent", HOFFSET(hObj, parent), PredType::NATIVE_INT);
//	hObj_t.insertMember( "level", HOFFSET(hObj, level), PredType::NATIVE_INT8);

//	vertex_t.insertMember("coords", HOFFSET(Vertex1, coords), coords_t);

//	edge_t.insertMember("sons", HOFFSET(Edge2, _sons), sons3_t);

//	side4_t.insertMember("_first", HOFFSET(Face3, neigh), ArrayType(PredType::NATIVE_INT,1, dims+2 ));
//    side4_t.insertMember("_sons", HOFFSET(Face3, _sons), ArrayType(PredType::NATIVE_INT,1, dims+Face3::sons_no ));

//	elemT4_t.insertMember("_vertices",  HOFFSET(ElemT4, _vertices), ArrayType(PredType::NATIVE_INT,1 ,dims + ElemT4::verts_no ));
//	elemT4_t.insertMember("_faces",     HOFFSET(ElemT4, _faces),    ArrayType(PredType::NATIVE_INT,1 ,dims + ElemT4::faces_no));
//	elemT4_t.insertMember("_edges",     HOFFSET(ElemT4, _edges),    ArrayType(PredType::NATIVE_INT,1 ,dims + ElemT4::edges_no));
//	elemT4_t.insertMember("_sides",     HOFFSET(ElemT4, _sides),    ArrayType(PredType::NATIVE_INT,1 ,dims + ElemT4::sides_no));
//	elemT4_t.insertMember("_sons",      HOFFSET(ElemT4, _sons),     ArrayType(PredType::NATIVE_INT,1 ,dims + ElemT4::sons_no));

//	elemPrism_t.insertMember("_vertices",  HOFFSET(ElemT4, _vertices), ArrayType(PredType::NATIVE_INT,1 ,dims + ElemPrism::verts_no));
//	elemPrism_t.insertMember("_faces",     HOFFSET(ElemT4, _faces),    ArrayType(PredType::NATIVE_INT,1 ,dims + ElemPrism::faces_no));
//	elemPrism_t.insertMember("_edges",     HOFFSET(ElemT4, _edges),    ArrayType(PredType::NATIVE_INT,1 ,dims + ElemPrism::edges_no));
//	elemPrism_t.insertMember("_sides",     HOFFSET(ElemT4, _sides),    ArrayType(PredType::NATIVE_INT,1 ,dims + ElemPrism::sides_no));
//	elemPrism_t.insertMember("_sons",      HOFFSET(ElemT4, _sons),     ArrayType(PredType::NATIVE_INT,1 ,dims + ElemPrism::sons_no));

}

bool	HDF5Exporter::Init()
{
	return true;
}

void    HDF5Exporter::Free()
{
}

void HDF5Exporter::WriteVerticesCount(const int noVerts){};

void HDF5Exporter::WriteVertex(const double coords[],const uTind type){};

void HDF5Exporter::WriteEdgesCount(const int noEdges){};

void HDF5Exporter::WriteEdge(const uTind verts[], const uTind type){};

void HDF5Exporter::WriteFacesCount(const int noFaces){};

void HDF5Exporter::WriteFace(const uTind edges[],const uTind type,const int8_t bc, const uTind neigh[]){};

void HDF5Exporter::WriteElementCount(const int noElems){};

void HDF5Exporter::WriteElement(const uTind faces[], const uTind neighbours[],const uTind type, const uTind father, const int8_t ref, const int8_t material){};
