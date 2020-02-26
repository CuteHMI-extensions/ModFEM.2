#include "BinaryFileWriter.h"
#include "../MeshModule/hHybridMesh.h"

#include <fstream>
#include <string>

namespace MeshWrite {
  
  BinaryFileWriter::BinaryFileWriter() : fileName_("not initialized!")
  {
  }

  BinaryFileWriter::BinaryFileWriter(const std::string & fileName) : fileName_(fileName)
  {
  }

  BinaryFileWriter::~BinaryFileWriter()
  {
	Free();
  }
  
  bool BinaryFileWriter::Init()
  {
	return Init(fileName_);
  }
  
bool BinaryFileWriter::Init(const std::string & filename)
{
    fileName_=filename;
    oFile_.close();
    oFile_.open(fileName_.c_str(),std::ios::in | std::ios::binary | std::ios::trunc);
	oFile_.clear();
	assert(oFile_.is_open() == true);
	return oFile_.is_open();
}

void BinaryFileWriter::Free()
{
    oFile_.close();
}
  
  // We use " " (space) as separator
bool BinaryFileWriter::doWrite(const hHybridMesh * m)
{
    if(oFile_.good()) {
	  //oFile_ << m->meshId_;
	  if(m->name_.empty()) {
		oFile_ << "no_mesh_name";
	  }
	  oFile_ << m->name_ << " ";
	  m->vertices_.write(oFile_);
	  oFile_ << " ";
	  m->edges_.write(oFile_);
	  oFile_ << " ";
	  m->faces_.write(oFile_);
	  oFile_ << " ";
	  m->elements_.write(oFile_);
	  oFile_ << " " << m->gen_
			 << " " << m->maxGen_
			 << " " << m->maxGenDiff_ << " ";
    }
	return true;
}
  
		/** /param noVert Number of all vertices (global) in mesh.
	*/
  void BinaryFileWriter::WriteVerticesCount(const int noVerts){}

	/** /param coords array of vertex coordinates
	*/
  void BinaryFileWriter::WriteVertex(const double coords[],const uTind type) {}

    /** /param noVert Number of all vertices (global) in mesh.
	*/
  void BinaryFileWriter::WriteEdgesCount(const int noEdges) {}

	/** /param coords array of vertex coordinates
	*/
  void BinaryFileWriter::WriteEdge(const uTind verts[], const uTind type) {} ;

	/** /param noVert Number of all vertices (global) in mesh.
	*/
  void BinaryFileWriter::WriteFacesCount(const int noFaces) {};

	/** /param coords array of vertex coordinates
	*/
  void BinaryFileWriter::WriteFace(const uTind edges[],const uTind type,const int8_t bc, const uTind neigh[]) {} ;

	/** /param noElems Number of all elements in mesh.
	*/
  void BinaryFileWriter::WriteElementCount(const int noElems) {};

	/** /param vertices array of indexes of vertices making this element
		/param neighbours array of indexes of neighbour elements
	*/
  void BinaryFileWriter::WriteElement(const uTind faces[], const uTind neighbours[],const uTind type, const uTind father, const int8_t ref, const int8_t material) {}
  
};
