#ifndef BINARY_FILE_WRITER_H_
#define BINARY_FILE_WRITER_H_

#include "IMeshWriter.h" 

#include <fstream>

namespace MeshWrite {

class BinaryFileWriter : public MeshWrite::IMeshWriter
	{
	public:
	BinaryFileWriter();
	BinaryFileWriter(const std::string & fileName);
	~BinaryFileWriter();

	bool Init();
	bool Init(const std::string & fileName);
    void Free();
    
    bool doWrite(const hHybridMesh *);
    
    template <typename T>
    friend BinaryFileWriter& operator<<(BinaryFileWriter& w, const T & obj);
    template <typename T>
    friend BinaryFileWriter& operator<<(BinaryFileWriter& w, const T * obj);
	
		/** /param noVert Number of all vertices (global) in mesh.
	*/
	void WriteVerticesCount(const int noVerts);

	/** /param coords array of vertex coordinates
	*/
	void WriteVertex(const double coords[],const uTind type) ;

    /** /param noVert Number of all vertices (global) in mesh.
	*/
	virtual void WriteEdgesCount(const int noEdges);

	/** /param coords array of vertex coordinates
	*/
	 void WriteEdge(const uTind verts[], const uTind type) ;

	/** /param noVert Number of all vertices (global) in mesh.
	*/
	 void WriteFacesCount(const int noFaces) ;

	/** /param coords array of vertex coordinates
	*/
	 void WriteFace(const uTind edges[],const uTind type,const int8_t bc, const uTind neigh[]) ;

	/** /param noElems Number of all elements in mesh.
	*/
	 void WriteElementCount(const int noElems) ;

	/** /param vertices array of indexes of vertices making this element
		/param neighbours array of indexes of neighbour elements
	*/
	 void WriteElement(const uTind faces[], const uTind neighbours[],const uTind type, const uTind father, const int8_t ref, const int8_t material) ;

	
private:
    std::string fileName_;
    std::ofstream oFile_;
};

template <typename T>
BinaryFileWriter& operator<<(BinaryFileWriter& w, const T & obj)
{
    w.oFile_.write(reinterpret_cast<char *>(&obj),sizeof(obj));
    return w;
}
    
template <typename T>
BinaryFileWriter& operator<<(BinaryFileWriter& w, const T * obj)
{
    w.oFile_.write(reinterpret_cast<char *>(obj),sizeof(T));
    return w;
}
 
};
 
#endif // BINARY_FILE_WRITER_H_
