#ifndef _HDF5EXPORTER_H_
#define _HDF5EXPORTER_H_

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
#endif  // H5_NO_STD
#endif

//#include "../HDF5/include/H5Cpp.h"
#include <H5Cpp.h>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

#include "IMeshWriter.h"


class HDF5Exporter : public MeshWrite::IMeshWriter
{
public:
	HDF5Exporter(const std::string & file_name) ;
	~HDF5Exporter(void) {} ;

    /** Initializes exporter. After this function call exporter should be ready for exporting the mesh.
		/return true if exporter initializes correctly, otherwise false.
				If this function doesn't succeeded exporting process is gently stopped and error is thrown.
	*/
	bool Init();
	bool Init(const std::string & name){ return false; }

	/** Frees all resources allocated by class during mesh export.
	*/
	void Free();

	/** /param noVert Number of all vertices (global) in mesh.
	*/
	void WriteVerticesCount(const int noVerts);

	/** /param coords array of vertex coordinates
	*/
	void WriteVertex(const double coords[],const uTind type);

    /** /param noVert Number of all vertices (global) in mesh.
	*/
	void WriteEdgesCount(const int noEdges);

	/** /param coords array of vertex coordinates
	*/
	void WriteEdge(const uTind verts[], const uTind type);

	/** /param noVert Number of all vertices (global) in mesh.
	*/
	void WriteFacesCount(const int noFaces);

	/** /param coords array of vertex coordinates
	*/
	void WriteFace(const uTind edges[],const uTind type,const int8_t bc, const uTind neigh[]);

	/** /param noElems Number of all elements in mesh.
	*/
	void WriteElementCount(const int noElems);

	/** /param vertices array of indexes of vertices making this element
		/param neighbours array of indexes of neighbour elements
	*/
	void WriteElement(const uTind faces[], const uTind neighbours[],const uTind type, const uTind father, const int8_t ref, const int8_t material);

protected:


	H5File		file;

	ArrayType	coords_t, sons3_t, sonsFace_t;
	CompType	hObj_t,
				vertex_t,
				edge_t,
				side4_t,

				elemT4_t,
				elemPrism_t;

    static const hsize_t     dims[25];

	//const int       NX = 5;                    // dataset dimensions
	//const int       NY = 6;
	//const int       RANK = 2;

};
#endif // _HDF5EXPORTER_H_
