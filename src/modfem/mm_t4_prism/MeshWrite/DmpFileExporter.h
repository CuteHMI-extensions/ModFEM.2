#ifndef DMPFILEEXPORTER_H
#define DMPFILEEXPORTER_H

#include "MeshFileWriter.h"

namespace MeshWrite
{

class DmpFileExporter : public MeshFileWriter
{
    public:
        DmpFileExporter();
        DmpFileExporter(const std::string & file_name);
        virtual ~DmpFileExporter();

    /** Initializes exporter. After this function call exporter should be ready for exporting the mesh.
		/return true if exporter initializes correctly, otherwise false.
				If this function doesn't succeeded exporting process is gently stopped and error is thrown.
	*/
	bool Init();

	bool Init(const std::string & file_name);

	/** Frees all resources allocated by class during mesh export.
	*/
	void Free();

	bool	doWrite(const hHybridMesh * mesh);

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
	void WriteEdge(const ID verts[], const uTind type);

	/** /param noVert Number of all vertices (global) in mesh.
	*/
	void WriteFacesCount(const int noFaces);

	/** /param coords array of vertex coordinates
	*/
	void WriteFace(const ID edges[],const uTind type,const int8_t bc, const ID neigh[]);

	/** /param noElems Number of all elements in mesh.
	*/
	void WriteElementCount(const int noElems);

	/** /param vertices array of indexes of vertices making this element
		/param neighbours array of indexes of neighbour elements
	*/
	void WriteElement(const ID faces[], const ID neighbours[],const uTind type, const ID father, const int8_t ref,const int8_t material);
};

};

#endif // DMPFILEEXPORTER_H
