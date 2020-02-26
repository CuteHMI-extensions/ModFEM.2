#ifndef IMESHWRITER_H_GUARD
#define IMESHWRITER_H_GUARD

/** /brief This interface should be implemented in all classes that want to export mesh from hpFEM.
	Exporting sequence is simple:
	A. Invoking hpMesh::Export(exporter) method, where "exporter" is an object implementing IMeshWriter interface.
	B. Export method calls on exporter:
		1. Init - which initializes exporter
		2. doWrite - which makes writing
		3. Free() which frees all resources used while exporting and finalize export.
*/

#include "../Common.h"
#include <string>

class hHybridMesh;

namespace MeshWrite
{

class IMeshWriter
{
public:
	IMeshWriter(void) : writeSequence(NULL) {};
	virtual ~IMeshWriter(void){};

	virtual bool doWrite(const hHybridMesh * mesh) = 0;

	static const char V_COUNT = 'v', VERTICES = 'V',
	E_COUNT = 'e', EDGES = 'E',
	F_COUNT ='f', FACES ='F',
	EL_COUNT = 'l', ELEMENTS = 'L',
	AUX='a', END = 0;

	/** In writeSequence is stored sequendce of actions used to correctly write output.
	*/
	char * writeSequence;

	/** Initializes exporter. After this function call exporter should be ready for exporting the mesh.
		/return true if exporter initializes correctly, otherwise false.
				If this function doesn't succeeded exporting process is gently stopped and error is thrown.
	*/
	virtual bool Init() = 0;
	virtual bool Init(const std::string & name) = 0;

	/** Frees all resources allocated by class during mesh export.
	*/
	virtual void Free() = 0;

	/** /param noVert Number of all vertices (global) in mesh.
	*/
	virtual void WriteVerticesCount(const int noVerts) = 0;

	/** /param coords array of vertex coordinates
	*/
	virtual void WriteVertex(const double coords[],const uTind type) = 0 ;

    /** /param noVert Number of all vertices (global) in mesh.
	*/
	virtual void WriteEdgesCount(const int noEdges) = 0;

	/** /param coords array of vertex coordinates
	*/
	virtual void WriteEdge(const uTind verts[], const uTind type) = 0 ;

	/** /param noVert Number of all vertices (global) in mesh.
	*/
	virtual void WriteFacesCount(const int noFaces) = 0;

	/** /param coords array of vertex coordinates
	*/
	virtual void WriteFace(const uTind edges[],const uTind type,const int8_t bc, const uTind neigh[]) = 0 ;

	/** /param noElems Number of all elements in mesh.
	*/
	virtual void WriteElementCount(const int noElems) = 0;

	/** /param vertices array of indexes of vertices making this element
		/param neighbours array of indexes of neighbour elements
	*/
	virtual void WriteElement(const uTind faces[], const uTind neighbours[],const uTind type, const uTind father, const int8_t ref, const int8_t material) = 0;



};

} // namespace

#endif

