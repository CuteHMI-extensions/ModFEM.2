#pragma once

#include "IMeshWriter.h"

#include <string>
#include <fstream>

class AMVisualExporter :
	public MeshWrite::IMeshWriter
{
public:
	AMVisualExporter(const std::string & file_name);
	~AMVisualExporter(void);

	/** Initializes exporter. After this function call exporter should be ready for exporting the mesh.
		/return true if exporter initializes correctly, otherwise false.
				If this function doesn't succeeded exporting process is gently stopped and error is thrown.
	*/
	bool Init() ;

	/** Frees all resources allocated by class during mesh export.
	*/
	void Free() ;

	/** /param noVert Number of all vertices (global) in mesh.
	*/
	void WriteVerticesCount(const int noVert) ;

	/** /param coords array of vertex coordinates
	*/
	void WriteVertex(double coords[])  ;

	/** /param noElems Number of all elements in mesh.
	*/
	void WriteElementCount(const int noElems) ;

	/** /param vertices array of indexes of vertices making this element
		/param neighbours array of indexes of neighbour elements
	*/
	void WriteElement(int vertices[], int neighbours[]) ;

protected:

	std::string	_file_name;
	std::fstream	_file;
};
