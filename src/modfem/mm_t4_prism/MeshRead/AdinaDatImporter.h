#ifndef ADINA_DAT_IMPORTER_H
#define ADINA_DAT_IMPORTER_H

#include "MeshFileImporter.h"

namespace MeshRead
{

class AdinaDatImporter :
	public MeshRead::MeshFileImporter
{
public:
	AdinaDatImporter(const std::string & file_name);
	~AdinaDatImporter(void);

	/** /return true if MeshBuilder initialized successfully, otherwise returns false.
	*/
	bool Init(const std::string & file_name) ;

	/** Frees resources used by builder.
	*/
	//void Free();

	/** returns number of mesh dimensions: 1D/2D/3D
	*/
	int GetCoordinatesDimension() const ;

	/** /return no. of vertices in initial mesh.
	*/
	int GetVerticesCount() ;

	/**	/param  coords array for coordinates of vertex
	/return true if there is next vertex and false if there isn't.
	*/
	bool GetNextVertex(double coords[]) ;

	/** /return no. of elements in inintial mesh
	*/
	int GetElementCount() ;

	/** /param  vertices array for current element vertices numbers
	/param  neighbours array for current element neighoubrs(other elements) numbers.
	/param	bc number of boundary condition at this element
	/return true if there is next element and false if this is last element.
	*/
	bool GetNextElement(int vertices[], int neighbours[]) ;

	/** /param bc array for boundary condition parameters
	bc[0] - Dirichlet BC
	bc[1] - Neumann BC
	bc[2] - Cauchy BC par1
	bc[3] - Cauchy BC par2
	*/
	bool GetBoundaryConditions(double ** bc, int & bcCount) ;
};

}// namespace

#endif //ADINA_DAT_IMPORTER_H
