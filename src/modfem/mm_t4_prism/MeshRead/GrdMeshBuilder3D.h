// ***************************************************************
//  GrdMeshBuilder3D   version:  1.0   ·  date: 05/27/2008
//  -------------------------------------------------------------
//  Kazimierz Michalik
//  -------------------------------------------------------------
//  Copyright (C) 2008 - All Rights Reserved
// ***************************************************************
//
// ***************************************************************
#ifndef GrdMeshBuilder3D_h__
#define GrdMeshBuilder3D_h__

#include "IMeshBuilder.h"

#include <fstream>
#include <string>

class GrdMeshBuilder3D :
	public MeshRead::IMeshReader
{
public:
	GrdMeshBuilder3D( const char  grdFileName[] );
	~GrdMeshBuilder3D(void);

	/** /return true if MeshBuilder initialized successfully, otherwise returns false.
	*/
	bool Init() ;

	/** Frees resources used by builder.
	*/
	virtual void Free();

	/** /return returns number of mesh dimensions: 1D/2D/3D
	*/
	 int GetCoordinatesDimension() const ;

	/** /return no. of vertices in initial mesh.
	*/
	 int GetVerticesCount() ;

	/**	/param  coords array for coordinates of vertex
		/return true if there is next vertex and false if there isn't.
	*/
	 bool GetNextVertex(double coords[]) ;

	/** /return no. of elements in inintial mesh.
	*/
	 int GetElementCount() ;

	/** /param  vertices array for current element vertices numbers
		/param  neighbours array for current element neighoubrs(other elements) numbers.
		/return true if there is next element and false if this is last element.
	*/
	 bool GetNextElement(int vertices[], int neighbours[]);

	 /** /param bc array for boundary condition parameters
		 bc[0] - Dirichlet BC
		 bc[1] - Neumann BC
		 bc[2] - Cauchy BC par1
		 bc[3] - Cauchy BC par2
	 */
	 bool GetBoundaryConditions(double ** bc,int & bcCount);

protected:
	int	getDim(const int n) const;

	std::ifstream	grid_file;
	std::string	    fileName;

	int	dimension;

	int	verticesCount;
	int	readedVertices;

	int	elementsCount;
	int	readedElements;
};
#endif // GrdMeshBuilder3D_h__
