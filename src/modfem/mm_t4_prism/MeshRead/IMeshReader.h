// ***************************************************************
//  IMeshReader   version:  1.0   �  date: 05/27/2008
//                          1.1            11/09/2009
//  -------------------------------------------------------------
//  Kazimierz Michalik
//  -------------------------------------------------------------
//  Copyright (C) 2008 - All Rights Reserved
// ***************************************************************
//
// ***************************************************************
#ifndef IMeshReader_h__
#define IMeshReader_h__

#include <string>
#include "../Common.h"

class hHybridMesh;

namespace MeshRead{

class IMeshReader
{
public:

	IMeshReader();
	virtual ~IMeshReader();

	virtual bool	doRead(hHybridMesh * mesh)=0;


	static const char V_COUNT = 'v', VERTICES = 'V',
	E_COUNT = 'e', EDGES = 'E',
	F_COUNT ='f', FACES ='F',
	EL_COUNT = 'l', ELEMENTS = 'L',
	AUX='a', END = 0;

	/** In readSequence is stored sequendce of actions used to correctly read input.
	*/
	char*	readSequence;

	/** /return true if MeshBuilder initialized successfully, otherwise returns false.
	*/
	virtual bool Init() = 0;
	virtual bool Init(const std::string & name) = 0;
	virtual std::string Name() const = 0;
	/** Frees resources used by builder.
	*/
	virtual void Free() = 0;

	/** returns number of mesh dimensions: 1D/2D/3D
	*/
	virtual int GetCoordinatesDimension() const = 0;

	/** /return no. of vertices in initial mesh.
	*/
	virtual int GetVerticesCount() = 0;

	/**	/param  coords array for coordinates of vertex
		/return true if there is next vertex and false if there isn't.
	*/
	virtual bool GetNextVertex(double coords[]) = 0;

	/** /param element_type no.of vertices in element
        /return no. of elements in inintial mesh
	*/
	virtual void GetElementCount(int type_count[]) = 0;

	/** /param  vertices OUT array for current element vertices numbers
		/param  neighbours OUT array for current element neighoubrs(other elements) numbers.
		/param	element_type OUT type of element (no. of vertices)
		/return true if there is next element and false if this is last element.
	*/
	virtual bool GetNextElement(Tind vertices[], Tind neighbours[], Tind faces[], Tind & element_type , Tind & father, int8_t  & material, int8_t & ref) = 0;

	/** /param bc array for boundary condition parameters
		bc[0] - Dirichlet BC
		bc[1] - Neumann BC
		bc[2] - Cauchy BC par1
		bc[3] - Cauchy BC par2
	*/
	virtual bool GetBoundaryConditions(double ** bc, int & bcCount) = 0;

	/**
	*/
	virtual int GetEdgesCount()
	{
	    return 0;
	}

	/** /param edge_type OUT type of edge (no. of vertices)
        /param vertices OUT
	*/
	virtual bool GetNextEdge(Tind vertices[],Tind & edge_type)
	{
	    throw("You shouldn't use this function!");
	    return false;
	}

	virtual void GetFacesCount(int face_type[])
	{
	}

	virtual bool GetNextFace(Tind edges[],Tind & face_type, int8_t & bc, Tind neigh[])
	{
	    return false;
	}
	
	virtual bool needsSetup() const {return true;}
};

}//!namespace hpFEm
#endif // IMeshReader_h__
