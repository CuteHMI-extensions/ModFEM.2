// ***************************************************************
//  GrdMeshBuilder3D   version:  1.0   ·  date: 05/27/2008
//  -------------------------------------------------------------
//  Kazimierz Michalik
//  -------------------------------------------------------------
//  Copyright (C) 2008 - All Rights Reserved
// ***************************************************************
//
// ***************************************************************

#include "GrdMeshBuilder3D.h"
#include <stdlib.h>
#include <string>
#include <iostream>
#include <stdexcept>

GrdMeshBuilder3D::GrdMeshBuilder3D(const char  grdFileName[] )
	: fileName(grdFileName), dimension(-1), verticesCount(0), readedVertices(0) {}

GrdMeshBuilder3D::~GrdMeshBuilder3D() {
	Free();
}

int	GrdMeshBuilder3D::getDim(const int n) const {
	if(n <= 0)
		return -1;
	else if(n < 20)
		return 2;
	else if(n < 100)
		return 3;
	else
		return 1;
}

bool	GrdMeshBuilder3D::Init() {
	readedElements	= 0;
	readedVertices	= 0;
	verticesCount	= 0;
	elementsCount	= 0;

	grid_file.open(fileName.c_str());
	if(grid_file.rdbuf()->is_open()) {
		using namespace std;
		char tmp[96];

		string str, pattern = "GRID_TYPE";
		grid_file >> str >> tmp;
		const int gtype = atoi(tmp);
		if(grid_file.fail() || gtype < 0 || (str != pattern)) {
			out_stream << "Error: bad grid_file format [GRID_TYPE]!\n";
			return false;
		}

		dimension = getDim(gtype);	// geometry dimension
		if(dimension != 3)
			return false;

		pattern = "DIMENSIONS";	// physical dimensions of domain
		grid_file >> str;
		if(grid_file.fail() || (str != pattern)) {
			out_stream << "Error: bad grid_file format [DIMENSIONS]!\n";
			return false;
		}
		grid_file >> str;
		grid_file >> str;
		grid_file >> str;

		pattern = "GRID_LEVELS";
		grid_file >> str;
		if(grid_file.fail() || (str != pattern)) {
			out_stream << "Error: bad grid_file format [LEVELS]!\n";
			return false;
		}
		grid_file >> str >> str;		// GRID LEVELS
	} else
		throw	std::runtime_error("Unable to open file.");
	return true;
}

void	GrdMeshBuilder3D::Free(){
	grid_file.close();
	readedVertices	= 0;
	verticesCount	= 0;
	readedElements	= 0;
	elementsCount	= 0;
	dimension		= 0;
}

int	GrdMeshBuilder3D::GetCoordinatesDimension() const {
	return dimension;
}

int	GrdMeshBuilder3D::GetVerticesCount() {
	if(readedVertices == 0){
		using namespace std;
		string	str, pattern = "POINTS";
		char	tmp[64];

		grid_file >> str >> tmp;
		verticesCount = atol(tmp);
		if(grid_file.fail() || verticesCount < 3 || (str != pattern)) {
			cout << "Error: bad grid_file format [POINTS]!\n";
			return -1;
		}
	}
	return verticesCount;
}

bool	GrdMeshBuilder3D::GetNextVertex(double coord[]){
	if(coord == 0)
		throw std::runtime_error("GetNextVertex: coord param is NULL pointer");
	if(verticesCount > readedVertices) {

		char	x[32], y[32], z[32];

		grid_file >> x;
		coord[0] = atof(x);

		grid_file >> y;
		coord[1] = atof(y);

		grid_file >> z;
		coord[2] = atof(z);

		++readedVertices;
		return true;
	}
	return false;
}

int	GrdMeshBuilder3D::GetElementCount() {
	if(readedElements == 0 ){
		std::string str, pattern = "ELEMENTS";
		char	tmp[64];

		grid_file >> str >> tmp;
		elementsCount =atol(tmp);
		grid_file >> tmp;
		if(grid_file.fail() || elementsCount < 1 || (str != pattern)) {
			out_stream << "Error: bad grid_file format [ELEMENTS]!\n";
			return -1;
		}
	}
	return elementsCount;
}

bool	 GrdMeshBuilder3D::GetNextElement(int vertices[], int neighbours[]){
	if(vertices == 0)
		throw std::runtime_error("GetNextElement: vertices param is NULL pointer");
	if(neighbours == 0)
		throw std::runtime_error("GetNextElement: neighbours param is NULL pointer");
	if( elementsCount > readedElements ) {
		char	tmp[64];
		grid_file >> tmp;		vertices[0] = atol(tmp);	// reading vertex ids
		grid_file >> tmp;		vertices[1] = atol(tmp);
		grid_file >> tmp;		vertices[2] = atol(tmp);
		grid_file >> tmp;		vertices[3] = atol(tmp);

		grid_file >> tmp;		neighbours[0] = atol(tmp);	// reading neighbor ids
		grid_file >> tmp;		neighbours[1] = atol(tmp);
		grid_file >> tmp;		neighbours[2] = atol(tmp);
		grid_file >> tmp;		neighbours[3] = atol(tmp);

		++readedElements;
		return true;
	}
	return false;
}

bool GrdMeshBuilder3D::GetBoundaryConditions(double ** bc, int & bcCount){

	if( elementsCount == readedElements ){
		char	tmp[64];
		grid_file >> tmp;		// BOUNDARIES
		grid_file >> tmp;		// 3

		const int newBcCount = atol(tmp)+1;
		if((bc == NULL) || (bcCount < newBcCount)) {
			// TODO: boudary conditions array
			//SAFE_DELETE_ARRAY(bc);
			//bcCount = newBcCount;
			//bc = new hpFEM::BC[bcCount];
		}

//		for(int i(1); i < bcCount; ++i) {
//			grid_file >> tmp;		// 1
//			bc[i].type = hpFEM::BC_TYPE(_wtol(tmp));
//
//			if(bc[i].type != hpFEM::BC_TYPE::BC_NON){
//				grid_file >> tmp;
//				bc[i].val[0] = atol(tmp);	// reading bc val
//				if(bc[i].type == hpFEM::BC_TYPE::BC_CAU) {
//					grid_file >> tmp;
//					bc[i].val[1] = atol(tmp);	// reading bc val
//				}
//			}
//		}
		return true;
	}
	return false;
}
