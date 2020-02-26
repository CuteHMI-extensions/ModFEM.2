
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "NasFileImporter.h"
#include "../MeshModule/hHybridMesh.h"
#include "../MeshModule/Face3.h"
#include "../MeshModule/Face4.h"
#include "../MeshModule/ElemPrism.h"
#include "../MeshModule/ElemT4.h"

namespace MeshRead
{

NasFileImporter::NasFileImporter()
	: MeshFileImporter(), _face_count(0), _face_readed(0), _edge_count(0), 
	  _edge_readed(0),_quads(0),_prisms(0)
{
  /*	readSequence = new char[7];
	readSequence[0] = V_COUNT;
	readSequence[1] = VERTICES;
	readSequence[2] = F_COUNT;
	readSequence[3] = FACES;
	readSequence[4] = EL_COUNT;
	readSequence[5] = ELEMENTS;
	readSequence[6] = END;
*/
}

NasFileImporter::NasFileImporter(const std::string & file_name)
	: MeshFileImporter(file_name), _face_count(0), 
	  _face_readed(0), _edge_count(0), _edge_readed(0),_quads(0),_prisms(0)
{
  /*readSequence = new char[7];
	readSequence[0] = V_COUNT;
	readSequence[1] = VERTICES;
	readSequence[2] = F_COUNT;
	readSequence[3] = FACES;
	readSequence[4] = EL_COUNT;
	readSequence[5] = ELEMENTS;
	readSequence[6] = END;
*/
	_file_name = file_name;
}

/** /return true if MeshBuilder initialized successfully, otherwise returns false.
*/
bool NasFileImporter::Init()
{
	MeshFileImporter::Free();

	if(MeshFileImporter::Init())
	{
		char tmp[256]={};
		std::string s_tmp;
		//std::vector<double> vec(1024);
		// ignore everything before grid data
		gridPos_= 0;
		elemPos_=0;
		facePos_=0;
		//size_t readed_no = 0;
		_file.clear();
		_file >> s_tmp;
		while(_file.good()) {
		    
		    if(s_tmp.compare("GRID")==0) {
			++_vert_count;
		    }
		    else if(s_tmp.compare("CTRIAX")==0) {
			if(_face_count==0) {
			    facePos_ = getPos(); //static_cast<int>(_file.tellg());
			    facePos_-= static_cast<int>(s_tmp.length());
			}
			++_face_count;
		    }
		    else if(s_tmp.compare("CQUADX")==0) {
			if(_face_count==0 && _quads==0) {
			    facePos_ = getPos(); //static_cast<int>(_file.tellg());
			    facePos_-= static_cast<int>(s_tmp.length());
			}
			++_quads;
		    }
		    else if(s_tmp.compare("CTETRA")==0) {
			if(_prisms==0 && _elem_count==0) {
			    elemPos_ = getPos(); //static_cast<int>(_file.tellg());
			    elemPos_-= static_cast<int>(s_tmp.length());
			}
			++_elem_count;
		    }
		    else if(s_tmp.compare("CPRIZM")==0) {
			if(_prisms==0 && _elem_count==0) {
			    elemPos_ = getPos(); //static_cast<int>(_file.tellg());
			    elemPos_-= static_cast<int>(s_tmp.length());
			}
			++_prisms;
		    }
		    _file.getline(tmp,255);
		    _file >> s_tmp;
		}
		_file.clear();
		_file.seekg(gridPos_);
	}
	return _file.good();
}

bool    NasFileImporter::Init(const std::string & name)
{
	_file_name = name;
	return MeshFileImporter::Init();
}


bool	NasFileImporter::doRead(hHybridMesh * mesh)
{
	if(_file.good())
	{
		// GRID

		mmv_out_stream<< "\n Reading vertices...";

		_file.seekg(gridPos_);
		mesh->vertices_.reserve(_vert_count);
        Vertex * vPtr(NULL);
		std::string	s_tmp;
		while(_vert_readed < _vert_count && _file.good())
		{
			double	val(0);
			//char tmp[12]={0};
			_file >> s_tmp;
            vPtr = mesh->vertices_.newObj<Vertex>(mesh,NULL);
			_file >> val >> vPtr->coords_[0] >> vPtr->coords_[1] >> vPtr->coords_[2] >> val;
			++_vert_readed;
		}
		assert(_vert_count == _vert_readed);

		mmv_out_stream << " Vertices readed: " << _vert_readed << "\n Reading faces... ";

		// FACES
		int quad_readed(0);
		_file.seekg(0/*facePos_*/);
		mesh->faces_.reserve(_face_count+_quads, sizeof(Face3) * _face_count+sizeof(Face4)*_quads);
		mesh->edges_.reserve(_elem_count*6+_prisms*9,sizeof(Edge) * (_elem_count*6+_prisms*9));
        hObj * fPtr(NULL);
		uTind	vertices[6]={UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN};
		while(_face_readed+quad_readed < (_face_count+_quads) && _file.good())
		{
			double	val(0);
			//char tmp[12]={0};
			//_file.getline(tmp,11,' ');
			//tmp[6]=0;
			_file >> s_tmp;
			int face_type(0);
			if(0 == s_tmp.compare("CTRIAX")) //if(0 == strcmp("CTRIAX",tmp ))
			{
				face_type = 3;
				if(face_type < 2 && face_type > 4)
				{
					throw "\nNasFileImporter::GetNextFace: wrong face_type.\n Make sure that file has platform-specyfic file endings (Win/Linux).";
				}
				_file >> val;	// number - ignore it
				int bc;
				_file >> bc;  // element group no. = bc no.
				_file >> vertices[0];// fPtr->verts(0) = tmp2; // vertices no. !!! (not edges)
				_file >> vertices[1]; //fPtr->verts(1) = tmp2; 
				_file >> vertices[2]; //fPtr->verts(2) = tmp2; 
				std::sort(vertices, vertices+3);

                fPtr = mesh->faces_.newObj<Face3>(mesh,vertices);
				fPtr->flags(B_COND)= bc;
				fPtr->components(0) = mesh->edge(vertices[0],vertices[1]).id_;
				fPtr->components(1) = mesh->edge(vertices[0],vertices[2]).id_;
				fPtr->components(2) = mesh->edge(vertices[1],vertices[2]).id_;
				fPtr->neighs(0) = UNKNOWN;
				fPtr->neighs(1) = UNKNOWN;
				++_face_readed;
			}
			else if(0 == s_tmp.compare("CQUADX")) // quads
			{
				face_type = 4;
				if(face_type < 2 && face_type > 4)
				{
					throw "\nNasFileImporter::GetNextFace: wrong face_type.\n Make sure that file has platform-specyfic file endings (Win/Linux).";
				}
				_file >> val;	// number - ignore it
				int bc;
				_file >> bc;  // element group no. = bc no.
				_file >> vertices[0];// fPtr->verts(0) = tmp2; // vertices no. !!! (not edges)
				_file >> vertices[1]; //fPtr->verts(1) = tmp2; 
				_file >> vertices[3]; //fPtr->verts(2) = tmp2; 
				_file >> vertices[2];
				//std::sort(vertices, vertices+3);

                fPtr = mesh->faces_.newObj<Face4>(mesh,vertices);
				fPtr->flags(B_COND)= bc;
				fPtr->components(0) = mesh->edge(vertices[0],vertices[1]).id_;
				fPtr->components(1) = mesh->edge(vertices[0],vertices[2]).id_;
				fPtr->components(2) = mesh->edge(vertices[1],vertices[3]).id_;
				fPtr->components(3) = mesh->edge(vertices[2],vertices[3]).id_;
				fPtr->neighs(0) = UNKNOWN;
				fPtr->neighs(1) = UNKNOWN;
				++quad_readed;
			}
		}
		assert(_face_count ==_face_readed);
		assert(_quads == quad_readed);
		//assert(mesh->checkUniqueness());


		mmv_out_stream << " Faces readed: " << _face_readed << "\n Reading elements... ";

		// 3D ELEMENTS
		int prism_readed(0);
		_file.seekg(0/*elemPos_*/);
		mesh->elements_.reserve(_elem_count+_prisms, sizeof(ElemT4) * _elem_count+sizeof(ElemPrism)*_prisms);
		mesh->faces_.reserve(_elem_count*4+_prisms*5, sizeof(Face3) * (_elem_count*4 + _prisms*2) + sizeof(Face4)*_prisms*3);	// allocate space for internal spaces
        ElemT4 * el(NULL);
        ElemPrism * pr(NULL);
		//int elem_type(7);
		//eFaceFlag f3Orient(F_OUT);
        //hObj *f0(NULL);
        //hObj *f1(NULL);
        //hObj *f2(NULL);
        //hObj* f3(NULL);
		while(_elem_readed+prism_readed < (_elem_count+_prisms) && _file.good())
		{
			double	val(0);
			//char tmp[12]={0};
			//_file.getline(tmp,11);
			//tmp[6]=0;
		
			_file >> s_tmp;
			if(0 == s_tmp.compare("CTETRA"))
			{
				_file >> val;	// number - ignore it
				int tmp2(0);
                _file >> tmp2; // element group no. = material
				_file >> vertices[0]; // vertices no. !!! (not edges)
				_file >> vertices[1]; 
				_file >> vertices[2]; 
				_file >> vertices[3]; 
				std::sort(vertices,vertices+4);
				//type=mesh->normalize(vertices);

                el = mesh->elements_.newObj<ElemT4>(mesh,vertices);
				// below is auto-done in constructor				
				//				el->components(3) = (f3=&mesh->face(vertices[1],vertices[2],vertices[3]))->id_;
				//				el->components(0) = (f0=&mesh->face(vertices[0],vertices[1],vertices[2]))->id_;
				//				el->components(1) = (f1=&mesh->face(vertices[0],vertices[1],vertices[3]))->id_;
				//				el->components(2) = (f2=&mesh->face(vertices[0],vertices[2],vertices[3]))->id_;
                el->flags(GROUP_ID) = tmp2;
				el->flags(REFINEMENT) = UNKNOWN;
				el->parent_=0; //NO_FATHER;
				  
				++_elem_readed;
			}
			else if(0 == s_tmp.compare("CPRIZM"))
			{
				_file >> val;	// number - ignore it
				int tmp2(0);
				_file >> tmp2; // element group no. = bc no.
				_file >> vertices[3]; // vertices no. !!! (not edges)
				_file >> vertices[4]; 
				_file >> vertices[5]; 
				_file >> vertices[0]; 
				_file >> vertices[1];
				_file >> vertices[2];
				//std::sort(vertices,vertices+4);
				//type=mesh->normalize(vertices);

                pr = mesh->elements_.newObj<ElemPrism>(mesh,vertices);
				
			
				//pr->components(0) = mesh->face(vertices[0],vertices[1],vertices[2]).id_;
				//pr->components(1) = mesh->face(vertices[3],vertices[4],vertices[5]).id_;
				//pr->components(2) = mesh->face(vertices[0],vertices[1],vertices[3],vertices[4]).id_;
				//pr->components(3) = mesh->face(vertices[0],vertices[2],vertices[3],vertices[5]).id_;
				//pr->components(4) = mesh->face(vertices[1],vertices[2],vertices[4],vertices[5]).id_;
				
				pr->flags(GROUP_ID) = 1;
				pr->flags(REFINEMENT) = UNKNOWN;
				pr->parent_=0; // NO_FATHER
				  
				++prism_readed;
				
			}

		}
		assert(_elem_count== _elem_readed);
		assert(_prisms == prism_readed);
		//assert(mesh->checkUniqueness());
		mmv_out_stream <<  " Elements readed: " << _elem_readed << ".\n Checking face neighbours.. ";
		// EDGES
		
		
		// checking faces neighbours
		{
			Memory::Field<2,ID> * faceNeighs( new Memory::Field<2,ID>[mesh->faces_.last()+1] );
				
			for(hHybridMesh::ElemPool::Iterator<hHybridMesh::ElemPool::allObj> it(& mesh->elements_); 
				!it.done(); ++it)
			{
				for(uTind f(0); f < it->typeSpecyfic_.nFaces_; ++f)
				{
					if(faceNeighs[hObj::posFromId(it->components(f))][0]==0)
					{
						faceNeighs[hObj::posFromId(it->components(f))][0]=it->id_;
					}
					else
					{
						faceNeighs[hObj::posFromId(it->components(f))][1]=it->id_;
					}
				}
			}


			hObj * fPtr(& mesh->faces_.at(FIRST));
			for(unsigned int i(FIRST); i <= mesh->faces_.last(); ++i)
			{
				if(fPtr->neighs(0) != faceNeighs[i][0])
				{
					fPtr->neighs(0) = faceNeighs[i][0];
				}
				if(fPtr->neighs(1) != faceNeighs[i][1])
				{
					fPtr->neighs(1) = faceNeighs[i][1];
				}
				assert(fPtr->neighs(0) != fPtr->neighs(1));
				fPtr=fPtr->next();
			}
				
			delete [] faceNeighs;
		}
		// end of checking faces neighs

		mmv_out_stream << " OK!\n";
	}
	else
	{
		throw "Error while reading";
	}
	return true;
}

/** returns number of mesh dimensions: 1D/2D/3D
*/
int  NasFileImporter::GetCoordinatesDimension() const
{
	return 0;
}

/** /return no. of vertices in initial mesh.
*/
int  NasFileImporter::GetVerticesCount()
{
	return _vert_count;
}

/**	/param  coords array for coordinates of vertex
	/return true if there is next vertex and false if there isn't.
*/
bool  NasFileImporter::GetNextVertex(double coords[])
{
	if(_file.good())
	{
		double	val(0);
		char tmp[12]={};
		_file.read(tmp,11);
		_file >> val >> coords[0] >> coords[1] >> coords[2] >> val;
	}
	return _file.good();
}

void     NasFileImporter::GetFacesCount(int type_count[])
{
	type_count[3]=_face_count;
}

bool    NasFileImporter::GetNextFace(Tind edges[],Tind & face_type, int8_t & bc, Tind neigh[])
{
	const bool nextExist((_face_readed < _face_count) && _file.good());
	//std::streampos	pos = _file.tellg();
	if(nextExist)
	{
        assert(edges != NULL);
        assert(neigh != NULL);
		
		double	val(0);
		char tmp[12]={};
		_file.read(tmp,11);
		_file >> val;
		int tmp2;
		_file >> tmp2; bc=tmp2;
		_file >> edges[0] >> edges[1] >> edges[2];
		
		face_type = 3;
		if(face_type < 2 && face_type > 4)
		{
			throw "\nNasFileImporter::GetNextFace: wrong face_type.\n Make sure that file has platform-specyfic file endings (Win/Linux).";
		}
		if(face_type == 4)
		{
			_file >> edges[3];
		}
		// check negative values
		for(Tind i(0); i < face_type; ++i)
		{
			 --edges[i]; // are numerated from 1, not 0!!
		}
		neigh[1] = -1;

		++_face_readed;
	}
	return nextExist;
}

/** /return no. of elements in inintial mesh
*/
void  NasFileImporter::GetElementCount(int type_count[])
{
	type_count[7]=_elem_count;
}

/** /param  vertices array for current element vertices numbers
	/param  neighbours array for current element neighoubrs(other elements) numbers.
	/param	bc number of boundary condition at this element
	/return true if there is next element and false if this is last element.
*/
bool    NasFileImporter::GetNextElement(Tind vertices[], Tind neighbours[], Tind faces[], Tind & element_type , Tind & father, int8_t & material, int8_t & ref)
{
	const bool nextExist(_elem_readed < _elem_count);
	if(nextExist)
	{
	  //double	val(0);
		char tmpBuf[12]={};
		_file.read(tmpBuf,11);
		int tmp;
		_file >> tmp;
	   element_type = 7;
		assert(element_type > 4 && element_type < 10);
		int no_faces(4);
		_file >> faces[0]>> faces[1]>> faces[2]>> faces[3];

		if(element_type == 5)   // prism
		{
			no_faces=5;
			_file >> faces[4];
		}

		for(int i(0); i < no_faces; ++i)
		{
			if(faces[i] < 0)
			{
				faces[i]*=-1;
			}
			--faces[i];
		}
		++_elem_readed;
	}
	return nextExist;
}

/** /param bc array for boundary condition parameters
	bc[0] - Dirichlet BC
	bc[1] - Neumann BC
	bc[2] - Cauchy BC par1
	bc[3] - Cauchy BC par2
*/
bool  NasFileImporter::GetBoundaryConditions(double ** bc, int & bcCount)
{
	return false;
}

} // namespace MeshRead
