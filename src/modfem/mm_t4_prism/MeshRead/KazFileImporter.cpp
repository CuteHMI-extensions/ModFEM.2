
#include <cassert>

#include "KazFileImporter.h"
#include "../MeshModule/hHybridMesh.h"
#include "../MeshModule/hObj.h"
#include "../MeshModule/Face3.h"
#include "../MeshModule/Face4.h"
#include "../MeshModule/ElemPrism.h"
#include "../MeshModule/ElemT4.h"

#include <algorithm>
#include <functional>

using MeshRead::KazFileImporter;
using MeshRead::MeshFileImporter;

KazFileImporter::KazFileImporter(const std::string & fileName)
: MeshFileImporter(fileName), _face_count(0), _face_readed(0), _edge_count(0), _edge_readed(0)
{
	readSequence = new char[9];
	readSequence[0] = V_COUNT;
	readSequence[1] = VERTICES;
	readSequence[2] = E_COUNT;
	readSequence[3] = EDGES;
	readSequence[4] = F_COUNT;
	readSequence[5] = FACES;
	readSequence[6] = EL_COUNT;
	readSequence[7] = ELEMENTS;
	readSequence[8] = END;
	//readSequence = "vVeEfFlL";
}

KazFileImporter::KazFileImporter()
: MeshFileImporter(), _face_count(0), _face_readed(0), _edge_count(0), _edge_readed(0)
{
	readSequence = new char[9];
	readSequence[0] = V_COUNT;
	readSequence[1] = VERTICES;
	readSequence[2] = E_COUNT;
	readSequence[3] = EDGES;
	readSequence[4] = F_COUNT;
	readSequence[5] = FACES;
	readSequence[6] = EL_COUNT;
	readSequence[7] = ELEMENTS;
	readSequence[8] = END;
	//readSequence = "vVeEfFlL";
}

KazFileImporter::~KazFileImporter()
{
	//dtor
}

bool    KazFileImporter::Init()
{
	MeshFileImporter::Free();

	if(MeshFileImporter::Init())
	{
		int tmp,next_index(0);
		_file >> tmp >> tmp >> tmp >> tmp >> _vert_count >> next_index;
	}

	return _file.is_open();
}

bool	KazFileImporter::doRead(hHybridMesh * mesh)
{
	mesh->wasNormalized_=true;
	// vertices count
	mesh->vertices_.reserve(_vert_count);
	// vertices
    Vertex * v( NULL );
	while(_vert_readed < _vert_count)
	{
        v = mesh->vertices_.newObj<Vertex>(mesh,NULL);
		_file >> v->coords_[0]; _file >> v->coords_[1]; _file >> v->coords_[2];
		++_vert_readed;
	}

	// edges count
	int tmp;
	_file >> _edge_count >>  tmp;
	mesh->edges_.reserve(_edge_count, _edge_count * sizeof(Edge));
	// edges
	  int edge_type,pos,parent,level;
	uTind vertices[6]={UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN};
    Edge * e(NULL);
	while(_edge_readed < _edge_count)
	{
	  _file >>  edge_type
			>> pos >> parent >> level; 
		_file >> vertices[0] >> vertices[1];
        e = mesh->edges_.newObj<Edge>(mesh,vertices);
		  e->pos_ = pos;
		  e->type_ = edge_type;
		  e->parent_ = parent;
		  e->level_ = level;
		e->components(0) = hObj::makeId(vertices[0],Vertex::myType);
		e->components(1) = hObj::makeId(vertices[1],Vertex::myType);
		++_edge_readed;
	}
	//assert(mesh->checkUniqueness());

	// faces count
	int face_type(0);
	if(_face_count == 0)
	{
		_file >> _face_count >> tmp;   // all - faces count
	}
	std::streampos position(getPos());
	int type_count[10]={0};
	for(int i(0); i < _face_count; ++i)
	{
	  _file >> face_type
			>> pos
			>> parent
			>> level;
		assert(face_type > 2 && face_type < 5);
		++type_count[face_type];
		_file.ignore(std::numeric_limits< int >::max(), '\n' );
		_file.ignore(std::numeric_limits< int >::max(), '\n' );  
	}
	_file.seekg(position);

	// faces
    //std::streampos	posit = getPos();
    hObj * f(NULL);
	int	edges[5]={0};
	int  typeFlag=0,params[3]={0};
      //Edge * edgePtrs[4]={NULL};
	mesh->faces_.reserve(type_count[3]+type_count[4], sizeof(Face4)*type_count[4] + sizeof(Face3)*type_count[3]);
	while((_face_readed < _face_count) && _file.good())
	{
	  _file >> face_type
			>> pos
			>> parent
			>> level;
		if(face_type < 2 && face_type > 4)
		{
			throw "\nKazFileImporter::GetNextFace: wrong face_type.\n Make sure that file has platform-specyfic file endings (Win/Linux).";
		}
		_file >> tmp; 
		_file >> params[0] >> params[1] >> typeFlag;

		// check negative values
		for(Tind i(0); i < face_type; ++i)
		{
			_file >> edges[i];
			  assert(edges[i] >= static_cast<Tind>(FIRST));
		}
		
		for(Tind i(0); i < face_type; ++i) {
		    _file >> vertices[i];
		}

		switch(face_type)
		{ 
		case Face3::myType:{
		vertices[3]=UNKNOWN;
        f = mesh->faces_.newObj<Face3>(mesh,vertices) ;
		  f->pos_ = pos;
		  f->components(0) = hObj::makeId(edges[0],Edge::myType);
		f->components(1) = hObj::makeId(edges[1],Edge::myType);
		f->components(2) = hObj::makeId(edges[2],Edge::myType);
		assert(f->components(0)==mesh->edge(f->verts(0),f->verts(1)).id_);
		assert(f->components(1)==mesh->edge(f->verts(0),f->verts(2)).id_);
		assert(f->components(2)==mesh->edge(f->verts(1),f->verts(2)).id_);
		} break;
		case Face4::myType:	{
        f = mesh->faces_.newObj<Face4>(mesh,vertices);
		  f->pos_ = pos;
		  f->components(0) = hObj::makeId(edges[0],Edge::myType);
		f->components(1) = hObj::makeId(edges[1],Edge::myType);
		f->components(2) = hObj::makeId(edges[2],Edge::myType);
		f->components(3) = hObj::makeId(edges[3],Edge::myType);
		assert(f->components(0)==mesh->edge(f->verts(0),f->verts(1)).id_);
		assert(f->components(1)==mesh->edge(f->verts(0),f->verts(2)).id_);
		assert(f->components(2)==mesh->edge(f->verts(1),f->verts(3)).id_);
		assert(f->components(3)==mesh->edge(f->verts(2),f->verts(3)).id_);
		} break;
		default: throw "Unknown face type!";
		}
		f->type_ = face_type;
		  f->parent_ = parent;
		  f->level_ = level;
		  f->flags(B_COND)=tmp;
		f->flags(F_TYPE)=typeFlag;
		f->neighs(0) = params[0];
		f->neighs(1) = params[1];
		// f->test();// neigs are not properly  set, so it will crash
		++_face_readed;
	}
	//assert(mesh->checkUniqueness());

	// element count
	int type(0);
	if(_elem_count == 0)
	{
		_file >> _elem_count >> tmp;
	}
	position = getPos();
	for(int i(0); i < _elem_count; ++i)
	{
		_file >> type >> tmp >> tmp >> tmp;
		++type_count[type];
		_file.ignore(std::numeric_limits< int >::max(), '\n' );
		_file.ignore(std::numeric_limits< int >::max(), '\n' );
	}
	_file.seekg(position);

	mesh->elements_.reserve(type_count[ElemT4::myType] + type_count[ElemPrism::myType],
		type_count[ElemT4::myType]*sizeof(ElemT4) + type_count[ElemPrism::myType]*sizeof(ElemPrism));

	// elements
	int element_type(0);
    hObj * elem(NULL);
	while(_elem_readed < _elem_count)
	{
		int mt,rt,father,faces[5];
		_file >> element_type;
		_file >> mt >> father >> rt >> typeFlag;

		for(unsigned int i(0); i < hObj::eTable[element_type].nComponents_; ++i)
		{
			_file >> faces[i];
		}

		for(unsigned int i(0); i < hObj::eTable[element_type].nVerts_; ++i)
		{
			_file >> vertices[i];
		}

		switch(element_type)
		{
		case ElemT4::myType:
		{
            elem = mesh->elements_.newObj<ElemT4>(mesh,vertices);
		    elem->components(0) = hObj::makeId(faces[0],Face3::myType);
		    elem->components(1) = hObj::makeId(faces[1],Face3::myType);
		    elem->components(2) = hObj::makeId(faces[2],Face3::myType);
		    elem->components(3) = hObj::makeId(faces[3],Face3::myType);
		    
		} break;
		case ElemPrism::myType: {
            elem = mesh->elements_.newObj<ElemPrism>(mesh,vertices);
		    elem->components(0) = hObj::makeId(faces[0],Face3::myType);
		    elem->components(1) = hObj::makeId(faces[1],Face3::myType);
		    elem->components(2) = hObj::makeId(faces[2],Face4::myType);
		    elem->components(3) = hObj::makeId(faces[3],Face4::myType); // CHANGE!!
		    elem->components(4) = hObj::makeId(faces[4],Face4::myType); // !!
		} break;
		default: { 
		    assert(!"KazFileImporter::doRead unknown element type");
		    throw "KazFileImporter::doRead unknown element type";
		} break;
		}//!switch(element_type)
		elem->type_ = element_type;
		elem->parent_= (father > 0 ? hObj::makeId(father,elem->type_) : 0);
		elem->flags(GROUP_ID) = mt;
		elem->flags(REFINEMENT) = rt;
		elem->flags(EL_TYPE) = typeFlag;
		//assert(elem->test());// see comment above near face test
		
		++_elem_readed;
	}
	//assert(mesh->checkUniqueness());

	for(hHybridMesh::FacePool::Iterator<hHybridMesh::FacePool::allObj> it(&mesh->faces_); !it.done(); ++it)
	{
		if(it->neighs(0) > 0)
			it->neighs(0) = mesh->elements_.at(it->neighs(0)).id_;

		if(it->neighs(1) > 0)
			it->neighs(1) = mesh->elements_.at(it->neighs(1)).id_;

        assert(it->test(mesh));
	}

	return mesh->test();
}

int    KazFileImporter::GetVerticesCount()
{
	return _vert_count;
}

bool    KazFileImporter::GetNextVertex(double vert[])
{
	if(_vert_readed < _vert_count)
	{
		_file >> vert[0]; _file >> vert[1]; _file >> vert[2];
		++_vert_readed;
	}
	return bool(_vert_readed < _vert_count);
}

int    KazFileImporter::GetEdgesCount()
{
	int tmp;
	_file >> _edge_count >>  tmp;
	return _edge_count;
}

bool KazFileImporter::GetNextEdge(Tind vertices[],Tind & edge_type)
{
	const bool nextExist(_edge_readed < _edge_count);
	if(nextExist)
	{
        assert(vertices != NULL);
		_file >>  edge_type >> vertices[0] >> vertices[1];
		++_edge_readed;
	}
	return nextExist;
}

void     KazFileImporter::GetFacesCount(int type_count[])
{

	int tmp, face_type(0);
	if(_face_count == 0)
	{
		_file >> _face_count >> tmp;   // all - faces count
	}
	const std::streampos position(getPos());
	for(int i(0); i < _face_count; ++i)
	{
		_file >> face_type;
		assert(face_type > 2 && face_type < 5);
		++type_count[face_type];
		_file.ignore(std::numeric_limits< int >::max(), '\n' );
	_file.ignore(std::numeric_limits< int >::max(), '\n' );  
	}
	_file.seekg(position);
}

bool    KazFileImporter::GetNextFace(Tind edges[],Tind & face_type, int8_t & bc, Tind neigh[])
{
	const bool nextExist((_face_readed < _face_count) && _file.good());
	  //std::streampos	pos = getPos();
	if(nextExist)
	{
        assert(edges != NULL);
        assert(neigh != NULL);
		_file >> face_type;
		if(face_type < 2 && face_type > 4)
		{
			throw "\nKazFileImporter::GetNextFace: wrong face_type.\n Make sure that file has platform-specyfic file endings (Win/Linux).";
		}
		int tmp;
		_file >> tmp; bc=tmp;
		_file >> neigh[0] >> neigh[1];
		_file >> edges[0] >> edges[1] >> edges[2];
		if(face_type == 4)
		{
			_file >> edges[3];
		}
		// check negative values
		for(Tind i(0); i < face_type; ++i)
		{
			if(edges[i] < 0)
			{
				edges[i]*=-1;
			}
		}

		++_face_readed;
	}
	return nextExist;
}

void     KazFileImporter::GetElementCount(int type_count[])
{
	int tmp, type(0);
	if(_elem_count == 0)
	{
		_file >> _elem_count >> tmp;
	}
	const std::streampos position = getPos();
	for(int i(0); i < _elem_count; ++i)
	{
		_file >> type >> tmp >> tmp >> tmp;
		++type_count[type];
		_file.ignore(std::numeric_limits< int >::max(), '\n' );
		_file.ignore(std::numeric_limits< int >::max(), '\n' );
	}
	_file.seekg(position);
}

bool    KazFileImporter::GetNextElement(Tind vertices[], Tind neighbours[], Tind faces[], Tind & element_type , Tind & father, int8_t & material, int8_t & ref)
{
	const bool nextExist(_elem_readed < _elem_count);
	if(nextExist)
	{
		int mt,rt;
		_file >> element_type >> mt >> father >> rt ;   material=mt; ref=rt;
		--father;
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
		}
		++_elem_readed;
	}
	return nextExist;
}
