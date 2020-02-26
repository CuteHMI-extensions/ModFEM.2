#include "KazFileExporter.h"

#include "../MeshModule/hHybridMesh.h"
#include "../MeshModule/ElemPrism.h"

using MeshWrite::KazFileExporter;

KazFileExporter::KazFileExporter() : MeshFileWriter()
{
	writeSequence = new char[9];
	writeSequence[0]=V_COUNT;
	writeSequence[1]=VERTICES;
	writeSequence[2]=E_COUNT;
	writeSequence[3]=EDGES;
	writeSequence[4]=F_COUNT;
	writeSequence[5]=FACES;
	writeSequence[6]=EL_COUNT;
	writeSequence[7]=ELEMENTS;
	writeSequence[8]=END;
}

KazFileExporter::KazFileExporter(const std::string & file_name) : MeshFileWriter(file_name)
{
	writeSequence = new char[9];
	writeSequence[0]=V_COUNT;
	writeSequence[1]=VERTICES;
	writeSequence[2]=E_COUNT;
	writeSequence[3]=EDGES;
	writeSequence[4]=F_COUNT;
	writeSequence[5]=FACES;
	writeSequence[6]=EL_COUNT;
	writeSequence[7]=ELEMENTS;
    writeSequence[8]=END;
}

KazFileExporter::~KazFileExporter()
{
	delete [] writeSequence;
}

bool    KazFileExporter::Init()
{
	return Init(_file_name);
}

bool	KazFileExporter::Init(const std::string & file_name)
{
	if(MeshFileWriter::Init(file_name))
	{
		_file << 300000 << " " << 1300000 << " " << 1600000 << " " << 600000 << std::endl ;
	}
	return _file.good();
}

void    KazFileExporter::Free()
{
	MeshFileWriter::Free();
}

bool	KazFileExporter::doWrite(const hHybridMesh * mesh)
{
	// vertices count
	const int noVerts( mesh->vertices_.size() );
	_file << noVerts << " " << (noVerts+1)<< std::endl;
	// vertices
	const Vertex * v( &mesh->vertices_.first() );
	for(int i(0); i < noVerts; ++i)
	{
		_file << v->coords_[0] << " " << v->coords_[1] << " " << v->coords_[2]<< std::endl;
		v = reinterpret_cast<const Vertex*>(v->next());
	}
	//edges count
	const int noEdges( mesh->edges_.size() );
	_file << noEdges << " " << (noEdges+1)<< std::endl;
	//edges
    //const Edge * e(NULL);
    //for(int i(FIRST); i <= noEdges; ++i)
    for(hHybridMesh::EdgePool::constIterator<> e(& mesh->edges_); !e.done(); ++e)
    {
//		e = & mesh->edges_[i];
	  _file << e->type_ << " "
			<< e->pos_ << " "
			<< static_cast<int>(e->parent_) << " "
			<< static_cast<int>(e->level_) << " "
			<< (e->verts(0)) << " " 
			<< (e->verts(1)) << std::endl;
	}
	// face count
	const int noFaces( mesh->faces_.size() );
	_file << noFaces <<" " << (noFaces+1) << std::endl;
	// faces
    //const hObj * f( NULL );
    //for(int i(FIRST); i <= noFaces; ++i)
    for(hHybridMesh::FacePool::constIterator<> f(& mesh->faces_); !f.done(); ++f)
	  {
//		f = & mesh->faces_[i];
	  _file << f->type_
			<< " " << static_cast<int>(f->pos_)
			<< " " << static_cast<int>(f->parent_)
			<< " " << static_cast<int>(f->level_)
			<< " " << static_cast<int>(f->flags(B_COND)) 
			<< " " << static_cast<int>(hObj::posFromId(f->neighs(0))) 
			<< " " << static_cast<int>(hObj::posFromId(f->neighs(1))) 
			<< " " << static_cast<int>(f->flags(F_TYPE))<< std::endl;
		
		for(uTind e(0); e < f->type_; ++e) {
		    _file << static_cast<int>(hObj::posFromId(f->components(e))) << " ";
		}
		
		for(uTind e(0); e < f->type_; ++e) {
		    _file << static_cast<int>(f->verts(e)) << " ";
		}
		_file << std::endl;
	}
	// elem count
	const int noElems( mesh->elements_.size() );
	_file << noElems << " " << (noElems+1)<< std::endl;
	// elements
    //const hObj * elem( &mesh->elements_.first() );
    //for(int i(0); i < noElems; ++i)
    for(hHybridMesh::ElemPool::constIterator<> it(& mesh->elements_); !it.done(); ++it)
	{
        _file << it->type_ << " " << static_cast<int>(it->flags(GROUP_ID))
                << " " << static_cast<int>(it->parent_)
                << " " << static_cast<int>(it->flags(REFINEMENT))
                << " " << static_cast<int>(it->flags(EL_TYPE)) << std::endl;
        for(uTind c(0); c < it->typeSpecyfic_.nComponents_; ++c) {
        _file << static_cast<int>(hObj::posFromId(it->components(c))) << " ";
	    }
	    
        for(uTind c(0); c < it->typeSpecyfic_.nVerts_; ++c) {
            _file << static_cast<int>(it->verts(c)) << " ";
	    }
	    
	    _file << std::endl;
	}
	// end
	return true;
}

void KazFileExporter::WriteVerticesCount(const int noVerts)
{
	_file << noVerts << " " << (noVerts+1)<< std::endl;
}

void KazFileExporter::WriteVertex(const double coords[],const uTind type)
{
	_file << coords[0] << " " << coords[1] << " " << coords[2]<< std::endl;
}

void KazFileExporter::WriteEdgesCount(const int noEdges)
{
	_file << noEdges << " " << (noEdges+1)<< std::endl;
}

void KazFileExporter::WriteEdge(const ID verts[], const uTind type)
{
	_file << type << " " << (verts[0]+1) << " " << (verts[1]+1) << std::endl;
}

void KazFileExporter::WriteFacesCount(const int noFaces)
{
	_file << noFaces <<" " << (noFaces+1) << std::endl;
}

void KazFileExporter::WriteFace(const ID edges[],const uTind type,const int8_t bc, const ID neigh[])
{
	_file << type <<" " << static_cast<int>(bc) 
		<< " " << static_cast<int>(neigh[0]+1) 
		<< " " << static_cast<int>(neigh[1]+1) << std::endl;
	for(uTind i(0);i < type; ++i)
	{
		_file << (edges[i]+1) << " ";
	}
	_file << std::endl;
}

void KazFileExporter::WriteElementCount(const int noElems)
{
	_file << noElems << " " << (noElems+1)<< std::endl;
}

void KazFileExporter::WriteElement(const ID faces[], const ID neighbours[],const uTind type, const ID father, const int8_t ref, const int8_t material)
{
	_file << type << " " << static_cast<int>(material) << " " << static_cast<int>(father) << " " << static_cast<int>(ref) << std::endl;
	const uTind nFaces = (type == 7 ? 4 : type);
	for(uTind i(0); i < nFaces; ++i)
	{
		_file << (faces[i]+1) << " ";
	}
	_file << std::endl;
}
