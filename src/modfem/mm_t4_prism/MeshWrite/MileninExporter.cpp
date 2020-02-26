#include "MileninExporter.h"

#include "../MeshModule/hHybridMesh.h"
#include "../MeshModule/ElemPrism.h"

using MeshWrite::MileninExporter;

MileninExporter::MileninExporter() : MeshFileWriter()
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

MileninExporter::MileninExporter(const std::string & file_name) : MeshFileWriter(file_name)
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

MileninExporter::~MileninExporter()
{
	delete [] writeSequence;
}

bool    MileninExporter::Init()
{
	return Init(_file_name);
}

bool	MileninExporter::Init(const std::string & file_name)
{
	if(MeshFileWriter::Init(file_name))
	{
		_file << 300000 << " " << 1300000 << " " << 1600000 << " " << 600000 << std::endl ;
	}
	return _file.good();
}

void    MileninExporter::Free()
{
	MeshFileWriter::Free();
}

bool	MileninExporter::doWrite(const hHybridMesh * mesh)
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
	const Edge * e( &mesh->edges_.first() );
	for(int i(0); i < noEdges; ++i)
	{
		_file << e->type_ << " " 
			<< (e->verts(0)) << " " 
			<< (e->verts(1)) << std::endl;
		e = reinterpret_cast<const Edge*>(e->next());
	}
	// face count
	const int noFaces( mesh->faces_.size() );
	_file << noFaces <<" " << (noFaces+1) << std::endl;
	// faces
	const hObj * f( &mesh->faces_.first() );
	for(int i(0); i < noFaces; ++i)
	{
		_file << f->type_ 
			<< " " << static_cast<int>(f->flags(B_COND)) 
			<< " " << static_cast<int>(hObj::posFromId(f->neighs(0))) 
			<< " " << static_cast<int>(hObj::posFromId(f->neighs(1))) << std::endl;
		
		for(uTind e(0); e < f->type_; ++e) {
		    _file << static_cast<int>(hObj::posFromId(f->components(e))) << " ";
		}
		
		//for(int e(0); e < f->type_; ++e) {
		//    _file << static_cast<int>(f->verts(e)) << " ";
		//}
		_file << std::endl;
		f = f->next();
	}
	// elem count
	const int noElems( mesh->elements_.size() );
	_file << noElems << " " << (noElems+1)<< std::endl;
	// elements
	const hObj * elem( &mesh->elements_.first() );
	for(int i(0); i < noElems; ++i)
	{
	    _file << elem->type_ << " " << static_cast<int>(elem->flags(GROUP_ID))
				<< " " << static_cast<int>(elem->parent_) 
				<< " " << static_cast<int>(elem->flags(REFINEMENT)) << std::endl;
	    for(uTind c(0); c < elem->typeSpecyfic_.nComponents_; ++c) {
		_file << static_cast<int>(hObj::posFromId(elem->components(c))) << " ";
	    }
	    
	    //for(int c(0); c < elem->typeSpecyfic_.nVerts_; ++c) {
	    //	_file << static_cast<int>(elem->verts(c)) << " ";
	    //}
	    
	    _file << std::endl;
	    elem = elem->next();
	}
	// end
	return true;
}

void MileninExporter::WriteVerticesCount(const int noVerts)
{
	_file << noVerts << " " << (noVerts+1)<< std::endl;
}

void MileninExporter::WriteVertex(const double coords[],const uTind type)
{
	_file << coords[0] << " " << coords[1] << " " << coords[2]<< std::endl;
}

void MileninExporter::WriteEdgesCount(const int noEdges)
{
	_file << noEdges << " " << (noEdges+1)<< std::endl;
}

void MileninExporter::WriteEdge(const ID verts[], const uTind type)
{
	_file << type << " " << (verts[0]+1) << " " << (verts[1]+1) << std::endl;
}

void MileninExporter::WriteFacesCount(const int noFaces)
{
	_file << noFaces <<" " << (noFaces+1) << std::endl;
}

void MileninExporter::WriteFace(const ID edges[],const uTind type,const int8_t bc, const ID neigh[])
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

void MileninExporter::WriteElementCount(const int noElems)
{
	_file << noElems << " " << (noElems+1)<< std::endl;
}

void MileninExporter::WriteElement(const ID faces[], const ID neighbours[],const uTind type, const ID father, const int8_t ref, const int8_t material)
{
	_file << type << " " << static_cast<int>(material) << " " << static_cast<int>(father) << " " << static_cast<int>(ref) << std::endl;
	const uTind nFaces = (type == 7 ? 4 : type);
	for(uTind i(0); i < nFaces; ++i)
	{
		_file << (faces[i]+1) << " ";
	}
	_file << std::endl;
}
