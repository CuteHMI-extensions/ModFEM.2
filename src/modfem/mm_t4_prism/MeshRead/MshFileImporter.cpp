#include "MshFileImporter.h"
#include "../MeshModule/hHybridMesh.h"
#include "../MeshModule/ElemT4.h"
// #include <algorithm>


namespace MeshRead
{
  const std::string MshFileImporter::nodesBegin="NBLOCK";
  const std::string MshFileImporter::elemsBegin="EBLOCK";
  const std::string MshFileImporter::cmBlockBegin="CMBLOCK";
  const std::string MshFileImporter::nodesEnd="N";
  
  MshFileImporter::MshFileImporter()
	:MeshFileImporter()
  {
	
  }

  MshFileImporter::MshFileImporter(const std::string & fn)
	: MeshFileImporter(fn)
  {
  }
  
  MshFileImporter::~MshFileImporter()
  {
	Free();
  }
  
  bool MshFileImporter::Init(const std::string & filename)
  {
	_file_name = filename;
	return MshFileImporter::Init();
  }
  
  bool MshFileImporter::Init()
  {
	if(MeshFileImporter::Init()) {	  
	  /*  const int sBuf=1000;
	  char buf[sBuf]={0};
	  std::string txt;
	  out_stream << "\nInFileImporter: analyzing verts";
	  // VERTS
	  _file.seekg(0);
	  do {
		_file >> txt;
		_file.ignore(sBuf,eol);
	  } while(txt != nodesBegin);

	  gridPos_ = getPos();
	  
	  do {
		_file >> txt;
		_file.ignore(sBuf,eol);
		++_vert_count;
	  } while(txt != nodesEnd);
	  _vert_count-=2;
	  // !VERTS
	  
	  // ELEMENTS 
	  do{
		out_stream << " ...found " << _vert_count  << "\nInFileImporter: analyzing elements ";
		do {
		  _file >> txt;
		  _file.ignore(sBuf,eol);
		} while(txt != elemsBegin);

		if(elemPos_ == 0) {
		  elemPos_ = getPos();
		}
		
		do {
		  _file >> txt;
		  _file.ignore(sBuf,eol);
		  ++_elem_count;
		  //out_stream << txt;
		} while((txt != "N") && (!_file.eof()));
		--_elem_count;
		out_stream << " ...found " << _elem_count;
		
	  } while(!_file.eof());
	  // !ELEMENTS

	  // CMBLOCK
	  // !CMBLOCK
	  */
	  }
	out_stream << "\nMshFileImporter: initialization end." << std::endl;
	//setPos(gridPos_);
	_file.clear();
	return _file.is_open();
  }

  void MshFileImporter::Free()
  {
	MeshFileImporter::Free();
  }

  bool MshFileImporter::doRead(hHybridMesh * mesh)
  {
	if(Init(_file_name)) {
	  assert(_file.good());
	  hHybridMesh & m(*mesh);
	  std::string txt;
	  const int bufS=1024;
	  char buf[bufS]={0};
	  
	  assert(_vert_count > 0);
	  assert(_elem_count > 0);
	  _file >> std::skipws;
	  
	  // VERTS
	  out_stream << "\nMshFileImporter: reading verts:\n";
	  _file.ignore(1000,eol); // ignore first line
	  _file >> _vert_count; // read number of vertices
	  m.vertices_.reserve(_vert_count);
	  
	  assert(_file.good());
	  while(_vert_readed < _vert_count) {
		int tmp0(0);
		Vertex & v = *m.vertices_.newObj(NULL);
		_file >> tmp0; // vert number	
		_file >> v.coords_[0] >> v.coords_[1] >> v.coords_[2];
		++_vert_readed;
		out_stream << "\r" << _vert_readed;
		assert(_file.good());
	  }
	  // !VERTS
	  
	  // ELEMENTS 
	  _file.ignore(1000,eol); // ignore line $ENDNOD
	  _file.ignore(1000,eol); // ignore line $ELM
	  _file >> _elem_count;

	  m.elements_.reserve(_elem_count, _elem_count * sizeof(ElemPrism));
	  m.faces_.reserve(_elem_count*5, _elem_count*5* sizeof(Face4));
	  m.edges_.reserve(_elem_count*9, _elem_count*9*sizeof(Edge));

	  assert(!_file.eof());
	  assert(!_file.bad());
	  if(_file.fail()) {
		_file.clear();
	  }
	  assert(_file.good());
	  do { // for each element block
		do { // find block begining
		  _file >> txt;
		  _file.ignore(bufS,eol);
		  //out_stream << txt << std::endl;
		} while((txt != "EBLOCK") && !_file.eof());
		
		_file.getline(buf,bufS);
		//out_stream << buf << std::endl;
		out_stream << "\nInFileImporter: reading elements:" << std::endl;
		register EBlockLine block;
		char nextChar(0);
		_file >> nextChar;
		while( (nextChar != '-') && !_file.eof()) { // read block elems
		  assert(_file.good());
		  _file.putback(nextChar);
		  _file >> block;
		  if(block.good) {
			hObj * obj(NULL);
			switch(block.N) {
			case 2: // 2 -this is number of nodes when an edge is in .in file
			  {
				//if(block.nNodes == Edge::nVerts) {
				//  obj = m.edges_.newObj<Edge>(block.nodes);
				//}
			  }
			case 4: // 4 - this is number of nodes when a face is in .in file
			  {
				//switch(block.nNodes) {
				//case Face3::nVerts:
				//  obj = m.faces_.newObj<Face3>(block.nodes);
				//  break;
				//case Face4::nVerts:
				//  obj = m.faces_.newObj<Face4>(block.nodes);
				//  break;
				//
				//}
			  } //// !switch(block.N) //faces
			  break;  
			case 8: // 8 - this is number of nodes when an element is in .in file.
			  switch(block.nNodes) {
			  case ElemT4::nVerts:
				std::swap(block.nodes[0],block.nodes[1]);
				obj = m.elements_.newObj<ElemT4>(block.nodes);
				break;
			  case ElemPrism::nVerts:
				obj = m.elements_.newObj<ElemPrism>(block.nodes);
				break;
			  defualt:
				throw "Unsupported element type readed from .in file";
				break;
			  }// !switch(block.nNodes) //elements
			  break;			  
			}// !switch(block.N)
			if(obj != NULL) {
			  for(int i(0);i < obj->typeSpecyfic_.nFaces_; ++i) {
				m.faces_.getById(obj->components(i)).flags(B_COND)= block.field[0];
				}
			}
			++_elem_readed;
			out_stream << "\r" << _elem_readed;
		  }
		  _file >> nextChar;
		}  // we and at "-1" at end of elements block
		out_stream << "/" << _elem_count;
	  }while( !_file.eof()); // elements blocks
	  // !ELEMENTS

	  // CMBLOCK - boundary definitions
	  _file.clear();
	  _file.seekg(elemPos_);
	  do { // find block begining
		_file >> txt;
		_file.ignore(bufS,eol);
		//out_stream << txt << std::endl;
	  } while((txt != cmBlockBegin) && !_file.eof());
	  
	  // !CMBLOCK

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


			Face * fPtr(& mesh->faces_.at(FIRST));
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
				// If both neighs are zero, than Face pointed by fPtr (teoretycally) isn't belonging to any element!
				assert(fPtr->neighs(0) != fPtr->neighs(1));
				fPtr=fPtr->next();
			}
				
			delete [] faceNeighs;
		}
		// end of checking faces neighs

	}
	return true;
  }

  std::istream& operator>>(std::istream & ifs,MshFileImporter::EBlockLine & b)
	{
		b.good=false;
		assert(ifs.good());
		if(ifs.good()){
		  ifs >> b.N;
		  ifs >> b.field[0];
		  ifs >> b.field[1];
		  ifs >> b.field[2];
		  ifs >> b.field[3];
		  ifs >> b.field[4];
		  ifs >> b.field[5];
		  ifs >> b.field[6];
		  ifs >> b.N;
		  ifs >> b.field[7];
		  ifs >> b.field[8];
		  assert(b.N > 0);
		  assert(b.N <= b.maxNodes);
		  for(int i(0); i < b.N; ++i) {
			ifs >> b.nodes[i];
		  }

		  b.nNodes = std::unique(b.nodes,b.nodes+b.N) - b.nodes;

		  assert(b.nNodes > 0);
		  assert(b.nNodes <= b.N);
		  
		  b.good=true;  
		}
		return ifs;
	  }


  
}; // namespace MeshRead
