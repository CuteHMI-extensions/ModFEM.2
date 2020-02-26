#include "Edge.h"
#include "hHybridMesh.h"

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


template<>
void Edge::init(hHybridMesh* myMesh)
{
  components(0)=hObj::makeId(verts(0),Vertex::myType);
  components(1)=hObj::makeId(verts(1),Vertex::myType);
  assert(components(0) != 0);
  assert(components(1) != 0);
  //out_stream << "\n Creating edge";

  for(int i=0; i < typeSpecyfic_.nComponents_; ++i) {
      myMesh->vertices_.getById(components(i)).incRef();
  }
}



void    EdgeSpace::mark2Ref(hHybridMesh* myMesh,hObj& edge,const int)
{
	assert(edge.type_ == Edge::myType);
	if(edge.type_ == Edge::myType)
	{
		//Edge   * const ptr(static_cast<Edge*>(obj));
		++edge.flags(0);
		if(!edge.isBroken())
		{
			assert(!edge.isBroken());
			edge.nMyClassSons_ = hObj::fullRefMark;
			// NOTE: we know that isBroken() == false
            myMesh->edges_.requestChange(2, 2*sizeof(Edge) + sizeof(EdgeD)-sizeof(Edge) );
            myMesh->vertices_.requestChange(1,sizeof(Vertex));

		}
	}
	else {
		throw "markEdge(hObj* obj): arg obj is not an edge.";
	}
}

int	EdgeSpace::refine(hHybridMesh *myMesh, hObj& edge, const int)
{
	//  0            --->            1
	//  *----------------------------*
	//
	//  0     ->      2     ->       1
	//  *-------------*--------------*
	//
	//EdgeD & edge(*static_cast<EdgeD*>(obj));
	edge.updatePointers();
	assert(edge.type_ == EdgeD::myType);
	if(!edge.isBroken())
	{
		edge.mySize_ = sizeof( EdgeD );
		edge.nMyClassSons_ = 2;
		// create 1 new vertex
        Vertex & v = * (myMesh->vertices_.newObj<Vertex>(myMesh,NULL));
		edge.sons(0) = v.id_;
		v.parent_ = edge.id_;
		v.level_ = edge.level_+1;
//		v.status(ACTIVE);
		//NOTE: call geometry module for new point coords
        const Vertex& v1( myMesh->vertices_[edge.verts(0)] );
        const Vertex& v2( myMesh->vertices_[edge.verts(1)] );
		v.coords_[0] = (v1.coords_[0] + v2.coords_[0])/2.0;
		v.coords_[1] = (v1.coords_[1] + v2.coords_[1])/2.0;
		v.coords_[2] = (v1.coords_[2] + v2.coords_[2])/2.0;
		  // this below should be called only at boundary
		  if(v1.isAtBoundary() || v2.isAtBoundary()) {
			GeometryModule::Instance().GetPointAt(v.coords_);
		  }
		// create 2 new edges_
		uTind	vertices[2]={edge.verts(0),v.pos_};
		BYTE*   adr = reinterpret_cast<BYTE*>(&edge) + edge.mySize_;
        Edge & child = *( myMesh->edges_.newObj<Edge>(myMesh,vertices,adr) );
		child.parent_ = edge.id_;
		child.level_ = edge.level_+1;
		child.components(0) = edge.components(0);
		child.components(1) = v.id_;

		vertices[0]= hObj::posFromId(v.id_);
		vertices[1]= hObj::posFromId(edge.components(1));
		adr += child.mySize_;
        Edge & child2 = *( myMesh->edges_.newObj<Edge>(myMesh,vertices,adr) );
		child2.parent_ = edge.id_;
		child2.level_ = edge.level_+1;
		child2.components(0) = v.id_;
		child2.components(1) = edge.components(1);

        ++(myMesh->edges_.dividedObjs_);
	}
	return ( sizeof(EdgeD) - sizeof(Edge) )+ 2*sizeof( Edge );
}

void EdgeSpace::derefine(hHybridMesh *myMesh, hObj& obj)
{
	assert(obj.isBroken());
	EdgeD & edge(*static_cast<EdgeD*>(&obj));
	//hHybridMesh & m(*edge.myMesh);
	edge.sons(0)=0;
	edge.nMyClassSons_=0;
	//return ( sizeof(EdgeD) - sizeof(Edge) )+ 2*sizeof( Edge );
}

bool	EdgeSpace::test(const hHybridMesh *myMesh, const hObj & edge)
{
	assert(edge.type_ == Edge::myType);
	//const EdgeD & edge(*static_cast<const EdgeD*>(obj));
	//edge.checkPointers();
	assert( edge.components(0) != edge.components(1));
	assert( edge.verts(0) != edge.verts(1));

    assert(myMesh->vertices_.getById(edge.components(0)).type_ == Vertex::myType);
    assert(myMesh->vertices_.getById(edge.components(1)).type_ == Vertex::myType);

	if(edge.isBroken())
	{
        assert(myMesh->vertices_.getById(edge.sons(0)).type_ == Vertex::myType);
        assert(myMesh->vertices_.getById(edge.sons(0)).parent_ == edge.id_);

		assert(edge.getMyClassChild<Edge>(0)->type_ == Edge::myType);
		assert(edge.getMyClassChild<Edge>(0)->parent_ == edge.id_);

		assert(edge.getMyClassChild<Edge>(1)->type_== Edge::myType);
		assert(edge.getMyClassChild<Edge>(1)->parent_ == edge.id_);
	}
	return true;
}

//ID	 EdgeSpace::components(const hObj & obj, const int i) { return obj.components(i); }
void EdgeSpace::mark2Deref(hHybridMesh *myMesh, hObj& obj)
{
  assert(obj.type_ == EdgeD::myType);
  EdgeD & edge(*static_cast<EdgeD*>(&obj));
  hHybridMesh & m(*myMesh);
  if(edge.isBroken())
	if(edge.getMyClassChild<Edge>(0)->isBroken()==false)
	  if(edge.getMyClassChild<Edge>(1)->isBroken()==false)
	{
	// derefinement is possible only if children don't have own children
	  static const int reqChangeBytes = -2* static_cast<int>(sizeof(Edge))
		  - ( static_cast<int>(sizeof(EdgeD)) - static_cast<int>(sizeof(Edge)) );
	  m.edges_.requestChange(-2,reqChangeBytes);
	  edge.getMyClassChild<Edge>(0)->nMyClassSons_ = hObj::delMark;
	  edge.getMyClassChild<Edge>(1)->nMyClassSons_ = hObj::delMark;
	  
      m.vertices_.requestChange(-1,- sizeof(Vertex));
	  edge.nMyClassSons_ = hObj::derefMark;
  }
}
void EdgeSpace::mark2Delete(hHybridMesh *myMesh, hObj& obj){
    Edge & edge = *static_cast<Edge*>(&obj);
    if(edge.nMyClassSons_ == hObj::nothingMark) {
        edge.nMyClassSons_ = hObj::delMark;
        myMesh->edges_.requestChange(-1,-edge.mySize_);

        for(int i=0; i < obj.typeSpecyfic_.nComponents_; ++i) {
            myMesh->vertices_.getById(obj.components(i)).decRef(myMesh);
        }
    }
}

/**  @} */
