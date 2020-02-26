#include "Face3.h"
#include "hHybridMesh.h"


template<>
void Face3::init(hHybridMesh* myMesh)
{
  if(myMesh != NULL)
	{
	  components(0)=myMesh->edge(verts(0),verts(1)).id_;
	  components(1)=myMesh->edge(verts(0),verts(2)).id_;
	  components(2)=myMesh->edge(verts(1),verts(2)).id_;  
	}

  for(int i=0; i < typeSpecyfic_.nComponents_; ++i) {
      myMesh->edges_.getById(components(i)).incRef();
  }
}



void    Face3Space::mark2Ref(hHybridMesh* myMesh,hObj& face,const int)
{
  assert(face.type_ == Face3::myType);
	
  if(!face.isBroken()) {
		face.nMyClassSons_ = hObj::fullRefMark;
		//Face3   & face(*static_cast<Face3*>(obj));
		// NOTE: we know that isBroken() == false
		face.nOtherSons_=3; // was 2 - why?
		for(uTind i(0); i < face.typeSpecyfic_.nComponents_; ++i)
		{
            myMesh->edges_.getById(face.components(i)).mark2Ref(myMesh);
			//ptr->myMesh->edges_.mark2ref(hObj::posFromId(ptr->components(i)));
		}
        myMesh->faces_.requestChange(4, 4*sizeof(Face3) + sizeof(Face3D)-sizeof(Face3));
        myMesh->edges_.requestChange(3, 3*sizeof(Edge));
  }
}

int Face3Space::refine(hHybridMesh* myMesh,hObj& f,const int)
{
	/// full break
	/// 0-3 - subfaces_, 4-6 - subedges_
	/// /*
	/// |\.
	/// |2 \.
	/// |    \.
	///  ---6--\.
	/// |\    3| \.
	/// |  4   5   \.
	/// |0___\_|____1\.   */

	Face3D & face( *static_cast<Face3D*>(&f) );
	face.updatePointers();
	assert(face.type_ == Face3D::myType);
	if(!face.isBroken())
	{
        hHybridMesh & m(*myMesh);
		face.mySize_= sizeof(Face3D);
		face.nMyClassSons_=4;

		assert(m.edges_.getById(face.components(0)).isBroken());
		assert(m.edges_.getById(face.components(1)).isBroken());
		assert(m.edges_.getById(face.components(2)).isBroken());
 
		const Vertex*	v[6]={ 
		    &m.vertices_[face.verts(0)],
		    &m.vertices_[face.verts(1)], 
		    &m.vertices_[face.verts(2)],
		    &m.vertices_.getById(m.edge(face.verts(0),face.verts(1)).sons(0)),
		    &m.vertices_.getById(m.edge(face.verts(0),face.verts(2)).sons(0)),
		    &m.vertices_.getById(m.edge(face.verts(1),face.verts(2)).sons(0)) };

		// new edges_
		// indexes of verts in v[6] for new edges
		static const int edgesVerts[3][2]={{3,4},{3,5},{4,5}};
		for(uTind i(0); i < 3; ++i)
		{
			const uTind vertices[2] ={ 
			v[ edgesVerts[i][0] ]->pos_,
			v[ edgesVerts[i][1] ]->pos_};
            Edge & child = *m.edges_.newObj<Edge>(myMesh,vertices);
			child.parent_= face.id_;
			child.level_= face.level_+1;
			face.sons(i)=child.id_;
//			child.status(INACTIVE);
			child.components(0)=v[ edgesVerts[i][0] ]->id_;
			child.components(1)=v[ edgesVerts[i][1] ]->id_;
		}

		// new faces3
		// define vertices  building childrens
		const uTind verts[4][4]=
		{
			{v[0]->pos_,v[3]->pos_,v[4]->pos_,UNKNOWN},
			{v[3]->pos_,v[1]->pos_,v[5]->pos_,UNKNOWN},
			{v[4]->pos_,v[5]->pos_,v[2]->pos_,UNKNOWN},
			{// this one is flipped
			v[3]->pos_,v[5]->pos_,v[4]->pos_,UNKNOWN}
		};

		#ifdef _DEBUG		
				double parentNormVec[3]={0.0},
				  childNormVec[3]={0.0};
		  
				m.faceNormal(face,parentNormVec,NULL);
		#endif
		
		BYTE * adr( reinterpret_cast<BYTE*>(&face)+ face.mySize_);
		for(uTind i(0); i < 4; ++i)
		{
          Face3 & faceChild =  *m.faces_.newObj<Face3>(myMesh,verts[i],adr);
			faceChild.parent_= face.id_;
			faceChild.level_= face.level_+ 1;
//			faceChild.status(INACTIVE);
			adr += faceChild.mySize_;
			// code below is no more needed, as this is coverd in Face3 constructor
			//faceChild.components(0)=verts_n_comps[i][3  ];
			//faceChild.components(1)=verts_n_comps[i][3+1];
			//faceChild.components(2)=verts_n_comps[i][3+2];
			
			faceChild.neighs(0) = 0;
			faceChild.neighs(1) = 0;
			faceChild.flags(B_COND) = face.flags(B_COND);
			//faceChild.neighs(0) = face.neighs(0);
			//faceChild.neighs(1) = face.neighs(1);
			#ifdef _DEBUG
			m.faceNormal(faceChild, childNormVec,NULL);
			assert(abs(parentNormVec[0]-childNormVec[0]) <= SMALL);
			assert(abs(parentNormVec[1]-childNormVec[1]) <= SMALL);
			assert(abs(parentNormVec[2]-childNormVec[2]) <= SMALL);
			#endif
			
			if(i==3) { // fourth is flipped
			  faceChild.flags(F_TYPE)=F_FLIPPED;
			}
		}


/* #ifdef _DEBUG
		out_stream <<"\n refine Face3("<<face.pos_<<") "
			<<face.verts(0)<<" "
			<<face.verts(1)<<" "
			<<face.verts(2);
		for(int i(0); i < 3; ++i)
		{
			Face3 & faceChild = *face.getMyClassChild<Face3>(i);
			out_stream <<"\n child Face3("<<faceChild.pos_<<") : "
			<<faceChild.verts(0)<<" "
			<<faceChild.verts(1)<<" "
			<<faceChild.verts(2);
		}
#endif */
	    ++(m.faces_.dividedObjs_);
	}
	return ( sizeof(Face3D) - sizeof(Face3) ) + 4*sizeof(Face3);
}

void Face3Space::derefine(hHybridMesh* myMesh,hObj& obj)
{
	assert(obj.isBroken());
	obj.sons(0)=0;
	obj.sons(1)=0;
	obj.sons(2)=0;
	obj.nMyClassSons_=0;
}

bool	Face3Space::test(const hHybridMesh* myMesh,const hObj & face)
{
	//const Face3D & face( *static_cast<const Face3D*>(obj) );

    const hHybridMesh & m(*myMesh);
	for(int i(0); i < 3; ++i)
	{
		assert( m.edges_.getById(face.components(i)).type_ == Edge::myType);
	}
	assert( face.verts(0) != face.verts(1));
	assert( face.verts(1) != face.verts(2) );
	assert( face.verts(2) != face.verts(0) );
	
	//if(face.flags(F_TYPE)==F_FLIPPED) {
	//assert(m.edges_.getById(face.components(0)).verts(0)==
	//	m.edges_.getById(face.components(1)).verts(0));
	//assert(m.edges_.getById(face.components(0)).verts(1)==
	//	m.edges_.getById(face.components(2)).verts(1));
	//assert(m.edges_.getById(face.components(1)).verts(1)==
	//	m.edges_.getById(face.components(2)).verts(0));
	//} 
	//else {
	//assert(m.edges_.getById(face.components(0)).verts(0)==
	//	m.edges_.getById(face.components(1)).verts(0));
	//assert(m.edges_.getById(face.components(0)).verts(1)==
	//	m.edges_.getById(face.components(2)).verts(0));
	//assert(m.edges_.getById(face.components(1)).verts(1)==
	//	m.edges_.getById(face.components(2)).verts(1));
	//}

	//assert( *m.edges_.getById(face.components(0)) < *m.edges_.getById(face.components(1)));
	//assert( *m.edges_.getById(face.components(1)) < *m.edges_.getById(face.components(2)));
	//assert( m.faceVertex(face, 0) < m.faceVertex(face, 1) );
	//assert( m.faceVertex(face, 1) < m.faceVertex(face, 2) );

	//if(face.neighs(0) != B_COND) {
    assert(face.neighs(0) != UNKNOWN);
    assert(face.neighs(1) != UNKNOWN);
    assert(face.neighs(0) != face.neighs(1) || (m.isRefining() && face.nMyClassSons_==hObj::delMark));

    if(face.neighs(0) != B_COND && face.neighs(0) != -2) { //subdomainn boundary
            const Elem & elem = m.elements_.getById(face.neighs(0));
            const int faceNo = m.whichFaceAmI(face,elem);
            assert(faceNo >= 0 && faceNo < 4);
            assert( face.id_ == elem.components(faceNo));

            for(int v=0; v < face.typeSpecyfic_.nVerts_; ++v) {
                bool found = false;
                for(int v2=0; v2 < elem.typeSpecyfic_.nVerts_ && !found; ++v2) {
                    if(face.verts(v) == elem.verts(v2)) {
                        found = true;
                    }
                }
                assert(found == true);
            }

        }
        if(face.neighs(1) != B_COND && face.neighs(1) != -2 ) { //subdomain boundary
            const Elem & elem = m.elements_.getById(face.neighs(1));
            const int faceNo = m.whichFaceAmI(face,elem);
            assert(faceNo >= 0 && faceNo < 4);
            assert( face.id_ == elem.components(faceNo));

            for(int v=0; v < face.typeSpecyfic_.nVerts_; ++v) {
                bool found = false;
                for(int v2=0; v2 < elem.typeSpecyfic_.nVerts_ && !found; ++v2) {
                    if(face.verts(v) == elem.verts(v2)) {
                        found = true;
                    }
                }
                assert(found == true);
            }
	    }
	//}
	if(face.isBroken()) {
		for(int i(0); i < 3; ++i)
		{
			assert( m.edges_.getById(face.sons(i)).type_ == Edge::myType);
			assert( m.edges_.getById(face.components(i)).type_ == Edge::myType);
		}
		for(int i(0); i < 4; ++i)
		{
			assert( face.getMyClassChild<Face3>(i)->type_ == Face3::myType);
		}
	}
	return true;
}

//ID	 Face3Space::components(const hObj& obj,const int i) { return obj.components(i);};
void Face3Space::mark2Deref(hHybridMesh* myMesh,hObj& obj)
{
  assert(obj.type_ == Face3D::myType);
  Face3D & face(*static_cast<Face3D*>(&obj));
  hHybridMesh & m(*myMesh);
  // we assume that element-level derefinement take care of permission to derefinement
  
  for(uTind i(0); i < face.typeSpecyfic_.nComponents_; ++i) {
    m.edges_.getById(face.components(i)).mark2Dref(myMesh);
  }
  m.faces_.requestChange(-4, -4*static_cast<Tind>(sizeof(Face3)) - static_cast<Tind>(sizeof(Face3D))+static_cast<Tind>(sizeof(Face3)));
  m.edges_.requestChange(-3, -3*static_cast<Tind>(sizeof(Edge)));
  
  face.nMyClassSons_ = hObj::derefMark;
}

void Face3Space::mark2Delete(hHybridMesh* myMesh,hObj& obj){
    Face3 & face = *static_cast<Face3*>(&obj);
    if(face.nMyClassSons_ == hObj::nothingMark) {
        face.nMyClassSons_ = hObj::delMark;
        myMesh->faces_.requestChange(-1,-face.mySize_);

        for(int i=0; i < obj.typeSpecyfic_.nComponents_; ++i) {
            myMesh->edges_.getById(obj.components(i)).decRef(myMesh);
        }
    }
}

