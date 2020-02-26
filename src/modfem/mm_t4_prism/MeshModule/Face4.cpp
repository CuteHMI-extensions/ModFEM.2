#include "Face4.h"
#include "hHybridMesh.h"

template<>
void Face4::init(hHybridMesh* myMesh)
{
  if(myMesh != NULL) {
	components(0)=myMesh->edge(verts(0),verts(1)).id_;
	components(1)=myMesh->edge(verts(0),verts(2)).id_;
	components(2)=myMesh->edge(verts(1),verts(3)).id_;
	components(3)=myMesh->edge(verts(2),verts(3)).id_;
  }

  for(int i=0; i < typeSpecyfic_.nComponents_; ++i) {
      myMesh->edges_.getById(components(i)).incRef();
  }
}


void    Face4Space::mark2Ref(hHybridMesh* myMesh,hObj& face,const int)
{
	if(face.type_ == Face4::myType)
	{
        hHybridMesh & m(*myMesh);
		if(!face.isBroken()) {
		  face.nMyClassSons_ = hObj::fullRefMark;
		  //Face4   & face(*static_cast<Face4*>(obj));
		  // NOTE: we know that isBroken() == false
		  if(face.nMyClassSons_==hObj::fullRefMark)
			{
			  face.nOtherSons_=4;
			  for(uTind i(0); i < face.typeSpecyfic_.nComponents_; ++i)
				{
                  m.edges_.getById(face.components(i)).mark2Ref(myMesh);
				  //m.edges_.mark2ref(hObj::posFromId(face.components(i)));
				}
			  m.faces_.requestChange(4, 4*sizeof(Face4) + sizeof(Face4D)-sizeof(Face4));
			  m.edges_.requestChange(4, 4*sizeof(Edge) );
              m.vertices_.requestChange(1, sizeof(Vertex));
			}
		}
	}
	else {
		throw "markFace4(hObj * obj): arh obj is not a face4!";
	}
}

int Face4Space::refine(hHybridMesh* myMesh,hObj& face,const int sonsCount)
{
	///
	/// (2)            (3)
	///  *------3----->*
	///  ^             ^
	///  |             |
	///  |             |
	///  1             2
	///  |             |
	///  |             |
	///  |             |
	///  *------0----->*
	/// (0)            (1)
	/// full break
	///  1 vert + 4 edges_ + 4 face4
	/// (2)            (3)
	///  *----->*----->*
	///  ^      ^      ^
	///  |  [2] 4  [3] |
	///  |      |      |
	///  *--2-->0---3->*
	///  ^      ^      ^
	///  |  [0] 1  [1] |
	///  |      |      |
	///  *----->*----->*
	/// (0)            (1)
	/// vert ordering()
	/// (2)    (7)     (3)
	///  *----->*----->*
	///  ^      ^      ^
	///  |  [2] 4  [3] |
	///  |      |      |
	///(5)--2--(8)--3->(6)
	///  ^      ^      ^
	///  |  [0] 1  [1] |
	///  |      |      |
	///  *----->*----->*
	/// (0)    (4)     (1)

	// TODO: CONTINUE
	int offset( sizeof(Face4D) - sizeof(Face4) );
	//Face4D & face(*static_cast<Face4D*>(obj));
	face.updatePointers();
	assert(face.type_ == Face4D::myType);
	if(!face.isBroken())
	{
        hHybridMesh & m(*myMesh);
		face.mySize_= sizeof( Face4D );
		face.nMyClassSons_ = 4;

		///NOTE: at this point we know that this element is NOT broken, becacouse of break()

		// 1 vertex
        Vertex & vert( *m.vertices_.newObj<Vertex>(myMesh,NULL) );
		vert.parent_= face.id_;
		vert.level_= face.level_+1;
//		vert.status(ACTIVE);
        for(int i(0); i < 4; ++i) {
            const Vertex & v(m.vertices_[face.verts(i)]);
            vert.coords_[X] += v.coords_[X];
            vert.coords_[Y] += v.coords_[Y];
            vert.coords_[Z] += v.coords_[Z];
        }
        vert.coords_[X]/=4;
        vert.coords_[Y]/=4;
        vert.coords_[Z]/=4;
		if(face.neighs(1)==B_COND) {
		  GeometryModule::Instance().GetPointAt(vert.coords_);
  
		}
		face.sons(0)=vert.id_;

		// vertices
		const uTind v[9]={
		face.verts(0),face.verts(1),face.verts(2),face.verts(3),
		hObj::posFromId(m.edges_(face.components(0)).sons(0)),
		hObj::posFromId(m.edges_(face.components(1)).sons(0)),
		hObj::posFromId(m.edges_(face.components(2)).sons(0)),
		hObj::posFromId(m.edges_(face.components(3)).sons(0)),
		vert.pos_
		};

		// 4 edges_
		assert(m.edges_.getById(face.components(0)).isBroken());
		assert(m.edges_.getById(face.components(1)).isBroken());
		assert(m.edges_.getById(face.components(2)).isBroken());
		assert(m.edges_.getById(face.components(3)).isBroken());

		// edges creation
		const uTind edges_verts[4][2]={{4,8},{5,8},{8,6},{8,7}};
		for(uTind i(1); i < 5; ++i)
		{
			const uTind vertices[]={
			v[edges_verts[i-1][0]],
			v[edges_verts[i-1][1]]
			};

            Edge & child( *m.edges_.newObj<Edge>(myMesh,vertices) );
			child.parent_= face.id_;
			child.level_= face.level_+1;
//			child.status(INACTIVE);
			child.components(0)=hObj::makeId(vertices[0],Vertex::myType);
			child.components(1)=hObj::makeId(vertices[1],Vertex::myType);
			face.sons(i)=child.id_;
		}

		// 4 face4

		const uTind verts_n_comps[4][8]=
{
		{ v[0],v[4],v[5],v[8],
		  m.edge(v[0],v[4]).id_,
		  m.edge(v[0],v[5]).id_,
		  m.edge(v[4],v[8]).id_,
		  m.edge(v[5],v[8]).id_
		},
		{ v[4],v[1],v[8],v[6],
		  m.edge(v[4],v[1]).id_,
		  m.edge(v[4],v[8]).id_,
		  m.edge(v[1],v[6]).id_,
		  m.edge(v[8],v[6]).id_
		},
		{ v[5],v[8],v[2],v[7],
		  m.edge(v[5],v[8]).id_,
		  m.edge(v[5],v[2]).id_,
		  m.edge(v[8],v[7]).id_,
		  m.edge(v[2],v[7]).id_
		},

		{ v[8],v[6],v[7],v[3],
		  m.edge(v[8],v[6]).id_,
		  m.edge(v[8],v[7]).id_,
		  m.edge(v[6],v[3]).id_,
		  m.edge(v[7],v[3]).id_
		}
		};

		BYTE * adr( reinterpret_cast<BYTE*>(&face)+face.mySize_);
		for(uTind i(0); i < 4; ++i)
		{
          Face4 & child = *m.faces_.newObj<Face4>(myMesh,verts_n_comps[i],adr);
			child.parent_= face.id_;
			child.level_= face.level_+1;
//			child.status(INACTIVE);
			adr += child.mySize_;

			child.components(0) = verts_n_comps[i][4+0];
			child.components(1) = verts_n_comps[i][4+1];
			child.components(2) = verts_n_comps[i][4+2];
			child.components(3) = verts_n_comps[i][4+3];

			child.flags(B_COND)=face.flags(B_COND);

			child.neighs(0) = 0; //face.neighs(0);
			child.neighs(1) = 0; //face.neighs(1);

            //assert(child.test());
//neighs are wrong
		}

		offset += 4*sizeof(Face4);
		++(m.faces_.dividedObjs_);
	}
	return offset;
}

void Face4Space::derefine(hHybridMesh* myMesh,hObj& obj)
{
	assert(obj.isBroken());
	obj.sons(0)=0;
	obj.sons(1)=0;
	obj.sons(2)=0;
	obj.sons(3)=0;
	obj.sons(4)=0;
	obj.nMyClassSons_=0;
}

//ID	 Face4Space::components(const hObj& obj,const int i) { return obj.components(i); };
void Face4Space::mark2Deref(hHybridMesh* myMesh,hObj& obj)
{
  assert(obj.type_ == Face4D::myType);
  Face4D & face(*static_cast<Face4D*>(&obj));
  hHybridMesh & m(*myMesh);
  for(uTind i(0); i < face.typeSpecyfic_.nComponents_; ++i)  {
        m.edges_.getById(face.components(i)).mark2Ref(myMesh);
		//m.edges_.mark2ref(hObj::posFromId(face.components(i)));
	}
	m.faces_.requestChange(4, 4*sizeof(Face4) + sizeof(Face4D)-sizeof(Face4));
	m.edges_.requestChange(4, 4*sizeof(Edge) );
    m.vertices_.requestChange(1, sizeof(Vertex) );
 
  face.nMyClassSons_ = hObj::derefMark;
  
}

void Face4Space::mark2Delete(hHybridMesh* myMesh,hObj& obj){
    Face4 & face = *static_cast<Face4*>(&obj);
    if(face.nMyClassSons_ == hObj::nothingMark) {
        face.nMyClassSons_ = hObj::delMark;
        myMesh->faces_.requestChange(-1,-face.mySize_);

        for(int i=0; i < obj.typeSpecyfic_.nComponents_; ++i) {
            myMesh->edges_.getById(obj.components(i)).decRef(myMesh);
        }
    }
}

bool Face4Space::test(const hHybridMesh* myMesh,const hObj& face)
{
    const hHybridMesh & m(*myMesh);

	//assert(m.edges_(face.components(0)).verts(0)==m.edges_(face.components(1)).verts(0));
	//assert(m.edges_(face.components(0)).verts(1)==m.edges_(face.components(2)).verts(0));
	//assert(m.edges_(face.components(1)).verts(1)==m.edges_(face.components(3)).verts(0));
	//assert(m.edges_(face.components(2)).verts(1)==m.edges_(face.components(3)).verts(1));

    assert( face.verts(0) != face.verts(1));
    assert( face.verts(1) != face.verts(2) );
    assert( face.verts(2) != face.verts(0) );

    assert(face.neighs(0) != UNKNOWN);
    assert(face.neighs(1) != UNKNOWN);
    assert(face.neighs(0) != face.neighs(1) || (m.isRefining() && face.nMyClassSons_==hObj::delMark));

//    for(int i=0;i < face.typeSpecyfic_.nVerts_; ++i){
//        assert( m.edges_(face.components(i)).verts(0) == face.verts(0));
//    }

	//const Face4 & face( *static_cast<const Face4*>(obj) );
	//face.checkPointers();
    if(face.neighs(0) != B_COND && face.neighs(0) != -2) { //subdomain boundary
            const Elem & elem = m.elements_.getById(face.neighs(0));
            const int faceNo = m.whichFaceAmI(face,elem);
            assert(faceNo >= 0 && faceNo < elem.typeSpecyfic_.nComponents_);
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
        if(face.neighs(1) != B_COND && face.neighs(1) != -2) { //subdomain boundary
            const Elem & elem = m.elements_.getById(face.neighs(1));
            const int faceNo = m.whichFaceAmI(face,elem);
            assert(faceNo >= 0 && faceNo < elem.typeSpecyfic_.nComponents_);
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
	return true;
}
