#include "ElemT4.h"
#include "hHybridMesh.h"
#include "mmh_vec3.h"
#include <modfem/uth_log.h>

template<>
void ElemT4::init(hHybridMesh * myMesh)
{
	if (myMesh != NULL) {
		const Vertex * v[4] = {&myMesh->vertices_[verts(0)],
						&myMesh->vertices_[verts(1)],
						&myMesh->vertices_[verts(2)],
						&myMesh->vertices_[verts(3)]
				};
		const double v0[3] = {v[1]->coords_[0] - v[0]->coords_[0],
						v[1]->coords_[1] - v[0]->coords_[1],
						v[1]->coords_[2] - v[0]->coords_[2]
				};
		const double v1[3] = {v[2]->coords_[0] - v[0]->coords_[0],
						v[2]->coords_[1] - v[0]->coords_[1],
						v[2]->coords_[2] - v[0]->coords_[2]
				};
		const double v2[3] = {v[3]->coords_[0] - v[0]->coords_[0],
						v[3]->coords_[1] - v[0]->coords_[1],
						v[3]->coords_[2] - v[0]->coords_[2]
				};
		const double volume = mmr_vec3_mxpr(v0, v1, v2) / 6.0; //for t4
		if (volume < 0.0) {
			swapVerts(1, 2);
		}

		components(0) = myMesh->face(verts(0), verts(1), verts(2)).id_;
		components(1) = myMesh->face(verts(0), verts(1), verts(3)).id_;
		components(2) = myMesh->face(verts(0), verts(2), verts(3)).id_;
		components(3) = myMesh->face(verts(1), verts(2), verts(3)).id_;

		for (int i = 0; i < typeSpecyfic_.nComponents_; ++i) {
			Face & f = myMesh->faces_.getById(components(i));
			f.incRef();

			if (f.neighs(0) == UNKNOWN
					|| f.neighs(0) == B_COND) {
				f.neighs(0) = id_;
			}
			else {
				f.neighs(1) = id_;
			}
		}

	}
}

void    ElemT4Space::mark2Ref(hHybridMesh * myMesh, hObj & el, const int) {
	if (el.type_ == ElemT4::myType) {
		ElemT4D & elem(*static_cast<ElemT4D *>(&el));
		hHybridMesh & m(*myMesh);
		assert(!elem.isBroken());
		elem.nMyClassSons_ = hObj::fullRefMark;

		// Checking if that forces neighbours to also break
		for (uTind i(0); i < elem.typeSpecyfic_.nComponents_; ++i) {
			const ID neigId(m.elemNeigh(elem, i));
			if (neigId != B_COND) {
				if (m.elements_(neigId).level_ < elem.level_) {
					m.elements_(neigId).mark2Ref(myMesh);
				}
			}
		}

		//ElemT4   * const ptr(static_cast<ElemT4*>(obj));
		// NOTE: we know that isBroken() == false
		if (elem.nMyClassSons_ == hObj::fullRefMark) {
			elem.nOtherSons_ = 8;
			for (uTind i(0); i < elem.typeSpecyfic_.nComponents_; ++i) {
				myMesh->faces_.getById(elem.components(i)).mark2Ref(myMesh);
				//ptr->myMesh->faces_.mark2ref(hObj::posFromId(ptr->components(i)));
			}
			myMesh->elements_.requestChange(8, 8 * sizeof(ElemT4) + (sizeof(ElemT4D) - sizeof(ElemT4)));
			myMesh->faces_.requestChange(4 + 4, 8 * sizeof(Face3) );
			myMesh->edges_.requestChange(1, sizeof(Edge) );
		}
	} else
		throw "markElemT4(hObj * obj): arg obj is not an ElemT4.";
}

int	ElemT4Space::refine(hHybridMesh * myMesh, hObj & el, const int sonsCount) {
#include "./ElemT4tables.hpp"

	int offset( sizeof(ElemT4D) - sizeof(ElemT4) );
	ElemT4D & elem(*static_cast<ElemT4D *>(&el));
	elem.updatePointers();
	assert(elem.type_ == ElemT4D::myType);

	if (!elem.isBroken()) {
		hHybridMesh & m(*myMesh);
		elem.mySize_ = sizeof(ElemT4D);
		elem.nMyClassSons_ = sonsCount;
		Face3D * const component[4] = {
			& m.faces_.getById<Face3D>(elem.components(0)),
			& m.faces_.getById<Face3D>(elem.components(1)),
			& m.faces_.getById<Face3D>(elem.components(2)),
			& m.faces_.getById<Face3D>(elem.components(3))
		};
		assert(component[0]->isBroken());
		assert(component[1]->isBroken());
		assert(component[2]->isBroken());
		assert(component[3]->isBroken());

		// 0-7 - (8)subintrs, 8-15 - (8)subfaces_, 16 - (1)subedge
		// Info subfaces_ 8-15:
		// 8-11 - subfaces_ outside octaedr, where:
		// 8 - near vertex 0 [457]
		// 9 - near vertex 1 [468]
		// 10- near vertex 2 [965]
		// 11- near vertex 3 [987]
		//
		// 12-15 - subfaces_ inside octaedr, where:
		// refKind	6-7	|	4-9	|	5-8
		// 12-		[647]		[945]		[458]
		// 13-		[657]		[649]		[658]
		// 14-		[678]		[749]		[758]
		// 15-		[976]		[948]		[958]
		//
		// edge
		// configure
		const Vertex * const	v[10] = {
			& m.vertices_.at(elem.verts(0)),
			& m.vertices_.at(elem.verts(1)),
			& m.vertices_.at(elem.verts(2)),
			& m.vertices_.at(elem.verts(3)),
			& m.vertices_.getById(m.edge(elem.verts(0), elem.verts(1)).sons(0)),
			& m.vertices_.getById(m.edge(elem.verts(0), elem.verts(2)).sons(0)),
			& m.vertices_.getById(m.edge(elem.verts(1), elem.verts(2)).sons(0)),
			& m.vertices_.getById(m.edge(elem.verts(0), elem.verts(3)).sons(0)),
			& m.vertices_.getById(m.edge(elem.verts(1), elem.verts(3)).sons(0)),
			& m.vertices_.getById(m.edge(elem.verts(2), elem.verts(3)).sons(0))
		};
		// computing distances to choose best refinement option
		const Tval	dist[3] = {	// 6-7
			(v[7]->coords_[0] - v[6]->coords_[0]) * (v[7]->coords_[0] - v[6]->coords_[0])
			+ (v[7]->coords_[1] - v[6]->coords_[1]) * (v[7]->coords_[1] - v[6]->coords_[1])
			+ (v[7]->coords_[2] - v[6]->coords_[2]) * (v[7]->coords_[2] - v[6]->coords_[2]),
			// 4-9
			(v[9]->coords_[0] - v[4]->coords_[0]) * (v[9]->coords_[0] - v[4]->coords_[0])
			+ (v[9]->coords_[1] - v[4]->coords_[1]) * (v[9]->coords_[1] - v[4]->coords_[1])
			+ (v[9]->coords_[2] - v[4]->coords_[2]) * (v[9]->coords_[2] - v[4]->coords_[2]),
			// 5-8
			(v[8]->coords_[0] - v[5]->coords_[0]) * (v[8]->coords_[0] - v[5]->coords_[0])
			+ (v[8]->coords_[1] - v[5]->coords_[1]) * (v[8]->coords_[1] - v[5]->coords_[1])
			+ (v[8]->coords_[2] - v[5]->coords_[2]) * (v[8]->coords_[2] - v[5]->coords_[2])
		};
		// choosing best ref option
		eRefKind	subRefKind;
		uTind	vertices[4] = {UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN};
		if (dist[0] < dist[1]) {
			if (dist[0] < dist[2]) {
				subRefKind = eRef_67;
				vertices[0] = std::min(v[6]->pos_, v[7]->pos_);
				vertices[1] = std::max(v[6]->pos_, v[7]->pos_);
			} else {
				subRefKind = eRef_58;
				vertices[0] = std::min(v[5]->pos_, v[8]->pos_);
				vertices[1] = std::max(v[5]->pos_, v[8]->pos_);
			}
		} else {
			if (dist[1] < dist[2]) {
				subRefKind = eRef_49;
				vertices[0] = std::min(v[4]->pos_, v[9]->pos_);
				vertices[1] = std::max(v[4]->pos_, v[9]->pos_);
			} else {
				subRefKind = eRef_58;
				vertices[0] = std::min(v[5]->pos_, v[8]->pos_);
				vertices[1] = std::max(v[5]->pos_, v[8]->pos_);
			}
		}
		// creating appropriate edge for chosen ref kind
		uTind i(0);
		Edge & e( *m.edges_.newObj<Edge>(myMesh, vertices) );
		e.parent_ = elem.id_;
		e.level_ = elem.level_ + 1;
		elem.sons(i) = e.id_;
		e.components(0) = hObj::makeId(vertices[0], Vertex::myType);
		e.components(1) = hObj::makeId(vertices[1], Vertex::myType);
//        e.status(INACTIVE);
		assert(e.test(myMesh));
		// !edge

		// faces_
		for (++i; i < elem.typeSpecyfic_.nSons_; ++i) {
			vertices[0] = v[ faceVrtsIdx[i][subRefKind][0] ]->pos_;
			vertices[1] = v[ faceVrtsIdx[i][subRefKind][1] ]->pos_;
			vertices[2] = v[ faceVrtsIdx[i][subRefKind][2] ]->pos_;
			Face3 & f3( *m.faces_.newObj<Face3>(myMesh, vertices) );
			f3.parent_ = elem.id_;
			f3.level_ = elem.level_ + 1;
			elem.sons(i) = f3.id_;
			// code below is already done in Face3 constructor
			//f3.components(0)=m.edge(f3.verts(0),f3.verts(1)).id_;
			//f3.components(1)=m.edge(f3.verts(0),f3.verts(2)).id_;
			//f3.components(2)=m.edge(f3.verts(1),f3.verts(2)).id_;
//            f3.status(INACTIVE);
		}
		//!faces
		/*
		// debug print
		printf("\n Elem type: %d, ref-kind: %d \n",elem.flags(EL_TYPE),subRefKind);
		*/
		// elem creation
		BYTE * adr( reinterpret_cast<BYTE *>(&elem) + elem.mySize_);
		for (uTind t(0); t < 8; ++t) {
			vertices[0] = v[ elemVrtsIdx[t][subRefKind][0] ]->pos_;
			vertices[1] = v[ elemVrtsIdx[t][subRefKind][1] ]->pos_;
			vertices[2] = v[ elemVrtsIdx[t][subRefKind][2] ]->pos_;
			vertices[3] = v[ elemVrtsIdx[t][subRefKind][3] ]->pos_;
			ElemT4 & t4( *m.elements_.newObj<ElemT4 >(myMesh, vertices, adr) );
			t4.parent_ = elem.id_;
			t4.level_ = elem.level_ + 1;
			adr += t4.mySize_;
			offset += t4.mySize_;
			t4.flags(GROUP_ID) = elem.flags(GROUP_ID);

			hObj * const f[4] = {
				& m.face(vertices[0], vertices[1], vertices[2]),
				& m.face(vertices[0], vertices[1], vertices[3]),
				& m.face(vertices[0], vertices[2], vertices[3]),
				& m.face(vertices[1], vertices[2], vertices[3])
			};
// THis is done inside ElemT4::init()
//            for(uTind i(0); i < elem.typeSpecyfic_.nFaces_; ++i) {
//			  //t4.components(i)=f[i]->id_;
//    	        const int freeNeig = ( f[i]->neighs(0)==0 ? 0:1);
//    	        f[i]->neighs(freeNeig)=t4.id_;
//            }

//            t4.status(INACTIVE);
			t4.flags(EL_TYPE) = elem.flags(EL_TYPE);
			//m.normalizeElem(t4);
			int & mask = t4.flags(EL_TYPE);
			mask = 0;
			mask |= (m.faceDirection(*f[0], t4) == -1 ? 1 : 0);
			mask |= (m.faceDirection(*f[1], t4) == -1 ? 1 : 0) << 1;
			mask |= (m.faceDirection(*f[2], t4) == -1 ? 1 : 0) << 2;
			mask |= (m.faceDirection(*f[3], t4) == -1 ? 1 : 0) << 3;

			// debug print
			/*
			printf("{ %d,%d,%d,%d }(%d),\n",m.faceDirection(*f[3],t4)==-1?1:0,
				   m.faceDirection(*f[2],t4)==-1?1:0,m.faceDirection(*f[1],t4)==-1?1:0,
				   m.faceDirection(*f[0],t4)==-1?1:0,mask);
			*/
			// elem types for [parent type][son no]
			//static const int elemType[][]={{0},
			//EL_TYPE0
			//EL_TYPE1
			// nothing
			//EL_TYPE2
			//EL_TYPE3
			/*	};

			// adjusting child elem types
			// check if element is boundary-special-care (faces are flipped)

			switch(mask) {
			case EL_TYPE0:
			case EL_TYPE1:
			case EL_TYPE2:

			  break;
			case EL_TYPE3:
			  break;
			default:
			  {
				// check if element is type II
				switch(t) {
				case 4: //4th element
				case 6: //6th element
				  mask^=EL_TYPE0; // toogle face 0 bit
				  mask^=EL_TYPE3; // toogle face 3 bit
				  break;
				}
			  }

			assert(m.faceDirection(*f[0],t4)==F_OUT || (mask&EL_TYPE0));
			assert(m.faceDirection(*f[1],t4)==F_OUT || (mask&EL_TYPE1));
			assert(m.faceDirection(*f[2],t4)==F_OUT || (mask&EL_TYPE2));
			assert(m.faceDirection(*f[3],t4)==F_OUT || (mask&EL_TYPE3));
			*/
		}
		// At this point all new mesh entities are created(allocated). Now setup.

		// Setting neighbors.
		//First: setting neighbours for internal faces
		for (i = 1; i < elem.typeSpecyfic_.nSons_; ++i) {
			Face3  & f3( m.faces_.getById<Face3>(elem.sons(i)));
			f3.neighs(0) = elem.getMyClassChild<ElemT4>(innerFaceNeigs[i][subRefKind][0])->id_;
			f3.neighs(1) = elem.getMyClassChild<ElemT4>(innerFaceNeigs[i][subRefKind][1])->id_;
			assert(f3.test(myMesh));
		}

		/* // the functionality from below was moved to elem creation
		// second: external faces: select cross boundary neighbours
				const ID neighsID[4]={m.elemNeigh(elem,0),
		// 'BIG' neighs IDs
									  m.elemNeigh(elem,1),
									  m.elemNeigh(elem,2),
									  m.elemNeigh(elem,3)  };
				ElemT4 *neighs[4]= { // BIG neighs pointers
					(neighsID[0] != 0) ? & m.elements_.getById<ElemT4>(neighsID[0]) : NULL,
					(neighsID[1] != 0) ? & m.elements_.getById<ElemT4>(neighsID[1]) : NULL,
					(neighsID[2] != 0) ? & m.elements_.getById<ElemT4>(neighsID[2]) : NULL,
					(neighsID[3] != 0) ? & m.elements_.getById<ElemT4>(neighsID[3]) : NULL,
				};

				ID smallNeighs[4][4]={ // basic setup of small neig elems
					{(neighs[0] != NULL) ? neighs[0]->id_ : 0},
					{(neighs[1] != NULL) ? neighs[1]->id_ : 0},
					{(neighs[2] != NULL) ? neighs[2]->id_ : 0},
					{(neighs[3] != NULL) ? neighs[3]->id_ : 0}
				};

				for (int i(0); i < 4; ++i) { // setup of small neigs elems if exist (==not a BC)
					if (neighs[i] != NULL && neighs[i]->isBroken()) {
						int faceNo=m.whichNeighAmI(elem,*neighs[i]);
						assert(faceNo >= 0 && faceNo < 4);
						for (int smallT4(0); smallT4 < 4; ++smallT4) {
							smallNeighs[i][smallT4] =
							neighs[i]->getMyClassChild<ElemT4>(offspringAtFace[subRefKind][faceNo][smallT4])->id_;
						}
					}
					else if(neighs[i]!=NULL) { // but neighs[i]->isBroken()==false
					}
				}

				// setting neighbours for external faces outside
				// (not sons of element, but sons of faces creating the element)
				for (int f(0); f < 4; ++f) {
					int neighNo = 0;
					if(neighs[f]!=NULL){
					neighNo = neighs[f]->isBroken() ? 1 : 0;
					}
					for (int child(0); child < 4; ++child) {
							component[f]->getMyClassChild<Face3>(child)->neighs(neighNo)=
							elem.getMyClassChild<ElemT4>(offspringAtFace[subRefKind][f][child])->id_;
					}
					assert(component[f]->test());
				}
		*/
		assert(elem.getMyClassChild<ElemT4>(0)->test(myMesh));
		assert(elem.getMyClassChild<ElemT4>(1)->test(myMesh));
		assert(elem.getMyClassChild<ElemT4>(2)->test(myMesh));
		assert(elem.getMyClassChild<ElemT4>(3)->test(myMesh));
		assert(elem.getMyClassChild<ElemT4>(4)->test(myMesh));
		assert(elem.getMyClassChild<ElemT4>(5)->test(myMesh));
		assert(elem.getMyClassChild<ElemT4>(6)->test(myMesh));
		assert(elem.getMyClassChild<ElemT4>(7)->test(myMesh));
		++(m.elements_.dividedObjs_);
	}
	return offset;
}

//ID	 ElemT4Space::components(const hObj & obj, const int i) { return obj.components(i); }
void ElemT4Space::mark2Deref(hHybridMesh * myMesh, hObj & obj) {
	if (obj.type_ == ElemT4::myType) {
		ElemT4D & elem(*static_cast<ElemT4D *>(&obj));
		hHybridMesh & m(*myMesh);
		assert(elem.isBroken());
		assert(elem.nMyClassSons_ > 1);
		// deref is possible if all neighbouring elems has (their level <= this elem level)
		bool allowed(true);
		for (uTind i(0); i < elem.typeSpecyfic_.nComponents_ && allowed; ++i) {
			if (m.elements_(m.elemNeigh(elem, i)).level_ > elem.level_) {
				allowed = false;
			}
		}
		if (allowed) {
			elem.nMyClassSons_ = hObj::derefMark;// here we lose explicit information abount sons no
			// but this information should be stored outside class objects - it's shared by many of them
		}
	}
}
void ElemT4Space::mark2Delete(hHybridMesh * myMesh, hObj & obj) {

	ElemT4D & elem(*static_cast<ElemT4D *>(&obj));
	if (elem.nMyClassSons_ == hObj::nothingMark) {
		elem.nMyClassSons_ = hObj::delMark;
		myMesh->elements_.requestChange(-1, -elem.mySize_);

		for (int i = 0; i < obj.typeSpecyfic_.nComponents_; ++i) {
			myMesh->faces_.getById(obj.components(i)).decRef(myMesh);
		}

	} else {
		mf_log_err("Deleting parent element with non-deleted children!");
	}

}
void ElemT4Space::derefine(hHybridMesh * myMesh, hObj & obj) {
	mmv_out_stream << "\nError: ElemT4::mark2Delete done something terribly wrong!\n";
}
bool ElemT4Space::test(const hHybridMesh * myMesh, const hObj & el) {
	//const ElemT4 * elem( reinterpret_cast<const ElemT4 *>(obj) );
	//elem.checkPointers();
	//assert( *m.faces_.getById<Face3>(elem.components(0)) < *m.faces_.getById<Face3>(elem.components(1)));
	//assert( *m.faces_.getById<Face3>(elem.components(1)) < *m.faces_.getById<Face3>(elem.components(2)));
	//assert( *m.faces_.getById<Face3>(elem.components(2)) < *m.faces_.getById<Face3>(elem.components(3)));
	const ElemT4 & elem(*static_cast<const ElemT4 *>(&el));
	const hHybridMesh & m(*myMesh);
	if (! (el.nMyClassSons_ == hObj::delMark && m.isRefining()) ) {
		assert(elem.verts(0) != elem.verts(1));
		assert(elem.verts(0) != elem.verts(2));
		assert(elem.verts(0) != elem.verts(3));
		assert(elem.verts(1) != elem.verts(2));
		assert(elem.verts(1) != elem.verts(3));
		assert(elem.verts(2) != elem.verts(3));

		assert( m.faces_.getById<Face3>(elem.components(0)).test(myMesh) );
		assert( m.faces_.getById<Face3>(elem.components(1)).test(myMesh) );
		assert( m.faces_.getById<Face3>(elem.components(2)).test(myMesh) );
		assert( m.faces_.getById<Face3>(elem.components(3)).test(myMesh) );

		assert( elem.components(0) == m.face(elem.verts(0), elem.verts(1), elem.verts(2)).id_);
		assert( elem.components(1) == m.face(elem.verts(0), elem.verts(1), elem.verts(3)).id_);
		assert( elem.components(2) == m.face(elem.verts(0), elem.verts(2), elem.verts(3)).id_);
		assert( elem.components(3) == m.face(elem.verts(1), elem.verts(2), elem.verts(3)).id_);


	}

	assert(m.elemNeigh(elem, 0) != UNKNOWN);
	assert(m.elemNeigh(elem, 1) != UNKNOWN);
	assert(m.elemNeigh(elem, 2) != UNKNOWN);
	assert(m.elemNeigh(elem, 3) != UNKNOWN);

	return true;
}

//int    testElemT4(hObj * obj) const
//{
//    int err(0);
//
//    // face existing
//    for(int i(0); i < faces__no; ++i)
//    {
//        if(_faces_[i] > 0) ++err;
//        err += myMesh->faces_.get<Face3>(_faces_[i])->test();
//    }
//
//    // edge consistency
//    if(err == 0)
//    {
//        err+=myMesh->edges_.get(edge(0))->test();
//        err+=myMesh->edges_.get(edge(1))->test();
//        err+=myMesh->edges_.get(edge(2))->test();
//        err+=myMesh->edges_.get(edge(3))->test();
//        err+=myMesh->edges_.get(edge(4))->test();
//        err+=myMesh->edges_.get(edge(5))->test();
//
//        if( myMesh->faces_.get<Face3>(_faces_[0])->_edges_[0] != myMesh->faces_.get<Face3>(_faces_[1])->_edges_[0] ) ++err; // e0
//        if( myMesh->faces_.get<Face3>(_faces_[0])->_edges_[1] != myMesh->faces_.get<Face3>(_faces_[2])->_edges_[0] ) ++err; // e1
//        if( myMesh->faces_.get<Face3>(_faces_[0])->_edges_[2] != myMesh->faces_.get<Face3>(_faces_[3])->_edges_[0] ) ++err; // e2
//        if( myMesh->faces_.get<Face3>(_faces_[1])->_edges_[1] != myMesh->faces_.get<Face3>(_faces_[2])->_edges_[2] ) ++err; // e3
//        if( myMesh->faces_.get<Face3>(_faces_[1])->_edges_[2] != myMesh->faces_.get<Face3>(_faces_[3])->_edges_[2] ) ++err; // e4
//        if( myMesh->faces_.get<Face3>(_faces_[2])->_edges_[2] != myMesh->faces_.get<Face3>(_faces_[3])->_edges_[2] ) ++err; // e5
//
//        if( myMesh->faces_.get<Face3>(_faces_[0])->_edges_[0] != edge(0) ) ++err; // e0
//        if( myMesh->faces_.get<Face3>(_faces_[0])->_edges_[1] != edge(1) ) ++err; // e1
//        if( myMesh->faces_.get<Face3>(_faces_[0])->_edges_[2] != edge(2) ) ++err; // e2
//        if( myMesh->faces_.get<Face3>(_faces_[1])->_edges_[1] != edge(3) ) ++err; // e3
//        if( myMesh->faces_.get<Face3>(_faces_[1])->_edges_[2] != edge(4) ) ++err; // e4
//        if( myMesh->faces_.get<Face3>(_faces_[2])->_edges_[2] != edge(5) ) ++err; // e5
//    }
//
//
//    // vertex consistency
//    if(err == 0)
//    {
//        if(myMesh->edges_.get(edge(0))->vert[0] != myMesh->edges_.get(edge(1))->vert[0]) ++err; // v0
//        if(myMesh->edges_.get(edge(0))->vert[1] != myMesh->edges_.get(edge(2))->vert[0]) ++err; // v1
//        if(myMesh->edges_.get(edge(1))->vert[1] != myMesh->edges_.get(edge(2))->vert[1]) ++err; // v2
//        if(myMesh->edges_.get(edge(4))->vert[1] != myMesh->edges_.get(edge(5))->vert[1]) ++err; // v3
//
//        for(int i(0),j(0); i < verts_no; ++i)
//        {
//            for(;j < verts_no; ++j)
//            {
//                if((i!=j) && (vertex(i) == vertex(j) )) ++err;
//            }
//        }
//    }
//
//    return err;
//}
