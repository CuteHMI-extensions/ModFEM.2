#include "ElemPrism.h"
#include "hHybridMesh.h"
#include "mmh_vec3.h"


template<>
void ElemPrism::init(hHybridMesh* myMesh){
  if(myMesh != NULL) {
	const Vertex* v[4]={&myMesh->vertices_[verts(0)],
						&myMesh->vertices_[verts(1)],
						&myMesh->vertices_[verts(2)],
						&myMesh->vertices_[verts(3)]};
	const double v0[3]={v[1]->coords_[0]-v[0]->coords_[0],
						v[1]->coords_[1]-v[0]->coords_[1],
						v[1]->coords_[2]-v[0]->coords_[2]};
	const double v1[3]={v[2]->coords_[0]-v[0]->coords_[0],
						v[2]->coords_[1]-v[0]->coords_[1],
						v[2]->coords_[2]-v[0]->coords_[2]};
	const double v2[3]={v[3]->coords_[0]-v[0]->coords_[0],
						v[3]->coords_[1]-v[0]->coords_[1],
						v[3]->coords_[2]-v[0]->coords_[2]};
	const double volume = mmr_vec3_mxpr(v0,v1,v2)/2.0; // /2.0 for prism

	if(volume <= 0.0) { // swapping top and bottom faces
	  swapVerts(0,3);
	  swapVerts(1,4);
	  swapVerts(2,5);
	}
	
	components(0)=myMesh->face(verts(0),verts(1),verts(2)).id_;
	components(1)=myMesh->face(verts(3),verts(4),verts(5)).id_;
	components(2)=myMesh->face(verts(0),verts(1),verts(3),verts(4)).id_;
	components(3)=myMesh->face(verts(0),verts(2),verts(3),verts(5)).id_;
	components(4)=myMesh->face(verts(1),verts(2),verts(4),verts(5)).id_;

    for(int i=0; i < typeSpecyfic_.nComponents_; ++i) {
        Face& f = myMesh->faces_.getById(components(i));
        f.incRef();

        if(f.neighs(0) == UNKNOWN
        || f.neighs(0) == B_COND) {
            f.neighs(0) = id_;
        }
        else {
           f.neighs(1) = id_;
        }
    }
  }
}

void    ElemPrismSpace::mark2Ref(hHybridMesh* myMesh,hObj& elem,const int)
{
	if(elem.type_ == ElemPrism::myType)
	{
	  if(!elem.isBroken()) {
		elem.nMyClassSons_ = hObj::fullRefMark;
		// NOTE: we know that isBroken() == false
		//ElemPrism   * const ptr(static_cast<ElemPrism*>(obj));

        hHybridMesh & m(*myMesh);
		
		// Cheking if that forces neighbours to also break
		for(uTind i(0); i < elem.typeSpecyfic_.nComponents_; ++i) {
		  if(m.elemNeigh(elem,i) != 0) {
			if(m.elements_(m.elemNeigh(elem,i)).level_ < elem.level_) {
              m.elements_(m.elemNeigh(elem,i)).mark2Ref(myMesh);
			}
		  }
		}

		
		if(elem.nMyClassSons_==hObj::fullRefMark)
		  {
			elem.nOtherSons_=13;
			for(uTind i(0); i < elem.typeSpecyfic_.nComponents_; ++i)
			  {
                myMesh->faces_.getById(elem.components(i)).mark2Ref(myMesh);
				//ptr->myMesh->faces_.mark2ref(hObj::posFromId(ptr->components(i)));
			  }
            myMesh->elements_.requestChange(8, 8*sizeof(ElemPrism)+ (sizeof(ElemPrismD)-sizeof(ElemPrism)));
            myMesh->faces_.requestChange(6+4, 6*sizeof(Face4)+ 4*sizeof(Face3) );
            myMesh->edges_.requestChange(3, 3*sizeof(Edge));
		  }
	  }
	}
	else
		throw "markElemPrism(hObj * obj): arg obj is not an ElemPrism.";
}

bool    ElemPrismSpace::test(const hHybridMesh* myMesh,const hObj & elem)
{
    const hHybridMesh & m(*myMesh);
    if( !(elem.nMyClassSons_==hObj::delMark && m.isRefining()) ) {
        //const ElemPrismD * const elem(static_cast<const ElemPrismD*>(obj));
        //elem.checkPointers();
        assert( myMesh->faces_.getById<Face3>(elem.components(0)).type_ == Face3::myType );
        assert( myMesh->faces_.getById<Face3>(elem.components(1)).type_ == Face3::myType );
        assert( myMesh->faces_.getById<Face4>(elem.components(2)).type_ == Face4::myType );
        assert( myMesh->faces_.getById<Face4>(elem.components(3)).type_ == Face4::myType );
        assert( myMesh->faces_.getById<Face4>(elem.components(4)).type_ == Face4::myType );

#ifdef PRINT_LISTS
        out_stream << "\n Elem "<<elem.pos_<<" vertices: ";
        for(int i(0); i < 6; ++i)
        {
            out_stream << elem.myMesh->elemVertex(elem,i) << ",";
        }

#endif

        assert( elem.verts(0) !=  elem.verts(1) );
        assert( elem.verts(1) !=  elem.verts(2) );
        assert( elem.verts(2) !=  elem.verts(3) );
        assert( elem.verts(3) !=  elem.verts(4) );
        assert( elem.verts(4) !=  elem.verts(5) );
        //assert( elem.myMesh->elemVertex(elem,1) !=  elem.myMesh->elemVertex(elem,2) );
        //assert( elem.myMesh->elemVertex(elem,2) !=  elem.myMesh->elemVertex(elem,3) );
        //assert( elem.myMesh->elemVertex(elem,3) !=  elem.myMesh->elemVertex(elem,4) );
        //assert( elem.myMesh->elemVertex(elem,4) !=  elem.myMesh->elemVertex(elem,5) );

        // main assumptions about vertices organization
        assert( elem.components(0) == m.face(elem.verts(0),elem.verts(1),elem.verts(2),UNKNOWN).id_ );
        assert( elem.components(1) == m.face(elem.verts(3),elem.verts(4),elem.verts(5),UNKNOWN).id_ );
        assert( elem.components(2) == m.face(elem.verts(0),elem.verts(1),elem.verts(3),elem.verts(4)).id_);
        assert( elem.components(3) == m.face(elem.verts(0),elem.verts(2),elem.verts(3),elem.verts(5)).id_);
        assert( elem.components(4) == m.face(elem.verts(1),elem.verts(2),elem.verts(4),elem.verts(5)).id_);
    }
	return true;
}

///
///    3---------5
///   /    __---^|
///  /___---     |
/// 4  |         |
/// |  |         |
/// |  |         |
/// |  |         |
/// |  |         |
/// |  0---------2
/// | /    __---^
/// |/___---
/// 1
///
int	ElemPrismSpace::refine(hHybridMesh* myMesh,hObj& elem,const int elem_sons)
{
	int offset(sizeof(ElemPrismD) - sizeof(ElemPrism));
	//ElemPrismD * const elem(static_cast<ElemPrismD*>(obj));
	elem.updatePointers();
	assert(elem.type_ == ElemPrismD::myType);

	if(elem.isBroken() == false)
	{
        hHybridMesh & m(*myMesh);
		assert(elem_sons == 2 || elem_sons == 8);   // slice and full breake
		static int nFailed(0);
		try
		  {
		  elem.mySize_= sizeof(ElemPrismD);
		  elem.nMyClassSons_ = elem_sons;

		  // slice break: 2 new elems + 1 new face
		  /// regular break: 8 new elems + (3*2+4=)10 new faces_ + 3 edges_

		  assert( elem.nOtherSons_ == 13);
		  elem.nMyClassSons_ = 8;
		  
		  // collecting data
		  const uTind v[18]={
		    // main ELEM verts
		    elem.verts(0),elem.verts(1),elem.verts(2),
		    elem.verts(3),elem.verts(4),elem.verts(5),
		    // middle height  triangle verts
		    hObj::posFromId(m.edge(elem.verts(0),elem.verts(3)).sons(0)),
		    hObj::posFromId(m.edge(elem.verts(1),elem.verts(4)).sons(0)),
		    hObj::posFromId(m.edge(elem.verts(2),elem.verts(5)).sons(0)),
		    // face0 (base) mid-edge verts
		    hObj::posFromId(m.edge(elem.verts(0),elem.verts(1)).sons(0)),
		    hObj::posFromId(m.edge(elem.verts(0),elem.verts(2)).sons(0)),
		    hObj::posFromId(m.edge(elem.verts(1),elem.verts(2)).sons(0)),
		    // face1 (top) mid-edge verts
		    hObj::posFromId(m.edge(elem.verts(3),elem.verts(4)).sons(0)),
		    hObj::posFromId(m.edge(elem.verts(3),elem.verts(5)).sons(0)),
		    hObj::posFromId(m.edge(elem.verts(4),elem.verts(5)).sons(0)),
		    // new triangle mid-edge verts
		    hObj::posFromId(m.faces_.getById<Face4D>(elem.components(2)).sons(0)),
		    hObj::posFromId(m.faces_.getById<Face4D>(elem.components(3)).sons(0)),
		    hObj::posFromId(m.faces_.getById<Face4D>(elem.components(4)).sons(0))
		  };

		  // edges
		  const static uTind comps_for_edge[3][2]={{2,3},{2,4},{3,4}};
		  uTind	vertices[3]={UNKNOWN},components[3]={UNKNOWN};
		  uTind i(0);
		  for(; i < 3; ++i)
			{
              components[0]=myMesh->faces_.getById<Face4D>(elem.components(comps_for_edge[i][0])).sons(0);
              components[1]=myMesh->faces_.getById<Face4D>(elem.components(comps_for_edge[i][1])).sons(0);
			  vertices[0]= hObj::posFromId(components[0]);
			  vertices[1]= hObj::posFromId(components[1]);

              Edge & edge = *myMesh->edges_.newObj<Edge>(myMesh,vertices);
			  edge.parent_= elem.id_;
			  edge.level_= elem.level_+1;
//			  edge.status(INACTIVE);
			  elem.sons(i) = edge.id_;
			  edge.components(0)=components[0];
			  edge.components(1)=components[1];
              assert(edge.test(myMesh));
			}

		  // FACES3
		  const uTind verts_f3[4][4]=
			{
			  {v[6],v[15],v[16],UNKNOWN},
			  {v[15],v[7],v[17],UNKNOWN},
			  {v[16],v[17],v[8],UNKNOWN},
			  {v[15],v[16],v[17],UNKNOWN}
			};
		  
		  for(int f(0);i<7; ++i,++f)
			{
              Face3 & child3( *myMesh->faces_.newObj<Face3>(myMesh,verts_f3[f]) );	// triangle faces_
			  child3.parent_= elem.id_;
			  child3.level_= elem.level_+1;
//			  child3.status(INACTIVE);
			  elem.sons(i) = child3.id_;
			  child3.components(0)=m.edge(verts_f3[f][0],verts_f3[f][1]).id_;
			  child3.components(1)=m.edge(verts_f3[f][0],verts_f3[f][2]).id_;
			  child3.components(2)=m.edge(verts_f3[f][1],verts_f3[f][2]).id_;
			  //assert(child3.test());
			}

		  // 2*3 face4
		  const uTind vert_f4[6][4]={
			// bottom faces
			{v[9],v[10],v[15],v[16]},
			{v[9],v[11],v[15],v[17]},
			{v[10],v[11],v[16],v[17]},
			// upper faces
			{v[12],v[13],v[15],v[16]},
			{v[12],v[14],v[15],v[17]},
			{v[13],v[14],v[16],v[17]}};
		  
		  Face4 * child4;
		  for(int f4(0); i < elem.typeSpecyfic_.nSons_; ++i,++f4)
			{
              child4 = myMesh->faces_.newObj<Face4>(myMesh,vert_f4[f4]);
			  child4->parent_= elem.id_;
			  child4->level_= elem.level_+1;
//			  child4->status(INACTIVE);
			  child4->components(0)=m.edge(vert_f4[f4][0],vert_f4[f4][1]).id_;
			  child4->components(1)=m.edge(vert_f4[f4][0],vert_f4[f4][2]).id_;
			  child4->components(2)=m.edge(vert_f4[f4][1],vert_f4[f4][3]).id_;
			  child4->components(3)=m.edge(vert_f4[f4][2],vert_f4[f4][3]).id_;
			  elem.sons(i) = child4->id_;
			}
		  
		  
		  // my class sons
		  
		  
		  const uTind vert_prism[8][6]={
			// lower part
			{v[0],v[9],v[10],v[6],v[15],v[16]},
			{v[9],v[1],v[11],v[15],v[7],v[17]},
			{v[10],v[11],v[2],v[16],v[17],v[8]},
			{v[9],v[11],v[10],v[15],v[17],v[16]},
			// upper part
			{v[6],v[15],v[16],v[3],v[12],v[13]},
			{v[15],v[7],v[17],v[12],v[4],v[14]},
			{v[16],v[17],v[8],v[13],v[14],v[5]},
			{v[15],v[17],v[16],v[12],v[14],v[13]}
		  };
		  
		  ElemPrism  * childPrism(NULL);
		  BYTE * adr( reinterpret_cast<BYTE *>(&elem)+elem.mySize_);
		  for(uTind ch(0); ch < 8; ++ch)
			{
              childPrism = myMesh->elements_.newObj<ElemPrism>(myMesh,vert_prism[ch],adr);
			  childPrism->parent_= elem.id_;
			  childPrism->level_= elem.level_+1;
			  childPrism->flags(GROUP_ID) = elem.flags(GROUP_ID);
			  childPrism->components(0)=m.face(vert_prism[ch][0],vert_prism[ch][1],vert_prism[ch][2]).id_;
			  childPrism->components(1)=m.face(vert_prism[ch][3],vert_prism[ch][4],vert_prism[ch][5]).id_;
			  childPrism->components(2)=m.face(vert_prism[ch][0],vert_prism[ch][1],vert_prism[ch][3],vert_prism[ch][4]).id_;
			  childPrism->components(3)=m.face(vert_prism[ch][0],vert_prism[ch][2],vert_prism[ch][3],vert_prism[ch][5]).id_;
			  childPrism->components(4)=m.face(vert_prism[ch][1],vert_prism[ch][2],vert_prism[ch][4],vert_prism[ch][5]).id_;
			  assert(childPrism->mySize_ == sizeof(ElemPrism));
			  m.normalizeElem(*childPrism);
			  adr += childPrism->mySize_;
			}

		  // setup neighbours for faces (internal)
		  static const int f_sons_neigs[10][2]={
			// face3 neighs
			{0,4},{1,5},{2,6},{3,7},
			// face 4 (bottom)
			{0,3},{1,3},{2,3},
			// face 4 (top)
			{4,7},{5,7},{6,7}};
		  for(int i(0); i < 10; ++i) {
		    hObj & iFace(m.faces_(elem.sons(3+i)));
		    iFace.neighs(0)=elem.getMyClassChild<ElemPrism>(f_sons_neigs[i][0])->id_;
		    iFace.neighs(1)=elem.getMyClassChild<ElemPrism>(f_sons_neigs[i][1])->id_;
            assert(iFace.test(myMesh));
		  }
		  
//		  // setup neigbours for faces (external)
//		  // 4 (face verts), 1 prism_neigh_id
//		  const int nExternFaces=20;
//		  const uTind externFaces[nExternFaces][5]={
//		    {v[0],v[9],v[10],UNKNOWN,0},
//		    {v[9],v[1],v[11],UNKNOWN,1},
//		    {v[10],v[11],v[2],UNKNOWN,2},
//		    {v[9],v[10],v[11],UNKNOWN,3},
		    
//		    {v[3],v[12],v[13],UNKNOWN,4},
//		    {v[12],v[4],v[14],UNKNOWN,5},
//		    {v[13],v[14],v[5],UNKNOWN,6},
//		    {v[12],v[13],v[14],UNKNOWN,7},
		    
//		    {v[0],v[9],v[6],v[15],0},
//		    {v[9],v[1],v[15],v[7],1},
//		    {v[6],v[15],v[3],v[12],4},
//		    {v[15],v[7],v[12],v[4],5},
		    
//		    {v[0],v[10],v[6],v[16],0},
//		    {v[10],v[2],v[16],v[8],2},
//		    {v[6],v[16],v[3],v[13],4},
//		    {v[16],v[8],v[13],v[5],6},
		    
//		    {v[1],v[11],v[7],v[17],1},
//		    {v[11],v[2],v[17],v[8],2},
//		    {v[7],v[17],v[4],v[14],5},
//		    {v[17],v[8],v[14],v[5],6}
//		  };
// This is done inside ElemPrism::init()
//		  for(int f(0); f < nExternFaces; ++f) {
//			hObj & face(m.face(externFaces[f]));
//			if(face.neighs(0) == 0) {
//			  face.neighs(0) = elem.getMyClassChild<ElemPrism>(externFaces[f][4])->id_;
//			}
//			else if(face.neighs(1) == 0) {
//			  face.neighs(1) = elem.getMyClassChild<ElemPrism>(externFaces[f][4])->id_;
//			}
//			else {
//			  assert(!"This should not happed!");
//			}
//			assert(face.neighs(0) != 0);
//		  }
		  
		  offset += 8*sizeof(ElemPrism);
          ++(myMesh->elements_.dividedObjs_);
		}
		catch(char const * e) {		  
		  ++nFailed;
		  elem.mySize_= sizeof(ElemPrism);
		  elem.nMyClassSons_ = 0;
		  elem.nOtherSons_ = 0;
		  offset=0;
		  if(nFailed > 0) {
			std::cerr << "\r" << nFailed  << " ElemPrisms break failed!" << e;
		  }
		}
	}
	return offset;
}

//ID	 ElemPrismSpace::components(const hObj & obj, const int i) { return obj.components(i); }
void ElemPrismSpace::mark2Deref(hHybridMesh* myMesh,hObj& obj) {
    if(obj.type_ == ElemPrism::myType) {
      ElemPrismD & elem(*static_cast<ElemPrismD*>(&obj));
      hHybridMesh & m(*myMesh);
      assert(elem.isBroken());
      assert(elem.nMyClassSons_ > 1);
      // deref is possible if all neighbouring elems has (their level <= this elem level)
      bool allowed(true);
      for(uTind i(0); i < elem.typeSpecyfic_.nComponents_ && allowed; ++i) {
        if(m.elements_(m.elemNeigh(elem,i)).level_ > elem.level_) {
          allowed=false;
        }
      }
      if(allowed) {
        elem.nMyClassSons_ = hObj::derefMark;// here we lose explicit information abount sons no
        // but this information should be stored outside class objects - it's shared by many of them
      }
    }
  }

void ElemPrismSpace::mark2Delete(hHybridMesh* myMesh,hObj& obj) {
    ElemPrismD & elem(*static_cast<ElemPrismD*>(&obj));
    if(elem.nMyClassSons_ == hObj::nothingMark) {
        elem.nMyClassSons_ = hObj::delMark;
        myMesh->elements_.requestChange(-1,-elem.mySize_);

        for(int i=0; i < obj.typeSpecyfic_.nComponents_; ++i) {
            myMesh->faces_.getById(obj.components(i)).decRef(myMesh);
        }

    } else {
        mf_log_err("Deleting parent element with non-deleted children!");
    }

}

void ElemPrismSpace::derefine(hHybridMesh* myMesh,hObj& obj) {
    mmv_out_stream << "\nError: ElemPrism::mark2Delete done something terribly wrong!\n";
  }
