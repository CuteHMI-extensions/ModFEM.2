#ifndef _HPARENT_HPP_
#define _HPARENT_HPP_

#include "../Common.h"
#include "hObj.h"
#include "../Field.hpp"
#include <algorithm>

class hHybridMesh;


using Memory::Field;

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <int I, int nTVerts, int nTComponents, class TComponents = BYTE, 
int nTFlags=0,int nTSons=0,int nTCoords=0, int nTNeighs=0>
class hParent : public hObj
{
public:
	static const uTind myType = I;
	static const uTind nVerts = nTVerts;
	static const uTind nComponents = nTComponents;
	static const uTind nSons = nTSons;
	static const uTind nFlags = nTFlags;
	static const uTind nNeighs = nTNeighs;
	static const uTind nNodes = nVerts+nComponents+nFlags+nNeighs+nSons;
	static const uTind nCoords = nTCoords;

  hParent(hHybridMesh* myMesh, const uTind vertices[],const uTind posID)
  : hObj(posID, eTable[myType])
	{ 
		// init & clear
		nodes_=nodesArr_; 
		mySize_ = sizeof( *this ); 
		clear(); 
		// set vertices
		if(vertices != NULL)
		{
            mf_check_debug(std::count(vertices,vertices+nVerts,0)==0, "Incorrect vertex number during elem creation!");
			memcpy(nodesArr_,vertices,sizeof(uTind)*typeSpecyfic_.nVerts_);
			/* updateHash(vertices);
			switch(type_)
			{
				case 1:
				case 2: myMesh->registerEdge(vertices,pos_); break;
				case 3:
				case 4: myMesh->registerFace(vertices,pos_); break;
				}*/
		}
        init(myMesh);
	}

  
  ~hParent(){}
  
  void init(hHybridMesh *myMesh);
  
	//void	updatePointers()
	//{
	//	nodes_ = nodesArr_;
	//}

	void    setComponentsPerPos(const Tind pos[])
	{
		for(uTind i(0); i < nComponents; ++i)
		{
			components(i) = hObj::makeId(pos[i],TComponents::myType);	
		}
	}

	bool	checkPointers() const
	{
		assert(nodes_ == nodesArr_);
		return (nodes_ == nodesArr_);
	}

	//void* operator new (std::size_t size) throw (std::bad_alloc);
    void print() const
    {
        hObj::print();
    }

private:
	ID	nodesArr_[nNodes];

	void    clear()
	{
		memset(nodesArr_,0,sizeof(ID)*nNodes);
	}

	//void updateHash(const uTind vertices[]);
	// forbidden operations 
	hParent(const hParent & other);
	hParent & operator=(const hParent & other);



public:
	Field<nCoords,double>	coords_;
};

//#include "hHybridMesh.h"

//template <int I, int nTVerts, int nTComponents, typename TComponents, 
//int nTFlags,int nTSons,int nTCoords, int nTNeighs>
//void	hParent<I, nTVerts, nTComponents,TComponents,nTFlags,nTSons,nTCoords, nTNeighs>::updateHash(const uTind vertices[])
//{
//	switch(type_)
//	{
//		case 1:
//		case 2: myMesh->registerEdge(vertices,pos_); break;
//		case 3:
//		case 4: myMesh->registerFace(vertices,pos_); break;
//	}
//}

typedef hParent<0,0,0> EmptyHParent;


namespace EmptyHParentSpace
{
    void mark2Ref(hHybridMesh* myMesh,hObj* ,const int );
    void mark2Deref(hHybridMesh* myMesh,hObj* );
    void mark2Delete(hHybridMesh* myMesh,hObj* );
    int  refine(hHybridMesh* myMesh,hObj* ,const int );
    void derefine(hHybridMesh* myMesh,hObj* );
    bool test(const hHybridMesh* myMesh,const hObj* );
    ID	 components(hHybridMesh* myMesh,const hObj* obj,const int i);

	const uTind nVertices=0;
	const uTind nEdges=0;
	const uTind nFaces=0;
};

#define MMT_H_MESH_TYPE(NAME,TYPE_ID,N_COMP,COMP_TYPE,N_VERT,N_FLAGS,N_NEIGS,N_SONS,N_EDGES,N_FACES,N_COORDS) \
const EntityAttributes shared(TYPE_ID,N_COMP,N_VERT,N_FLAGS,N_NEIGS,N_SONS,N_EDGES,N_FACES,N_COORDS, \
&mark2Ref,&mark2Deref,&mark2Delete,&refine,&derefine,&test); \
typedef hParent<TYPE_ID,N_VERT, N_COMP,COMP_TYPE,N_FLAGS,0,N_COORDS,N_NEIGS>     NAME; \
typedef hParent<TYPE_ID,N_VERT, N_COMP,COMP_TYPE,N_FLAGS,N_SONS,N_COORDS,N_NEIGS>     NAME##D;

/**  @} */

#endif // _HPARENT_HPP_
