#ifndef _ELEMPRISM_H_
#define _ELEMPRISM_H_

#include "Face3.h"
#include "Face4.h"

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


namespace ElemPrismSpace{
    void mark2Ref(hHybridMesh* myMesh,hObj& obj,const int i) ;
    void mark2Deref(hHybridMesh* myMesh,hObj& obj);
    void mark2Delete(hHybridMesh* myMesh,hObj& obj);
    int  refine(hHybridMesh* myMesh,hObj& obj,const int i);
    void derefine(hHybridMesh* myMesh,hObj& obj);
    bool test(const hHybridMesh* myMesh,const hObj& obj);

	MMT_H_MESH_TYPE(ElemPrism,5,5,Face3,6,3,0,13,9,5,0)

	//const EntityAttributes shared(5,5,6,2,0,13,9,5,0,
	//		&mark2Ref,&mark2Deref,&mark2Delete,&refine,&derefine,&test,&uniqueId);

const int faceNeigByEdge[5][4]={{2,3,4,0},{2,3,4,0},{0,2,3,1},{0,1,3,1},{0,1,2,1}};

};

typedef ElemPrismSpace::ElemPrism     ElemPrism;
typedef ElemPrismSpace::ElemPrismD  ElemPrismD;
/**  @} */
#endif // _ELEMPRISM_H_
