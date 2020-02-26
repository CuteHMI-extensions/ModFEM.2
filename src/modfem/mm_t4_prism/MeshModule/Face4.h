#ifndef _FACE4_H_
#define _FACE4_H_

#include "Edge.h"
/// faces_
/// components = edges_
/// params_[0-1]	= neighbours IDs
/// params2_[B_COND] = BC

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


namespace Face4Space{
    void mark2Ref(hHybridMesh* myMesh,hObj& obj,const int i) ;
    void mark2Deref(hHybridMesh* myMesh,hObj& obj);
    void mark2Delete(hHybridMesh* myMesh,hObj& obj);
    int  refine(hHybridMesh* myMesh,hObj& obj,const int i);
    void derefine(hHybridMesh* myMesh,hObj& obj);
    bool test(const hHybridMesh* myMesh,const hObj& obj);

	MMT_H_MESH_TYPE(Face4,4,4,Edge,4,2,2,5,4,0,0)

	//const EntityAttributes shared(4,4,4,2,2,5,4,0,0,
	//		&mark2Ref,&mark2Deref,&mark2Delete,&refine,&derefine,&test,&uniqueId);
};

typedef Face4Space::Face4       Face4;
typedef Face4Space::Face4D     Face4D;
/**  @} */
#endif // _FACE4_H_
