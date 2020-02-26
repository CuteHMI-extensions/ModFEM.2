#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "EntityAttributes.hpp"
#include "hParent.hpp"

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


namespace VertexSpace{

    void mark2Ref(hHybridMesh* myMesh,hObj& obj,const int i) ;
    void mark2Deref(hHybridMesh* myMesh,hObj& obj);
    void mark2Delete(hHybridMesh* myMesh,hObj& obj);
    int  refine(hHybridMesh* myMesh,hObj& obj,const int i);
    void derefine(hHybridMesh* myMesh,hObj& obj);
    bool test(const hHybridMesh* myMesh,const hObj& obj);

	MMT_H_MESH_TYPE(Vertex1,0,0,BYTE,0,1,0,0,0,0,3)

	//const EntityAttributes shared(0,0,0,1,0,0,0,0,3,
	//	&mark2Ref,&mark2Deref,&mark2Delete,&refine,&derefine,&test,&uniqueId);
};

//typedef hParent<0,0,0,BYTE,1,0,3>           Vertex;
typedef VertexSpace::Vertex1 Vertex;

/**  @} */
#endif // _VERTEX_H_
