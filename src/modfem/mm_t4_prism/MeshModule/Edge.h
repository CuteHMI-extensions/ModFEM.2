#ifndef _EDGE_H_
#define _EDGE_H_

#include "Vertex.h"

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


namespace EdgeSpace{
    void mark2Ref(hHybridMesh *myMesh, hObj& obj, const int i) ;
    void mark2Deref(hHybridMesh *myMesh,hObj& obj);
    void mark2Delete(hHybridMesh *myMesh,hObj& obj);
    int  refine(hHybridMesh *myMesh,hObj& obj,const int i);
    void derefine(hHybridMesh *myMesh,hObj& obj);
    bool test(const hHybridMesh *myMesh,const hObj& obj);

	MMT_H_MESH_TYPE(Edge2,1,2,Vertex,2,1,0,1,0,0,0)
	//const EntityAttributes shared(1,2,2,1,0,1,0,0,0,
	//	&mark2Ref,&mark2Deref,&mark2Delete,&refine,&derefine,&test,&uniqueId);
};



//typedef hParent<1,2,2,Vertex,1>     Edge;
//typedef hParent<1,2,2,Vertex,1,1>   EdgeD;	// in 3D divided edge stores info about how many times it was requested to split

typedef EdgeSpace::Edge2 Edge;
typedef EdgeSpace::Edge2D EdgeD;
/*
class Edge : public EdgeSpace::Edge2 {
public:
    Edge(const uTind vertices[]): EdgeSpace::Edge2(vertices){ }
};

class EdgeD : public EdgeSpace::Edge2D {
public:
    EdgeD(const uTind vertices[]) : EdgeSpace::Edge2D(vertices){}
};
*/
/**  @} */
#endif // _EDGE_H_
