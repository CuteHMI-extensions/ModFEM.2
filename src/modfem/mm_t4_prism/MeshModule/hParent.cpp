#include "hParent.hpp"
#include "hHybridMesh.h"

/*template<int I, int nTVerts, int nTComponents, typename TComponents,
int nTFlags, int nTSons, int nTCoords, int nTNeighs>
void	hParent<I,nTVerts,nTComponents,TComponents,nTFlags,nTSons,nTCoords,nTNeighs>::updateHash(const uTind vertices[])
{
    switch(type_)
    {
	case 1:
	case 2: myMesh->registerEdge(vertices,pos_); break;
	case 3:
	case 4: myMesh->registerFace(vertices,pos_); break;
    }
}*/

void EmptyHParentSpace::mark2Ref(hHybridMesh *myMesh, hObj* , const int ){ throw "Not implemented";};
void EmptyHParentSpace::mark2Deref(hHybridMesh* myMesh,hObj* ){ throw "Not implemented";};
void EmptyHParentSpace::mark2Delete(hHybridMesh* myMesh,hObj* ){ throw "Not implemented";};
int  EmptyHParentSpace::refine(hHybridMesh* myMesh,hObj* ,const int ){ throw "Not implemented";return 0;};
void EmptyHParentSpace::derefine(hHybridMesh* myMesh,hObj* ){ throw "Not implemented";};
bool EmptyHParentSpace::test(const hHybridMesh* myMesh,const hObj* ){ throw "Not implemented";return true;};
ID	 EmptyHParentSpace::components(hHybridMesh* myMesh, const hObj* obj,const int i) { throw "Not implemented";return 0;};
