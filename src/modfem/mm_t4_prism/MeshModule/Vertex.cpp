#include "Vertex.h"
#include "hHybridMesh.h"

template<>
void Vertex::init(hHybridMesh * myMesh){}

template<>
void       Vertex::print() const
{
    hObj::print();
    std::cout << "(";
       for(int i=0; i < typeSpecyfic_.nCoords_;++i) {
           std::cout << coords_[i] << ",";
       }
       std::cout << ")\n";
}


void VertexSpace::mark2Ref(hHybridMesh* myMesh,hObj& ,const int ) {}
void VertexSpace::mark2Deref(hHybridMesh* myMesh,hObj& ){}
void VertexSpace::mark2Delete(hHybridMesh* myMesh,hObj& obj){
    if(obj.nMyClassSons_ != hObj::delMark) {
        obj.nMyClassSons_ = hObj::delMark;
        myMesh->vertices_.requestChange(-1, -sizeof(Vertex));
    }
}
int  VertexSpace::refine(hHybridMesh* myMesh,hObj& ,const int ){ return 0;}
void VertexSpace::derefine(hHybridMesh* myMesh,hObj& ){}
bool VertexSpace::test(const hHybridMesh *myMesh, const hObj& ){ return false;}
//ID	 VertexSpace::components(const hObj *ptr,const int i){assert(!"Not implemented");return 0;}

