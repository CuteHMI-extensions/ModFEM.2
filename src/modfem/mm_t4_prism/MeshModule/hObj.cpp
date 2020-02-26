#include "hObj.h"
#include "hHybridMesh.h"
#include <algorithm>
//hHybridMesh *  hObj::myMesh = NULL;


hObj::hObj(const uTind posID,const EntityAttributes & typeSpecyfic)
  : pos_(posID),type_(typeSpecyfic.myType),parent_(0), 
	mySize_(sizeof(*this)),level_(0),
    nMyClassSons_(0),nOtherSons_(0), nRefs(0), // /*, status_(1)/*==MMC_ACTIVE*/
#ifdef TURBULENTFLOW
	dist2bound_(0.0),nearstFace_(0),
#endif
	typeSpecyfic_(typeSpecyfic),nodes_(NULL)
{
#ifdef TURBULENTFLOW
  std::fill(planeCoords_,planeCoords_+4,0.0);
#endif
}

hObj::~hObj()
{}

//	This function below is error prone and de facto useless.
//hObj*		hObj::getMyClassChild(const int number)
//{
//	assert( static_cast<int>(nMyClassSons_) > number);
//	return  reinterpret_cast<hObj*>( reinterpret_cast<BYTE*>(this)+mySize_ + sizeof(hObj)*number );
//}
//
//const hObj* hObj::getMyClassChild(const int number) const
//{
//	assert( static_cast<int>(nMyClassSons_) > number);
//	return  reinterpret_cast<const hObj*>( reinterpret_cast<const BYTE*>(this)+mySize_ + sizeof(hObj)*number );
//}

void    hObj::mark2Ref(hHybridMesh *myMesh, const int refType)
{
    return typeSpecyfic_.mark2Ref_(myMesh,*this,refType);
}

void    hObj::mark2Dref(hHybridMesh* myMesh)
{
    typeSpecyfic_.mark2Deref_(myMesh,*this);
}

void    hObj::mark2Del(hHybridMesh* myMesh)
{
    if( nMyClassSons_ == hObj::nothingMark) {
//        mfp_debug("Deleting type %d number %d",hObj::type_,hObj::pos_);
        typeSpecyfic_.mark2Delete_(myMesh,*this);
    }
}

int	hObj::hBreak(hHybridMesh* myMesh,const Tind sonsCount)
{
	assert(this->nMyClassSons_ == hObj::fullRefMark);
    return typeSpecyfic_.refine_(myMesh,*this,sonsCount);
}

void	hObj::derefine(hHybridMesh* myMesh)
{
	assert(isBroken());
	//return typeSpecyfic_.derefine_(*this);
	for(uTind i(0); i < typeSpecyfic_.nSons_; ++i) {
		sons(i)=0;
	}
	nMyClassSons_=0;
}

bool    hObj::test(const hHybridMesh* myMesh) const
{
    bool result=true;
    if(this->nMyClassSons_ != hObj::delMark) {
        for(uTind i(0); i < typeSpecyfic_.nComponents_; ++i) {
            assert(components(i) != 0);
        }

        for(uTind i(0); i < typeSpecyfic_.nVerts_; ++i) {
            assert(verts(i) != 0);
        }
        result = typeSpecyfic_.test_(myMesh,*this);
    }
    return result;

}

void	hObj::updateHash(hHybridMesh * myMesh) const
{
	switch(type_)
	  {
		case 2:
			myMesh->edge(verts(0),verts(1));
		break;
		case 3:
			myMesh->face(verts(0),verts(1),verts(2));
			break;
		case 4:
			myMesh->face(verts(0),verts(1),verts(2),verts(3));
			break;
	}
}

double* hObj::geoCenter(const hHybridMesh * myMesh, IN double center[3]) const
{
    if(center != NULL) {
	center[0]=0.0;
	center[1]=0.0;
	center[2]=0.0;
	for(uTind v(0); v < typeSpecyfic_.nVerts_; ++v) {
	    center[0]+=myMesh->vertices_[verts(v)].coords_[0];
	    center[1]+=myMesh->vertices_[verts(v)].coords_[1];
	    center[2]+=myMesh->vertices_[verts(v)].coords_[2];
	}
	center[0]/=typeSpecyfic_.nVerts_;
	center[1]/=typeSpecyfic_.nVerts_;
	center[2]/=typeSpecyfic_.nVerts_;
    }
    return center;
}


void hObj::swapVerts(const int v1,const int v2) 
{ 
    assert(v1>=0 && v1 < static_cast<int>(typeSpecyfic_.nVerts_));
    assert(v2>=0 && v2 < static_cast<int>(typeSpecyfic_.nVerts_));
    std::swap(nodes_[v1],nodes_[v2]);
}

bool hObj::equals(const hObj& other,bool) const 
{
    int diff(0);
    if(type_ != other.type_) ++diff;
    assert(diff==0);
    if(pos_ != other.pos_) ++diff;
    assert(diff==0);
    if(parent_ != other.parent_) ++diff;
    assert(diff==0);
    if(&typeSpecyfic_ != &other.typeSpecyfic_) ++diff;
    assert(diff==0);
    for(uTind i(0); i < typeSpecyfic_.nVerts_; ++i) {
	if(verts(i) != other.verts(i)) ++diff;
    }
    assert(diff==0);
    for(uTind i(0); i < typeSpecyfic_.nComponents_; ++i) {
	if(components(i) != other.components(i)) ++diff;
    }
    assert(diff==0);
    for(uTind i(0); i < typeSpecyfic_.nNeighs_; ++i) {
	if(neighs(i) != other.neighs(i)) ++diff;
    }
    assert(diff==0);

    
    return diff==0;
}

hObj::operator const EntityAttributes() const
{
  return this->typeSpecyfic_;
}

std::ostream& operator<<(std::ostream & os, const hObj & obj)
	  {
		os << "pos:" << obj.pos_ << " type:" << obj.type_ << " verts:";
		  for(uTind i(0);i < obj.typeSpecyfic_.nVerts_; ++i) {
			os << obj.verts(0) << " ";
		  }
		  os << "\n";
		return os;
	  }

void hObj::print() const {
    std::cout << "ID=[" << type_ << "," << pos_ << "] V[";
    for(int i=0; i < typeSpecyfic_.nVerts_;++i) {
        std::cout << verts(i) << ",";
    }
    std::cout << "] C[";

    for(int i=0; i < typeSpecyfic_.nComponents_;++i) {
        std::cout << "(" << hObj::typeFromId(components(i)) << "," << hObj::posFromId(components(i)) << "),";
    }

    std::cout << "]\n";
}
