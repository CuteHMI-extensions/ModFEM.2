#ifndef H_HYBRID_MESH_H
#define H_HYBRID_MESH_H

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


#include <string.h>
#include <vector>
#include <fstream>
//#include "mmh_vec3.h"

#include "StaticPool.hpp"
//#include "ArrayPool.hpp"
//#include "FixedSizeAllocator.hpp"
//#include "IndexedMemPool.hpp"


#include "../MeshRead/IMeshReader.h"
#include "../MeshWrite/IMeshWriter.h"
#include "../GeometryModule/GeometryModule.hpp"

#include "hObj.h"
#include "Vertex.h"
#include "Edge.h"

#include "ElemPrism.h"

//enum Entities {ELEM=0,FACE,EDGE,VERT,N_ENTITIES};

//struct Representation {

//    bool entityConnectivity[N_ENTITIES][N_ENTITIES] ;

//};

//template<class TRepresentation, class TContainers, class TAllocators>
//class Mesh
//{
//public:
//    typedef Containers::ElementContainer ElementPool;
//    ElementPool elements;
//};

#ifdef _MSC_VER

template<>
void StaticPool<Vertex, 0>::rebuildVts2Pos();

template<>
void StaticPool<Vertex, 0>::checkIfInVts2Pos(const uTind verts[]) const ;

template<>
void StaticPool<Vertex, 0>::insertIntoVts2Pos(const uTind verts[],const int pos);

#endif

class hHybridMesh
{
	static int lastMesh;
public:
	typedef StaticPool<Elem, ElemPrism::nVerts>			ElemPool;
	typedef StaticPool<Face, Face4::nVerts>			FacePool;
	typedef StaticPool<Edge, Edge::nVerts>            EdgePool;
	typedef StaticPool<Vertex, Vertex::nVerts>			VertexPool;
	typedef ElemPool::Iterator<ElemPool::notBroken>          Iterator;
	typedef ElemPool::constIterator<ElemPool::notBroken>    constIterator;
	typedef ElemPool::Iterator<>          allIterator;
	typedef ElemPool::constIterator<>    allConstIterator;
	
	const int    meshId_;

	hHybridMesh(MeshWrite::IMeshWriter * defaultWriter = NULL,
				MeshRead::IMeshReader * defaultReader = NULL);

	~hHybridMesh();

	bool read(MeshRead::IMeshReader & reader);
	bool read(MeshRead::IMeshReader* readers[],const int noReaders);
	bool write(MeshWrite::IMeshWriter & writer) const ;

	bool	free();

    virtual bool normalizationProcessor();
    bool printSetup();
	bool checkUniqueness() const;
	bool checkGeometry() const;
	bool checkTypes() const;
	bool checkAll() const;

	bool createBoundaryLayer(const double thicknessProc=0.1, const int nLayers=1,
	const bool quadraticDistribution=false, const double * vecIgnoreNormal=NULL);

	bool initRefine() ;
	bool refine();
	bool refineElem(hObj & el);
	bool refineElem(const int elemID);
	bool rRefine(const int BC_id, void (*reallocation_func)(double * x,double * y, double * z));
	bool derefine();
	bool derefineElem(hObj & elx);
	bool derefineElem(const int elemID);
	bool finalRef();
    bool isRefining() const { return isRefining_; }

	Edge&	edge(const uTind v0, const uTind v1);
	hObj&	face(const uTind v0, const uTind v1, const uTind v2, const uTind v3=UNKNOWN);
    const hObj& face(const uTind v0, const uTind v1, const uTind v2, const uTind v3=UNKNOWN) const ;
	hObj&	element(const uTind v0, const uTind v1, const uTind v2, const uTind v3, const uTind v4=UNKNOWN, const uTind v5=UNKNOWN);

	Edge&	edge(const uTind * v);
	hObj&	face(const uTind v[4]);
    const hObj&	face(const uTind v[4]) const;

	ID	faceVertex(const hObj & face, const Tind vertNo) const;
	ID	elemVertex(const hObj & elem, const Tind vertNo) const;
	ID  elemNeigh(const hObj & elem, const Tind neighNo) const;
	int whichNeighAmI(const hObj& elem, const hObj& otherElem) const;
	int whichFaceAmI(const hObj & face, const hObj& elem) const;
	int whichNodeAmI(const hObj & node, const hObj& elem) const;
	eElemFlag normalizeElem(hObj & elem);
	eFaceFlag normalizeFace(hObj & face,hObj & elem);
	double elemHSize(const hObj & elem) ;
	// Returns area of the face.
	void faceNormal(IN const hObj& face,OUT double vecNorm[3], OUT double * area) const;
	// Returns indicator whether face normal vectors points inside or outside given element.
	eFaceFlag faceDirection(const hObj& face, const hObj & elem) const;

	Iterator begin() { return elements_.begin(); }
	constIterator begin() const  { return elements_.begin(); }


	size_t	totalSize() const {return sizeof(*this)+vertices_.totalMemory()+edges_.totalMemory()+faces_.totalMemory()+elements_.totalMemory();}
	bool    test() const;
	void    findEdgeElems();
	void computeDist2Bound(const int BCs[],const int nBCs);
	
	friend std::ofstream &	operator << (std::ofstream & stream, const hHybridMesh & mesh);
	
//  FIELDS
	VertexPool  vertices_;
	EdgePool    edges_;
	FacePool    faces_;
	ElemPool    elements_;

	std::string     name_;	// Mesh name.
	int 		gen_,		// current mesh generation (level_od division
				maxGen_,    // maximal generation
				maxGenDiff_,
				maxEdgesPerVertex_,
				maxFacesPerVertex_;

	MeshWrite::IMeshWriter*		defaultWriter_;
	MeshRead::IMeshReader*		defaultReader_;

	std::vector< std::vector<uTind> >	edge_elems;
	bool wasNormalized_;

    void    print() const;
	
    //std::vector<uTind> vts4check2del;

protected:
	hHybridMesh(const hHybridMesh & other);
	hHybridMesh&    operator=(const hHybridMesh & other);

	bool actualRefine();
	bool actualDerefine();
	bool actualDelete();
	bool isRefining_;

#ifdef TURBULENTFLOW
	void computePlaneCoords(Face & face);
	void findNearestBoundFace(Vertex & v,const int BCs[],const int nBCs);
	double dist2Face(const Vertex & v,const Face & f);
#endif
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**  @} */
#endif // H_HYBRID_MESH_H
