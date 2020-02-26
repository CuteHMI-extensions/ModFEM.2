#include "hObj.h"
#include "ElemPrism.h"
#include "ElemT4.h"

//#define ENTITY(NO,CLASS,SPACE) EntityAttributes(CLASS::myType, 
//CLASS::nComponents, CLASS::nSons, 
//CLASS::nParams, CLASS::nParams2, 
//SPACE::nVertices, SPACE::nEdges, SPACE::nFaces, 
//&SPACE::mark2Ref, &SPACE::mark2Deref, 
//&SPACE::mark2Delete, &SPACE::refine, 
//	&SPACE::derefine, &SPACE::test, &SPACE::uniqueId, &SPACE::components )

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


const EntityAttributes hObj::eTable[]={
	VertexSpace::shared,
	EdgeSpace::shared,
	EntityAttributes(),
	Face3Space::shared,
	Face4Space::shared,
	ElemPrismSpace::shared,
	EntityAttributes(),
	ElemT4Space::shared
};

/**  @} */
