#include "../Common.h"
#include "hHybridMesh.h"
#include "hObj.h"
#include "ElemPrism.h"
#include "ElemT4.h"
#include "mmh_vec3.h"
#include "../MeshRead/SimplePointReader.h"

#include <modfem/uth_system.h>

#include <stdio.h>
#include <vector>
#include <stack>
#include <algorithm>
#include <functional>
#include <cmath>
#include <fstream>
#include <sstream>

#ifndef M_PI_2
const double M_PI_2 = 1.57079632679489662;
#endif

//#define PRINT_BC_GEN_INFO 1
#ifdef _MSC_VER
template<>
void StaticPool<Vertex, 0>::rebuildVts2Pos() {}

template<>
void StaticPool<Vertex, 0>::checkIfInVts2Pos(const uTind verts[]) const {}

template<>
void StaticPool<Vertex, 0>::insertIntoVts2Pos(const uTind verts[], const int pos) {}

#endif
int hHybridMesh::lastMesh = 0; // no. mesh created

hHybridMesh::hHybridMesh(MeshWrite::IMeshWriter * defaultWriter,
		MeshRead::IMeshReader * defaultReader)
	: meshId_(++lastMesh), name_("noname_hybrid_mesh"), gen_(0), maxGen_(99), maxGenDiff_(1),
	  defaultWriter_(defaultWriter), defaultReader_(defaultReader),
	  wasNormalized_(false), isRefining_(false), elements_(this), faces_(this), edges_(this), vertices_(this)
{
	//	    for(int i(0); i < hashSize; ++i) {
	//edgesHash[i]=NULL;
	//facesHash[i]=NULL;
	//}
	free();
}

hHybridMesh::~hHybridMesh()
{
	//utr_io_write_img_to_pnm(NULL,NULL,NULL,0,0,0,0,NULL,NULL);
	free();
}



// Iterator and constIterator operators:
bool	hHybridMesh::checkGeometry() const
{
	// TODO:
	return true;
}

//// Warnign: function depracted as we start to use std::map in static pool
//// and std::map assure that all entries are unique
bool	hHybridMesh::checkUniqueness() const
{
	bool success(true);
	mmv_out_stream << "\n Checking uniqueness of mesh entities..."; // omp_num_threads="<<omp_get_num_threads()<< " ...";
	//#pragma omp parallel sections
	{
		//#pragma omp section
		{
			for (EdgePool::constIterator<> it(&edges_); !it.done(); ++it)
			{
				for (EdgePool::constIterator<> it2(&edges_, it->pos_ + 1); !it2.done(); ++it2)
				{
					if ( ((it->verts(0) == it2->verts(0)) && (it->verts(1) == it2->verts(1)))
							|| ((it->components(0) == it2->components(0)) && (it->components(1) == it2->components(1))) )
						throw "Test error: the same edge twice!";
				}
			}
		}

		//#pragma omp section
		{
			for (FacePool::constIterator<> it(&faces_); !it.done(); ++it)
			{
				for (FacePool::constIterator<> it2(&faces_, it->pos_ + 1); !it2.done(); ++it2)
				{
					if ( ((it->verts(0) == it2->verts(0)) && (it->verts(1) == it2->verts(1)) && (it->verts(2) == it2->verts(2)))
							|| ((it->components(0) == it2->components(0)) && (it->components(1) == it2->components(1)) && (it->components(2) == it2->components(2)) ) )
						throw "Test error: the same face twice!";
				}
			}
		}

	}

	mmv_out_stream << "OK!";
	return success;
}

bool	hHybridMesh::checkTypes() const
{
	//ElemT4 & element = *(_elementsT4.alloc());
	//
	//element.vertex(0,_vertices_.alloc());
	//
	//hObj * obj = dynamic_cast<hObj*>(element.vertex(0));

	return true;
}

bool    hHybridMesh::checkAll() const
{
	elements_.checkHash();
	faces_.checkHash();
	edges_.checkHash();
	//vertices_.checkHash();
//#ifdef _DEBUG
//	assert(checkUniqueness());
//#endif
	return checkGeometry() && checkTypes() ;
}

bool hHybridMesh::read(MeshRead::IMeshReader & reader)
{
	assert(mmv_out_stream << "\n >>> Warning: assertions are enabled! <<<\n" || true);
	const std::string name(reader.Name() + "> ");
	/*assert(out_stream << name <<"Creating hHybridMesh with sizeof(Vertex)="<<sizeof(Vertex)<<
		 ",\n" << name <<" sizeof(Edge)="<<sizeof(Edge)<<
	    ", sizeof(EdgeD)="<<sizeof(EdgeD)<<
	    ",\n"<< name <<" sizeof(Face3)="<<sizeof(Face3)<<
	    ", sizeof(Face3D)="<<sizeof(Face3D)<<
	    ",\n" << name << " sizeof(Face4)="<<sizeof(Face4)<<
	    ", sizeof(Face4D)="<<sizeof(Face4D)<<
	    ",\n" << name << " sizeof(ElemT4)="<<sizeof(ElemT4)<<
	    ", sizeof(ElemT4D)="<<sizeof(ElemT4D)<<
	    ",\n"<< name <<" sizeof(ElemPrism)="<<sizeof(ElemPrism)<<
	    ", sizeof(ElemPrismD)="<<sizeof(ElemPrismD)  ); */

	if (reader.Init()) {
		//hObj::myMesh = this;
		if (reader.doRead(this)) {
			wasNormalized_ = false;
			normalizationProcessor();
			printSetup();
		}
		else
		{
			std::cerr << "\n" << name << "MeshFileImporter: Unable to initialize mesh reader correctly.\n Check filename and file existence.\n";
			//assert(!"Unable to initialize mesh reader correctly.");
			return false;
		}
	}
	return this->test();

}

bool hHybridMesh::normalizationProcessor()
{
	// checking faces neighbours
	if (wasNormalized_ == false)
	{
		mf_log_info("Normalizing mesh.");
		Memory::Field<2, ID> * faceNeighs( new Memory::Field<2, ID>[faces_.last() + 1] );

		for (hHybridMesh::ElemPool::Iterator<hHybridMesh::ElemPool::allObj> it(& elements_);
				!it.done(); ++it)
		{
			for (uTind f(0); f < it->typeSpecyfic_.nFaces_; ++f)
			{
				if (faceNeighs[hObj::posFromId(it->components(f))][0] == 0)
				{
					faceNeighs[hObj::posFromId(it->components(f))][0] = it->id_;
				}
				else
				{
					assert(faceNeighs[hObj::posFromId(it->components(f))][1] == 0);
					faceNeighs[hObj::posFromId(it->components(f))][1] = it->id_;
				}
			}
		}


		//hObj * fPtr(& faces_.at(FIRST));
		for (FacePool::Iterator<> it(& faces_); !it.done(); ++it)
			//for(unsigned int i(FIRST); i <= faces_.last(); ++i)
		{
			if (it->neighs(0) != faceNeighs[it->pos_][0]) {
				assert(faceNeighs[it->pos_][0] != B_COND);
				it->neighs(0) = faceNeighs[it->pos_][0];
			}
			if (it->neighs(1) != faceNeighs[it->pos_][1])
			{
				it->neighs(1) = faceNeighs[it->pos_][1];
			}
			assert(it->neighs(0) != B_COND);
			assert(it->neighs(0) != UNKNOWN);
			assert(it->neighs(0) != it->neighs(1));
			assert(it->neighs(1) != UNKNOWN);
		}

		delete [] faceNeighs;

		///////////////////////////////////////

		for (allIterator el(&elements_); !el.done(); ++el) {
			normalizeElem(*el);
		}

		// Faces at BOUNDARY must have NORMAL vector pointing OUTSIDE.
		// All elems must have defined type (I or II).

		// To check this, and mark if needed we do as follows:
		hObj * elem(NULL);
		//bool flip(false);
		std::vector<int>	vertBCs(vertices_.last() + FIRST);
		for (FacePool::Iterator<> it(&faces_); !it.done(); ++it) {
			if (it->neighs(1) == B_COND) { // if boundary face
				elem = & elements_.getById(it->neighs(0));
#ifdef TURBULENTFLOW
				computePlaneCoords(*it);
#endif
				for (uTind i(0); i < it->typeSpecyfic_.nVerts_; ++i) {
					vertices_[it->verts(i)].setAtBoundary(true);
#ifdef TURBULENTFLOW
					vertices_[it->verts(i)].nearstFace_ = it->id_;
#endif
					vertBCs.at(it->verts(i)) = it->flags(B_COND);
				}

				normalizeFace(*it, *elem);
				//flip=false;
				assert(faceDirection(*it, *elem) == F_OUT);
			}
		}
		// edge elems
		findEdgeElems();

		wasNormalized_ = true;
	}
	return true;
}

bool hHybridMesh::printSetup()
{
	if (!vertices_.empty()) {
		assert(sizeof(Vertex) == vertices_.first().mySize_);
	}
	//assert(sizeof(Edge) == edges_.first().mySize_);

	// Reading detailed boundary description (if exist).
	//std::string detFile(reader.Name()+".detail");
	//out_stream << name << "Trying to read detailed boundary description from file " << detFile << "... ";
	//	MeshRead::SimplePointReader reader2(detFile);

	//reader.Free();
	//	if(reader2.Init()) {
	//	  GeometryModule::Instance().ReadPoints(reader2);
	//	  out_stream << " succeeded.\n";
	//	}
	//	else {
	//	  out_stream << " failed. File not found.\n";
	//	}

	mmv_out_stream << "Analizing mesh...\n";
	//wasNormalized_ = false;//! reader.needsSetup();
	/// Marking element types
	int nElemT4(0), nElemPrism(0);
	double hsize(0.0), avgHSize(0.0),
		   minHSize(std::numeric_limits<double>::max() ),
		   maxHSize(std::numeric_limits<double>::min() );
	for (allIterator el(&elements_); !el.done(); ++el) {

		hsize = elemHSize(*el);
		avgHSize += hsize;
		minHSize = std::min(minHSize, hsize);
		maxHSize = std::max(maxHSize, hsize);
		switch (el->type_) {
			case ElemPrism::myType :
				++nElemPrism;
				break;
			case ElemT4::myType:
				++nElemT4;
				break;
		}
	}
	avgHSize /= elements_.size();

	mmv_out_stream << "\n" << "Found " << nElemT4 << " tetrahedrons and " << nElemPrism << " prismatic elements.";
	mmv_out_stream << "\n" << "Avg element h_size=" << avgHSize
			<< "\n" << " h_size span <" << minHSize << ", " << maxHSize << ">.";

	if (nElemT4 == 0 && nElemPrism == 0) {
		mmv_out_stream << "\n Empty mesh.";
		return true;
	}

	// Faces at BOUNDARY must have NORMAL vector pointing OUTSIDE.
	// All elems must have defined type (I or II).

	// To check this, and mark if needed we do as follows:
	const int bcSpan(1000);
	std::vector<int>	bcCount(bcSpan * 2);
	std::vector<Memory::Field<3, double>> bcNormal(bcCount.size());
	std::vector<Memory::Field<3, double>> bcCenter(bcCount.size());
	std::vector<Memory::Field<2, double>> bcXspan(bcCount.size());
	std::vector<Memory::Field<2, double>> bcYspan(bcCount.size());
	std::vector<Memory::Field<2, double>> bcZspan(bcCount.size());
	for (unsigned int i(0); i < bcCount.size(); ++i) {
		bcXspan[i][0] = std::numeric_limits<double>::max();
		bcXspan[i][1] = -std::numeric_limits<double>::max();
		bcYspan[i][0] = std::numeric_limits<double>::max();
		bcYspan[i][1] = -std::numeric_limits<double>::max();
		bcZspan[i][0] = std::numeric_limits<double>::max();
		bcZspan[i][1] = -std::numeric_limits<double>::max();
	}
	double	vecNorm[3] = {0.0};
	double  fCenter[3] = {0.0};
	hObj * elem(NULL);
	//bool flip(false);
	std::vector<int>	vertBCs(vertices_.last() + FIRST);
	for (FacePool::Iterator<> it(&faces_); !it.done(); ++it) {
		if (it->neighs(1) == B_COND) { // if boundary face

			const int bcId(bcSpan + it->flags(B_COND));
			++bcCount[bcId];
			elem = & elements_.getById(it->neighs(0));

			faceNormal(*it, vecNorm, NULL);

			bcNormal[bcId][0] += vecNorm[0];
			bcNormal[bcId][1] += vecNorm[1];
			bcNormal[bcId][2] += vecNorm[2];
			it->geoCenter(this, fCenter);
			bcCenter[bcId][0] += fCenter[0];
			bcCenter[bcId][1] += fCenter[1];
			bcCenter[bcId][2] += fCenter[2];

			for (uTind i(0); i < it->typeSpecyfic_.nVerts_; ++i) {
				const Memory::Field<3, double> & v(vertices_[it->verts(i)].coords_);
				bcXspan[bcId][0] = std::min(bcXspan[bcId][0], v[X]);
				bcXspan[bcId][1] = std::max(bcXspan[bcId][1], v[X]);

				bcYspan[bcId][0] = std::min(bcYspan[bcId][0], v[Y]);
				bcYspan[bcId][1] = std::max(bcYspan[bcId][1], v[Y]);

				bcZspan[bcId][0] = std::min(bcZspan[bcId][0], v[Z]);
				bcZspan[bcId][1] = std::max(bcZspan[bcId][1], v[Z]);
			}

		}

	}

	for (uTind i(0); i < bcCount.size(); ++i) {
		if (bcCount[i] > 0) {
			bcNormal[i][0] /= bcCount[i]; // computing avg
			bcNormal[i][1] /= bcCount[i];
			bcNormal[i][2] /= bcCount[i];
			bcCenter[i][0] /= bcCount[i];
			bcCenter[i][1] /= bcCount[i];
			bcCenter[i][2] /= bcCount[i];
			double spanCenter[3] = {(bcXspan[i][1] + bcXspan[i][0]) * 0.5,
							(bcYspan[i][1] + bcYspan[i][0]) * 0.5,
							(bcZspan[i][1] + bcZspan[i][0]) * 0.5
					};
			mmv_out_stream << "\n" << "Detected " << bcCount[i] << " faces with BC = " << i - bcSpan
					<< "\n" << " X span <" << bcXspan[i][0] << ", " << bcXspan[i][1] << ">"
					<< "\n" << " Y span <" << bcYspan[i][0] << ", " << bcYspan[i][1] << ">"
					<< "\n" << " Z span <" << bcZspan[i][0] << ", " << bcZspan[i][1] << ">"
					<< "\n" << " span center point (" << spanCenter[0] << ", " << spanCenter[1] << ", " << spanCenter[2]
					<< "\n" << " avg center point (" << bcCenter[i][0] << ", " << bcCenter[i][1] << ", " << bcCenter[i][2]
					<< ")\n" << " avg normal vec [" << bcNormal[i][0] << ", " << bcNormal[i][1] << ", " << bcNormal[i][2] << "].";
		}
	}
	// dump out bc conditions info
	//out_stream << "\n Dumping field_bc.dmp...";
	//std::ofstream	bcFile("field_bc.dmp");
	//if(bcFile.is_open() && bcFile.good()) {
	//    bcFile << "1 1\n" << vertices_.size() << "\n";
	//    for(uTind i(FIRST); i <= vertices_.size(); ++i) {
	//	bcFile << i << " " << vertBCs[i] << "\n";
	//    }
	//    bcFile.close();
	//}
	mmv_out_stream << "\n" << " Collecting edge elems...";

	//out_stream << "OK!";

	const double mem(totalSize() / 10e6);
	mmv_out_stream << "\n" << "Memory used by hHybridMesh: ~" << mem << " MB";
	//std::cout << mmv_out_stream.str() << std::endl;



	return true;
}

void    hHybridMesh::findEdgeElems()
{
	edge_elems.clear();
	//edge_elems.resize(edges_.size()+FIRST+1);
	edge_elems.resize(elements_.last() * 9 + FIRST);

	// new algo
	for (allConstIterator el(&elements_); !el.done(); ++el) {
		edge_elems.at(edge(el->verts(0), el->verts(1)).pos_).push_back(el->pos_);
		edge_elems.at(edge(el->verts(0), el->verts(2)).pos_).push_back(el->pos_);
		edge_elems.at(edge(el->verts(0), el->verts(3)).pos_).push_back(el->pos_);
		edge_elems.at(edge(el->verts(1), el->verts(2)).pos_).push_back(el->pos_);
		switch (el->type_) {
			case ElemT4::myType:
				edge_elems.at(edge(el->verts(1), el->verts(3)).pos_).push_back(el->pos_);
				edge_elems.at(edge(el->verts(2), el->verts(3)).pos_).push_back(el->pos_);
				break;
			case ElemPrism::myType:
				edge_elems.at(edge(el->verts(1), el->verts(4)).pos_).push_back(el->pos_);
				edge_elems.at(edge(el->verts(2), el->verts(5)).pos_).push_back(el->pos_);
				edge_elems.at(edge(el->verts(3), el->verts(4)).pos_).push_back(el->pos_);
				edge_elems.at(edge(el->verts(3), el->verts(5)).pos_).push_back(el->pos_);
				edge_elems.at(edge(el->verts(4), el->verts(5)).pos_).push_back(el->pos_);
				break;
		}
	}
}

bool hHybridMesh::read(MeshRead::IMeshReader * readers[], const int noReaders)
{
	bool ok(true);
	for (int i(0); i < noReaders; ++i)
	{
		read(*readers[i]);
	}
	return ok;
}

bool hHybridMesh::write(MeshWrite::IMeshWriter & writer) const
{
	writer.Free();
	if (writer.Init())
	{
		if (writer.doWrite(this))
		{
			return true;
		}

		for (const char * what_now(writer.writeSequence); *what_now != MeshWrite::IMeshWriter::END ; ++what_now)
			switch (*what_now)
			{
				case MeshWrite::IMeshWriter::AUX:
					break;

				case MeshWrite::IMeshWriter::V_COUNT:
					writer.WriteVerticesCount(static_cast<int>(vertices_.size()));
					break;

				case MeshWrite::IMeshWriter::VERTICES:
				{
					const Vertex * it(&vertices_.first());
					for (const Vertex * const end(it + vertices_.last()); it <= end; ++it)
					{
						writer.WriteVertex(it->coords_, Vertex::myType);
					}
				}
				break;

				case MeshWrite::IMeshWriter::E_COUNT:
					writer.WriteEdgesCount(edges_.size());
					break;

				case MeshWrite::IMeshWriter::EDGES:
					for (EdgePool::constIterator<> it(&edges_); !it.done(); ++it)
					{
						writer.WriteEdge(it->components(), Edge::myType);
					}
					break;

				case MeshWrite::IMeshWriter::F_COUNT:
					writer.WriteFacesCount(faces_.size());
					break;

				case MeshWrite::IMeshWriter::FACES:
					for (FacePool::constIterator<> it(&faces_); !it.done(); ++it)
					{
						assert(it->type_ <= 4);
						/*uTind   comp[4];
						uTind	params[2];
						for (uTind i(0); i < 4; ++i)
						{
							comp[i] = hObj::posFromId(it->components(i));
							if(i < 2)
							{
								params[i] = hObj::posFromId(it->params_[i]);
							}
						}*/

						writer.WriteFace(it->components(),
								it->type_,
								it->flags(B_COND),
								it->neighs());
					}
					break;

				case MeshWrite::IMeshWriter::EL_COUNT:
					writer.WriteElementCount(elements_.size());
					break;

				case MeshWrite::IMeshWriter::ELEMENTS:
					for (ElemPool::constIterator<> it(& elements_); !it.done(); ++it)
					{
						assert(it->type_ >= 5);
						uTind   comp[6];
						uTind  size = 0;
						size = it->typeSpecyfic_.nComponents_;
						for (uTind i(0); i < size; ++i)
						{
							comp[i] = hObj::posFromId(it->components(i));
						}

						writer.WriteElement(comp,
								NULL,
								it->type_,
								hObj::posFromId(it->parent_),
								it->flags(REFINEMENT),
								it->flags(GROUP_ID));

					}
					break;
			}
	}
	return false;
}

bool	hHybridMesh::free()
{
	vertices_.clear();
	edges_.clear();
	faces_.clear();
	elements_.clear();

	//for(int i(0); i < hashSize; ++i) {
	//    if(edgesHash[i] != NULL) {
	//	safeDeleteArray(edgesHash[i]);
	//   }
	//   if(facesHash[i] != NULL) {
// 		safeDeleteArray(facesHash[i]);
// 	    }
// 	}

// 	memset(edgesHash,NULL,hashSize*sizeof(Tind));
// 	memset(facesHash,NULL,hashSize*sizeof(Tind));

	name_.clear();

	return true;
}

bool    hHybridMesh::initRefine()
{
	bool canRefine(true);

	if (isRefining_)
	{
		canRefine = false;
	}
	else
	{
		isRefining_ = true;
	}
	return canRefine;
}

bool hHybridMesh::refine()
{
	if (isRefining_)
	{
		for (hHybridMesh::Iterator   it(begin()); !it.done(); ++it)
		{
			refineElem(*it);
		}
	}
	return isRefining_;
}

bool hHybridMesh::refineElem(const int elemPos)
{
	assert(static_cast<uTind>(elemPos) <= elements_.last());
	return refineElem(elements_[elemPos]);
}

bool hHybridMesh::refineElem(Elem & el)
{
	if (isRefining_)
	{

		// in some cases it is possible to try mark same element again.
		assert(!el.isBroken());
		if (!el.isBroken() && !el.isMarkedBreak())
		{
			el.mark2Ref(this);
			wasNormalized_ = false;
		}
		//else{
		//    throw "Refinement of already refined elem!";
		//}
	}
	return isRefining_;
}

bool hHybridMesh::rRefine(const int BC_id, void (*reallocation_func)(double * x, double * y, double * z))
{
	assert(reallocation_func != NULL);
	if (reallocation_func != NULL)
	{
		// find faces_ with given bc
		// collect vertices_ at given boundary
		std::vector<ID>	bound_faces;
		bound_faces.reserve(faces_.size() / 3);

		std::vector<ID>	bound_points;
		bound_points.reserve(bound_faces.size());

		for (FacePool::Iterator<> iter(&faces_); !iter.done(); ++iter)
		{
			Face3 & it = *iter.as<Face3>();
			if (it.flags(B_COND) == BC_id)
			{
				bound_faces.push_back(it.pos_);
				if (vertices_.getById(edges_.getById(it.components(0)).components(0)).isAtBoundary() == false)
				{
					vertices_.getById(edges_.getById(it.components(0)).components(0)).setAtBoundary(true);
					bound_points.push_back(vertices_.getById(edges_.getById(it.components(0)).components(0)).pos_);
				}

				if (vertices_.getById(edges_.getById(it.components(0)).components(1)).isAtBoundary() == false)
				{
					vertices_.getById(edges_.getById(it.components(0)).components(1)).setAtBoundary(true);
					bound_points.push_back(vertices_.getById(edges_.getById(it.components(0)).components(1)).pos_);
				}

				if (vertices_.getById(edges_.getById(it.components(1)).components(1)).isAtBoundary() == false)
				{
					vertices_.getById(edges_.getById(it.components(1)).components(1)).setAtBoundary(true);
					bound_points.push_back(vertices_.getById(edges_.getById(it.components(1)).components(1)).pos_);
				}
			}
		}

		for (size_t i(0); i < bound_points.size(); ++i)
		{
			reallocation_func(& vertices_.getById(bound_points[i]).coords_[0],
					& vertices_.getById(bound_points[i]).coords_[1],
					& vertices_.getById(bound_points[i]).coords_[2]);
		}
	}
	else
		throw "hHybridMesh::r_refine reallocation_func=NULL!";
	return true;
}

bool hHybridMesh::derefine()
{
	if (isRefining_)
	{
		for (hHybridMesh::Iterator   it(begin()); !it.done(); ++it)
		{
			derefineElem(*it);
		}
	}
	return isRefining_;
}

bool hHybridMesh::derefineElem(const int elemPos)
{
	assert(static_cast<uTind>(elemPos) <= elements_.last());
	return derefineElem(elements_[elemPos]);
}

bool hHybridMesh::derefineElem(hObj & el)
{
	if (isRefining_)
	{

		if (el.isBroken()) {
			el.mark2Dref(this);
			//mark2drf(*elements_.get(elemID)); // -1 for temporary parent
		}
	}
	return isRefining_;
}

//bool hHybridMesh::derefineElem(const int type, const int pos_)
//{
//    if(isRefining_)
//    {
//        mark2drf(hObj(pos,type));
//    }
//    return isRefining_;
//}

/** \brief Assume that: no more than one-level_adaptation is done in one REFINEMENT.
*/
bool    hHybridMesh::finalRef() // Actual refine.
{
	if (isRefining_)
	{
		mf_log_info("hHybridMesh::finalRef begin");
		//computing new space requirements is done on-the-fly while markings
		//save old mesh to backupfile
		std::string file_name(name_);
		file_name.append("_");
		{
			char tmp[16] = {};
			sprintf(tmp, "%d", gen_);
			file_name.append(tmp);
		}

		if (defaultWriter_ != NULL)
		{
			defaultWriter_->Init(file_name);
			write(*defaultWriter_);
		}

		// derefine mesh
		//actualDerefine();

		//try to reallocate memory
		// if expanding, try dynamically
		try
		{
			time_init();
//#pragma omp parallel
			{
//#pragma omp sections
				{
//#pragma omp section
					vertices_.adjust();
//#pragma omp section
					{
						edges_.adjust();

						faces_.adjust();

						elements_.adjust();
					}
				}
			}
			time_print();
			assert(test());
		}
		catch (const std::bad_alloc & err)	// memory error occures
		{
			// if not succeeded, then free mesh and try to allocate enugh space for mesh after REFINEMENT, then read mesh from backup
			free();	// TODO: NOTE: but do not clear div/del marks!

			try
			{

				vertices_.adjust();
				edges_.adjust();
				faces_.adjust();
				elements_.adjust();

				if (defaultReader_ != NULL)
				{
					defaultReader_->Init(file_name);
					read(*defaultReader_);
				}
				else
					throw err;
			}
			catch (const std::bad_alloc & e)
			{
				throw e;
			}

		}
		// here we still have old mesh, but with memory  adjusted for new elements

		// derefine first
		//         actualDerefine();
		//         actualDelete();
		//         // refine mesh
		//         actualRefine();

		// End of de- and re-finement


		isRefining_ = false;
		edge_elems.clear();
		normalizationProcessor();
		//checkSetup(); // NOTE: bad idea?
		assert(test());	// self-check in DEBUG mode
		mf_log_info("hHybridMesh::finalRef end");
	}
	return !isRefining_;
}

// bool	hHybridMesh::actualRefine()
// {
//     for (ElemPool::Iterator  it(&elements); !it.done(); ++it)
//     {
//         it->hBreak();
//     }
//     return true;
// }

// bool	hHybridMesh::actualDerefine()
// {
//     for (ElemPool::Iterator  it(&elements); !it.done(); ++it)
//     {
//         it->derefine();
//     }
//     return true;
// }
//
//
// bool	hHybridMesh::actualDelete()
// {
//     return actualDerefine();
// }

bool    hHybridMesh::test() const
{
	checkAll();
	int	iter_count(0);
//#ifdef _DEBUG
//	out_stream << "\nPrint elements? [t/n]";
//	char a(0);
//	std::cin >> a;
//#endif

	for (hHybridMesh::constIterator it(begin()); !it.done(); ++it)
	{
		assert(it->test(this));
		if (it->isBroken())
		{
			throw "hHybridMesh::test failed: iteration fail.";
		}
		if (++iter_count > 1000000000)
		{
			throw "hHybridMesh::test failed: iteration condition (i > 1000000000) fails (neverending loop)";
		}
//#ifdef _DEBUG
//		if('t'==a)
//		{
//			out_stream << "\nhHybridMesh::test: Element type: " << it->type_ << ", pos: " << it->pos_ << ", vertices: ";
//			for(unsigned int i(0); i < it->typeSpecyfic_.nVerts_; ++i)
//			{
//				out_stream << it->verts(i)/*elemVertex(*it,i)*/ << " ";
//			}
//			out_stream << ", type: " << it->flags(EL_TYPE);
//		}
//		//if(iter_count%1000 == 0)
//		//{
//		//	out_stream << " Iter:" << iter_count;
//		//}
//#endif
	}
	//	std::cerr << "\nhHybridMesh::test: Passed.\n";
	return true;
}

ID  hHybridMesh::elemNeigh(const hObj & elem, const Tind neighNo) const
{
	const ID * nbrs( & faces_.getById(elem.components(neighNo)).neighs(0));
	return nbrs[0] == elem.id_ ? nbrs[1] : nbrs[0] ;
}

int hHybridMesh::whichNeighAmI(const hObj & elem, const hObj & other) const
{
	assert(elem.id_ != other.id_);
	int result(-1);
	for (register unsigned int n(0); (n < elem.typeSpecyfic_.nComponents_) && (result == -1); ++n)
	{
		if (faces_.getById(other.components(n)).neighs(0) == elem.id_
				|| faces_.getById(other.components(n)).neighs(1) == elem.id_) {
			result = n;
		}
	}
	return result;
}

int hHybridMesh::whichFaceAmI(const hObj & face, const hObj & elem) const
{
	assert(face.type_ < ElemPrism::myType);
	assert(elem.type_ >= ElemPrism::myType);
	const ID * faceIds(elem.components());
	int faceCount(elem.typeSpecyfic_.nComponents_);
	register int i(0);
	while (i < faceCount)
	{
		if (faceIds[i] == face.id_)
			break;
		++i;
	}
	assert((i == faceCount) || (elem.components(i) == face.id_));
	return i == faceCount ? -1 : i;
}

int hHybridMesh::whichNodeAmI(const hObj & node, const hObj & elem) const
{
	register int i(0), nodeCount(elem.typeSpecyfic_.nVerts_);
	while (i < nodeCount)
	{
		if (elem.components(i) == node.id_)
			break;
		++i;
	}
	return i;
}


ID	hHybridMesh::faceVertex(const hObj & face, const Tind vertNo) const
{
	assert(vertNo >= 0 && vertNo < 4);
	assert(face.type_ > 2);
	assert(face.type_ < 5);
	ID vertID(UNKNOWN);
	switch (vertNo)
	{
		case 0:
		case 1:
			vertID = edges_.getById(face.components(0)).components(vertNo);
			break;
		case 2:
			vertID = edges_.getById(face.components(1)).components(1);
			break;
		case 3:
			vertID = edges_.getById(face.components(3)).components(1);
			break;
		case 4:
			throw "Not implemented!";
	}
	return vertID;
}

ID	hHybridMesh::elemVertex(const hObj & elem, const Tind vertNo) const
{
	assert(vertNo >= 0 && vertNo < 6);
	assert(elem.type_ > 4);
	ID	vertID(elem.verts(vertNo));
	/*
	switch(vertNo)
	{
	case 0:
	case 1:
	case 2: vertID= faceVertex(faces_.getById(elem.components(0)),vertNo);break;
	case 3: vertID= faceVertex(faces_.getById(elem.components(2)),2);break;
	case 4: vertID= faceVertex(faces_.getById(elem.components(1)),1);break;
	case 5: vertID= faceVertex(faces_.getById(elem.components(1)),2);break;
	}*/
	return vertID;
}


Edge	& hHybridMesh::edge(const uTind V0, const uTind V1)
{
	//	assert( V0 < hashSize);
	assert(V0 >= FIRST);
	assert(V1 >= FIRST);
	uTind v[2] = {V0, V1};
	EdgePool::Vts2posIter it(edges_.vts2pos_.find(v));

	Edge * ptr(NULL);
	if (it == edges_.vts2pos_.end()) {
		// create new
		ptr = edges_.newObj<Edge>(this, v);
	} else {
		ptr = & edges_[it->second];
	}

// 	if(edgesHash[v0] == NULL)
// 	{
// 		edgesHash[v0]=new Tind[edgesPerVertex];
// 		memset(edgesHash[v0],0,edgesPerVertex*sizeof(Tind));
// 	}
// 	register Edge*	ptr(NULL);
// 	Tind * hash(edgesHash[v0]);
// 	int i(0);
// 	for(;(hash[i] > 0) && (ptr == NULL);i+=2)
// 	{
// 		assert( i < edgesPerVertex);
// 		if(hash[i]==v1)
// 		{
// 			ptr = &edges_[hash[i+1]];
// 		}
// 	}

// 	if(i >= edgesPerVertex)
// 	{
// 	  //assert(!"hHybridMesh::edge: edgeesPerVertex exceeded!");
// 		throw "hHybridMesh::edge: edgeesPerVertex exceedeed!";
// 	}

// 	if(ptr == NULL)
// 	{
// 		const uTind	vertices[2]={V0,V1};
// 		ptr = edges_.newObj<Edge>(vertices);
// 		ptr->components(0)=hObj::makeId(V0,Vertex::myType);
// 		ptr->components(1)=hObj::makeId(V1,Vertex::myType);
// 	}

	assert(ptr != NULL);
	return *ptr;
}

hObj	& hHybridMesh::face(const uTind v[4])
{
	assert( v[0] >= FIRST);
	assert( v[1] >= FIRST);
	assert( v[2] >= FIRST);
	assert( v[3] >= FIRST);
	FacePool::Vts2posIter it(faces_.vts2pos_.find(v));

	hObj * ptr(NULL);
	if (it == faces_.vts2pos_.end()) {
		// not found in hash, so create and register new
		if (v[3] == UNKNOWN) {
			ptr = faces_.newObj<Face3>(this, v);
		} else {
			ptr = faces_.newObj<Face4>(this, v);
		}
	} else {
		// found pos, so get it
		ptr = & faces_[it->second];
	}
	assert(ptr != NULL);
	return *ptr;
}

const hObj	& hHybridMesh::face(const uTind v[4]) const
{
	assert( v[0] >= FIRST);
	assert( v[1] >= FIRST);
	assert( v[2] >= FIRST);
	assert( v[3] >= FIRST);
	FacePool::Vts2posCIter it(faces_.vts2pos_.find(v));

	const hObj * ptr(& faces_[it->second]);

	assert(ptr != NULL);
	return *ptr;
}

const hObj & hHybridMesh::face(const uTind v0, const uTind v1, const uTind v2, const uTind v3) const
{
	const uTind v[4] = {v0, v1, v2, v3};
	return face(v);
}

hObj	& hHybridMesh::face(const uTind v0, const uTind v1, const uTind v2, const uTind v3)
{
	const uTind v[4] = {v0, v1, v2, v3};
	return face(v);
	// 	assert( v[0] < hashSize);
// 	if(facesHash[v[0]] == NULL)
// 	{
// 		facesHash[v[0]]=new Tind[facesPerVertex];
// 		memset(facesHash[v[0]],0,facesPerVertex*sizeof(Tind));
// 	}
// 	Tind * hash(facesHash[v[0]]);
// 	register hObj*	ptr(NULL);
// 	int i(0);
// 	for(;(hash[i] > 0) && (ptr == NULL);i+=4)
// 	{
// 		assert( i < facesPerVertex);
// 		if(hash[i]==v[1]
// 		&& hash[i+1]==v[2])
// 		{
// 			if(v3 == UNKNOWN || hash[i+2]==v[3])
// 			{
// 				ptr = &faces_[hash[i+3]];
// 			}
// 		}
// 	}

// 	if(i >= facesPerVertex)
// 	{
// 	  //assert(!"hHybridMesh::face: facesPerVertex exceedeed!");
// 		throw "hHybridMesh::face: facesPerVertex exceedeed!";
// 	}


// 	if(ptr == NULL)
// {
// 	    const uTind vertices[4]={v0,v1,v2,v3};
// 		if(v3==UNKNOWN)
// 		{
// 			ptr = faces_.newObj<Face3>(vertices);
// 			ptr->components(0)=edge(v0,v1).id_;
// 			ptr->components(1)=edge(v0,v2).id_;
// 			ptr->components(2)=edge(v1,v2).id_;
// 		}
// 		else
// 		{
// 			assert(v2 != v3);
// 			ptr = faces_.newObj<Face4>(vertices);
// 			ptr->components(0)=edge(v0,v1).id_;
// 			ptr->components(1)=edge(v0,v2).id_;
// 			ptr->components(2)=edge(v1,v3).id_;
// 			ptr->components(3)=edge(v2,v3).id_;
// 		}
// 	}
//	assert(ptr != NULL);
//	return *ptr;
}



hObj	& hHybridMesh::element(const uTind v0, const uTind v1, const uTind v2, const uTind v3, const uTind v4, const uTind v5)
{
	register hObj	* ptr(NULL);
	//register int		i(0);
	//while(i < v1 && ptr == NULL)
	//{
	//  if(v0index[v0+i]->verts(1)==v1)
	//  {
	//		ptr = v0index[v0+i];
	//  }
	//	++i;
	//}
	return *ptr;
}

//void	hHybridMesh::registerEdge(const uTind V[],const uTind pos)
//{
//	assert( V[0] < hashSize);
//	uTind	v[2]={V[0],V[1]};
//	if(v[0] > v[1]) {
//		std::swap(v[0],v[1]);
//	}
//	register Edge*	ptr(NULL);
//	if(edgesHash[v[0]] == NULL)
//	{
//		edgesHash[v[0]]=new Tind[edgesPerVertex];
//		memset(edgesHash[v[0]],0,edgesPerVertex*sizeof(Tind));
//	}
//	Tind * hash(edgesHash[v[0]]);
//	int i(0);
//	for(;(hash[i] > 0) && (ptr == NULL); i+=2)
//	{
//		assert( i < edgesPerVertex);
//		if(hash[i]==v[1])
//		{
//			ptr = &edges_[hash[i+1]];
//		}
//	}
//
//	if(i >= edgesPerVertex)
//	{
//	  //assert(!"hHybridMesh::registerEdge: edgesPerVertex exceedeed!");
//		throw "hHybridMesh::registerEdge: edgesPerVertex exceedeed!";
//	}
//
//
//	if(ptr != NULL)
//	{
//		throw "hHybridMesh::hash: trying to add existing face!";
//	}
//
//	hash[i]=v[1];
//	hash[i+1]=pos;
//}
//
//void	hHybridMesh::registerFace(const uTind V[],const uTind pos)
//{
// 	assert( V[0] >= FIRST);
// 	assert( V[1] >= FIRST);
// 	assert( V[2] >= FIRST);
// 	assert( V[3] >= FIRST);
// 	uTind v[4]={V[0],V[1],V[2],V[3]};
// 	std::sort(v,v+ (V[3]==UNKNOWN?3:4));
// 	assert( v[0] < hashSize);
// 	if(facesHash[v[0]] == NULL)
// 	{
// 		facesHash[v[0]]=new Tind[facesPerVertex];
// 		memset(facesHash[v[0]],0,facesPerVertex*sizeof(Tind));
// 	}
// 	assert(facesHash[v[0]] != NULL);
// 	Tind * hash(facesHash[v[0]]);
// 	int i(0);
// 	register hObj * ptr(NULL);
// 	for(;hash[i] > 0 && ptr == NULL && i < facesPerVertex;i+=4)
// 	{
// 		assert( i < facesPerVertex);
// 		if(hash[i]==v[1]
// 		&& hash[i+1]==v[2])
// 		{
// 			if((V[3] == UNKNOWN) || (hash[i+2]==v[3]))
// 			{
// 				ptr = &faces_[hash[i+3]];
// 			}
// 		}
// 	}

// 	if(i >= facesPerVertex)
// 	{
// 	  //assert(!"hHybridMesh::registerFace: facesPerVertex exceedeed!");
// 		throw "hHybridMesh::registerFace: facesPerVertex exceedeed!";
// 	}

// 	if(ptr != NULL)
// 	{
// 	  //assert(!"hHybridMesh::hash: trying to add existing face!");
// 		throw "hHybridMesh::hash: trying to add existing face!";
// 	}
// 	else
// 	{
// 		assert( i+3 < facesPerVertex);
// 		hash[i]=v[1];
// 		hash[i+1]=v[2];
// 		hash[i+2]=(V[3]==UNKNOWN ? 0 : v[3]);
// 		hash[i+3]=pos;
// 	}
// }

void hHybridMesh::faceNormal(IN const hObj & face, OUT double vecNorm[3], OUT double * area) const
{
	switch (face.type_) {
		case Face3::myType: {
			const double vec_a[3] = {
				vertices_[face.verts(1)].coords_[0] - vertices_[face.verts(0)].coords_[0],
				vertices_[face.verts(1)].coords_[1] - vertices_[face.verts(0)].coords_[1],
				vertices_[face.verts(1)].coords_[2] - vertices_[face.verts(0)].coords_[2]
			};
			const double vec_b[3] = {
				vertices_[face.verts(2)].coords_[0] - vertices_[face.verts(1)].coords_[0],
				vertices_[face.verts(2)].coords_[1] - vertices_[face.verts(1)].coords_[1],
				vertices_[face.verts(2)].coords_[2] - vertices_[face.verts(1)].coords_[2]
			};

			mmr_vec3_prod(vec_a, vec_b, vecNorm);
			const double len(mmr_vec3_length(vecNorm));
			vecNorm[0] /= len;
			vecNorm[1] /= len;
			vecNorm[2] /= len;
			if (area != NULL) {
				*area = 0.5 * len;
			}
		}
		break;
		case Face4::myType: {
			const double vec_a[3] = {
				vertices_[face.verts(1)].coords_[0] - vertices_[face.verts(0)].coords_[0],
				vertices_[face.verts(1)].coords_[1] - vertices_[face.verts(0)].coords_[1],
				vertices_[face.verts(1)].coords_[2] - vertices_[face.verts(0)].coords_[2]
			};
			const double vec_b[3] = {
				vertices_[face.verts(2)].coords_[0] - vertices_[face.verts(1)].coords_[0],
				vertices_[face.verts(2)].coords_[1] - vertices_[face.verts(1)].coords_[1],
				vertices_[face.verts(2)].coords_[2] - vertices_[face.verts(1)].coords_[2]
			};
			const double vec_c[3] = {
				vertices_[face.verts(3)].coords_[0] - vertices_[face.verts(2)].coords_[0],
				vertices_[face.verts(3)].coords_[1] - vertices_[face.verts(2)].coords_[1],
				vertices_[face.verts(3)].coords_[2] - vertices_[face.verts(2)].coords_[2]
			};
			const double vec_d[3] = {
				vertices_[face.verts(0)].coords_[0] - vertices_[face.verts(3)].coords_[0],
				vertices_[face.verts(0)].coords_[1] - vertices_[face.verts(3)].coords_[1],
				vertices_[face.verts(0)].coords_[2] - vertices_[face.verts(3)].coords_[2]
			};

			mmr_vec3_prod(vec_a, vec_b, vecNorm);
			const double len(mmr_vec3_length(vecNorm));
			vecNorm[0] /= len;
			vecNorm[1] /= len;
			vecNorm[2] /= len;
			if (area != NULL) {
				double vec_f[3] = {0.0};
				mmr_vec3_prod(vec_c, vec_d, vec_f);
				*area = 0.5 * ( len + mmr_vec3_length(vec_f) );
			}

		}
		break;
	}//!switch
}

eFaceFlag hHybridMesh::faceDirection(const hObj & face, const hObj & el) const
{
	// 1. Compute normal vector.
	double vec_norm[3] = {0.0}, vec_ops[3] = {0.0};
	faceNormal(face, vec_norm, NULL);
	/// 2. Compute "center vector" from face center to existing neigh center.
	/// (pointing inside element)
	double fCenter[3] = {0.0};
	face.geoCenter(this, fCenter);
	double elCenter[3] = {0.0};
	el.geoCenter(this, elCenter);
	mmr_vec3_subst(elCenter, fCenter, vec_ops);
	const double len_ops(mmr_vec3_length(vec_ops));
	vec_ops[0] /= len_ops; // normalization of vector
	vec_ops[1] /= len_ops;
	vec_ops[2] /= len_ops;
	/// 3. Determine angle between normal vector and center vector.
	/// If angle < 90(deg) normal vector IS POINTING INSIDE
	/// and must points outside, so we mark appr. flag.
	// Check if vectors are not parallel (Yes, it is necessarry: it already happend!).
	double angle(0.0);
	if (!((fabs(vec_norm[0] - vec_ops[0]) < SMALL) &&
					(fabs(vec_norm[1] - vec_ops[1]) < SMALL) &&
					(fabs(vec_norm[2] - vec_ops[2]) < SMALL))) {
		angle = acos(mmr_vec3_dot(vec_norm, vec_ops));
	}
	return (angle <= M_PI_2 ? F_IN : F_OUT);
}

eElemFlag hHybridMesh::normalizeElem(hObj & el)
{
	eElemFlag type(EL_TYPEI);
	if (el.type_ == ElemT4::myType) {
		hObj * f[4] = {NULL};
		f[0] = &face(el.verts(0), el.verts(1), el.verts(2));
		f[1] = &face(el.verts(0), el.verts(1), el.verts(3));
		f[2] = &face(el.verts(0), el.verts(2), el.verts(3));
		f[3] = &face(el.verts(1), el.verts(2), el.verts(3));
		assert(f[0]->id_ == el.components(0));
		assert(f[1]->id_ == el.components(1));
		assert(f[2]->id_ == el.components(2));
		assert(f[3]->id_ == el.components(3));
		// if(faceDirection(*f[3],el)==F_IN) {
		//   //assert(faceDirection(*f[1],el)==F_IN);
		//   assert(faceDirection(*f[0],el)==F_OUT);
		//   //assert(faceDirection(*f[1],el)==F_IN);
		//   assert(faceDirection(*f[2],el)==F_OUT);
		//   el.flags(EL_TYPE) = EL_TYPEII; // just temporary
		//   el.swapVerts(1,2);// v:0,2,1,3 = type II
		//   // convert to type I  - rotation of face 2 (v[1] not moves)
		//   el.swapVerts(0,2);// v:1,2,0,3
		//   el.swapVerts(2,3);// v:1,2,3,0
		//   const ID comp3ID = el.components(3);
		//   el.components(3)=el.components(2); //c:
		//   el.components(2)=el.components(1);
		//   el.components(1)=el.components(0); //c:
		//   el.components(0)=comp3ID;

		//    // update pointers to changes
		//    f[0] = &face(el.verts(0),el.verts(1),el.verts(2));
		//    f[1] = &face(el.verts(0),el.verts(1),el.verts(3));
		//    f[2] = &face(el.verts(0),el.verts(2),el.verts(3));
		//    f[3] = &face(el.verts(1),el.verts(2),el.verts(3));
		// }

		//el.flags(EL_TYPE) = EL_TYPEI;
		int mask = 0;
		for (uTind i(0); i < el.typeSpecyfic_.nComponents_; ++i) {
			if (faceDirection(*f[i], el) == F_IN) {
				mask |= 1 << i;
			}
		}
		el.flags(EL_TYPE) = static_cast<eElemFlag>(mask);

		assert(f[0]->id_ == el.components(0));
		assert(f[1]->id_ == el.components(1));
		assert(f[2]->id_ == el.components(2));
		assert(f[3]->id_ == el.components(3));
		assert(faceDirection(*f[0], el) == F_OUT || (mask & EL_TYPE0));
		assert(faceDirection(*f[1], el) == F_OUT || (mask & EL_TYPE1));
		assert(faceDirection(*f[2], el) == F_OUT || (mask & EL_TYPE2));
		assert(faceDirection(*f[3], el) == F_OUT || (mask & EL_TYPE3));

		{	// volume tesing
			const Vertex * v[4] = {&vertices_[el.verts(0)],
							&vertices_[el.verts(1)],
							&vertices_[el.verts(2)],
							&vertices_[el.verts(3)]
					};
			const double v0[3] = {v[1]->coords_[0] - v[0]->coords_[0],
							v[1]->coords_[1] - v[0]->coords_[1],
							v[1]->coords_[2] - v[0]->coords_[2]
					};
			const double v1[3] = {v[2]->coords_[0] - v[0]->coords_[0],
							v[2]->coords_[1] - v[0]->coords_[1],
							v[2]->coords_[2] - v[0]->coords_[2]
					};
			const double v2[3] = {v[3]->coords_[0] - v[0]->coords_[0],
							v[3]->coords_[1] - v[0]->coords_[1],
							v[3]->coords_[2] - v[0]->coords_[2]
					};
			const double volume = mmr_vec3_mxpr(v0, v1, v2) / 6.0; //for t4
			if (volume <= 0.0) {
				std::stringstream ss;
				ss << "\nNormalizeElem: Error: tetra volume <= 0.0! Element pos: ";
				ss << el.pos_ << ". Vertices:\n"
						<< v[0]->coords_ << "\n"
						<< v[1]->coords_ << "\n"
						<< v[2]->coords_ << "\n"
						<< v[3]->coords_ << "\n";

				throw ss.str();
			}
		}//! volume testing
	}
	else if (el.type_ == ElemPrism::myType) {
		hObj * f[5] = {0};
		f[0] = &face(el.verts(0), el.verts(1), el.verts(2));
		f[1] = &face(el.verts(3), el.verts(4), el.verts(5));
		f[2] = &face(el.verts(0), el.verts(1), el.verts(4), el.verts(3));
		f[3] = &face(el.verts(0), el.verts(2), el.verts(5), el.verts(3));
		f[4] = &face(el.verts(1), el.verts(2), el.verts(5), el.verts(4));
		assert(f[0]->id_ == el.components(0));
		assert(f[1]->id_ == el.components(1));
		assert(f[2]->id_ == el.components(2));
		assert(f[3]->id_ == el.components(3));
		assert(f[4]->id_ == el.components(4));
		int mask(0);
		for (uTind i(0); i < el.typeSpecyfic_.nComponents_; ++i) {
			if (faceDirection(*f[i], el) == F_IN) {
				mask |= 1 << i;
			}
		}
		el.flags(EL_TYPE) = static_cast<eElemFlag>(mask);

		/*if(faceDirection(*f0,el)==F_OUT) {
		 el.flags(EL_TYPE)=type=EL_TYPE_PII;
		 assert(faceDirection(*f1,el)==F_IN);
		 assert(faceDirection(*f2,el)==F_IN);
		 assert(faceDirection(*f3,el)==F_OUT);
		 assert(faceDirection(*f4,el)==F_IN);
		 // cross flipping upside-down
		 //new: first rotate
		 el.swapVerts(0,1); el.swapVerts(1,2);
		 el.swapVerts(3,4); el.swapVerts(4,5);
		 const ID tmp=el.components(3);	//OUT face
		 el.components(2)=el.components(3);
		 el.components(3)=el.components(4);
		 el.components(4)=tmp; // OUT as last face
		 // swap upside <=> down
		 el.swapVerts(0,3);
		 el.swapVerts(1,4);
		 el.swapVerts(2,5);
		 std::swap(el.components(0),el.components(1));
		 std::swap(el.components(2),el.components(3));
		 f0 = &face(el.verts(0),el.verts(1),el.verts(2));
		 f1 = &face(el.verts(3),el.verts(4),el.verts(5));
		 f2 = &face(el.verts(0),el.verts(1),el.verts(4),el.verts(3));
		 f3 = &face(el.verts(0),el.verts(2),el.verts(5),el.verts(3));
		 f4 = &face(el.verts(1),el.verts(2),el.verts(5),el.verts(4));
		 assert(faceDirection(*f4,el)==F_OUT);
		}
		else {
		 if(faceDirection(*f4,el)==F_IN) {
		  el.flags(EL_TYPE)=type=EL_TYPE_PII;
		 }
		 else {
		  el.flags(EL_TYPE)=type=EL_TYPE_PI;
		 }
		}
		assert(f0->id_ == el.components(0));
		assert(f1->id_ == el.components(1));
		assert(f2->id_ == el.components(2));
		assert(f3->id_ == el.components(3));
		assert(f4->id_ == el.components(4));
		assert(faceDirection(*f0,el)==F_IN);
		assert(faceDirection(*f1,el)==F_OUT);
		assert(faceDirection(*f2,el)==F_OUT);
		assert(faceDirection(*f3,el)==F_IN);*/

		{	// volume tesing
			const Vertex * v[6] = {&vertices_[el.verts(0)],
							&vertices_[el.verts(1)],
							&vertices_[el.verts(2)],
							&vertices_[el.verts(3)],
							&vertices_[el.verts(4)],
							&vertices_[el.verts(5)]
					};
			const double v0[3] = {v[1]->coords_[0] - v[0]->coords_[0],
							v[1]->coords_[1] - v[0]->coords_[1],
							v[1]->coords_[2] - v[0]->coords_[2]
					};
			const double v1[3] = {v[2]->coords_[0] - v[0]->coords_[0],
							v[2]->coords_[1] - v[0]->coords_[1],
							v[2]->coords_[2] - v[0]->coords_[2]
					};
			const double v2[3] = {v[3]->coords_[0] - v[0]->coords_[0],
							v[3]->coords_[1] - v[0]->coords_[1],
							v[3]->coords_[2] - v[0]->coords_[2]
					};
			const double volume = mmr_vec3_mxpr(v0, v1, v2) / 2.0; // /2.0 for prism
			if (volume <= 0.0) {
				std::stringstream ss("NormalizeElem: Error: prism volume =");
				ss << volume << "! Element pos: "
						<< el.pos_ << ". Vertices: "
						<< v[0]->coords_
						<< v[1]->coords_
						<< v[2]->coords_
						<< v[3]->coords_
						<< v[4]->coords_
						<< v[5]->coords_;

				throw ss.str();
			}
		}//! volume testing
	}

	return type;
}

eFaceFlag hHybridMesh::normalizeFace(hObj & face, hObj & elem)
{
	if (face.neighs(1) == B_COND) {
		if (faceDirection(face, elem) == F_IN) {
			assert(false == wasNormalized_);
			face.flags(F_TYPE) = F_FLIPPED;
			face.swapVerts(1, 2);
			std::swap(face.components(0), face.components(1));
			if (face.type_ == Face4::myType) {
				std::swap(face.components(2), face.components(3)); // other two edges also
			}
			// update elem type flag (bitmask)
			int & ftype(elem.flags(EL_TYPE));
			switch (whichFaceAmI(face, elem)) {
				case 0:
					ftype &= (~EL_TYPE0);
					break;
				case 1:
					ftype &= (~EL_TYPE1) ;
					break;
				case 2:
					ftype &= (~EL_TYPE2) ;
					break;
				case 3:
					ftype &= (~EL_TYPE3) ;
					break;
				case 4:
					ftype &= (~EL_TYPE4) ;
					break;
			}//!switch
		}
		assert(faceDirection(face, elem) == F_OUT);
		assert(face.components(0) == edge(face.verts(0), face.verts(1)).id_);
		assert(face.components(1) == edge(face.verts(0), face.verts(2)).id_);
		switch (face.type_) {
			case Face3::myType:
				assert(face.components(2) == edge(face.verts(1), face.verts(2)).id_);
				break;
			case Face4::myType:
				assert(face.components(2) == edge(face.verts(1), face.verts(3)).id_);
				assert(face.components(3) == edge(face.verts(2), face.verts(3)).id_);
				break;
		}
	}
	return eFaceFlag(face.flags(F_TYPE));
}

double hHybridMesh::elemHSize(const hObj & elem)
{
	/* form three vectors from nodes coordinates */
	//double hsize(0.0);
	Vertex * v[4] = {
		&vertices_[elem.verts(0)],
		&vertices_[elem.verts(1)],
		&vertices_[elem.verts(2)],
		&vertices_[elem.verts(3)]
	};
	double v1[3] = {0.0}, v2[3] = {0.0}, v3[3] = {0.0};
	for (int i(0); i < 3; ++i) {
		v1[i] = v[1]->coords_[i] - v[0]->coords_[i];
		v2[i] = v[2]->coords_[i] - v[0]->coords_[i];
		v3[i] = v[3]->coords_[i] - v[0]->coords_[i];
	}
	/// /* compute hsize as third root of volume (computed as mixed produc/
	const double hsize = pow(fabs(mmr_vec3_mxpr(v1, v2, v3)), 1.0 / 3.0);
	//hsize = fabs(mmr_vec3_mxpr(v1,v2,v3));
	//hsize = pow(hsize,1.0/3.0);
	assert(hsize > 0);
	return hsize;
}

//void hHybridMesh::rotateVerts(IN OUT uTind verts[4],const int i,const int j,const int k) const {
//    register uTind tmp(verts[i]);
//    verts[i]=verts[j];
//    verts[j]=verts[k];
//    verts[k]=tmp;
//}
//
//// returns false if elem is type II
//eElemFlag hHybridMesh::normalize(IN OUT uTind verts[4]) const
//{
//    register Tind i_m = 0;
//    if(verts[1] < verts[i_m])	i_m = 1;
//    if(verts[2] < verts[i_m])	i_m = 2;
//    if(verts[3] < verts[i_m])	i_m = 3;
//
//    switch(i_m) {	// 1st rotation
//    case 1: rotateVerts(verts,0, 1, 2);break;
//    case 2: rotateVerts(verts,0, 2, 1);break;
//    case 3: rotateVerts(verts,0, 3, 1);break;
//    }
//    // searching for max vertex numbers
//    i_m = 1;
//    if(verts[2] > verts[i_m])	i_m = 2;
//    if(verts[3] > verts[i_m])	i_m = 3;
//
//    switch(i_m) {	// 2nd rotation
//    case 1: rotateVerts(verts,3, 1, 2);break;
//    case 2: rotateVerts(verts,3, 2, 1);break;
//    }
//    eElemFlag type = verts[1] < verts[2] ? EL_TYPEI : EL_TYPEII; // type I (true) or II (false)
//    if(type==EL_TYPEII) {
//		std::swap(verts[1],verts[2]);
//    }
//	return type;
//}
#ifdef TURBULENTFLOW
void hHybridMesh::computePlaneCoords(Face & face)
{
	double	vn[3] = {0.0};
	const double * v(vertices_[face.verts(0)].coords_);
	faceNormal(face, vn, NULL);
	face.planeCoords_[0] = vn[0];
	face.planeCoords_[1] = vn[1];
	face.planeCoords_[2] = vn[2];
	face.planeCoords_[3] = - vn[0] * v[X] - vn[1] * v[Y] - vn[2] * v[Z];

}

double hHybridMesh::dist2Face(const Vertex & v, const Face & f)
{
	double dist(0.0);
	if (!v.isAtBoundary()) {
		dist = (f.planeCoords_[0] * v.coords_[X]
						+ f.planeCoords_[1] * v.coords_[Y]
						+ f.planeCoords_[2] * v.coords_[Z]
						+ f.planeCoords_[3])
				/ sqrt(f.planeCoords_[0] * f.planeCoords_[0]
						+    f.planeCoords_[1] * f.planeCoords_[1]
						+    f.planeCoords_[2] * f.planeCoords_[2]);

	}
	return dist;
}

void hHybridMesh::findNearestBoundFace(Vertex & v, const int BCs[], const int nBCs)
{
	if (!v.isAtBoundary()) {
		FacePool::Iterator<> it(&faces_);
		double newdist(0.0), dist( std::numeric_limits<double>::max()  );

		for (; it.done(); ++it) {
			if (it->neighs(1) == B_COND) { // if at boundary
				if (std::find(BCs, BCs + nBCs, it->flags(B_COND)) != BCs + nBCs) { // if good BC number
					newdist = dist2Face(v, *it);
					if (newdist < dist) {
						dist = newdist;
						v.nearstFace_ = it->id_;
					}
				}
			}
		}
	}
}

void hHybridMesh::computeDist2Bound(const int BCs[], const int nBCs)
{
	for (VertexPool::Iterator it(&vertices_.first()); it < vertices_.end(); ++it) {
		if (!it->isAtBoundary()) {
			findNearestBoundFace(*it, BCs, nBCs);
		}
	}


}
#endif

bool	hHybridMesh::createBoundaryLayer(const double thicknessProc, const int nLayers,
		const bool quadraticDistribution,
		const double * vecIgnoreNormal)
{
	assert(thicknessProc > 0.0);
	assert(nLayers >= 1);
	const double correctThicknessProc(thicknessProc * 0.02); // 0.01 = 1%; *2 becouse we multiply vector to center (so it's halved).

	mmv_out_stream << "\nGenerating " << nLayers << " prism boundary layer" << ((nLayers > 1) ? "s" : "") << " (thickness=" << thicknessProc ;
	mmv_out_stream << ", " << (quadraticDistribution ? "quadratic" : "linear") << " distribution ) ...";

	try {

		bool ox(vecIgnoreNormal[0] != 1.0), oy(vecIgnoreNormal[1] != 1.0), oz(vecIgnoreNormal[2] != 1.0);

		// 1st PART - GATHER BOUNDARY SURFACE DATA
		// create list of vertices_ on boundary
		for (Vertex * it(&vertices_.first()), *end( &vertices_.first() + vertices_.size() ); it < end; ++it)
		{
			it->setAtBoundary(false);
		}

		std::vector<ID>	bound_faces;
		bound_faces.reserve(faces_.size() / 3);
		for (FacePool::Iterator<> iter(&faces_); !iter.done(); ++iter)
		{
			Face3 * it = iter.as<Face3>();
			if (it->neighs(1) == 0) {
				assert(it->type_ == Face3::myType);	// must be triangular space
				if (vecIgnoreNormal == NULL) {
					bound_faces.push_back(it->pos_);
					vertices_[it->verts(0)].setAtBoundary(true);
					vertices_[it->verts(1)].setAtBoundary(true);
					vertices_[it->verts(2)].setAtBoundary(true);
				}
				else {
					double vecNorm[3] = {0.0};
					faceNormal(*it, vecNorm, NULL);
					if ((fabs(vecIgnoreNormal[X] - fabs(vecNorm[X])) < SMALL) &&
							(fabs(vecIgnoreNormal[Y] - fabs(vecNorm[Y])) < SMALL) &&
							(fabs(vecIgnoreNormal[Z] - fabs(vecNorm[Z])) < SMALL) ) {
						//vertices_[it->verts(0)].setAtBoundary(false);
						//vertices_[it->verts(1)].setAtBoundary(false);
						//vertices_[it->verts(2)].setAtBoundary(false);
						mmv_out_stream << "ignoring " << it->pos_ << ", ";
					}
					else {
						bound_faces.push_back(it->pos_);
						vertices_[it->verts(0)].setAtBoundary(true);
						vertices_[it->verts(1)].setAtBoundary(true);
						vertices_[it->verts(2)].setAtBoundary(true);
					}
				}
			}
		}

		std::vector<ID>	bound_points;
		bound_points.reserve(bound_faces.size());
		for (const Vertex * it(&vertices_.first()), *end( &vertices_.first() + vertices_.size() ); it < end; ++it)
		{
			if (it->isAtBoundary())
			{
				bound_points.push_back(it->pos_);
			}
		}

		mmv_out_stream << "\nFound " << bound_points.size() << " surface vertices, in " << bound_faces.size() << " faces." ;
		// compute geometrical center of object
		std::vector<Field<3, double>>	origPts(bound_points.size());
		double center[3] = {0};
		int pts(0);
		for (std::vector<ID>::const_iterator it(bound_points.begin()); it != bound_points.end(); ++it, ++pts)
		{
			origPts[pts][X] = vertices_.at(*it).coords_[X];
			origPts[pts][Y] = vertices_.at(*it).coords_[Y];
			origPts[pts][Z] = vertices_.at(*it).coords_[Z];
			center[X] += origPts[pts][X];
			center[Y] += origPts[pts][Y];
			center[Z] += origPts[pts][Z];
		}
		center[0] /= bound_points.size();
		center[1] /= bound_points.size();
		center[2] /= bound_points.size();

		// 2nd PART - GENERATE NEW BOUNDARY PRISM LAYERS
		//const ID last_original_vertex = vertices_.last();
		const uTind newNVerts(vertices_.size() + bound_points.size()*nLayers);
		mmv_out_stream << "\nRising vertices count to " << newNVerts;;
		vertices_.reserve(newNVerts);

		const uTind newNEdges(edges_.size() + bound_faces.size() * 6 * nLayers);
		mmv_out_stream << "\nRising edges count to " << newNEdges;
		edges_.reserve(newNEdges);	// 3 for new external face + 3 for vertical edges

		const uTind newNFace3( faces_.size() + bound_faces.size()*nLayers),
			  newNFace4(bound_points.size() * 3 * nLayers + 1);
		mmv_out_stream << "\nRising faces count to " << (newNFace3 + newNFace4);
		faces_.reserve(newNFace3 + newNFace4, static_cast<uTind>(sizeof(Face3)*newNFace3 + sizeof(Face4)*newNFace4 ));

		const uTind newNElemT4(elements_.size()),
			  newNElemPrism(bound_faces.size()*nLayers);
		mmv_out_stream << "\nRising elements count to " << (newNElemT4 + newNElemPrism);
		elements_.reserve(newNElemT4 + newNElemPrism,
				static_cast<uTind>(sizeof(ElemT4)*newNElemT4 + sizeof(ElemPrism)*newNElemPrism));


		// for each layer execute loop once
		for (int layer(1); layer <= nLayers; ++layer)
		{
			std::vector<ID>	face_map(bound_faces.back() + 1, 0);
			std::vector<ID>	edge_mapI(edges_.size() + 1, 0);
			std::vector<ID>	new_edge_mapI(bound_points.back() + 1, 0);
			std::vector<ID>	vertex_mapI(bound_points.back() + 1, 0);
			std::vector<ID>	vertcial_faces_mapI(faces_.capacity() + 1, 0);
			const int nOldVerts(static_cast<int>(vertices_.size()));

			// 2.1! Changing positions of vertices.
			if (1 == layer) // Shift into geometrical center all vertices_, except lastly created.
			{	// Do it only once for first layer, to create space for all layers.
				Vertex * it(&vertices_.first());
				for (const Vertex * const end(it + nOldVerts); it < end; ++it)
				{
#ifdef PRINT_BC_GEN_INFO
					out_stream << "\n(Once)Moving point " << it->pos_ << " from: ("
							<< it->coords_[X] << " " << it->coords_[Y] << " " <<  it->coords_[Z];
#endif
					if (ox) {
						it->coords_[0] += (center[0] - it->coords_[0]) * correctThicknessProc;
					}
					if (oy) {
						it->coords_[1] += (center[1] - it->coords_[1]) * correctThicknessProc;
					}
					if (oz) {
						it->coords_[2] += (center[2] - it->coords_[2]) * correctThicknessProc;
					}
#ifdef PRINT_BC_GEN_INFO
					out_stream << ") to: ("
							<< it->coords_[X] << " " << it->coords_[Y] << " " <<  it->coords_[Z] << ")";
#endif
				}
			}

			// 2.2 Generate new vertices
			mmv_out_stream << "\nLayer " << layer << "/" << nLayers << " ";
			size_t	size = bound_points.size();
			const double thicnessCoof( quadraticDistribution ?
					(nLayers > 1 ? correctThicknessProc / static_cast<double>(1 << layer) : 0 )
					: ((1.0 - static_cast<double>(layer) / static_cast<double>(nLayers))*correctThicknessProc) );
			for (size_t v(0); v < size; ++v)
			{
				Vertex & oldVert(vertices_[bound_points[v]]);
				Vertex & newVert( *vertices_.newObj<Vertex>(this, NULL) );
				newVert.coords_[X] = origPts[v][X]; // transfer params from old vertex
				newVert.coords_[Y] = origPts[v][Y];
				newVert.coords_[Z] = origPts[v][Z];
				oldVert.setAtBoundary(false);
				newVert.setAtBoundary(true);
//			newVert.status(1); //-1); //==MMC_INACTIVE	// do this do persuade aproximation module to allocate new dofs...
				vertex_mapI[oldVert.pos_] = newVert.id_;

				// Shift.
#ifdef PRINT_BC_GEN_INFO
				out_stream << "\nMoving point " << newVert.pos_  << " from: ("
						<< newVert.coords_[X] << " " << newVert.coords_[Y] << " " <<  newVert.coords_[Z];
#endif
				if (ox) {
					newVert.coords_[0] += (center[0] - newVert.coords_[0]) * (thicnessCoof);
				}
				if (oy) {
					newVert.coords_[1] += (center[1] - newVert.coords_[1]) * (thicnessCoof);
				}
				if (oz) {
					newVert.coords_[2] += (center[2] - newVert.coords_[2]) * (thicnessCoof);
				}
#ifdef PRINT_BC_GEN_INFO
				out_stream << ") to: ("
						<< newVert.coords_[X] << " " << newVert.coords_[Y] << " " <<  newVert.coords_[Z] << ")";
#endif
			}

			// 2.3 Creating vertical edges.
			size = bound_points.size();
			for (size_t v(0); v < size; ++v)
			{
				const uTind	vertices[2] = {bound_points[v], hObj::posFromId(vertex_mapI[bound_points[v]])};
				Edge & newVerticalE( *edges_.newObj<Edge>(this, vertices) );
				newVerticalE.parent_ = 0; //NO_FATHER
				newVerticalE.level_ = 0;
				newVerticalE.components(0) = hObj::makeId(bound_points[v], Vertex::myType);
				newVerticalE.components(1) = vertex_mapI[bound_points[v]];
				new_edge_mapI[bound_points[v]] = newVerticalE.id_;
				assert(newVerticalE.test(this));
			}

			// 2.4 Creating new prism elements.
			for (size_t f(0); f < bound_faces.size(); ++f)
			{
				assert(bound_faces[f] > 0);
				Face3 & original_face( faces_.at<Face3>(bound_faces[f]) );
				hObj & originalEl(elements_.getById(original_face.neighs(0)));
				int whichFace(whichFaceAmI(original_face, originalEl));
				assert(whichFace >= 0 && whichFace < 4);
				if (original_face.flags(F_TYPE) == F_FLIPPED) {
					static const int flipEdgeSwap[3] = {1, 0, 2};
					whichFace = flipEdgeSwap[whichFace];
				}
				const uTind	vertices[7] = { original_face.verts(0),
								original_face.verts(1),
								original_face.verts(2),
								hObj::posFromId(vertex_mapI[original_face.verts(0)]),
								hObj::posFromId(vertex_mapI[original_face.verts(1)]),
								hObj::posFromId(vertex_mapI[original_face.verts(2)]),
								UNKNOWN
						}; // this is necessary to inform registerFace that this is Face3 not Face4
				ElemPrism & new_prism( *elements_.newObj<ElemPrism>(this, vertices) );
				new_prism.flags(GROUP_ID) = originalEl.flags(GROUP_ID);

				// 2.4.1 Creating new triangular face (external).
				Face3 & new_face( *faces_.newObj<Face3>(this, vertices + 3) );
				new_face.neighs(0) = new_prism.id_; // pointing interior
				new_face.neighs(1) = original_face.neighs(1); //pointing exterior (boundary)
				new_face.flags(B_COND) = original_face.flags(B_COND);	//boundary condition (bc number)
				new_face.flags(F_TYPE) = original_face.flags(F_TYPE);
				original_face.flags(B_COND) = -1;	// canceling bc at old external face
				original_face.neighs(1) = new_prism.id_;
				face_map[original_face.pos_] = new_face.pos_; // associate new face with old one
				new_prism.components(0) = original_face.id_;	// interior face
				new_prism.components(1) = new_face.id_; // setup prism exterior face

				// 2.4.2 Creating new faces' edges.
				for (int e(0); e < 3; ++e)
				{	// We need 3 edges, check if they exist, if not - create.
					int bc_number = 1;
					//int parallelFaceInOrigElem(-1);
					//if(originalEl.type_ == ElemT4::myType) {
					//    parallelFaceInOrigElem=ElemT4Space::faceNeigByEdge[whichFace][e];
					//}else {
					//    parallelFaceInOrigElem=ElemPrismSpace::faceNeigByEdge[whichFace][e];
					//}
					//bc_number=faces_.getById(originalEl.components(parallelFaceInOrigElem)).flags(B_COND);

					Edge & original_edge = edges_.getById(original_face.components(e));
					Edge * clone_edge(NULL);
					if (edge_mapI[original_edge.pos_] == 0)	// isn't new edge already created?
					{
						const uTind	edgeVerts[2] = {hObj::posFromId(vertex_mapI[original_edge.verts(0)]),
										hObj::posFromId(vertex_mapI[original_edge.verts(1)])
								};

						clone_edge = edges_.newObj<Edge>(this, edgeVerts);	// if not, create it
						clone_edge->level_ = 0;
						clone_edge->components(0) = vertex_mapI[original_edge.verts(0)];
						clone_edge->components(1) = vertex_mapI[original_edge.verts(1)];
						assert(clone_edge->id_ != 0);

						edge_mapI[ original_edge.pos_] = clone_edge->id_;	// mark that created

						if (vertcial_faces_mapI[hObj::posFromId(original_face.components(e))] == 0)
						{
							const uTind face4verts[4] = {original_edge.verts(0), original_edge.verts(1), clone_edge->verts(0), clone_edge->verts(1)};
							Face4 & verticalFace( *faces_.newObj<Face4>(this, face4verts) );
							verticalFace.components(0) = original_face.components(e);
							verticalFace.components(1) = new_edge_mapI[vertices_.getById(edges_.getById(original_face.components(e)).components(0)).pos_];
							verticalFace.components(2) = new_edge_mapI[vertices_.getById(edges_.getById(original_face.components(e)).components(1)).pos_];
							verticalFace.components(3) = clone_edge->id_;
							verticalFace.neighs(0) = new_prism.id_;
							verticalFace.neighs(1) = B_COND;
							verticalFace.flags(B_COND) = bc_number;
							vertcial_faces_mapI[hObj::posFromId(original_face.components(e))] = verticalFace.id_;
							//assert(verticalFace.test());
						}
						else
						{
							//assert(!"This should not happend!");
							throw "Generating boundary layer error: This should not happend!";
						}
					}
					else
					{
						hObj & face4 = faces_.getById<Face4>(vertcial_faces_mapI[hObj::posFromId(original_face.components(e))] );
						assert(face4.neighs(1) == B_COND);
						clone_edge = & edges_.getById(edge_mapI[hObj::posFromId(original_face.components(e))]);
						face4.neighs(1) = new_prism.id_;
						face4.flags(B_COND) = -1; //internal boundary
						assert(face4.neighs(1) != B_COND);
					}
					assert(clone_edge->test(this));
					new_face.components(e) = clone_edge->id_;
				}//! FACE LEVEL

				// 2.4 Rectangle vertical faces_
				new_prism.components(2) = vertcial_faces_mapI[ hObj::posFromId(original_face.components(0)) ];
				new_prism.components(3) = vertcial_faces_mapI[ hObj::posFromId(original_face.components(1)) ];
				new_prism.components(4) = vertcial_faces_mapI[ hObj::posFromId(original_face.components(2)) ];
				// Testing
				new_prism.test(this);
				normalizeElem(new_prism);
			}
			// 3. If this is NOT a last loop run, replace old external points with new. Do the same for faces.
			if (layer < nLayers)
			{
				const size_t nPoints(bound_points.size()),
					  nFaces(bound_faces.size());
				bound_points.erase(bound_points.begin(), bound_points.end());

				Vertex * v (& vertices_[static_cast<uTind>(vertices_.size() - nPoints + 1)]);
				for (size_t i(0); i < nPoints; ++i, ++v)
				{
					bound_points.push_back(v->pos_);
				}

				//hObj * f (& faces_[static_cast<uTind>(faces_.size()-nFaces+1)]);
				//FacePool::constIterator<>	f(&faces_,faces_.size()-nFaces+1);
				//bound_faces.erase(bound_faces.begin(), bound_faces.end());
				for (size_t i(0); i < nFaces; ++i)
				{
					bound_faces[i] = face_map[bound_faces[i]];
				}
			}
			checkUniqueness();
		}
		// check Face(4) at boundary and normalize if necessary
		for (FacePool::Iterator<> it(&faces_); !it.done(); ++it) {
			if (it->neighs(1) == B_COND) {
				hObj & elem(elements_.getById(it->neighs(0)));
				normalizeFace(*it, elem);
			}
		}
		// check all
		mmv_out_stream << "\nDone. Checking... ";
	}
	catch (char * err)
	{
		mmv_out_stream << "\nException (char * ): " << err;
	}
	return checkAll();
}


std::ofstream &	operator << (std::ofstream & stream, const hHybridMesh & mesh)
{
	if (stream.good())
	{
		//stream << mesh.name << " " << mesh.mesh_id << std::endl;
		//stream << mesh.vertices_;
		//stream << mesh.edges_;
		//stream << mesh.faces_;
		//stream << mesh.elements;
	}
	return stream;
}

std::ifstream &	operator >> (std::ifstream & stream, hHybridMesh & mesh)
{
	if (stream.good())
	{
		//stream >> mesh.name;
		//stream >> mesh.vertices_;
		//stream >> mesh.edges_;
		//stream >> mesh.faces_;
		//stream >> mesh.elements;
	}
	return stream;
}

void    hHybridMesh::print() const
{
	std::cout << "Elements:" << elements_.size() << "\n";
	for (ElemPool::constIterator<> iter(&elements_); !iter.done(); ++iter)
	{
		(*iter).print();
	}

	std::cout << "Faces:" << faces_.size() << "\n";
	for (FacePool::constIterator<> iter(&faces_); !iter.done(); ++iter)
	{
		(*iter).print();
	}

	std::cout << "Edges:" << edges_.size() << "\n";
	for (EdgePool::constIterator<> iter(&edges_); !iter.done(); ++iter)
	{
		(*iter).print();
	}

	std::cout << "Vertices:" << vertices_.size() << "\n";
	for (VertexPool::constIterator<> iter(&vertices_); !iter.done(); ++iter)
	{
		(*iter).print();
	}

}

///////////////////////////////////////////////////////////////////////////////

// ID	 VertexSpace::uniqueId()	{ assert(hObj::myMesh != NULL); return hObj::myMesh->vertices_.uniqueID();}
// ID	 EdgeSpace::uniqueId()		{ assert(hObj::myMesh != NULL); return hObj::myMesh->edges_.uniqueID();}
// ID	 Face3Space::uniqueId()		{ assert(hObj::myMesh != NULL); return hObj::myMesh->faces_.uniqueID();}
// ID	 Face4Space::uniqueId()		{ assert(hObj::myMesh != NULL); return hObj::myMesh->faces_.uniqueID();}
// ID	 ElemPrismSpace::uniqueId()	{ assert(hObj::myMesh != NULL); return hObj::myMesh->elements_.uniqueID();}
// ID	 ElemT4Space::uniqueId()	{ assert(hObj::myMesh != NULL); return hObj::myMesh->elements_.uniqueID();}

