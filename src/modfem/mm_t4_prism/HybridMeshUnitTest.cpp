//// UnitTests using Boost.Test

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */

#include "MeshModule/hHybridMesh.h"
#include "MeshModule/Face4.h"
#include "MeshModule/ElemT4.h"

#include <modfem/mmh_intf.h>
#include "MeshRead/NasFileImporter.h"
#include "MeshRead/BinaryFileReader.h"
#include "MeshRead/InFileImporter.h"
#include "MeshRead/MshFileImporter.h"
#include "MeshWrite/BinaryFileWriter.h"

#define BOOST_TEST_MODULE HybridMeshUnitTest
#include <boost/test/included/unit_test.hpp>

using namespace MeshRead;

BOOST_AUTO_TEST_SUITE(MeshModule)

const uTind verts0[4] = {1, 2, 3, 4};
const uTind verts1[6] = {5, 6, 7, 8, 9, 10};
const uTind verts2[4] = {11, 12, 13, 14};
const uTind verts3[6] = {15, 16, 17, 1, 2, 3};
const uTind verts4[4] = {18, 5, 6, 7};

const Tval pointsCoords[5][3] = {
	{1, 2, 3},
	{11, 12, 13},
	{21, 22, 23},
	{31, 32, 33},
	{41, 42, 43}
};

BOOST_AUTO_TEST_CASE( hObjTest_basic )
{
	BOOST_CHECK( hObj::myMesh == NULL );
}


BOOST_AUTO_TEST_CASE( VtsSqId_test  )
{
	VtsSqId<4> id(verts0), id2(verts4);
	VtsSqId<4> id1(id);
	//std::cout << id << "\n" << id1 << "\n" << id2 << "\n";


	BOOST_CHECK_EQUAL(id, id1);
	BOOST_CHECK(id >= id1);
	BOOST_CHECK(id <= id1);
	BOOST_CHECK_EQUAL(id != id1, false);

	BOOST_CHECK(id != id2);
	BOOST_CHECK(id < id2);
	BOOST_CHECK(id <= id2);
	BOOST_CHECK(id2 > id);
	BOOST_CHECK(id2 >= id);

}


BOOST_AUTO_TEST_CASE( ArrayPool_test )
{
	ArrayPool<Vertex>  ar;

	BOOST_CHECK_EQUAL(ar.totalMemory(), 0);
	BOOST_CHECK_EQUAL(ar.usedMemory(), 0);
	BOOST_CHECK_EQUAL(ar.size(), 0);
	BOOST_CHECK_EQUAL(ar.last(), 0);
	BOOST_CHECK_EQUAL(ar.capacity(), 0);
	BOOST_CHECK_EQUAL(ar.empty(), true);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 0);
	//BOOST_CHECK_EQUAL(ar.end(),(void*)NULL);

	ar.reserve(4);

	BOOST_CHECK_EQUAL(ar.totalMemory(), 4 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.usedMemory(), 0);
	BOOST_CHECK_EQUAL(ar.size(), 0);
	BOOST_CHECK_EQUAL(ar.capacity(), 4);
	BOOST_CHECK_EQUAL(ar.empty(), true);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 0);
	//BOOST_CHECK_EQUAL(ar.end(),NULL);

	const  Vertex & v0 = *ar.newObj(verts0);
	BOOST_CHECK(&v0 != NULL);
	BOOST_CHECK(ar.first() == v0);
	BOOST_CHECK(ar.at(v0.pos_) == v0);
	BOOST_CHECK(ar.getById(v0.id_) == v0);
	BOOST_CHECK(ar[v0.pos_] == v0);

	BOOST_CHECK_EQUAL(ar.totalMemory(), 4 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.usedMemory(), sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.size(), 1);
	BOOST_CHECK_EQUAL(ar.capacity(), 4);
	BOOST_CHECK_EQUAL(ar.empty(), false);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 1 * sizeof(Vertex));

	const Vertex & v1 = *ar.newObj(verts1);
	BOOST_CHECK(ar.at(v1.pos_) == v1);
	BOOST_CHECK(ar.getById(v1.id_) == v1);
	BOOST_CHECK(ar[v1.pos_] == v1);

	const Vertex & v2 = *ar.newObj(verts2);
	BOOST_CHECK(ar.at(v2.pos_) == v2);
	BOOST_CHECK(ar.getById(v2.id_) == v2);
	BOOST_CHECK(ar[v2.pos_] == v2);

	const  Vertex & v3 = *ar.newObj(verts3);
	BOOST_CHECK(ar.at(v3.pos_) == v3);
	BOOST_CHECK(ar.getById(v3.id_) == v3);
	BOOST_CHECK(ar[v3.pos_] == v3);

	BOOST_CHECK_EQUAL(ar.totalMemory(), 4 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.usedMemory(), 4 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.size(), 4);
	BOOST_CHECK_EQUAL(ar.capacity(), 4);
	BOOST_CHECK_EQUAL(ar.empty(), false);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 4 * sizeof(Vertex));

	ar.requestChange(1);

	BOOST_CHECK_EQUAL(ar.memoryNeeded(), (ar.size() + 1)*sizeof(Vertex));

	ar.adjust();

	BOOST_CHECK_EQUAL(ar.totalMemory(), 5 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.usedMemory(), 4 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.size(), 4);
	BOOST_CHECK_EQUAL(ar.capacity(), 5);
	BOOST_CHECK_EQUAL(ar.empty(), false);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 4 * sizeof(Vertex));

	const Vertex & v4 = *ar.newObj(verts4);

	BOOST_CHECK_EQUAL(ar.totalMemory(), 5 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.usedMemory(), 5 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.size(), 5);
	BOOST_CHECK_EQUAL(ar.capacity(), 5);
	BOOST_CHECK_EQUAL(ar.empty(), false);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 5 * sizeof(Vertex));

	for (int i(FIRST); i <= ar.size(); ++i) {
		ar[i].coords_[0] = pointsCoords[i - FIRST][0];
		ar[i].coords_[1] = pointsCoords[i - FIRST][1];
		ar[i].coords_[2] = pointsCoords[i - FIRST][2];
	}

	for (int i(FIRST); i <= ar.size(); ++i) {
		BOOST_CHECK_EQUAL(ar[i].coords_[0], pointsCoords[i - FIRST][0]);
		BOOST_CHECK_EQUAL(ar[i].coords_[1], pointsCoords[i - FIRST][1]);
		BOOST_CHECK_EQUAL(ar[i].coords_[2], pointsCoords[i - FIRST][2]);
	}

	const char fileName[] = {"ArrayPool.test"};
	std::ofstream file_(fileName);
	file_.clear();
	BOOST_CHECK(file_.good());
	ar.write(file_);
	file_.close();

	ar.clear();
	BOOST_CHECK_EQUAL(ar.totalMemory(), 0);
	BOOST_CHECK_EQUAL(ar.usedMemory(), 0);
	BOOST_CHECK_EQUAL(ar.size(), 0);
	BOOST_CHECK_EQUAL(ar.capacity(), 0);
	BOOST_CHECK_EQUAL(ar.empty(), true);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 0);

	std::ifstream file2_ (fileName);
	file2_.clear();
	BOOST_CHECK(file2_.good());
	ar.read(file2_);
	file2_.close();

	BOOST_CHECK_EQUAL(ar.totalMemory(), 5 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.usedMemory(), 5 * sizeof(Vertex));
	BOOST_CHECK_EQUAL(ar.size(), 5);
	BOOST_CHECK_EQUAL(ar.capacity(), 5);
	BOOST_CHECK_EQUAL(ar.empty(), false);
	BOOST_CHECK_EQUAL(ar.memoryNeeded(), 5 * sizeof(Vertex));

	for (int i(FIRST); i <= ar.size(); ++i) {
		BOOST_CHECK_EQUAL(ar[i].coords_[0], pointsCoords[i - FIRST][0]);
		BOOST_CHECK_EQUAL(ar[i].coords_[1], pointsCoords[i - FIRST][1]);
		BOOST_CHECK_EQUAL(ar[i].coords_[2], pointsCoords[i - FIRST][2]);
	}
}


BOOST_AUTO_TEST_CASE( StaticPool_test )
{
	StaticPool<hObj, Face4::nVerts> sp;
	BOOST_CHECK_EQUAL(sp.capacity(), 0);
	BOOST_CHECK_EQUAL(sp.size(), 0);
	BOOST_CHECK_EQUAL(sp.last(), 0);
	BOOST_CHECK_EQUAL(sp.empty(), true);
	BOOST_CHECK_EQUAL(sp.totalMemory(), 0);
	BOOST_CHECK_EQUAL(sp.usedMemory(), 0);
	BOOST_CHECK_EQUAL(sp.memoryNeeded(), 0);
	BOOST_CHECK_EQUAL(sp.nonDividedObjs(), 0);

	const int totalMem = 3 * sizeof(Face4) + 2 * sizeof(Face3);

	sp.reserve(5, totalMem);
	BOOST_CHECK_EQUAL(sp.capacity(), 5);
	BOOST_CHECK_EQUAL(sp.size(), 0);
	BOOST_CHECK_EQUAL(sp.last(), 0);
	BOOST_CHECK_EQUAL(sp.empty(), true);
	BOOST_CHECK_EQUAL(sp.totalMemory(), totalMem);
	BOOST_CHECK_EQUAL(sp.usedMemory(), 0);
	BOOST_CHECK_EQUAL(sp.memoryNeeded(), 0);
	BOOST_CHECK_EQUAL(sp.nonDividedObjs(), 0);

	const Face4 & f0 = *sp.newObj<Face4>(verts0);
	BOOST_CHECK_EQUAL(sp.capacity(), 5);
	BOOST_CHECK_EQUAL(sp.size(), 1);
	BOOST_CHECK_EQUAL(sp.last(), 1);
	BOOST_CHECK_EQUAL(sp.empty(), false);
	BOOST_CHECK_EQUAL(sp.totalMemory(), totalMem);
	BOOST_CHECK_EQUAL(sp.usedMemory(), sizeof(Face4));
	BOOST_CHECK_EQUAL(sp.memoryNeeded(), sizeof(Face4));
	BOOST_CHECK_EQUAL(sp.nonDividedObjs(), 1);

	const Face4 & f1 = *sp.newObj<Face4>(verts1);
	BOOST_CHECK_EQUAL(sp.capacity(), 5);
	BOOST_CHECK_EQUAL(sp.size(), 2);
	BOOST_CHECK_EQUAL(sp.last(), 2);
	BOOST_CHECK_EQUAL(sp.empty(), false);
	BOOST_CHECK_EQUAL(sp.totalMemory(), totalMem);
	BOOST_CHECK_EQUAL(sp.usedMemory(), 2 * sizeof(Face4));
	BOOST_CHECK_EQUAL(sp.memoryNeeded(), 2 * sizeof(Face4));
	BOOST_CHECK_EQUAL(sp.nonDividedObjs(), 2);
	BOOST_CHECK(f1 == sp[f1.pos_]);

	const Face3 & f2 = *sp.newObj<Face3>(verts2);
	BOOST_CHECK(f2 == sp[f2.pos_]);
	const Face3 & f3 = *sp.newObj<Face3>(verts3);
	BOOST_CHECK(f3 == sp[f3.pos_]);
	const Face4 & f4 = *sp.newObj<Face4>(verts4);
	BOOST_CHECK(f4 == sp[f4.pos_]);
	BOOST_CHECK_EQUAL(sp.capacity(), 5);
	BOOST_CHECK_EQUAL(sp.size(), 5);
	BOOST_CHECK_EQUAL(sp.last(), 5);
	BOOST_CHECK_EQUAL(sp.empty(), false);
	BOOST_CHECK_EQUAL(sp.totalMemory(), totalMem );
	BOOST_CHECK_EQUAL(sp.usedMemory(), totalMem );
	BOOST_CHECK_EQUAL(sp.memoryNeeded(), totalMem );
	BOOST_CHECK_EQUAL(sp.nonDividedObjs(), 5);

	// BYTE - write and read
	const char fileName[] = {"StaticPool.test"};
	std::ofstream file_(fileName);
	file_.clear();
	BOOST_CHECK(file_.good());
	sp.write(file_);
	file_.close();
	BOOST_CHECK(sp.at<Face4>(1).equals(f0));
	BOOST_CHECK(sp.at<Face4>(2).equals(f1));
	BOOST_CHECK(sp.at<Face3>(3).equals(f2));
	BOOST_CHECK(sp.at<Face3>(4).equals(f3));
	BOOST_CHECK(sp.at<Face4>(5).equals(f4));

	sp.clear();

	BOOST_CHECK_EQUAL(sp.capacity(), 0);
	BOOST_CHECK_EQUAL(sp.size(), 0);
	BOOST_CHECK_EQUAL(sp.last(), 0);
	BOOST_CHECK_EQUAL(sp.empty(), true);
	BOOST_CHECK_EQUAL(sp.totalMemory(), 0);
	BOOST_CHECK_EQUAL(sp.usedMemory(), 0);
	BOOST_CHECK_EQUAL(sp.memoryNeeded(), 0);
	BOOST_CHECK_EQUAL(sp.nonDividedObjs(), 0);

	std::ifstream file2_(fileName, std::ios::in | std::ios::binary);
	BOOST_CHECK(file2_.is_open());
	BOOST_CHECK(file2_.good());
	sp.read(file2_);
	file2_.close();

	BOOST_CHECK_EQUAL(sp.capacity(), 5);
	BOOST_CHECK_EQUAL(sp.size(), 5);
	BOOST_CHECK_EQUAL(sp.last(), 5);
	BOOST_CHECK_EQUAL(sp.empty(), false);
	BOOST_CHECK_EQUAL(sp.totalMemory(), totalMem  );
	BOOST_CHECK_EQUAL(sp.usedMemory(), totalMem  );
	BOOST_CHECK_EQUAL(sp.memoryNeeded(), totalMem );
	BOOST_CHECK_EQUAL(sp.nonDividedObjs(), 5);

	//  BOOST_CHECK(sp.at<Face4>(1).equals(f0));
	//BOOST_CHECK(sp.at<Face4>(2).equals(f1));
	//BOOST_CHECK(sp.at<Face4>(3).equals(f2));
	//BOOST_CHECK(sp.at<Face4>(4).equals(f3));
	//BOOST_CHECK(sp.at<Face4>(5).equals(f4));

}

BOOST_AUTO_TEST_CASE( hHybridMesh_test  )
{
	using namespace MeshRead;
	char filename[] = {"mesh.nas"};
	hHybridMesh hm;
	IMeshReader * r = new NasFileImporter(filename);
	hm.read(*r);
	safeDelete(r);

	const int nEdges = hm.edges_.size();
	for (hHybridMesh::ElemPool::constIterator<> el(& hm.elements_); !el.done(); ++el) {
		BOOST_CHECK_EQUAL(nEdges, hm.edges_.size());
		BOOST_CHECK(hm.edge(el->verts(0), el->verts(1)).pos_ <= nEdges);
		BOOST_CHECK(hm.edge(el->verts(0), el->verts(2)).pos_ <= nEdges);
		BOOST_CHECK(hm.edge(el->verts(0), el->verts(3)).pos_ <= nEdges);
		BOOST_CHECK(hm.edge(el->verts(1), el->verts(2)).pos_ <= nEdges);
		switch (el->type_) {
			case ElemT4::myType:
				BOOST_CHECK(hm.edge(el->verts(1), el->verts(3)).pos_ <= nEdges);
				BOOST_CHECK(hm.edge(el->verts(2), el->verts(3)).pos_ <= nEdges);
				break;
			case ElemPrism::myType:
				BOOST_CHECK(hm.edge(el->verts(1), el->verts(4)).pos_ <= nEdges);
				BOOST_CHECK(hm.edge(el->verts(2), el->verts(5)).pos_ <= nEdges);
				BOOST_CHECK(hm.edge(el->verts(3), el->verts(4)).pos_ <= nEdges);
				BOOST_CHECK(hm.edge(el->verts(3), el->verts(5)).pos_ <= nEdges);
				BOOST_CHECK(hm.edge(el->verts(4), el->verts(5)).pos_ <= nEdges);
				break;
		}
		BOOST_CHECK_EQUAL(nEdges, hm.edges_.size());

	}

	/// testing refinement
	hm.initRefine();
	hm.refineElem(*hm.begin());
	hm.finalRef();

	const char binFileName[] = {"test_mesh.bin"};
	MeshWrite::BinaryFileWriter binWriter(binFileName);
	hm.write(binWriter);
	binWriter.Free();
	hm.free();
	MeshRead::BinaryFileReader binReader(binFileName);
	hm.read(binReader);
	binReader.Free();
}

BOOST_AUTO_TEST_CASE( MshFileImporter )
{
	hHybridMesh hm;
	const std::string file("model_2011.09.26.msh");
	IMeshReader * r = new MeshRead::MshFileImporter(file);
	BOOST_CHECK(hm.read(*r));
	safeDelete(r);
}

BOOST_AUTO_TEST_CASE( InFileImporter )
{

	//hHybridMesh hm;
	//const std::string fn("mesh.in");
	//IMeshReader *r = new MeshRead::InFileImporter(fn);
	//BOOST_CHECK(hm.read(*r));
	//safeDelete(r);
}

BOOST_AUTO_TEST_CASE( mmd_t4_prism )
{
	char filename[] = {"mesh.nas"};

	const int mesh_id = mmr_init_mesh(MMC_MOD_FEM_HYBRID_DATA, filename);

	const int BCs[4] = {3, 4, 6, 9};
	const int nBCs = 4;

	mmr_init_dist2bound(mesh_id, BCs, nBCs);


	double coords[3] = {0, 0, 0};
	double dist2 = mmr_get_el_dist2boundary(mesh_id, 1, coords);
	BOOST_CHECK(dist2 == -1.0);

	//// adaptation

}

BOOST_AUTO_TEST_SUITE_END()

//BOOST_AUTO_TEST_SUITE(MeshWrite_MeshWrite)
//BOOST_AUTO_TEST_SUITE_END()

/** @}*/
