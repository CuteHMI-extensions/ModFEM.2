#include <cstdlib>
#include <modfem/mm_remesh/ticooMesh3D.h>

#include <ctime>


double frandd(double wielkosc) {
	return (double)rand() / (RAND_MAX) * (wielkosc);
}

int frandi(double wielkosc) {

	return int((double)rand() / (RAND_MAX) * (wielkosc));
}

bool TicooMesh3D::empty()
{
	return getNumberPoints() == 8 && getNumberElements() == 6;
}


void TicooMesh3D::zmienWarunkeSciana(double r_max, double r_min, int zWar, int naWar) {

	siatka.zmienWarunkeSciana(r_max, r_min, zWar, naWar);

}


void TicooMesh3D::mesh_quality_stats(FILE * out) {

	siatka.stosunekRdor(out);

}

void TicooMesh3D::count_bc_surface(FILE * out) {

	siatka.count_bc_surface(out);

}


int TicooMesh3D::mmr_grups_number() {

	return siatka.materialNumber();

}
int TicooMesh3D::mmr_grups_ids(int * tab, int l_tab) {

	return siatka.materialIDs(tab, l_tab);

}


//void TicooMesh3D::mmr_set_fluid(int el_id,char a){

//	siatka.elements.getElement(el_id).setFluid(a);


//}
//int TicooMesh3D::mmr_is_fluid(int el_id){

//	return siatka.elements.getElement(el_id).isFluid();

//}


void TicooMesh3D::utr_smoothing(double l_step, int repetition) {

	siatka.smoothing(l_step, repetition);

}

void TicooMesh3D::mmr_split_into_blocks_add_contact(const char * workdir, int * tabMat, int l_mat, int * tabBlock, int l_block, int * war, int l_war, double * tempBlock, int l_tempBlock, int * groups_fluid, int l_groups_fluid) {




	string a = workdir;



//std::cout<<"-----------------!!!!!!!!!!!!!!!!!!----------------------------";
//std::cout<<a;
	if ((workdir[a.size() - 1] == '.')) {
		a = "init_temp.dat";
	}
	else if (workdir[a.size() - 1] == '/') {
		a += "init_temp.dat";
	}
	else {
		a += "/init_temp.dat";
	}


//double tempB[4]={1.0,295.0, 2.0,1102.0},dl_tempB=4;
	/*
	int tabBlock[16]={1,1, 2,1, 3,1, 4,1, 5,2, 6,2, 7,2, 8,2},l_block=16;
	int tabMat[16]={1,1, 3,1, 5,1, 7,1, 2,2, 4,2, 6,3, 8,3},l_mat=16;
	int war[24]={2,2,3,1,    2,1,-1,2,   1,2,-1,3,    0,1,5,4,     0,3,7,4,     0,1,-1,5},l_war=24;
	*/
	/*
	double tempB[4]={2.0,295.0, 2.0,1102.0},dl_tempB=4;
	int tabBlock1[4]={1,1, 2,2},dl_block=4;
	int tabMat1[16]={1,1,2,1},dl_mat=4;
	int war1[24]={1,1,2,2},dl_war=4;
	*/

//zmienWarunkeSciana(0.235,0.213375,26,30);
//zmienWarunkeSciana(0.235,0.205,26,30);
//ZapiszDoPlikuNasZPrzesunieciem_BC_Face("new_mesh.dat");

	for (int i = 0; i < l_groups_fluid; ++i) {

		siatka.fluidGroup.setElement(groups_fluid[i]);

	}

	for (int i = 0; i < l_block; i += 2) {

		if (siatka.fluidGroup.sprawdzCzyJest(tabBlock[i])) {
			siatka.fluidBlock.setElement(tabBlock[i + 1]);
		}

	}
// Na razie wpisane ręcznie
	siatka.uppermostBlock.setElement(4287);

	siatka.podzielNaObjekty(tabBlock, l_block);
	siatka.sasiednieElementy();
//siatka.usunPunktyNiepotrzebne();
	siatka.zapiszWarunkiPoczatkoweNaPodstawieDanych(a.c_str(), tempBlock, l_tempBlock, tabBlock, l_block);

	siatka.sasiednieElementy();
	siatka.creatFaceAndEdgeRememberBC();
	siatka.tablicaFC_connect();
	siatka.sasiednieElementy();
	siatka.warunkiNaPodstawieDanych(tabMat, l_mat, tabBlock, l_block, war, l_war);

	siatka.wyszukajNormalnePowierzchniDlaPunktowGranicznych();


}

int TicooMesh3D::mmr_el_groupID(int idEl) {
	return siatka.elements.getElement(idEl).getRodzajZiarna();
}

int TicooMesh3D::mmr_get_fa_el_bc_connect(int face_id, int * el_id) {

	int fa_connect = siatka.get_facebc_connect(face_id), zmien = 0;;

	if (el_id[0] == -1) {
		zmien = 1;
	}
	if (fa_connect != -1) {
		el_id[0] = siatka.faces.getElement(face_id).getInEl() + 1;
	}
	else {
		el_id[0] = 0;
	}


	//if(siatka.elements.getElement(siatka.faces.getElement(face_id).getInEl()).getRodzajZiarna()==3){return 0;}

	if (fa_connect != -1)siatka.mmr_dopasuj_punkty(face_id, fa_connect, el_id, zmien);

	return fa_connect + 1;

}

void TicooMesh3D::mmr_get_coor_from_motion_element(int idEl, int idLP, double * coor, int flagaSiatki) {

	switch (flagaSiatki) {

		case 0:
			siatka.polozeniePunktowWybranegoEl(idEl, idLP, coor);
			break;

		case 1:
			siatka1.polozeniePunktowWybranegoEl(idEl, idLP, coor);
			break;

	}

}

void TicooMesh3D::copyMesh(int flaga) {
// 0 \\ siatka = siatka1
// 1 \\ siatka1 = siatka

	IntList a;
	a.setElement(4287);

	switch (flaga) {

		case 0:
			siatka1.kopiujSiatke(siatka);
			break;

		case 1:
			siatka.kopiujSiatke(siatka1);
			break;
		case 2:
			siatka.info_bc_connect();
			break;
		case 3:
			siatka.kopiujSiatke(siatka1);

			//siatka.cutlowestGroupID(temp1);
			siatka.usunWybraneElementy(true, false, siatka.fluidGroup);
			break;
		case 4:
			siatka.kopiujSiatke(siatka1);


			//siatka.cutlowestGroupID(temp1);
			siatka.usunWybraneElementy(false, false, a);
			break;

	}

}
/*
void TicooMesh3D::mmr_set_new_and_get_old_uk_val(int idEl,int idGP,double *new_uk_val,double *old_uk_val){


if(old_uk_val!=NULL){siatka.elements.getElement(idEl).get_uk_val(idGP,old_uk_val);}
if(new_uk_val!=NULL){siatka.elements.getElement(idEl).set_uk_val(idGP,new_uk_val);}

}
*/

void TicooMesh3D::getElwithFace(int Fa, int * El, int * Neig_sides) {

	if (El != NULL) {

		El[0] = siatka.faces.getElement(Fa).getInEl();

		char nf = siatka.faces.getElement(Fa).getFaceInEl();

		switch (nf) {
			case 'a':
				El[1] = siatka.elements.getElement(El[0]).getE1();

				//if(El[1]<-1){El[1]=-El[1]-2;}
				if (El[1] < 0) {
					El[1] = -1;
				}

				if (Neig_sides != NULL) {
					Neig_sides[0] = 0;
					Neig_sides[1] = siatka.elements.getNumerSasSacianE1(El[0]);

				}
				break;
			case 'b':
				El[1] = siatka.elements.getElement(El[0]).getE2();

				//if(El[1]<-1){El[1]=-El[1]-2;}
				if (El[1] < 0) {
					El[1] = -1;
				}
				if (Neig_sides != NULL) {

					Neig_sides[0] = 1;
					Neig_sides[1] = siatka.elements.getNumerSasSacianE2(El[0]);

				}
				break;
			case 'c':
				El[1] = siatka.elements.getElement(El[0]).getE3();

				//if(El[1]<-1){El[1]=-El[1]-2;}
				if (El[1] < 0) {
					El[1] = -1;
				}
				if (Neig_sides != NULL) {
					Neig_sides[0] = 2;
					Neig_sides[1] = siatka.elements.getNumerSasSacianE3(El[0]);

				}
				break;
			case 'd':
				El[1] = siatka.elements.getElement(El[0]).getE4();

				//if(El[1]<-1){El[1]=-El[1]-2;}
				if (El[1] < 0) {
					El[1] = -1;
				}
				if (Neig_sides != NULL) {
					Neig_sides[0] = 3;
					Neig_sides[1] = siatka.elements.getNumerSasSacianE4(El[0]);

				}
				break;
				//switch
		}

		El[0]++;
		El[1]++;
	}

}

void TicooMesh3D::mmr_fa_edges(int Fa, int * Fa_edges, int * Ed_orient) {

	if (Fa_edges != NULL) {

		Fa_edges[0] = 3;

		Fa_edges[1] = siatka.faces.getElement(Fa).getEd1();
		Fa_edges[2] = siatka.faces.getElement(Fa).getEd2();
		Fa_edges[3] = siatka.faces.getElement(Fa).getEd3();
	}

	if (Ed_orient != NULL) {

		Ed_orient[0] = 1;
		Ed_orient[0] = 2;
		Ed_orient[0] = 3;

	}

}

void TicooMesh3D::mmr_edge_nodes(int Ed, int * Edge_nodes) {

	if (Edge_nodes != NULL) {
		Edge_nodes[0] = siatka.edges.getElement(Ed).getP1();
		Edge_nodes[1] = siatka.edges.getElement(Ed).getP2();
	}
}

int TicooMesh3D::mmr_el_edges(int El, int * Edges) {

	int b = 0;

	if (Edges != NULL) {
		b = siatka.wyszukajEdgesInEl(El, Edges);
	}

	return b;

}
int TicooMesh3D::mmr_edge_elems(int Ed, int * Edge_elems) {

	int b = 0;
	if (Edge_elems != NULL) {


		//b = siatka.wyszukajElementsWithEdge(Ed,Edge_elems);
		siatka.wyszukajElementsWithEdge(Ed, Edge_elems);

	}

	return b;
}

void TicooMesh3D::creatEdgeFace() {

	siatka.creatFaceAndEdge();

}

void TicooMesh3D::czyscEdgeFace() {

	siatka.czyscFaceAndEdge();

}

void TicooMesh3D::czyscEdge() {

	siatka.czyscEdge();

}

void TicooMesh3D::czyscFace() {

	siatka.czyscFace();

}

int TicooMesh3D::mmr_el_faces(int El, int * Face, int * Orient) {

	if (Face != NULL) {

		Face[0] = 4;
		Face[1] = siatka.mapFace[El].getElement(0);
		Face[2] = siatka.mapFace[El].getElement(1);
		Face[3] = siatka.mapFace[El].getElement(2);
		Face[4] = siatka.mapFace[El].getElement(3);

		if (Orient != NULL) {
			if (siatka.faces.getElement(Face[1]).getInEl() == El) {
				Orient[0] = 1;
			}
			else {
				Orient[0] = -1;
			}

			if (siatka.faces.getElement(Face[2]).getInEl() == El) {
				Orient[1] = 1;
			}
			else {
				Orient[1] = -1;
			}

			if (siatka.faces.getElement(Face[3]).getInEl() == El) {
				Orient[2] = 1;
			}
			else {
				Orient[2] = -1;
			}

			if (siatka.faces.getElement(Face[4]).getInEl() == El) {
				Orient[3] = 1;
			}
			else {
				Orient[3] = -1;
			}
		}

	}
	else if (Orient != NULL) {

		int Fac[5];
		Fac[0] = 4;
		Fac[1] = siatka.mapFace[El].getElement(0);
		Fac[2] = siatka.mapFace[El].getElement(1);
		Fac[3] = siatka.mapFace[El].getElement(2);
		Fac[4] = siatka.mapFace[El].getElement(3);


		if (siatka.faces.getElement(Fac[1]).getInEl() == El) {
			Orient[0] = 1;
		}
		else {
			Orient[0] = -1;
		}

		if (siatka.faces.getElement(Fac[2]).getInEl() == El) {
			Orient[1] = 1;
		}
		else {
			Orient[1] = -1;
		}

		if (siatka.faces.getElement(Fac[3]).getInEl() == El) {
			Orient[2] = 1;
		}
		else {
			Orient[2] = -1;
		}

		if (siatka.faces.getElement(Fac[4]).getInEl() == El) {
			Orient[3] = 1;
		}
		else {
			Orient[3] = -1;
		}
	}

	return 4;


}

void TicooMesh3D::mmr_node_coor(int Node, double * Coor) {

	Punkt & p = siatka.points.getElement(Node);

	if (Coor != NULL) {
		Coor[0] = p.getX();
		Coor[1] = p.getY();
		Coor[2] = p.getZ();
	}
}

int TicooMesh3D::mmr_node_status(int Node) {

	if (siatka.points.getElement(Node).getGranica() == 'z') {
		return 0;
	}
	return 1;

}


void TicooMesh3D::mmr_fa_node_coor(int Fa, int * Nodes, double * Node_coor) {

	Face & fa1 = siatka.faces.getElement(Fa);
	PunktList & p = siatka.zwrocRefPunktow();
	siatka.faces.getElement(Fa).getInEl();

	if (Nodes != NULL) {


		Nodes[0] = 3;
		int nr_e = fa1.getInEl();
		Element & el = siatka.elements.getElement(nr_e);
		char nf = fa1.getFaceInEl();

		switch (nf) {
			case 'a':
				Nodes[1] = el.getP4();
				Nodes[2] = el.getP2();
				Nodes[3] = el.getP1();
				break;
			case 'b':
				Nodes[1] = el.getP4();
				Nodes[2] = el.getP3();
				Nodes[3] = el.getP2();
				break;
			case 'c':
				Nodes[1] = el.getP4();
				Nodes[2] = el.getP1();
				Nodes[3] = el.getP3();
				break;
			case 'd':
				Nodes[1] = el.getP1();
				Nodes[2] = el.getP2();
				Nodes[3] = el.getP3();
				break;
				//switch
		}

		if (Node_coor != NULL) {

			for (int i = 1, it = 0, ile = Nodes[0] + 1; i < ile; ++i) {

				Node_coor[it++] = p.getElement(Nodes[i]).getX();
				Node_coor[it++] = p.getElement(Nodes[i]).getY();
				Node_coor[it++] = p.getElement(Nodes[i]).getZ();

			}
		}

	}
	else if (Node_coor != NULL) {

		int Node[4];
		Node[0] = 3;
		int nr_e = fa1.getInEl();
		Element & el = siatka.elements.getElement(nr_e);
		char nf = fa1.getFaceInEl();

		switch (nf) {
			case 'a':
				Node[1] = el.getP4();
				Node[2] = el.getP2();
				Node[3] = el.getP1();
				break;
			case 'b':
				Node[1] = el.getP4();
				Node[2] = el.getP3();
				Node[3] = el.getP2();
				break;
			case 'c':
				Node[1] = el.getP4();
				Node[2] = el.getP1();
				Node[3] = el.getP3();
				break;
			case 'd':
				Node[1] = el.getP1();
				Node[2] = el.getP2();
				Node[3] = el.getP3();
				break;
				//switch
		}

		for (int i = 1, it = 0, ile = Node[0] + 1; i < ile; ++i) {

			Node_coor[it++] = p.getElement(Node[i]).getX();
			Node_coor[it++] = p.getElement(Node[i]).getY();
			Node_coor[it++] = p.getElement(Node[i]).getZ();

		}
	}

}

void TicooMesh3D::mmr_get_vec_norm_at_node(int node, double * vec_norm) {

	PunktList & p = siatka.zwrocRefPunktow();

	p.getElement(node).getVecNorm(vec_norm);



}

void TicooMesh3D::mmr_el_node_coor(int Nel, int * Nodes, double * Node_coor) {

	Element & el = siatka.elements.getElement(Nel);
	PunktList & p = siatka.zwrocRefPunktow();

	if (Nodes != NULL) {
		Nodes[0] = 4;
		Nodes[1] = el.getP1();
		Nodes[2] = el.getP2();
		Nodes[3] = el.getP3();
		Nodes[4] = el.getP4();

		if (Node_coor != NULL) {
			for (int i = 1, it = 0, ile = Nodes[0] + 1; i < ile; ++i) {

				Node_coor[it++] = p.getElement(Nodes[i]).getX();
				Node_coor[it++] = p.getElement(Nodes[i]).getY();
				Node_coor[it++] = p.getElement(Nodes[i]).getZ();

			}
		}
	}
	else if (Node_coor != NULL) {

		int Node[5];
		Node[0]	= 4;
		Node[1] = el.getP1();
		Node[2] = el.getP2();
		Node[3] = el.getP3();
		Node[4] = el.getP4();

		for (int i = 1, it = 0, ile = Node[0] + 1; i < ile; ++i) {

			Node_coor[it++] = p.getElement(Node[i]).getX();
			Node_coor[it++] = p.getElement(Node[i]).getY();
			Node_coor[it++] = p.getElement(Node[i]).getZ();

		}

	}


}

void TicooMesh3D::mmr_fa_area(int Fa, double * Area, double * Vec_norm) {


	Element & el = siatka.elements.getElement(siatka.faces.getElement(Fa).getInEl());
	char nf = siatka.faces.getElement(Fa).getFaceInEl();
	int p1 = 0, p2 = 0, p3 = 0;
	bool flag = true;

	switch (nf) {
		case 'a':
			p1 = el.getP4();
			p2 = el.getP2();
			p3 = el.getP1();
			break;
		case 'b':
			p1 = el.getP4();
			p2 = el.getP3();
			p3 = el.getP2();
			break;
		case 'c':
			p1 = el.getP4();
			p2 = el.getP1();
			p3 = el.getP3();
			break;
		case 'd':
			p1 = el.getP1();
			p2 = el.getP2();
			p3 = el.getP3();
			break;
		default :
			flag = false;
			//switch
	}

	double * vec_n; //vx,vy,vz,vL [4]

	if (flag) {

		vec_n = siatka.zworcNormalnaDoFace(p1, p2, p3);

		if (Vec_norm != NULL) {

			Vec_norm[0] = vec_n[0];
			Vec_norm[1] = vec_n[1];
			Vec_norm[2] = vec_n[2];

		}

		if (Area != NULL) {
			*Area = 0.5 * vec_n[3];
		}

		delete []vec_n;

	}



}

double TicooMesh3D::mmr_el_hsize(int El) {

	double hsize = pow(siatka.V_objetoscT(El) * 6, 1.0 / 3.0);
	return hsize;

}


int TicooMesh3D::mmr_el_eq_neig(int El, int * Neig, int * Neig_sides) {

	Element & el = siatka.elements.getElement(El);
	int it = 3;

	if (Neig != NULL) {
		Neig[0] = el.getE1() + 1;
		Neig[1] = el.getE2() + 1;
		Neig[2] = el.getE3() + 1;
		Neig[3] = el.getE4() + 1;

		for (int i = 0; i < 4; ++i) {
			if (Neig[i] == 0) {
				--it;
			}
		}

	}

	if (Neig_sides != NULL) {
		Neig_sides[0] = el.getE1() + 1;
		Neig_sides[1] = el.getE2() + 1;
		Neig_sides[2] = el.getE3() + 1;
		Neig_sides[3] = el.getE4() + 1;
	}

	return it;
}


TicooMesh3D::TicooMesh3D(int wielkoscTab) : plansza(wielkoscTab, wielkoscTab, wielkoscTab, false),
	rozrost(&plansza),
	siatka(wielkoscTab, wielkoscTab, wielkoscTab),
	siatka1(wielkoscTab, wielkoscTab, wielkoscTab)
{


	srand((unsigned)time(&losowa));
	wielkoscTabX = wielkoscTab;
	wielkoscTabY = wielkoscTab;
	wielkoscTabZ = wielkoscTab;
	wielkoscTabMax = wielkoscTab;
	aktualneIDElementu = 0;

}


// Rozrost ////////////////////////////////////////////////
// Rozrost ////////////////////////////////////////////////
// Rozrost ////////////////////////////////////////////////

void TicooMesh3D::rozrostLosuj(int iloscZiaren) {

	rozrost.rozrostLosuj(iloscZiaren);

}

void TicooMesh3D::rozrostLosujRange(int iloscZiaren, int range) {

	rozrost.rozrostLosujRange(iloscZiaren, range);

}

void TicooMesh3D::rozrostZiarenGranicy(int wspPrawdop) {

	rozrost.rozrostZiarenGranicy(wspPrawdop);

}

void TicooMesh3D::pomniejszenieZiarnaGranicy(int wspPomniejszaniaZiaren) {

	rozrost.pomniejszenieZiarnaGranicy(wspPomniejszaniaZiaren);

}

void TicooMesh3D::rozrostZiarnaZOtoczeniem(int wspPrawdop, int wspPomniejszaniaZiaren, int ileRozrost, int ileOtocz, int ilePetla) {

	for (int i = 0; i < ilePetla; i++) {

		for (int j = 0; j < ileRozrost; j++) {
			rozrost.rozrostZiarenGranicy(wspPomniejszaniaZiaren);
		}

		for (int j = 0; j < ileOtocz; j++) {
			rozrost.pomniejszenieZiarnaGranicy(wspPrawdop);
		}

	}

}

void TicooMesh3D::komplekosoweTworzenieZiaren(int wspPrawdop, int ilePetla, int ileRozrost, int ileOtocz, int wspPomniejszaniaZiaren) {

	rozrost.komplekosoweTworzenieZiaren(wspPrawdop, ilePetla, ileRozrost, ileOtocz, wspPomniejszaniaZiaren);

}

void TicooMesh3D::granica() {

	if (plansza.getWarPer()) {
		rozrost.graniceP1Per();
	}
	else {
		rozrost.graniceP1NiePer();
	}

}

// Plansza ////////////////////////////////////////////////
// Plansza ////////////////////////////////////////////////
// Plansza ////////////////////////////////////////////////

void TicooMesh3D::wczytajAC(const char * nazwa) {

	plansza.wczytajAC(nazwa);

}

void TicooMesh3D::setPlansza(double dx, double dy, double dz, bool periodyczne, bool wlaczAC) {

	wielkoscTabX = dx;
	wielkoscTabY = dy;
	wielkoscTabZ = dz;

	if (dx > dy) {
		if (dx > dz) {
			wielkoscTabMax = dx;
		}
		else {
			wielkoscTabMax = dz;
		}
	}
	else {
		if (dy > dz) {
			wielkoscTabMax = dy;
		}
		else {
			wielkoscTabMax = dz;
		}
	}

	plansza.setPlansza(dx, dy, dz, periodyczne, wlaczAC);
	siatka.ustawWIelkoscObszaru(dx, dy, dz);
	siatka1.ustawWIelkoscObszaru(dx, dy, dz);

}

void TicooMesh3D::elementyDoZiarnaAC() {

	plansza.elementyDoZiaren(siatka.zwrocRefElementow());

}
void TicooMesh3D::elementyDoZiarnaPunktyGraniczne() {

	int ileZ = siatka.dopasujElementDoZiarna(true) + 10000;
	plansza.resetZiaren(ileZ);

}

// Losuj ////////////////////////////////////////////////
// Losuj ////////////////////////////////////////////////
// Losuj ////////////////////////////////////////////////

void TicooMesh3D::losowanieDodatkowychPunktow(bool pobierzPunkty, bool losuj, int ileP, int rozdz, int odG, int odP, bool flagaPer, double ax, double b0, int rSzuk) {

	double wx = plansza.getWielkoscX();
	double wy = plansza.getWielkoscY();
	double wz = plansza.getWielkoscZ();

	if (pobierzPunkty) {
		losowane.setPG(rozrost.pGranicy);
	}
	losowane.losuj(losuj, wx, wy, wz, ileP, rozdz, odG, odP, flagaPer, ax, b0, rSzuk);

}

void TicooMesh3D::ustawPunktySpaw1() {

	losowane.dodajPunktyPodSpawanie(plansza.getWielkoscX(), plansza.getWielkoscY(), plansza.getWielkoscZ());

}

void TicooMesh3D::zmianaKolejnosciPunktow(int ileRazy, int wielkoscKomurki, int ileNaKomurke) {

	for (int i = 0; i < ileRazy; ++i) {
		losowane.punktyPoczatkoweOctree(wielkoscKomurki, ileNaKomurke);
	}

}

// Siatka ////////////////////////////////////////////////
// Siatka ////////////////////////////////////////////////
// Siatka ////////////////////////////////////////////////



void TicooMesh3D::init_all_change(int a) {

	switch (a) {


		case 1:
			siatka.maxZgranica();
			break;
		case 2:
			siatka.maxZgranica1();
			break;



	}

}

void TicooMesh3D::delanouy(bool dokladneWyszukanie) {


	double wymX = plansza.getWielkoscX();
	double wymY = plansza.getWielkoscY();
	double wymZ = plansza.getWielkoscZ();

	siatka.pryzmy.czysc(10000, 10000);
	siatka.delanoyTablicaOP(wymX, wymY, wymZ, false, losowane.zwrocPunktyZGranica(), 0.3, 3, dokladneWyszukanie);
	siatka.zapiszNajVTrojkat();

}

void TicooMesh3D::reMes() {

	siatka.reMes();

}


void TicooMesh3D::wygladzanieReMES(int ileRazy, double waga, bool wyglZWag) {

	for (int j = 0; j < ileRazy; ++j) {


		siatka.sasiedniePunkty();
		//Memo1->Lines->Add(j);
		siatka.stosunekRdorCout();
		for (int i = 0; i < 20; ++i) {

			if (wyglZWag) {
				siatka.wygladzanieLaplaceWaga(2, waga, 1);
			}
			else {
				siatka.wygladzanieLaplace();
			}

		}
		siatka.reMes();
		//if(test.rozbijajDuzeTrojkaty()){--j;}

	}

	siatka.reMes();
	siatka.sasiedniePunkty();
	siatka.stosunekRdorCout();
}



double TicooMesh3D::objetoscT(FILE * out) {

	return siatka.V_objetosc(out);

}

void TicooMesh3D::sasiedniePunkty() {

	siatka.sasiedniePunkty();

}

void TicooMesh3D::sasiednieElementy() {

	siatka.sasiednieElementy();

}

void TicooMesh3D::wygladzanieLaplaceWaga(int ileR, double waga) {

	for (int i = 0; i < ileR; ++i) {
		siatka.wygladzanieLaplaceWaga(2, waga, 1);
	}

}

void TicooMesh3D::wygladzanieLaplace(int ileR) {

	for (int i = 0; i < ileR; ++i) {
		siatka.wygladzanieLaplace();
	}

}

void TicooMesh3D::optymalizacjaStosunekRdor(int ileR, double dokladnosc) {

	for (int i = 0; i < ileR; ++i) {
		siatka.optymalizacjaStosunekRdor(dokladnosc);
	}

}



void TicooMesh3D::poprawWszystkieElFAST(double docelowyStosunekRdor, bool narozneTez) {

	siatka.poprawaElementWszystkieRdorPoprawaMiejscowaPrzysp(5, docelowyStosunekRdor, narozneTez);

}

void TicooMesh3D::conwersjaDoTicooMesh3DZapisu() {

	siatka.uzupelnienieElementow();
	siatka.wyszukanieSasidnichElementowE();
	siatka.oznaczWezlyGraniczneNaPodstawieScianG();

}

//operacje na siatkach //skalowanie przesuwanie itp.

void TicooMesh3D::siatkaNaSrodek(double powiekszPrzestrzen) {

	double * wielkosc;
	wielkosc = siatka.ustawaElementyNaSrodku(powiekszPrzestrzen);

	plansza.setPlansza(int(wielkosc[0]) + 1, int(wielkosc[1]) + 1, int(wielkosc[2]) + 1, false, false);
	siatka.ustawWIelkoscObszaru(int(wielkosc[0]) + 1, int(wielkosc[1]) + 1, int(wielkosc[2]) + 1);

	delete []wielkosc;
}

// Zapis ////////////////////////////////////////////////
// Zapis ////////////////////////////////////////////////
// Zapis ////////////////////////////////////////////////

void TicooMesh3D::wczytajPunktyGraniczne(const char * nazwa) {

	siatka.pryzmy.czysc(10000, 10000);
	rozrost.wczytajPunktyGraniczne(nazwa);
	siatka.points.wczytajZPliku(nazwa, false);

}

void TicooMesh3D::zapisParaView(const char * nazwa) {

	siatka.zapisParaView(nazwa);

}

void TicooMesh3D::zapiszAbqus(const char * nazwa) {

	siatka.ZapiszDoPlikuAbaqus(nazwa);

}


void TicooMesh3D::ZapiszDoPlikuNasZPrzesunieciem_BC_Face(const char * nazwa) {

	siatka.ZapiszDoPlikuNasZPrzesunieciem_BC_Face(nazwa, 1, 0, 0, 0, 1);

}

void TicooMesh3D::zapiszDoPlikuNas(const char * nazwa, double podzielZ, bool brick, double tX, double tY, double tZ, double podzielAll) {

	siatka.ZapiszDoPlikuNasZPrzesunieciem(nazwa, podzielZ, brick, tX, tY, tZ, podzielAll, NULL);

}

void TicooMesh3D::zapiszPunkty(const char * nazwa) {


	int liczba = (int)siatka.points.getElement(7).getX();
	string tmp; // brzydkie rozwi¹zanie
	sprintf((char *)tmp.c_str(), "%d", liczba);
	string str = nazwa;
	str += tmp.c_str();

	if (siatka.points.getIter() > 1) {
		siatka.points.zapisDoPliku((char *)str.c_str());
	}

}


void TicooMesh3D::wczytajPunktyPSS(const char * nazwa) {

	siatka.wczytajPunktyPSS(nazwa, 1);

}


void TicooMesh3D::wczytajPunktyPSSnast(const char * nazwa) {

	siatka.wczytajPunktyPSSnast(nazwa, 1);

}

void TicooMesh3D::zapiszPSS1DoPlikuNasWOparciuOPowTetra(const char * nazwa) {

	siatka.ZapiszDoPlikuNasWOparciuOPowTetra(nazwa, 1, true, true);

}

void TicooMesh3D::zapiszPSS1DoPlikuNasWOparciuOPowHybrydPSS(const char * nazwa) {

	siatka.ZapiszDoPlikuNasWOparciuOPowHybrydPSS(nazwa, 1, true, true);

}


void TicooMesh3D::wczytajAbaqus(const char * nazwa, int mnoznik, double px, double py, double pz) {

	siatka.pryzmy.czysc(10000, 10000);
	siatka.wczytajAnsys1(nazwa, mnoznik, px, py, pz);

	//siatka.elements.ustaw_ukval_zero();
	siatka.uzupelnienieElementow();
	siatka.wyszukanieSasidnichElementowE();

}

void TicooMesh3D::wczytajNAS(const char * nazwa, int mnoznik, double px, double py, double pz) {

	siatka.pryzmy.czysc(10000, 10000);
	siatka1.pryzmy.czysc(10000, 10000);
	siatka.wczytajPlikNas(nazwa, mnoznik, true, px, py, pz);

	//siatka.elements.ustaw_ukval_zero();
	//siatka.uzupelnienieElementow();
	//siatka.wyszukanieSasidnichElementowE();
	siatka.kopiujSiatke(siatka1);

}

double TicooMesh3D::ruchSpawanie(double obecnyKrok, double krok_start, double minPoprawy, double doX, double doY, double dl, double zmiejszaPrzes, double limit, double szerokosc) {

	return siatka1.ruchSpawanie(obecnyKrok, krok_start, minPoprawy, doX, doY, dl, zmiejszaPrzes, limit, szerokosc);


}

void TicooMesh3D::ruchPrzestrzen(int ileWarstw, int obecny_krok, int od_krok, int ileKrok, double minPoprawy, double px0, double py0, double pz0, double px1, double py1, double pz1, double endX, double endY, double endZ) {

	//siatka.kopiujSiatke(siatka1);
	siatka1.ruchPrzestrzen(ileWarstw, obecny_krok, od_krok, ileKrok, minPoprawy, px0, py0, pz0, px1, py1, pz1, endX, endY, endZ);
	//
}

void TicooMesh3D::createCube(const char * nazwa, int node_x, int node_y, int node_z, double size_x, double size_y, double size_z, double divide, int * warunki) {


	//int warunki[6]={1,1,1,1,1,1}; //0...5 z0,z1,x0,x1,y0,y1

	plansza.setPlansza(size_x, size_y, size_z, false, false);
	siatka1.ustawWIelkoscObszaru(size_x, size_y, size_z);
	siatka1.siatkaStatyczna(size_x, size_y, size_z, node_x, node_y, node_z);

	siatka1.uzupelnienieElementow();
	siatka1.wyszukanieSasidnichElementowE();
	siatka1.oznaczWezlyGraniczneNaPodstawieScianG();

	siatka1.ZapiszDoPlikuNasZPrzesunieciem(nazwa, 1, true, 0, 0, 0, divide, warunki);

}
