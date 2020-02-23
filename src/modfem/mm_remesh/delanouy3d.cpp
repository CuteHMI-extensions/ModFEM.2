
#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdlib>


#include <modfem/mm_remesh/delanouy3d.h>
#include <modfem/mm_remesh/texty.h>



Delanouy3D::~Delanouy3D() {

	delete []laplas;
	delete []mapFace;

}


void Delanouy3D::count_bc_surface(FILE * Interactive_output) {

	IntList warBC(1000, 1000);
	double * warBC_pow, zew = 0, wew = 0, pow = 0;

	for (int i = 0, ileF = faces.getIter(); i < ileF; i++) {

		if (faces.getElement(i).getBC() >= 1) {
			warBC.setElementUniqat(faces.getElement(i).getBC());
		}

	}

	warBC.sortowanie();
	warBC_pow = new double[warBC.getIter()];

	for (int i = 0, ileF = warBC.getIter(); i < ileF; i++) {

		warBC_pow[i] = 0;

	}

	int pon[3];
	for (int i = 0, id, bc_id, ileF = faces.getIter(); i < ileF; i++) {

		bc_id = faces.getElement(i).getBC();
		id = warBC.uniqatPozycja(bc_id);
		if (id != -1) {

			getFacePoints(i, pon);
			pow = powSciany(pon[0], pon[1], pon[2]);
			warBC_pow[id] += pow;



			if (bc_id == 2) {
				wew += pow;
			}
			if ((bc_id != 2) && (bc_id >= 1) && (bc_id <= 25)) {
				zew += pow;
			}
		}


	}




	for (int i = 0, ileF = warBC.getIter(); i < ileF; i++) {

		fprintf(Interactive_output, "Surface of bc_num_ %d value = %lf \n", warBC.getElement(i), warBC_pow[i]);

	}

	fprintf(Interactive_output, "Surface of external_inflow : %lf \n", wew);
	fprintf(Interactive_output, "Surface of internal_inflow : %lf \n", zew);


	delete []warBC_pow;

}

void Delanouy3D::zmienWarunkeSciana(double r_max, double r_min, int zWar, int naWar) {

	int pon[3];

	double x1, x2, x3, y1, y2, y3, r;

	for (int i = 0, ileF = faces.getIter(); i < ileF; ++i) {
		if (faces.getElement(i).getBC() == zWar) {


			getFacePoints(i, pon);

			x1 = points.getElement(pon[0]).getX();
			y1 = points.getElement(pon[0]).getY();

			x2 = points.getElement(pon[1]).getX();
			y2 = points.getElement(pon[1]).getY();

			x3 = points.getElement(pon[2]).getX();
			y3 = points.getElement(pon[2]).getY();

			x1 = (x1 + x2 + x3) / 3;
			y1 = (y1 + y2 + y3) / 3;
			r = sqrt(x1 * x1 + y1 * y1);

			if (r_min <= r && r <= r_max) {

				faces.getElement(i).setBC(naWar);

			}

		}
	}

}




void Delanouy3D::stosunekRdor(FILE * out) {



	fprintf(out, " Volume :  %lf \n", V_objetosc(out));

//bool flaga=false;
	double V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;

	double Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, nx, ny, nz;
//double maxX = dx;
	int z5 = 0, z6 = 0, z7 = 0, z8 = 0, z9 = 0, z10 = 0, z11 = 0, z12 = 0, z13 = 0, z15 = 0, z17 = 0, z19 = 0, z21 = 0, z23 = 0, zZ = 0;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {


		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		//R=sqrt(elements.getElement(i).getR());

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);

		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );

		pol = 0;
		Ax = x2 - x1;
		Bx = x3 - x1;
		Cx = x4 - x1;
		Ay = y2 - y1;
		By = y3 - y1;
		Cy = y4 - y1;
		Az = z2 - z1;
		Bz = z3 - z1;
		Cz = z4 - z1;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = Ay * Cz - Cy * Az;
		ny = Cx * Az - Ax * Cz;
		nz = Ax * Cy - Cx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = By * Cz - Cy * Bz;
		ny = Cx * Bz - Bx * Cz;
		nz = Bx * Cy - Cx * By;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		V = (nx * Ax + ny * Ay + nz * Az) * 0.5;

		Ax = x2 - x4;
		Bx = x3 - x4;
		Ay = y2 - y4;
		By = y3 - y4;
		Az = z2 - z4;
		Bz = z3 - z4;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		r = V / pol;
		R = R / r;

		if (R < 5) {
			++z5;
		}
		else if (R < 6) {
			++z6;
		}
		else if (R < 7) {
			++z7;
		}
		else if (R < 8) {
			++z8;
		}
		else if (R < 9) {
			++z9;
		}
		else if (R < 10) {
			++z10;
		}
		else if (R < 11) {
			++z11;
		}
		else if (R < 12) {
			++z12;
		}
		else if (R < 13) {
			++z13;
		}
		else if (R < 15) {
			++z15;
		}
		else if (R < 17) {
			++z17;
		}
		else if (R < 19) {
			++z19;
		}
		else if (R < 21) {
			++z21;
		}
		else if (R < 23) {
			++z23;
		}
		else {
			++zZ;

		}
	}



	fprintf(out, " R<5 :  %d \n", z5);
	fprintf(out, " R<6 :  %d \n", z6);
	fprintf(out, " R<7 :  %d \n", z7);
	fprintf(out, " R<8 :  %d \n", z8);
	fprintf(out, " R<9 :  %d \n", z9);
	fprintf(out, " R<10 :  %d \n", z10);
	fprintf(out, " R<11:  %d \n", z11);
	fprintf(out, " R<12 :  %d \n", z12);
	fprintf(out, " R<13 :  %d \n", z13);
	fprintf(out, " R<15 :  %d \n", z15);
	fprintf(out, " R<17 :  %d \n", z17);
	fprintf(out, " R<19 :  %d \n", z19);
	fprintf(out, " R<21 :  %d \n", z21);
	fprintf(out, " R<23 :  %d \n", z23);
	fprintf(out, " R>23 :  %d \n", zZ);




}

void Delanouy3D::smoothing(double l_step, int repetition) {

//uzupelnienieElementow();
//wyszukanieSasidnichElementowE();
	oznaczWezlyGraniczneNaPodstawieScianG();
	sasiednieElementy();

	for (int i = 0; i < repetition; ++i) {
		optymalizacjaStosunekRdorCalosc_zmienna_dl_kroku(l_step, true);
	}

}



int Delanouy3D::getlowestGroupID() {

	int lowestID = elements.getElement(0).getRodzajZiarna();

	for (int i = 1; i < elements.getIter(); ++i) {
		if (lowestID > elements.getElement(i).getRodzajZiarna()) {
			lowestID = elements.getElement(i).getRodzajZiarna();
		}
	}

	return lowestID;

}

void Delanouy3D::zapiszWarunkiBrzegowe(int * warTab) {
	//warTab ilosc elementow*4
	int el, bc, nrFCinEL;
	char nf;

	for (int i = 0; i < elements.getIter() * 4; ++i) {
		warTab[i] = 0;
	}

	for (int i = 0, ileF = faces.getIter(); i < ileF; ++i) {

		bc = faces.getElement(i).getBC();

		if (bc != 0) {

			el = faces.getElement(i).getInEl();
			nf = faces.getElement(i).getFaceInEl();

			switch (nf) {

				case 'a':
					nrFCinEL = 0;
					break;
				case 'b':
					nrFCinEL = 1;
					break;
				case 'c':
					nrFCinEL = 2;
					break;
				case 'd':
					nrFCinEL = 3;
					break;

			}

			warTab[el * 4 + nrFCinEL] = bc;

		}

	}

}

void Delanouy3D::wczytajWarunkiBrzegowe(int * warTab) {

	//warTab ilosc elementow*4
	int el, bc, nrFCinEL;
	char nf;



	for (int i = 0, ileF = faces.getIter(); i < ileF; ++i) {


		el = faces.getElement(i).getInEl();
		nf = faces.getElement(i).getFaceInEl();

		switch (nf) {

			case 'a':
				nrFCinEL = 0;
				break;
			case 'b':
				nrFCinEL = 1;
				break;
			case 'c':
				nrFCinEL = 2;
				break;
			case 'd':
				nrFCinEL = 3;
				break;

		}

		bc = warTab[el * 4 + nrFCinEL];

		if (bc != 0) {
			faces.getElement(i).setBC(bc);
		}

	}


}

void Delanouy3D::creatFaceAndEdgeRememberBC() {


	int * warTab;
	warTab = new int[elements.getIter() * 4];

	sasiednieElementy();
	zapiszWarunkiBrzegowe(warTab);
	creatFaceAndEdge();
	wczytajWarunkiBrzegowe(warTab);

	delete []warTab;

}

void Delanouy3D::zapiszWarunkiPoczatkoweNaPodstawieDanych(const char * nazwa, double * tempB, int dl_tempB, int * tabBlock, int dl_block) {


	DoubleList tempBlock(100, 100);
	for (int i = 0, num; i < dl_tempB; i += 2) {

		num = tempB[i];

		for (int j = 0; j < dl_block; j += 2) {

			if (num == int(tabBlock[j + 1])) {
				tempBlock.setElement(tabBlock[j]);
				tempBlock.setElement(tempB[i + 1]);
			}

		}

	}

	ofstream zapis(nazwa);

	int ileP = points.getIter();

	zapis << 1 << " " << 3 << endl << ileP << endl;
	float * temp, tmp;
	temp = new float[ileP];


	for (int i = 0; i < ileP; ++i) {
		temp[i] = 0.0;
	}

	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		for (int ii = 0, z, ileD = tempBlock.getIter(); ii < ileD; ii += 2) {
			z = elements.getElement(i).getRodzajZiarna();

			if (int(tempBlock.getElement(ii)) == z) {
				tmp = tempBlock.getElement(ii + 1);
				temp[p1] = tmp;
				temp[p2] = tmp;
				temp[p3] = tmp;
				temp[p4] = tmp;

				break;
			}

		}
	}

	for (int i = 0; i < ileP; i++) {

		zapis << i + 1 << " " << temp[i] << endl << temp[i] << endl << temp[i] << endl;

	}


	delete []temp;
	zapis.close();

}



int Delanouy3D::materialNumber() {

	IntList listaMat(100, 100);
	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		listaMat.setElementUniqat(elements.getElement(i).getRodzajZiarna());

	}

	return listaMat.getIter();

}

int Delanouy3D::materialIDs(int * tab, int l_tab) {

	IntList listaMat(100, 100);
	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		listaMat.setElementUniqat(elements.getElement(i).getRodzajZiarna());

	}

	for (int i = 0; i < l_tab; ++i) {

		tab[i] = listaMat.getElement(i);

	}

	return 0;
}

void Delanouy3D::warunkiNaPodstawieDanych(int * matMap_in, int dl_mat, int * blockMap_in, int dl_block, int * war, int dl_war) {

	IntList matMap(dl_mat + 10, dl_mat + 10);
	IntList blockMap(dl_block + 10, dl_block + 10);

	for (int i = 0; i < dl_mat; ++i) {
		matMap.setElement(matMap_in[i]);
	}
	for (int i = 0; i < dl_block; ++i) {
		blockMap.setElement(blockMap_in[i]);
	}

	IntList war_grup(100, 100);
	IntList war_mat(100, 100);
	IntList war_block(100, 100);

	for (int i = 0; i < dl_war; i += 4) {

		switch (war[i]) {

			case 0:
				war_grup.setElement(war[i + 1]);
				war_grup.setElement(war[i + 2]);
				war_grup.setElement(war[i + 3]);
				break;

			case 1:
				war_block.setElement(war[i + 1]);
				war_block.setElement(war[i + 2]);
				war_block.setElement(war[i + 3]);
				break;

			case 2:
				war_mat.setElement(war[i + 1]);
				war_mat.setElement(war[i + 2]);
				war_mat.setElement(war[i + 3]);
				break;

		}

	}

	int grupID1, grupID2, blockID1, blockID2, matID1, matID2, nr_war;
	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		for (int j = 0, e; j < 4; ++j) {

			e = elements.getElement(i).getE(j);

			if (e < 0) { //if war

				grupID1 = elements.getElement(i).getRodzajZiarna();
				blockID1 = blockMap.getElementKey(grupID1);
				matID1 = matMap.getElementKey(grupID1);

				if (e != -1) { // if contact
					e = -e - 2;
					grupID2 = elements.getElement(e).getRodzajZiarna();
					blockID2 = blockMap.getElementKey(grupID2);
					matID2 = matMap.getElementKey(grupID2);


					//m1->Lines->Add(blockID2);

					if (-10 != (nr_war = war_grup.getElement2Kay(grupID1, grupID2))) {
						faces.getElement(mapFace[i].getElement(j)).setBC(nr_war);
					}
					else if (-10 != (nr_war = war_block.getElement2Kay(blockID1, blockID2))) {
						faces.getElement(mapFace[i].getElement(j)).setBC(nr_war);
					}
					else if (-10 != (nr_war = war_mat.getElement2Kay(matID1, matID2))) {
						faces.getElement(mapFace[i].getElement(j)).setBC(nr_war);
					}
					//else{m1->Lines->Add("zle contact");}


				}
				else { //if normal boundary
					if (faces.getElement(mapFace[i].getElement(j)).getBC() == 0) {
						if (-10 != (nr_war = war_grup.getElement2Kay(grupID1, -1))) {
							faces.getElement(mapFace[i].getElement(j)).setBC(nr_war);
						}
						else if (-10 != (nr_war = war_block.getElement2Kay(blockID1, -1))) {
							faces.getElement(mapFace[i].getElement(j)).setBC(nr_war);
						}
						else if (-10 != (nr_war = war_mat.getElement2Kay(matID1, -1))) {
							faces.getElement(mapFace[i].getElement(j)).setBC(nr_war);
						}
						//else{m1->Lines->Add("zle normal boundary");}
					}

				}

			}

		}
	}


}

void Delanouy3D::rozdzielSasiednieMaterialy(bool block) {

	if (block) {



		for (int i = 0, rodzaj0, ileE = elements.getIter(); i < ileE; ++i) {
			rodzaj0 = elements.getElement(i).getTemp();


			for (int j = 0, e, rodzaj1; j < 4; ++j) {
				e = elements.getElement(i).getE(j);
				if (e >= 0) {
					rodzaj1 = elements.getElement(e).getTemp();

					if (rodzaj0 != rodzaj1) {

						elements.getElement(i).setE(-e - 2, j);
					}
				}

			}
		}



	}
	else {

		for (int i = 0, rodzaj0, ileE = elements.getIter(); i < ileE; ++i) {
			rodzaj0 = elements.getElement(i).getRodzajZiarna();


			for (int j = 0, e, rodzaj1; j < 4; ++j) {
				e = elements.getElement(i).getE(j);
				if (e >= 0) {
					rodzaj1 = elements.getElement(e).getRodzajZiarna();

					if (rodzaj0 != rodzaj1) {

						elements.getElement(i).setE(-e - 2, j);
					}
				}

			}
		}

	}


}

void Delanouy3D::podzielNaObjekty(int * tabBlock, int dl_tab) {



	if (dl_tab == 0) {
		elements.ustawTempnaZiarna();
	}
	else {

		elements.ustawTempZtabB(tabBlock, dl_tab);

	}

	sasiednieElementy();
	rozdzielSasiednieMaterialy(true);

	int ilePunktow = points.getIter(), punk[3];

	IntList listaMaterialow(100, 100);
	int * rodzajZiarna = new int[elements.getIter()];
	//int *pomocPunkty = new int[points.getIter()];
	//for(int i=0,ileP=points.getIter();i<ileP;++i){pomocPunkty[0] = -1;}

	for (int i = 0, z, ileE = elements.getIter(); i < ileE; ++i) {
		z = elements.getElement(i).getTemp();
		rodzajZiarna[i] = z;
		listaMaterialow.setElementUniqat(z);
	}



	listaMaterialow.sortowanie();
	for (int i = 0, ile = listaMaterialow.getIter(); i < ile; ++i) {
		if (fluidBlock.sprawdzCzyJest(listaMaterialow.getElement(i))) {

			listaMaterialow.zamienMiejscami(i, --ile);

		}
	}


	for (int i = 0, ile = listaMaterialow.getIter(); i < ile; ++i) {
		if (uppermostBlock.sprawdzCzyJest(listaMaterialow.getElement(i))) {

			listaMaterialow.zamienMiejscami(i, 0);

		}
	}

	for (int o = 0, rodzaj0, ileZ = listaMaterialow.getIter(); o < ileZ; ++o) {
		rodzaj0 = listaMaterialow.getElement(o);

		for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
			if (rodzaj0 == rodzajZiarna[i]) {
				//szukamy sciany do rozdzielenia
				for (int j = 0, e, rodzaj1; j < 4; ++j) {
					e = elements.getElement(i).getE(j);
					if (e < -1 && rodzajZiarna[-e - 2] > rodzaj0) {
						elements.getElement(i).getNumeryElSciany(j, punk);
						for (int ij = 0; ij < 3; ++ij) {
							if (punk[ij] < ilePunktow) {

								points.setElement(points.getElement(punk[ij]));
								for (int l = 0, nrLE, ileL = laplas[punk[ij]].getIter(); l < ileL; ++l) {
									nrLE = laplas[punk[ij]].getElement(l);
									if (rodzaj0 == rodzajZiarna[nrLE]) {

										elements.getElement(nrLE).podmienPunkt(punk[ij], points.getIter() - 1);
									}

								}


							}
						}
					}

				}

			}

		}

	}

	delete []rodzajZiarna;

	//wrocZSasiadami();

	//delete []pomocPunkty;
}



void Delanouy3D::getFacePoints(int face_id, int * points) {

	int el = faces.getElement(face_id).getInEl();
	char nf = faces.getElement(face_id).getFaceInEl();

	switch (nf) {
		case 'a':
			points[0] = elements.getElement(el).getP4();
			points[1] = elements.getElement(el).getP2();
			points[2] = elements.getElement(el).getP1();
			break;

		case 'b':
			points[0] = elements.getElement(el).getP4();
			points[1] = elements.getElement(el).getP3();
			points[2] = elements.getElement(el).getP2();
			break;

		case 'c':
			points[0] = elements.getElement(el).getP4();
			points[1] = elements.getElement(el).getP1();
			points[2] = elements.getElement(el).getP3();
			break;

		case 'd':
			points[0] = elements.getElement(el).getP1();
			points[1] = elements.getElement(el).getP2();
			points[2] = elements.getElement(el).getP3();
			break;
			//switch
	}

}

void Delanouy3D::zmien_do_kontaktu(int pozycja) {


	switch (do_kontaktu[pozycja]) {

		case 0:
			do_kontaktu[pozycja + 1] = 0;
			do_kontaktu[pozycja + 2] = 1;
			do_kontaktu[pozycja + 3] = 2;
			do_kontaktu[pozycja] = 1;
			break;
		case 1:
			do_kontaktu[pozycja + 1] = 0;
			do_kontaktu[pozycja + 2] = 2;
			do_kontaktu[pozycja + 3] = 1;
			do_kontaktu[pozycja] = 2;
			break;
		case 2:
			do_kontaktu[pozycja + 1] = 1;
			do_kontaktu[pozycja + 2] = 0;
			do_kontaktu[pozycja + 3] = 2;
			do_kontaktu[pozycja] = 3;
			break;
		case 3:
			do_kontaktu[pozycja + 1] = 1;
			do_kontaktu[pozycja + 2] = 2;
			do_kontaktu[pozycja + 3] = 0;
			do_kontaktu[pozycja] = 4;
			break;
		case 4:
			do_kontaktu[pozycja + 1] = 2;
			do_kontaktu[pozycja + 2] = 1;
			do_kontaktu[pozycja + 3] = 0;
			do_kontaktu[pozycja] = 5;
			break;
		case 5:
			do_kontaktu[pozycja + 1] = 2;
			do_kontaktu[pozycja + 2] = 0;
			do_kontaktu[pozycja + 3] = 1;
			do_kontaktu[pozycja] = 6;
			break;
		case 6:
			std::cout << " ZLEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE do_kontaktu" << endl;
			break;

	}


}
void Delanouy3D::info_bc_connect() {

	for (int i = 0; i < 16; ++i) {
		std::cout << endl;
		for (int ii = 0, ij = 0; ii < 3; ++ii, ij += 4) {

			std::cout << "case " << i * 10 + ij << ":" << "el_id[0] = " << do_kontaktu[i * 12 + ii * 4 + 1] << ";" << "el_id[1] = " << do_kontaktu[i * 12 + ii * 4 + 2] << ";" << "el_id[2] = " << do_kontaktu[i * 12 + ii * 4 + 3] << ";break;" << endl;



			//std::cout<<do_kontaktu[i*12+ii*4]<<"   "<<do_kontaktu[i*12+ii*4+1]<<"   "<<do_kontaktu[i*12+ii*4+2]<<"   "<<do_kontaktu[i*12+ii*4+3]<<endl;


		}
		std::cout << endl;
	}

}

void Delanouy3D::mmr_dopasuj_punkty(int face_id, int fa_connect, int * el_id, int zmien) {



	int f1tab[3], f2tab[3], e1, e2;

	char f1 = faces.getElement(face_id).getFaceInEl(), f2 = faces.getElement(fa_connect).getFaceInEl();
	e1 = faces.getElement(face_id).getInEl();
	e2 = faces.getElement(fa_connect).getInEl();
	elements.getElement(e1).getNumeryElSciany(f1, f1tab);
	elements.getElement(e2).getNumeryElSciany(f2, f2tab);
	el_id[3] = e1 + 1;
	el_id[4] = e2 + 1;


	el_id[0] = 0;
	el_id[1] = 1;
	el_id[2] = 2;

	int pozycja = ((int)f1 - 97) + 4 * ((int)f2 - 97), miejsce, sup_pozycja;

	if (points.porownajPunkty(f1tab[0], f2tab[0])) {
		miejsce = 0;
	}
	else if (points.porownajPunkty(f1tab[0], f2tab[1])) {
		miejsce = 4;
	}
	else {
		miejsce = 8;
	}

	sup_pozycja = pozycja * 10 + miejsce;

	/*
	if(zmien){zmien_do_kontaktu(12*pozycja+miejsce);}

	el_id[0] = do_kontaktu[12*pozycja+miejsce+1];
	el_id[1] = do_kontaktu[12*pozycja+miejsce+2];
	el_id[2] = do_kontaktu[12*pozycja+miejsce+3];
	*/







	switch (sup_pozycja) {

		case 0:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;
		case 4:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;
		case 8:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;


		case 10:
			el_id[0] = 1;
			el_id[1] = 2;
			el_id[2] = 0;
			break;
		case 14:
			el_id[0] = 0;
			el_id[1] = 1;
			el_id[2] = 2;
			break;
		case 18:
			el_id[0] = 2;
			el_id[1] = 0;
			el_id[2] = 1;
			break;


		case 20:
			el_id[0] = 0;
			el_id[1] = 1;
			el_id[2] = 2;
			break;
		case 24:
			el_id[0] = 2;
			el_id[1] = 0;
			el_id[2] = 1;
			break;
		case 28:
			el_id[0] = 1;
			el_id[1] = 2;
			el_id[2] = 0;
			break;


		case 30:
			el_id[0] = 2;
			el_id[1] = 0;
			el_id[2] = 1;
			break;
		case 34:
			el_id[0] = 1;
			el_id[1] = 2;
			el_id[2] = 0;
			break;
		case 38:
			el_id[0] = 0;
			el_id[1] = 1;
			el_id[2] = 2;
			break;


		case 40:
			el_id[0] = 2;
			el_id[1] = 0;
			el_id[2] = 1;
			break;
		case 44:
			el_id[0] = 0;
			el_id[1] = 1;
			el_id[2] = 2;
			break;
		case 48:
			el_id[0] = 1;
			el_id[1] = 2;
			el_id[2] = 0;
			break;


		case 50:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;
		case 54:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;
		case 58:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;


		case 60:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;
		case 64:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;
		case 68:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;


		case 70:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;
		case 74:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;
		case 78:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;


		case 80:
			el_id[0] = 0;
			el_id[1] = 1;
			el_id[2] = 2;
			break;
		case 84:
			el_id[0] = 1;
			el_id[1] = 2;
			el_id[2] = 0;
			break;
		case 88:
			el_id[0] = 2;
			el_id[1] = 0;
			el_id[2] = 1;
			break;


		case 90:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;
		case 94:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;
		case 98:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;


		case 100:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;
		case 104:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;
		case 108:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;


		case 110:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;
		case 114:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;
		case 118:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;


		case 120:
			el_id[0] = 1;
			el_id[1] = 2;
			el_id[2] = 0;
			break;
		case 124:
			el_id[0] = 2;
			el_id[1] = 0;
			el_id[2] = 1;
			break;
		case 128:
			el_id[0] = 0;
			el_id[1] = 1;
			el_id[2] = 2;
			break;


		case 130:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;
		case 134:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;
		case 138:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;


		case 140:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;
		case 144:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;
		case 148:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;


		case 150:
			el_id[0] = 0;
			el_id[1] = 2;
			el_id[2] = 1;
			break;
		case 154:
			el_id[0] = 1;
			el_id[1] = 0;
			el_id[2] = 2;
			break;
		case 158:
			el_id[0] = 2;
			el_id[1] = 1;
			el_id[2] = 0;
			break;


	}



	/*
	switch(pozycja){

		case 0:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 1:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 2:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 3:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 4:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 5:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 6:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 7:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 8:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 9:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 10:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 11:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 12:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 13:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 14:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;

		case 15:
		if(points.porownajPunkty(f1tab[0],f2tab[0])){
			el_id[0] = do_kontaktu[12*pozycja+1];
			if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;el_id[5]=0;}
			else{el_id[1] = do_kontaktu[12*pozycja+2];el_id[2] = do_kontaktu[12*pozycja+3];el_id[5]=1;}
		}
		else if(points.porownajPunkty(f1tab[0],f2tab[1])){
			el_id[0] = do_kontaktu[12*pozycja+5];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = do_kontaktu[12*pozycja+6];el_id[2] = do_kontaktu[12*pozycja+7];el_id[5]=2;}
			else{el_id[1] = 0;el_id[2] = 0;el_id[5]=3;}
		}
		else{
			el_id[0] = do_kontaktu[12*pozycja+9];
			if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;el_id[5]=4;}
			else{el_id[1] = do_kontaktu[12*pozycja+10];el_id[2] = do_kontaktu[12*pozycja+11];el_id[5]=5;}
		}
		break;


	}
	*/





	/*
	int punkty1[3],punkty2[3];

	getFacePoints(face_id,punkty1);
	getFacePoints(fa_connect,punkty2);

	std::cout<<"wspolzedne f0: "<<points.getElement(punkty1[0]).getX()<<"  "<<points.getElement(punkty1[0]).getY()<<"  "<<points.getElement(punkty1[0]).getZ()<<" - "<<points.getElement(punkty1[1]).getX()<<"  "<<points.getElement(punkty1[1]).getY()<<"  "<<points.getElement(punkty1[1]).getZ()<<" - "<<points.getElement(punkty1[2]).getX()<<"  "<<points.getElement(punkty1[2]).getY()<<"  "<<points.getElement(punkty1[2]).getZ()<<endl;
	std::cout<<"wspolzedne f1: "<<points.getElement(punkty2[0]).getX()<<"  "<<points.getElement(punkty2[0]).getY()<<"  "<<points.getElement(punkty2[0]).getZ()<<" - "<<points.getElement(punkty2[1]).getX()<<"  "<<points.getElement(punkty2[1]).getY()<<"  "<<points.getElement(punkty2[1]).getZ()<<" - "<<points.getElement(punkty2[2]).getX()<<"  "<<points.getElement(punkty2[2]).getY()<<"  "<<points.getElement(punkty2[2]).getZ()<<endl;
	*/
	/*
	string a;

	switch(f2){






	case 'a':
	//02 11 22
	a+=" a ";
	if(points.porownajPunkty(f1tab[0],f2tab[0])){
		a+=" 0 ";
		el_id[0] = 2;
		if(points.porownajPunkty(f1tab[1],f2tab[1])){el_id[1] = 1;el_id[2] = 0;a+=" 1 ";}
		else{el_id[1] = 0;el_id[2] = 1;a+=" 2 ";}
	}
	else if(points.porownajPunkty(f1tab[0],f2tab[1])){
		a+=" 1 ";
		el_id[0] = 1;
		if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 0;a+=" 1 ";}
		else{el_id[1] = 0;el_id[2] = 2;a+=" 2 ";}
	}
	else{
		a+=" 2 ";
		el_id[0] = 0;
		if(points.porownajPunkty(f1tab[1],f2tab[0])){el_id[1] = 2;el_id[2] = 1;a+=" 1 ";}
		else{el_id[1] = 1;el_id[2] = 2;a+=" 2 ";}

	}

	//std:cout<<a<<endl;




	case 'a':

	if(points.porownajPunkty(f1tab[0],f2tab[0])){el_id[0] = 2;el_id[1] = 0;el_id[2] = 1;}
	else if(points.porownajPunkty(f1tab[0],f2tab[1])){el_id[0] = 1;el_id[1] = 2;el_id[2] = 0;}
	else{el_id[0] = 0;el_id[1] = 1;el_id[2] = 2;}



	break;
	case 'b':

	if(points.porownajPunkty(f1tab[0],f2tab[0])){el_id[0] = 1;el_id[1] = 0;el_id[2] = 2;}
	else if(points.porownajPunkty(f1tab[0],f2tab[1])){el_id[0] = 2;el_id[1] = 1;el_id[2] = 0;}
	else{el_id[0] = 0;el_id[1] = 2;el_id[2] = 1;}

	break;
	case 'c':

	if(points.porownajPunkty(f1tab[0],f2tab[0])){el_id[0] = 2;el_id[1] = 1;el_id[2] = 0;}
	else if(points.porownajPunkty(f1tab[0],f2tab[1])){el_id[0] = 0;el_id[1] = 2;el_id[2] = 1;}
	else{el_id[0] = 1;el_id[1] = 0;el_id[2] = 2;}

	break;
	case 'd':

	if(points.porownajPunkty(f1tab[0],f2tab[0])){el_id[0] = 0;el_id[1] = 2;el_id[2] = 1;}
	else if(points.porownajPunkty(f1tab[0],f2tab[1])){el_id[0] = 1;el_id[1] = 0;el_id[2] = 2;}
	else{el_id[0] = 2;el_id[1] = 1;el_id[2] = 0;}
	break;

	}
	*/









}



void Delanouy3D::tablicaFC_connect() {

	delete []facebc_connect;
	int ileFB = faces.getIter();
	facebc_connect = new int[ileFB + 1];
	for (int i = 0; i < ileFB; ++i) {
		facebc_connect[i] = -1;
	}
	flag_bc_connect = true;

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
		for (int ii = 0, e, f1, f2; ii < 4; ++ii) {
			e = elements.getElement(i).getE(ii);
			if (e < -1) {


				f1 = mapFace[i].getElement(ii);
				e = -e - 2;

				for (int iii = 0, ee; iii < 4; ++iii) {
					ee = elements.getElement(e).getE(iii);
					ee = -ee - 2;
					if (ee == i) {
						f2 = mapFace[e].getElement(iii);
						facebc_connect[f1] = f2;
						//facebc_connect[f2]=f1;
					}
				}

				/*
				int punkty1[3],punkty2[3];

				getFacePoints(f1,punkty1);
				getFacePoints(f2,punkty2);

				std::cout<<"wspolzedne f0: "<<points.getElement(punkty1[0]).getX()<<"  "<<points.getElement(punkty1[0]).getY()<<"  "<<points.getElement(punkty1[0]).getZ()<<" - "<<points.getElement(punkty1[1]).getX()<<"  "<<points.getElement(punkty1[1]).getY()<<"  "<<points.getElement(punkty1[1]).getZ()<<" - "<<points.getElement(punkty1[2]).getX()<<"  "<<points.getElement(punkty1[2]).getY()<<"  "<<points.getElement(punkty1[2]).getZ()<<endl;
				std::cout<<"wspolzedne f1: "<<points.getElement(punkty2[0]).getX()<<"  "<<points.getElement(punkty2[0]).getY()<<"  "<<points.getElement(punkty2[0]).getZ()<<" - "<<points.getElement(punkty2[1]).getX()<<"  "<<points.getElement(punkty2[1]).getY()<<"  "<<points.getElement(punkty2[1]).getZ()<<" - "<<points.getElement(punkty2[2]).getX()<<"  "<<points.getElement(punkty2[2]).getY()<<"  "<<points.getElement(punkty2[2]).getZ()<<endl;
				*/

			}


		}

	}

	/*
	for(int i=0,f1;i<ileFB;++i){
			if(facebc_connect[i]!=-1){

					f1=facebc_connect[i];
					if(facebc_connect[f1]!=i){std:cout<<"ble";}

			}
	}
	*/


}



void Delanouy3D::siatkaStatyczna(double dl_x, double dl_y, double dl_z, double il_w_x, double il_w_y, double il_w_z) {

	points.czysc(il_w_x * il_w_y * il_w_z * 0.2, il_w_x * il_w_y * il_w_z + 100);
	elements.czysc(il_w_x * il_w_y * il_w_z, il_w_x * il_w_y * il_w_z * 6);


	double ii = dl_z / (il_w_z - 1);
	double jj = dl_y / (il_w_y - 1);
	double kk = dl_x / (il_w_x - 1);

	for (int i = 0; i < il_w_z; i++) {

		for (int j = 0; j < il_w_y; j++) {

			for (int k = 0; k < il_w_x; k++) {

				points.setElement(k * kk, j * jj, i * ii);

			}
		}
	}

	for (int i = 0; i < il_w_z - 1; i++) {

		for (int j = 0; j < il_w_y - 1; j++) {

			for (int k = 0; k < il_w_x - 1; k++) {


				elements.setElement((i + 1)*il_w_x * il_w_y + k + il_w_x * (j + 1), i * il_w_x * il_w_y + k + il_w_x * (j + 1) + 1, (i + 1)*il_w_x * il_w_y + k + il_w_x * (j + 1) + 1, (i + 1)*il_w_x * il_w_y + k + il_w_x * j);
				elements.setElement(i * il_w_x * il_w_y + k + il_w_x * (j + 1), (i + 1)*il_w_x * il_w_y + k + il_w_x * (j + 1), i * il_w_x * il_w_y + k + il_w_x * j, i * il_w_x * il_w_y + k + il_w_x * (j + 1) + 1);
				elements.setElement((i + 1)*il_w_x * il_w_y + k + il_w_x * j, i * il_w_x * il_w_y + k + il_w_x * j, (i + 1)*il_w_x * il_w_y + k + il_w_x * (j + 1), i * il_w_x * il_w_y + k + il_w_x * (j + 1) + 1);
				elements.setElement((i + 1)*il_w_x * il_w_y + k + 1 + il_w_x * j, (i + 1)*il_w_x * il_w_y + k + il_w_x * (j + 1) + 1, i * il_w_x * il_w_y + k + 1 + il_w_x * j, (i + 1)*il_w_x * il_w_y + k + il_w_x * j);
				elements.setElement(i * il_w_x * il_w_y + k + 1 + il_w_x * j, i * il_w_x * il_w_y + k + il_w_x * (j + 1) + 1, i * il_w_x * il_w_y + k + il_w_x * j, (i + 1)*il_w_x * il_w_y + k + il_w_x * j);
				elements.setElement((i + 1)*il_w_x * il_w_y + k + il_w_x * (j + 1) + 1, (i + 1)*il_w_x * il_w_y + k + il_w_x * j, i * il_w_x * il_w_y + k + il_w_x * (j + 1) + 1, i * il_w_x * il_w_y + k + 1 + il_w_x * j);


			}
		}
	}
//koniec 2

}

void Delanouy3D::wczytajAnsys1(const char * nazwa, double powieksz, double transX, double transY, double transZ) {

	ifstream wczytaj(nazwa);
	wczytaj.precision(10);
	string text;
	int war = 0;
	double xW, yW, zW, numer;
	double maxX = -999999, maxY = -999999, maxZ = -999999, minX = 999999, minY = 999999, minZ = 999999;
	int n1, n2, n3, n4, n5, n6, n7, n8;
	char granicaW;
	int ileP;
	string kontrola, tempStr;
	bool * kon;

	int rodzajZiarna = 1;
	int temp0;
	int flagaSwitch = 1;
	char znakZminay;
	StringList nazwyZiaren(10, 10);

	while (flagaSwitch) {

		switch (flagaSwitch) {
			case 1: {
				flagaSwitch = 0;
				for (int i = 0; i < 100; i++) {

					wczytaj >> kontrola;
					//m1->Lines->Add(kontrola.c_str());
					if (kontrola == "TYPE=C3D4,") {
						flagaSwitch = 3;
						break;
					}
					else if (kontrola == "TYPE=C3D6,") {
						flagaSwitch = 4;
						break;
					}
					else if (kontrola == "*ELSET,") {
						flagaSwitch = 0;
						break;
					}
					else if (kontrola == "*NODE,") {
						flagaSwitch = 2;
						break;
					}

				}
				//m1->Lines->Add("Zmina");
				break;
			}

			case 2: {

				//m1->Lines->Add("punkty");
				wczytaj >> kontrola;
				points.czysc(100000, 100000);
				elements.czysc(500000, 500000);
				pryzmy.czysc(100000, 100000);

				for (int i = 0;; ++i) {

					wczytaj >> numer;
					wczytaj.get();
					wczytaj >> xW;
					wczytaj.get();
					wczytaj >> yW;
					wczytaj.get();
					wczytaj >> zW;
					wczytaj.get();


					xW = xW * powieksz + transX;
					yW = yW * powieksz + transY;
					zW = zW * powieksz + transZ;

					points.setElement(xW, yW, zW, 'a');

					znakZminay = wczytaj.get();
					wczytaj.unget();
					if (znakZminay == '*') {
						flagaSwitch = 1;
						break;
					}
				}

				break;
			}

			case 3: {
				//m1->Lines->Add("elementy");
				wczytaj >> kontrola;

				rodzajZiarna = nazwyZiaren.setElementUniqatReturnPosition_plus1(kontrola);

				for (int i = 0;; ++i) {

					wczytaj >> numer;
					wczytaj.get();

					wczytaj >> n1;
					wczytaj.get();
					wczytaj >> n2;
					wczytaj.get();
					wczytaj >> n3;
					wczytaj.get();
					wczytaj >> n4;
					wczytaj.get();

					elements.setElement(n1 - 1, n2 - 1, n3 - 1, n4 - 1, rodzajZiarna);

					znakZminay = wczytaj.get();
					wczytaj.unget();
					if (znakZminay == '*') {
						flagaSwitch = 1;
						break;
					}
				}

				break;
			}
				/*
					case 4:{

						//m1->Lines->Add("elementy");
						wczytaj>>kontrola;

						rodzajZiarna=nazwyZiaren.setElementUniqatReturnPosition_plus1(kontrola);

						for(int i=0;;++i){

							wczytaj>>numer;wczytaj.get();

							wczytaj>>n1;wczytaj.get();
							wczytaj>>n2;wczytaj.get();
							wczytaj>>n3;wczytaj.get();
							wczytaj>>n4;wczytaj.get();
							wczytaj>>n5;wczytaj.get();
							wczytaj>>n6;wczytaj.get();

							pryzmy.setElement(n1-1,n2-1,n3-1,n4-1,n5-1,n6-1,rodzajZiarna);

							znakZminay = wczytaj.get(); wczytaj.unget();
							if(znakZminay=='*'){flagaSwitch=1;break;}
						}

						break;
					}
					*/


		}//SWITCH

	}//WHILE  ELSET=PID9

	//uzupelnienieElementow();
	wyszukanieSasidnichElementowE();
	oznaczWezlyGraniczneNaPodstawieScianG();

	sasiednieElementy();
	creatFaceAndEdge();
	sasiednieElementy();
	warunkiNaPodstawieZiaren();


	wczytaj.close();

}

void Delanouy3D::warunkiNaPodstawieZiaren_pop_in_out(int numer) {
	double zet = -100;
	for (int i = 0, ileF = faces.getIter(); i < ileF; ++i) {

		if (faces.getElement(i).getBC() == numer) {
			zet = points.getElement(edges.getElement(faces.getElement(i).getEd1()).getP1()).getZ();//break;
		}

	}

	for (int i = 0, ileF = faces.getIter(); i < ileF; ++i) {

		if (faces.getElement(i).getBC() == numer) {
			if (zet != points.getElement(edges.getElement(faces.getElement(i).getEd1()).getP1()).getZ()) {
				faces.getElement(i).setBC(1);
			}
		}

	}


}

void Delanouy3D::warunkiNaPodstawieZiaren_reset(int numer) {

	for (int i = 0, ileF = faces.getIter(); i < ileF; ++i) {

		if (faces.getElement(i).getBC() == numer) {
			faces.getElement(i).setBC(0);
		}

	}


}

void Delanouy3D::warunkiNaPodstawieZiaren() {

	int numerWar = 2;
	IntList przydzial(100, 100);

	for (int i = 0, z_i, ileE = elements.getIter(); i < ileE; ++i) {
		z_i = elements.getElement(i).getRodzajZiarna();

		for (int j = 0, e, z_e, war; j < 4; ++j) {


			e = elements.getElement(i).getE(j);
			if (e != -1) {
				z_e = elements.getElement(e).getRodzajZiarna();
				if (z_i != z_e) {

					war = przydzial.funPodWarunki(z_i, z_e, numerWar);
					faces.getElement(mapFace[i].getElement(j)).setBC(war);
				}
			}
			else {

				war = przydzial.funPodWarunki(z_i, e, numerWar);
				faces.getElement(mapFace[i].getElement(j)).setBC(war);

			}


		}
	}




}


void Delanouy3D::polozeniePunktowWybranegoEl(int el, int localP, double * coor) {

	int p1, p2, p3, p4, pp1, pp2, pp3, pp4;

	pp1 = elements.getElement(el).getP1();
	pp2 = elements.getElement(el).getP2();
	pp3 = elements.getElement(el).getP3();
	pp4 = elements.getElement(el).getP4();

	switch (localP) {

		case 0:
			p1 = pp4;
			p2 = pp1;
			p3 = pp2;
			p4 = pp3;
			break;
		case 1:
			p1 = pp3;
			p2 = pp1;
			p3 = pp4;
			p4 = pp2;
			break;
		case 2:
			p1 = pp2;
			p2 = pp4;
			p3 = pp1;
			p4 = pp3;
			break;
		case 3:
			p1 = pp1;
			p2 = pp4;
			p3 = pp3;
			p4 = pp2;
			break;


	}


	coor[0] = points.getElement(p1).getX();
	coor[1] = points.getElement(p1).getY();
	coor[2] = points.getElement(p1).getZ();

	coor[3] = points.getElement(p2).getX();
	coor[4] = points.getElement(p2).getY();
	coor[5] = points.getElement(p2).getZ();

	coor[6] = points.getElement(p3).getX();
	coor[7] = points.getElement(p3).getY();
	coor[8] = points.getElement(p3).getZ();

	coor[9] = points.getElement(p4).getX();
	coor[10] = points.getElement(p4).getY();
	coor[11] = points.getElement(p4).getZ();



}

void Delanouy3D::kopiujSiatke(Delanouy3D & a) {

	int ileE = elements.getIter();
	int ileP = points.getIter();
	a.dx = dx;
	a.dy = dy;
	a.dz = dz;
	a.ustawWIelkoscObszaru(dx, dy, dz);
	a.elements.czysc(ileE + 10, ileE + 10);
	a.points.czysc(ileP + 10, ileP + 10);
	for (int i = 0; i < ileE; ++i) {
		a.elements.setElement(elements.getElement(i));
	}
	for (int i = 0; i < ileP; ++i) {
		a.points.setElement(points.getElement(i));
	}

}

void Delanouy3D::wyszukajPunktuDoWygladzeniaNaPodstawiePunktow(IntList & punk, IntList & punktyZaz, int ilePok) {
	punk.ustawIter(0);
	sasiedniePunkty();
	int ile = 0, ile1;
	//start pierwsza warstwa
	for (int i = 0, p, ileP = punktyZaz.getIter(); i < ileP; ++i) {
		p = punktyZaz.getElement(i);

		for (int j = 0, ilePP = laplas[p].getIter(), pp; j < ilePP; ++j) {

			pp = laplas[p].getElement(j);
			if (points.getElement(pp).getGranica() == 'a') {
				points.getElement(pp).setGranica('m');
				punk.setElement(pp);
			}

		}

	}

	//kolejne warstwy

	while (--ilePok) {

		ile1 = punk.getIter();
		for (int i = ile, p, ileP = punk.getIter(); i < ileP; ++i) {

			p = punk.getElement(i);

			for (int j = 0, ilePP = laplas[p].getIter(), pp; j < ilePP; ++j) {

				pp = laplas[p].getElement(j);
				if (points.getElement(pp).getGranica() == 'a') {
					points.getElement(pp).setGranica('m');
					punk.setElement(pp);
				}

			}

		}
		ile = ile1;

	}

	for (int i = 0, p, ileP = punk.getIter(); i < ileP; ++i) {
		p = punk.getElement(i);
		if (points.getElement(p).getGranica() == 'm') {
			points.getElement(p).setGranica('a');
		}
	}

}


void Delanouy3D::tworzenieGranicy2D_XY(int rodzaj) {

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		points.getElement(i).setGranica('a');
	}

	if (rodzaj != 0) {
		double z1, z2, z3, z4;

		for (int i = 0, p1, p2, p3, p4, e1, e2, e3, e4, ileE = elements.getIter(); i < ileE; ++i) {


			p1 = elements.getElement(i).getP1();
			p2 = elements.getElement(i).getP2();
			p3 = elements.getElement(i).getP3();
			p4 = elements.getElement(i).getP4();

			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();
			z4 = points.getElement(p4).getZ() * 2;

			e1 = elements.getElement(i).getE1();
			e2 = elements.getElement(i).getE2();
			e3 = elements.getElement(i).getE3();
			e4 = elements.getElement(i).getE4();

			if (e1 <= -1 && (z1 + z2) != z4) {
				points.getElement(p1).setGranica('g');
				points.getElement(p2).setGranica('g');
				points.getElement(p4).setGranica('g');
			}
			if (e3 <= -1 && (z1 + z3) != z4) {
				points.getElement(p1).setGranica('g');
				points.getElement(p3).setGranica('g');
				points.getElement(p4).setGranica('g');
			}
			if (e2 <= -1 && (z2 + z3) != z4) {
				points.getElement(p3).setGranica('g');
				points.getElement(p2).setGranica('g');
				points.getElement(p4).setGranica('g');
			}
			if (e4 <= -1 && (z1 + z3) != (z2 * 2)) {
				points.getElement(p1).setGranica('g');
				points.getElement(p2).setGranica('g');
				points.getElement(p3).setGranica('g');
			}


		}
	}


}

void Delanouy3D::initPrzestrzenRuch(double px0, double py0, double pz0, double px1, double py1, double pz1) {

	double x, y, z;
	int ile = 0;
	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {


		x = points.getElement(i).getX();
		y = points.getElement(i).getY();
		z = points.getElement(i).getZ();

		if (px0 <= x && x <= px1 &&
				py0 <= y && y <= py1 &&
				pz0 <= z && z <= pz1) {

			points.getElement(i).setGranica('b');
			ile++;

		}

	}

	punktyPR.czysc(ile * 0.1, ile + 10);

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'b') {
			points.getElement(i).setGranica('g');
			punktyPR.setElement(i);
		}


	}



}


void Delanouy3D::ruchPrzestrzen(int ileWarstw, int obecny_krok, int od_krok, int ileKrok, double minPoprawy, double px0, double py0, double pz0, double px1, double py1, double pz1, double endX, double endY, double endZ) {

	int do_krok = od_krok;
	do_krok = 1000;
	int ileRazy = 20;

	if (od_krok <= obecny_krok && do_krok >= obecny_krok) {

		if (od_krok == obecny_krok) {
			tworzenieGranicy2D_XY(0);
			initPrzestrzenRuch(px0, py0, pz0, px1, py1, pz1);

		}
		int wart = 0; // = obecny_krok%3;
		//int wart=obecny_krok%3;// = obecny_krok%3;
		if (!wart) {
			//printf("ble \n");
			IntList punktyDoOptym(10000, 10000);
			wyszukajPunktuDoWygladzeniaNaPodstawiePunktow(punktyDoOptym, punktyPR, ileWarstw);

			int znak;
			znak = (obecny_krok / (ileKrok + od_krok)) % 2;
			if (znak) {
				znak = -1;
			}
			else {
				znak = 1;
			}

			ileKrok *= znak;


			double px = (endX - px0) / (ileKrok * ileRazy);
			double py = (endY - py0) / (ileKrok * ileRazy);
			double pz = (endZ - pz0) / (ileKrok * ileRazy);



			for (int il = 0; il < ileRazy; ++il) {

				for (int i = 0, p, ileP = punktyPR.getIter(); i < ileP; ++i) {

					p = punktyPR.getElement(i);
					points.getElement(p).przesunPunkt(px, py, pz);

				}

				///OPTYMALIZACJA
				sasiedniePunkty();
				for (int i = 0; i < 5; ++i) {
					wygladzanieLaplaceWagaWybranePunkty(2, 1, punktyDoOptym, 1);
				}


				sasiednieElementy();
				for (int i = 0; i < 5; ++i) {
					optymalizacjaStosunekRdorNajlepszayWynikCaloscDlaWybranychPunktow(punktyDoOptym, minPoprawy, true);
					//optymalizacjaStosunekRdorNajlepszayWynikCaloscDlaWybranychPunktow_volume(punktyDoOptym,minPoprawy,true);
					//optymalizacjaStosunekRdorNajlepszayWynikCalosc(minPoprawy,true);
				}


			}


		}
	}



}


void Delanouy3D::optymalizacjaStosunekRdorNajlepszayWynikCaloscDlaWybranychPunktow(IntList & wybranePunkty, double dlKroku, bool bezGranicy) {

	char znakG = 'g';
	if (!bezGranicy) {
		znakG = 'w';
	}
	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R = 0;
	double maxR = 0;
	int ktory;

	for (int it = 0, i, ileP = wybranePunkty.getIter(); it < ileP; ++it) {
		i = wybranePunkty.getElement(it);

		if (points.getElement(i).getGranica() != znakG) {
			bool flaga = true;
			int razy = 1000;

			//if(it==0){cout<<"  "<<stosunekRdorNajlepszayWynikDlaPunktu1(i);}


			while (flaga && razy > 0) {

				tX = points.getElement(i).getX();
				tY = points.getElement(i).getY();
				tZ = points.getElement(i).getZ();

				//if(tX==0 || tX==maxX || tY==0 || tY==maxY || tZ==0 || tZ==maxZ){break;}

				--razy;
				maxR = stosunekRdorNajlepszayWynikDlaPunktu(i);
				ktory = 0;


				if (0 < tX - dlKroku && tX != maxX) {

					points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);

					if (maxR > R) {
						maxR = R;
						ktory = 1;
					}
				}

				if (maxX > tX + dlKroku && tX != 0) {

					points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 2;
					}
				}

				if (0 < tY - dlKroku && tY != maxY) {

					points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 3;
					}
				}

				if (maxY > tY + dlKroku && tY != 0) {

					points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 4;
					}
				}


				if (0 < tZ - dlKroku && tZ != maxZ) {

					points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 5;
					}
				}

				if (maxZ > tZ + dlKroku && tZ != 0) {

					points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						ktory = 6;
					}
				}

				switch (ktory) {

					case 0:
						points.getElement(i).setPunkt(tX, tY, tZ);
						flaga = false;
						break;
					case 1:
						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						break;
					case 2:
						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						break;
					case 3:
						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						break;
					case 4:
						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						break;
					case 5:
						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						break;
					case 6:
						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						break;

				}

				//while
			}
		}
	}

}


void Delanouy3D::maxZgranica() {

	double maxZ = -900000;


	for (int i = 0, ile = points.getIter(); i < ile; ++i) {
		points.getElement(i).setGranica('a');
		if (points.getElement(i).getZ() > maxZ) {
			maxZ = points.getElement(i).getZ();
		}
	}

	for (int i = 0, ile = points.getIter(); i < ile; ++i) {
		if (points.getElement(i).getZ() == maxZ) {

			points.getElement(i).setGranica('g');

		}

	}


}

void Delanouy3D::maxZgranica1() {


	for (int i = 0, ile = points.getIter(); i < ile; ++i) {
		if (points.getElement(i).getGranica() == 'g') {

			printf(" nowyZ %f", points.getElement(i).getZ());

		}

	}

	cout << endl << endl << endl;
}

double Delanouy3D::ruchSpawanie(double obecnyKrok, double krok_start, double minPoprawy, double doX, double doY, double dl, double zmiejszaPrzes, double limit, double szerokosc) {

	if (krok_start <= obecnyKrok) {


		if (powSinus(doX, doY, dl, zmiejszaPrzes, limit, szerokosc)) {

			sasiednieElementy();

			for (int i = 0; i < 5; ++i) {
				optymalizacjaStosunekRdorNajlepszayWynikCalosc(minPoprawy * 0.1, true);
			}

		}

	}

	return minZwezglyGraniczne();

}

double Delanouy3D::minZwezglyGraniczne() {

	double minZ = 1000000;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'g') {
			if (points.getElement(i).getZ() < minZ) {
				minZ = points.getElement(i).getZ();
			}
		}
	}

	return minZ;
}

bool Delanouy3D::powSinus(double doX, double doY, double odl, double zmniejszPrzesuniecie, double limit, double szerokosc) {


	double zmiena = 0;
	double X, Y, R, Z, pi = 3.1415926535, Nz, maxZ = 1000, dP;
	for (int i = 0, ile = points.getIter(); i < ile; ++i) {

		if (points.getElement(i).getGranica() == 'g') {

			X = points.getElement(i).getX();
			Y = points.getElement(i).getY();
			Z = points.getElement(i).getZ();

			//printf(" Z %f", Z);

			if (sqrt((X - doX) * (X - doX) + (Y - doY) * (Y - doY)) <= odl * (1 / szerokosc)) {

				X -= doX;
				Y -= doY;
				dP =  (sqrt( (X * pi) * (X * pi) + (Y * pi) * (Y * pi) )) * szerokosc;
				Z = sin(dP) / dP ;
				Nz = points.getElement(i).getZ() - Z * zmniejszPrzesuniecie;

				// points.getElement(i).przesunPunkt(0,0,-Z*zmniejszPrzesuniecie);

				if (Nz < maxZ) {
					maxZ = Nz;
				}

			}
		}

	}

	//m1->Lines->Add(maxZ);

	// printf(" maxZ %f, limit %f", maxZ, limit);

	if (maxZ > limit) {
		for (int i = 0, ile = points.getIter(); i < ile; ++i) {

			if (points.getElement(i).getGranica() == 'g') {

				X = points.getElement(i).getX();
				Y = points.getElement(i).getY();
				if (sqrt((X - doX) * (X - doX) + (Y - doY) * (Y - doY)) <= odl * (1 / szerokosc)) {

					X -= doX;
					Y -= doY;
					dP =  (sqrt( (X * pi) * (X * pi) + (Y * pi) * (Y * pi) )) * szerokosc;
					Z = sin(dP) / dP ;
					Nz = points.getElement(i).getZ() - Z;

					points.getElement(i).przesunPunkt(0, 0, -Z * zmniejszPrzesuniecie);

					if (Nz > maxZ) {
						maxZ = Nz;
					}

				}
			}

		}
	}
	else {
		return false;
	}


	return true;
}


double * Delanouy3D::zworcNormalnaDoFace(int p1, int p2, int p3) {

	double * vec_norm;
	vec_norm = new double[4];

	double vA[3];
	double vB[3];

	vA[0] = points.getElement(p2).getX() - points.getElement(p1).getX();
	vA[1] = points.getElement(p2).getY() - points.getElement(p1).getY();
	vA[2] = points.getElement(p2).getZ() - points.getElement(p1).getZ();

	vB[0] = points.getElement(p2).getX() - points.getElement(p3).getX();
	vB[1] = points.getElement(p2).getY() - points.getElement(p3).getY();
	vB[2] = points.getElement(p2).getZ() - points.getElement(p3).getZ();

	vec_norm[0] = vA[1] * vB[2] - vA[2] * vB[1];
	vec_norm[1] = vA[2] * vB[0] - vA[0] * vB[2];
	vec_norm[2] = vA[0] * vB[1] - vA[1] * vB[0];

	vec_norm[3] = sqrt( vec_norm[0] * vec_norm[0] + vec_norm[1] * vec_norm[1] + vec_norm[2] * vec_norm[2] );
	vec_norm[0] /= vec_norm[3];
	vec_norm[1] /= vec_norm[3];
	vec_norm[2] /= vec_norm[3];

	return vec_norm;

}

void Delanouy3D::wyszukajElementsWithEdge(int nrEd, int * nrElements) {

	if (flagaLaplasa != 1) {
		sasiednieElementy();
	}

	int ed1, ed2;

	ed1 = edges.getElement(nrEd).getP1();
	ed2 = edges.getElement(nrEd).getP2();
	int ii = 0;


	for (int i = 0, ii = 1, ileI = laplas[ed1].getIter(), el; i < ileI; ++i) {

		el = laplas[ed1].getElement(i);
		if (elements.getElement(el).sprawdzCzyJestTakiPunkt(ed2)) {
			nrElements[++ii] = el;
		}

	}

	nrElements[0] = ii;

}


int Delanouy3D::wyszukajEdgesInEl(int nrEl, int * nrEdges) {

	int nrF = mapFace[nrEl].getElement(0);

	nrEdges[0] = 6;
	nrEdges[6] = -10;
	nrEdges[1] = faces.getElement(nrF).getEd1();
	nrEdges[2] = faces.getElement(nrF).getEd2();
	nrEdges[3] = faces.getElement(nrF).getEd3();

	nrF = mapFace[nrEl].getElement(1);

	if (faces.getElement(nrF).getEd1() == nrEdges[1]) {
		nrEdges[4] = faces.getElement(nrF).getEd2();
		nrEdges[5] = faces.getElement(nrF).getEd3();
	}
	else if (faces.getElement(nrF).getEd2() == nrEdges[1]) {
		nrEdges[4] = faces.getElement(nrF).getEd1();
		nrEdges[5] = faces.getElement(nrF).getEd3();
	}
	else if (faces.getElement(nrF).getEd3() == nrEdges[1]) {
		nrEdges[4] = faces.getElement(nrF).getEd1();
		nrEdges[5] = faces.getElement(nrF).getEd2();
	}

	nrF = mapFace[nrEl].getElement(2);

	if (faces.getElement(nrF).getEd1() == nrEdges[1]) {


		if (nrEdges[4] == faces.getElement(nrF).getEd2()) {
			nrEdges[6] = faces.getElement(nrF).getEd3();
		}
		else if (nrEdges[4] == faces.getElement(nrF).getEd3()) {
			nrEdges[6] = faces.getElement(nrF).getEd2();
		}
		else if (nrEdges[5] == faces.getElement(nrF).getEd2()) {
			nrEdges[6] = faces.getElement(nrF).getEd3();
		}
		else if (nrEdges[5] == faces.getElement(nrF).getEd3()) {
			nrEdges[6] = faces.getElement(nrF).getEd2();
		}

	}
	else if (faces.getElement(nrF).getEd2() == nrEdges[1]) {

		if (nrEdges[4] == faces.getElement(nrF).getEd1()) {
			nrEdges[6] = faces.getElement(nrF).getEd3();
		}
		else if (nrEdges[4] == faces.getElement(nrF).getEd3()) {
			nrEdges[6] = faces.getElement(nrF).getEd1();
		}
		else if (nrEdges[5] == faces.getElement(nrF).getEd1()) {
			nrEdges[6] = faces.getElement(nrF).getEd3();
		}
		else if (nrEdges[5] == faces.getElement(nrF).getEd3()) {
			nrEdges[6] = faces.getElement(nrF).getEd1();
		}

	}
	else if (faces.getElement(nrF).getEd3() == nrEdges[1]) {

		if (nrEdges[4] == faces.getElement(nrF).getEd1()) {
			nrEdges[6] = faces.getElement(nrF).getEd2();
		}
		else if (nrEdges[4] == faces.getElement(nrF).getEd2()) {
			nrEdges[6] = faces.getElement(nrF).getEd1();
		}
		else if (nrEdges[5] == faces.getElement(nrF).getEd1()) {
			nrEdges[6] = faces.getElement(nrF).getEd2();
		}
		else if (nrEdges[5] == faces.getElement(nrF).getEd2()) {
			nrEdges[6] = faces.getElement(nrF).getEd1();
		}

	}

	if (nrEdges[6] == -10) {
		return 0;
	}
	return 1;

}



double * Delanouy3D::ustawaElementyNaSrodku(double stosunekWielkoscPlanszy) {


	if (stosunekWielkoscPlanszy < 1) {
		stosunekWielkoscPlanszy = 1;
	}

	double xW, yW, zW; //numer
	double maxX = -999999, maxY = -999999, maxZ = -999999, minX = 999999, minY = 999999, minZ = 999999;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		xW = points.getElement(i).getX();
		yW = points.getElement(i).getY();
		zW = points.getElement(i).getZ();

		if (maxX < xW) {
			maxX = xW;
		}
		if (maxY < yW) {
			maxY = yW;
		}
		if (maxZ < zW) {
			maxZ = zW;
		}

		if (minX > xW) {
			minX = xW;
		}
		if (minY > yW) {
			minY = yW;
		}
		if (minZ > zW) {
			minZ = zW;
		}

	}

	double * wielkosc = new double[3];

	wielkosc[0] = (maxX - minX) * 0.5 * stosunekWielkoscPlanszy ;
	wielkosc[1] = (maxY - minY) * 0.5 * stosunekWielkoscPlanszy ;
	wielkosc[2] = (maxZ - minZ) * 0.5 * stosunekWielkoscPlanszy ;

	xW = minX + (maxX - minX) * 0.5 ;
	yW = minY + (maxY - minY) * 0.5 ;
	zW = minZ + (maxZ - minZ) * 0.5 ;
	xW = wielkosc[0] - xW ;
	yW = wielkosc[1] - yW ;
	zW = wielkosc[2] - zW ;

	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {
		points.getElement(i).przesunPunkt(xW, yW, zW);
	}

	wielkosc[0] *= 2 ;
	wielkosc[1] *= 2 ;
	wielkosc[2] *= 2 ;

	globalTranX = xW;
	globalTranY = yW;
	globalTranZ = zW;

	return wielkosc;

}

void Delanouy3D::wyszukajPowZew_zmienKoloryWezlow(char ** mapE_F) {

	int ileE = elements.getIter();
	char zero = 0, temp;
	for (int i = 0, p1, p2, p3, p4; i < ileE; ++i) {
		for (int ii = 0; ii < 4; ++ii) {

			if (mapE_F[i][ii] != zero) {
				temp = mapE_F[i][ii];
				p1 = elements.getElement(i).getP1();
				p2 = elements.getElement(i).getP2();
				p3 = elements.getElement(i).getP3();
				p4 = elements.getElement(i).getP4();

				switch (ii) {

					case 0:
						points.getElement(p4).setGranica(temp);
						points.getElement(p2).setGranica(temp);
						points.getElement(p1).setGranica(temp);
						break;
					case 1:
						points.getElement(p4).setGranica(temp);
						points.getElement(p3).setGranica(temp);
						points.getElement(p2).setGranica(temp);
						break;
					case 2:
						points.getElement(p4).setGranica(temp);
						points.getElement(p1).setGranica(temp);
						points.getElement(p3).setGranica(temp);
						break;
					case 3:
						points.getElement(p1).setGranica(temp);
						points.getElement(p2).setGranica(temp);
						points.getElement(p3).setGranica(temp);
						break;
				}

			}
		}
	}

}

void Delanouy3D::wyszukajPowZew_wyszkiwaniePow(char ** mapE_F, IntList & pkScian, bool * pointFlag, char obecnyNumer) {

	char scianaPusta = 1;

	for (int i = 0, max = pkScian.getIter(); i < max; ) {

		for (int nrP, ileAktualnie = pkScian.getIter(); i < ileAktualnie ; i++) {

			nrP = pkScian.getElement(i);

			for (int j = 0, nrE, p1, p2, p3, p4, ile = laplas[nrP].getIter() ; j < ile ; ++j) {
				nrE = laplas[nrP].getElement(j);

				p1 = elements.getElement(nrE).getP1();
				p2 = elements.getElement(nrE).getP2();
				p3 = elements.getElement(nrE).getP3();
				p4 = elements.getElement(nrE).getP4();

				if (mapE_F[nrE][0] == scianaPusta) {
					if (p4 == nrP || p2 == nrP || p1 == nrP) {
						mapE_F[nrE][0] = obecnyNumer;
						if (pointFlag[p4]) {
							pkScian.setElement(p4);
						}
						if (pointFlag[p2]) {
							pkScian.setElement(p2);
						}
						if (pointFlag[p1]) {
							pkScian.setElement(p1);
						}
					}
				}

				if (mapE_F[nrE][1] == scianaPusta) {
					if (p4 == nrP || p3 == nrP || p2 == nrP) {
						mapE_F[nrE][1] = obecnyNumer;
						if (pointFlag[p4]) {
							pkScian.setElement(p4);
						}
						if (pointFlag[p3]) {
							pkScian.setElement(p3);
						}
						if (pointFlag[p2]) {
							pkScian.setElement(p2);
						}
					}
				}

				if (mapE_F[nrE][2] == scianaPusta) {
					if (p4 == nrP || p1 == nrP || p3 == nrP) {
						mapE_F[nrE][2] = obecnyNumer;
						if (pointFlag[p4]) {
							pkScian.setElement(p4);
						}
						if (pointFlag[p1]) {
							pkScian.setElement(p1);
						}
						if (pointFlag[p3]) {
							pkScian.setElement(p3);
						}
					}
				}

				if (mapE_F[nrE][3] == scianaPusta) {
					if (p1 == nrP || p2 == nrP || p3 == nrP) {
						mapE_F[nrE][3] = obecnyNumer;
						if (pointFlag[p1]) {
							pkScian.setElement(p1);
						}
						if (pointFlag[p2]) {
							pkScian.setElement(p2);
						}
						if (pointFlag[p3]) {
							pkScian.setElement(p3);
						}
					}
				}
			}
		}

		max = pkScian.getIter();
	}

}

void Delanouy3D::wyszukajPowZew_wygladScianyGranicznej(int rodzajGranicy, char ** mapE_F) {

	char scianaPusta = 1;
	int ileE = elements.getIter();
//ustawienia
	switch (rodzajGranicy) {

		case 1:

			for (int i = 0, e1, e2, e3, e4; i < ileE; ++i) {

				e1 = elements.getElement(i).getE1();
				e2 = elements.getElement(i).getE2();
				e3 = elements.getElement(i).getE3();
				e4 = elements.getElement(i).getE4();

				if (e1 <= -1) {
					mapE_F[i][0] = scianaPusta;
				}
				if (e2 <= -1) {
					mapE_F[i][1] = scianaPusta;
				}
				if (e3 <= -1) {
					mapE_F[i][2] = scianaPusta;
				}
				if (e4 <= -1) {
					mapE_F[i][3] = scianaPusta;
				}

			}
			break;

		case 2:

			for (int i = 0, e1, e2, e3, e4, z; i < ileE; ++i) {

				z = elements.getElement(i).getRodzajZiarna();
				e1 = elements.getElement(i).getE1();
				e2 = elements.getElement(i).getE2();
				e3 = elements.getElement(i).getE3();
				e4 = elements.getElement(i).getE4();

				if (e1 != -1 && elements.getElement(e1).getRodzajZiarna() != z) {
					mapE_F[i][0] = scianaPusta;
				}
				if (e2 != -1 && elements.getElement(e2).getRodzajZiarna() != z) {
					mapE_F[i][1] = scianaPusta;
				}
				if (e3 != -1 && elements.getElement(e3).getRodzajZiarna() != z) {
					mapE_F[i][2] = scianaPusta;
				}
				if (e4 != -1 && elements.getElement(e4).getRodzajZiarna() != z) {
					mapE_F[i][3] = scianaPusta;
				}

			}
			break;

		case 3:

			for (int i = 0, e1, e2, e3, e4, z; i < ileE; ++i) {

				z = elements.getElement(i).getRodzajZiarna();
				e1 = elements.getElement(i).getE1();
				e2 = elements.getElement(i).getE2();
				e3 = elements.getElement(i).getE3();
				e4 = elements.getElement(i).getE4();

				if (e1 <= -1 || elements.getElement(e1).getRodzajZiarna() != z) {
					mapE_F[i][0] = scianaPusta;
				}
				if (e2 <= -1 || elements.getElement(e2).getRodzajZiarna() != z) {
					mapE_F[i][1] = scianaPusta;
				}
				if (e3 <= -1 || elements.getElement(e3).getRodzajZiarna() != z) {
					mapE_F[i][2] = scianaPusta;
				}
				if (e4 <= -1 || elements.getElement(e4).getRodzajZiarna() != z) {
					mapE_F[i][3] = scianaPusta;
				}

			}

			break;
	}

}

void Delanouy3D::wyszukajPowZew_start(int rodzajGranicy) {

//rodzajGranicy=1 powierznie
//rodzajGranicy=2 granica miêdzy ziarnami
//rodzajGranicy=3 granica i powierzchnia

	sasiednieElementy();

	int iterCharStart = 1;
	char element = 0;
	char scianaPusta = 1;
	int ileE = elements.getIter();
	int ileP = points.getIter();
	IntList punktyScian(100000, 100000);

//ustawienie MapE_F
	bool * pointFlag = new bool[ileP];
	for (int i = 0; i < ileP; ++i) {
		pointFlag[i] = true;
	}
	char ** mapE_F;
	mapE_F = new char * [ileE];
	for (int i = 0; i < ileE; ++i) {
		mapE_F[i] = new char[4];
	}
	for (int i = 0; i < ileE; ++i) {
		mapE_F[i][0] = element;
		mapE_F[i][1] = element;
		mapE_F[i][2] = element;
		mapE_F[i][3] = element;
	}

	wyszukajPowZew_wygladScianyGranicznej(rodzajGranicy, mapE_F);


	for (int i = 0; i < ileE; ++i) {

		for (int j = 0, p1, p2, p3, p4; j < 4; ++j) {

			if (mapE_F[i][j] == scianaPusta) {

				++iterCharStart;

				//m1->Lines->Add(int(iterCharStart));

				punktyScian.ustawIter(0);
				mapE_F[i][j] = iterCharStart;

				p1 = elements.getElement(i).getP1();
				p2 = elements.getElement(i).getP2();
				p3 = elements.getElement(i).getP3();
				p4 = elements.getElement(i).getP4();

				switch (j) {

					case 0:
						punktyScian.setElement(p4);
						punktyScian.setElement(p2);
						punktyScian.setElement(p1);
						pointFlag[p4] = false;
						pointFlag[p2] = false;
						pointFlag[p1] = false;
						break;
					case 1:
						punktyScian.setElement(p4);
						punktyScian.setElement(p3);
						punktyScian.setElement(p2);
						pointFlag[p4] = false;
						pointFlag[p3] = false;
						pointFlag[p2] = false;
						break;
					case 2:
						punktyScian.setElement(p4);
						punktyScian.setElement(p1);
						punktyScian.setElement(p3);
						pointFlag[p4] = false;
						pointFlag[p1] = false;
						pointFlag[p3] = false;
						break;
					case 3:
						punktyScian.setElement(p1);
						punktyScian.setElement(p2);
						punktyScian.setElement(p3);
						pointFlag[p1] = false;
						pointFlag[p2] = false;
						pointFlag[p3] = false;
						break;
				}

				wyszukajPowZew_wyszkiwaniePow(mapE_F, punktyScian, pointFlag, iterCharStart);

			}

		}
	}


	wyszukajPowZew_zmienKoloryWezlow(mapE_F);


}



void Delanouy3D::czyscFace() {
	delete []mapFace;
	mapFace = new IntList[1];
	faces.czysc(10, 10);
}
void Delanouy3D::czyscEdge() {
	edges.czysc(10, 10);
}

void Delanouy3D::czyscFaceAndEdge() {

	delete []mapFace;
	mapFace = new IntList[1];
	faces.czysc(10, 10);
	edges.czysc(10, 10);

}

void Delanouy3D::creatFaceAndEdge() {

	delete []mapFace;

	int ileP = points.getIter();
	int ileE = elements.getIter();
	IntList * creatEdge;
	mapFace = new IntList[ileE];

	faces.czysc(ileE * 1.1, ileE * 1.1);
	edges.czysc(ileE * 1.1, ileE * 1.1);

	creatEdge = new IntList[ileP];
	bool ** flagaFace;
	flagaFace = new bool*[ileE];
	for (int i = 0; i < ileE; ++i) {
		flagaFace[i] = new bool[4];
	}
	for (int i = 0; i < ileE; ++i) {
		for (int j = 0; j < 4; ++j) {
			flagaFace[i][j] = true;
		}
	}



	for (int i = 0; i < ileE; ++i) {
		mapFace[i].setElement(-1);
		mapFace[i].setElement(-1);
		mapFace[i].setElement(-1);
		mapFace[i].setElement(-1);
	}


//creat edges
	for (int i = 0, p1, p2, p3, p4, iterEdge = 0; i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		//42 21 14 43 32 31


		if (creatEdge[p4].setElementUniqatEvenNumbers(p2)) {
			creatEdge[p4].setElement(iterEdge);
			creatEdge[p2].setElement(p4);
			creatEdge[p2].setElement(iterEdge);
			edges.setElement(p4, p2, -1);
			++iterEdge;
		}
		if (creatEdge[p2].setElementUniqatEvenNumbers(p1)) {
			creatEdge[p2].setElement(iterEdge);
			creatEdge[p1].setElement(p2);
			creatEdge[p1].setElement(iterEdge);
			edges.setElement(p2, p1, -1);
			++iterEdge;
		}
		if (creatEdge[p1].setElementUniqatEvenNumbers(p4)) {
			creatEdge[p1].setElement(iterEdge);
			creatEdge[p4].setElement(p1);
			creatEdge[p4].setElement(iterEdge);
			edges.setElement(p1, p4, -1);
			++iterEdge;
		}
		if (creatEdge[p4].setElementUniqatEvenNumbers(p3)) {
			creatEdge[p4].setElement(iterEdge);
			creatEdge[p3].setElement(p4);
			creatEdge[p3].setElement(iterEdge);
			edges.setElement(p4, p3, -1);
			++iterEdge;
		}
		if (creatEdge[p3].setElementUniqatEvenNumbers(p2)) {
			creatEdge[p3].setElement(iterEdge);
			creatEdge[p2].setElement(p3);
			creatEdge[p2].setElement(iterEdge);
			edges.setElement(p3, p2, -1);
			++iterEdge;
		}
		if (creatEdge[p3].setElementUniqatEvenNumbers(p1)) {
			creatEdge[p3].setElement(iterEdge);
			creatEdge[p1].setElement(p3);
			creatEdge[p1].setElement(iterEdge);
			edges.setElement(p3, p1, -1);
			++iterEdge;
		}

	}


//creat faces
	for (int i = 0, p1, p2, p3, p4, e1, e2, e3, e4, temp, nrF = 0; i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		//421 432 413 123

		if (flagaFace[i][0]) {
			faces.setElement(-1, -1, -1, 'a', i);
			mapFace[i].podmienElement(0, nrF++);
			if (e1 >= 0) {
				temp = elements.getElement(e1).getNumerSciany(p4, p2, p1);
				flagaFace[e1][temp] = false;
				mapFace[e1].podmienElement(temp, nrF - 1);
			}
		}
		if (flagaFace[i][1]) {
			faces.setElement(-1, -1, -1, 'b', i);
			mapFace[i].podmienElement(1, nrF++);
			if (e2 >= 0) {
				temp = elements.getElement(e2).getNumerSciany(p4, p3, p2);
				flagaFace[e2][temp] = false;
				mapFace[e2].podmienElement(temp, nrF - 1);
			}
		}
		if (flagaFace[i][2]) {
			faces.setElement(-1, -1, -1, 'c', i);
			mapFace[i].podmienElement(2, nrF++);
			if (e3 >= 0) {
				temp = elements.getElement(e3).getNumerSciany(p4, p1, p3);
				flagaFace[e3][temp] = false;
				mapFace[e3].podmienElement(temp, nrF - 1);
			}
		}
		if (flagaFace[i][3]) {
			faces.setElement(-1, -1, -1, 'd', i);
			mapFace[i].podmienElement(3, nrF++);
			if (e4 >= 0) {
				temp = elements.getElement(e4).getNumerSciany(p1, p2, p3);
				flagaFace[e4][temp] = false;
				mapFace[e4].podmienElement(temp, nrF - 1);
			}
		}
	}
//numbered edges in faces


	for (int i = 0, el, p1 = 0, p2 = 0, p3 = 0, ed1 = 0, ed2 = 0, ed3 = 0, ileF = faces.getIter(); i < ileF; ++i) {

		el = faces.getElement(i).getInEl();

		switch (faces.getElement(i).getFaceInEl()) {

			case 'a':
				p1 = elements.getElement(el).getP4();
				p2 = elements.getElement(el).getP2();
				p3 = elements.getElement(el).getP1();
				break;

			case 'b':
				p1 = elements.getElement(el).getP4();
				p2 = elements.getElement(el).getP3();
				p3 = elements.getElement(el).getP2();
				break;

			case 'c':
				p1 = elements.getElement(el).getP4();
				p2 = elements.getElement(el).getP1();
				p3 = elements.getElement(el).getP3();
				break;

			case 'd':
				p1 = elements.getElement(el).getP1();
				p2 = elements.getElement(el).getP2();
				p3 = elements.getElement(el).getP3();
				break;

		}

		for (int j = 0, ile = creatEdge[p1].getIter(); j < ile; j += 2) {
			if (creatEdge[p1].getElement(j) == p2) {
				ed1 = creatEdge[p1].getElement(j + 1);
				break;
			}
		}
		for (int j = 0, ile = creatEdge[p1].getIter(); j < ile; j += 2) {
			if (creatEdge[p1].getElement(j) == p3) {
				ed2 = creatEdge[p1].getElement(j + 1);
				break;
			}
		}
		for (int j = 0, ile = creatEdge[p2].getIter(); j < ile; j += 2) {
			if (creatEdge[p2].getElement(j) == p3) {
				ed3 = creatEdge[p2].getElement(j + 1);
				break;
			}
		}


		faces.getElement(i).setAllE(ed1, ed2, ed3);

	}




//testy faces
	for (int i = 0, el, ile = faces.getIter(); i < ile; ++i) {

		el = faces.getElement(i).getInEl();
		if (-1 == faces.getElement(i).getEd1()) {}
		else {
			edges.getElement(faces.getElement(i).getEd1()).setEl(el);
		}
		if (-1 == faces.getElement(i).getEd2()) {}
		else {
			edges.getElement(faces.getElement(i).getEd2()).setEl(el);
		}
		if (-1 == faces.getElement(i).getEd3()) {}
		else {
			edges.getElement(faces.getElement(i).getEd3()).setEl(el);
		}
	}

//testy edges
	/*
	for(int i=0,ile=edges.getIter();i<ile;++i){
		if(-1==edges.getElement(i).getEl()){m1->Lines->Add("blad dopasowania edge do edge");}

	}
	*/

//test mapFeces
	/*
	for(int i=0;i<ileE;++i){
		if(mapFace[i].getElement(0)==-1){m1->Lines->Add("blad mapyFaces");}
		if(mapFace[i].getElement(1)==-1){m1->Lines->Add("blad mapyFaces");}
		if(mapFace[i].getElement(2)==-1){m1->Lines->Add("blad mapyFaces");}
		if(mapFace[i].getElement(3)==-1){m1->Lines->Add("blad mapyFaces");}
	}
	m1->Lines->Add("ile faces , edges");
	m1->Lines->Add(faces.getIter());
	m1->Lines->Add(edges.getIter());
	*/

	for (int i = 0; i < ileE; ++i) {
		delete []flagaFace[i];
	}
	delete []flagaFace;
	delete []creatEdge;

}

void Delanouy3D::wczytajPrzerob(const char * nazwa) {

	PunktList noweP(100000, 100000);
	noweP.wczytajZPlikuPrzerob(nazwa);


}




void Delanouy3D::pobierzWarunkiZFacow(IntList & warunkiBrzegowe, int * mapaPunktow1) {

	int ileFwBC = 0;
	for (int i = 0, ileI = faces.getIter(); i < ileI; ++i) {

		if (faces.getElement(i).getBC()) {
			++ileFwBC;
		}

	}

	warunkiBrzegowe.czysc(ileFwBC, ileFwBC * 4 + 100);


	char nf;
	for (int i = 0, p1, p2, p3, el, bc, ileI = faces.getIter(); i < ileI; ++i) {

		bc = faces.getElement(i).getBC();

		if (bc) {

			el = faces.getElement(i).getInEl();
			nf = faces.getElement(i).getFaceInEl();

			switch (nf) {
				case 'a':
					p1 = elements.getElement(el).getP4();
					p2 = elements.getElement(el).getP2();
					p3 = elements.getElement(el).getP1();

					break;
				case 'b':
					p1 = elements.getElement(el).getP4();
					p2 = elements.getElement(el).getP3();
					p3 = elements.getElement(el).getP2();
					break;
				case 'c':
					p1 = elements.getElement(el).getP4();
					p2 = elements.getElement(el).getP1();
					p3 = elements.getElement(el).getP3();
					break;
				case 'd':
					p1 = elements.getElement(el).getP1();
					p2 = elements.getElement(el).getP2();
					p3 = elements.getElement(el).getP3();
					break;
					//switch
			}

			p1 = mapaPunktow1[p1];
			p2 = mapaPunktow1[p2];
			p3 = mapaPunktow1[p3];

			warunkiBrzegowe.setElement(bc);
			warunkiBrzegowe.setElement(p1);
			warunkiBrzegowe.setElement(p2);
			warunkiBrzegowe.setElement(p3);

		}
	}

}

void Delanouy3D::pobierzWarunkiZFacow(IntList & warunkiBrzegowe) {

	int ileFwBC = 0;
	for (int i = 0, ileI = faces.getIter(); i < ileI; ++i) {
		if (faces.getElement(i).getBC()) {
			++ileFwBC;
		}
	}

	warunkiBrzegowe.czysc(ileFwBC, ileFwBC * 4 + 100);

	char nf;
	for (int i = 0, p1, p2, p3, el, bc, ileI = faces.getIter(); i < ileI; ++i) {

		bc = faces.getElement(i).getBC();

		if (bc) {

			el = faces.getElement(i).getInEl();
			nf = faces.getElement(i).getFaceInEl();

			switch (nf) {
				case 'a':
					p1 = elements.getElement(el).getP4();
					p2 = elements.getElement(el).getP2();
					p3 = elements.getElement(el).getP1();

					break;
				case 'b':
					p1 = elements.getElement(el).getP4();
					p2 = elements.getElement(el).getP3();
					p3 = elements.getElement(el).getP2();
					break;
				case 'c':
					p1 = elements.getElement(el).getP4();
					p2 = elements.getElement(el).getP1();
					p3 = elements.getElement(el).getP3();
					break;
				case 'd':
					p1 = elements.getElement(el).getP1();
					p2 = elements.getElement(el).getP2();
					p3 = elements.getElement(el).getP3();
					break;
					//switch
			}

			warunkiBrzegowe.setElement(bc);
			warunkiBrzegowe.setElement(p1);
			warunkiBrzegowe.setElement(p2);
			warunkiBrzegowe.setElement(p3);

		}
	}

}

void Delanouy3D::ZapiszDoPlikuNasZPrzesunieciem_BC_Face(const char * nazwa, double dzielnikZ, double transX, double transY, double transZ, double globalnyDzielnik) {

	const char * pause = "      ";
	double dzielnik;

	dzielnik = 1 / globalnyDzielnik;
	dzielnikZ = 1 / dzielnikZ;

	IntList facesTemp;
	pobierzWarunkiZFacow(facesTemp);


	ofstream zapis(nazwa);
	zapis.precision(12);


	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << (points.getElement(i).getX() - transX)*dzielnik << " " << pause << (points.getElement(i).getY() - transY)*dzielnik << " " << pause << (points.getElement(i).getZ() - transZ)*dzielnik * dzielnikZ << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {
		//if(elements.getElement(i).getRodzajZiarna()==3){t1=2;}
		//else{t1=1;}
		zapis << "CTETRA" << pause << i + 1 << " " << pause << elements.getElement(i).getRodzajZiarna() << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	zapis << endl << endl;


	for (int j = 0, rodzaj = -1, ileF = facesTemp.getIter(); j < ileF; ++i, j += 4) {

		rodzaj = facesTemp.getElement(j);
		if (rodzaj) {
			zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << facesTemp.getElement(j + 1) + 1 << " " << pause << facesTemp.getElement(j + 2) + 1 << " " << pause << facesTemp.getElement(j + 3) + 1 << endl;
		}
	}

	if (facebc_connect[0] != -2) {

		for (int i = 0, p1[3], p2[3], f, e1, e2, ile = faces.getIter(); i < ile; ++i) {

			if (facebc_connect[i] != -1) {
				f = facebc_connect[i];
				e1 = faces.getElement(i).getInEl();
				e2 = faces.getElement(f).getInEl();

				elements.getElement(e1).getNumeryElSciany(faces.getElement(i).getFaceInEl(), p1);
				elements.getElement(e2).getNumeryElSciany(faces.getElement(f).getFaceInEl(), p2);

				zapis << "CTRIAX_BC" << pause << e1 + 1 << pause << p1[0] + 1 << pause << p1[1] + 1 << pause << p1[2] + 1 << pause << e2 + 1 << pause << p2[0] + 1 << pause << p2[1] + 1 << pause << p2[2] + 1 << endl;
			}

		}

	}

	/*
	for(int i=0,ileE=elements.getIter();i<ileE;++i){
		for(int ii=0,e;ii<4;++ii){
				e=elements.getElement(i).getE(ii);
				if(e<-1){
						e = -e-2;
						zapis<<"CTRIAX_BC"<<pause<<i+1<<" "<<pause<<ii<<pause<<e+1<<endl;
				}
		}
	}
	 */


	zapis.close();
}



void Delanouy3D::ZapiszDoPlikuNasZPrzesunieciem(const char * nazwa, double dzielnikZ, bool scianyProstopadl, int transX, int transY, int transZ, double globalnyDzielnik, int * warunki) {

	int tabWarNULL[6] = {1, 2, 3, 4, 5, 6};
	int * war;

	if (warunki == NULL) {
		war = tabWarNULL;
	}
	else {
		war = warunki;
	}

	double maxX =  dx;
	double maxY =  dy;
	double maxZ =  dz;

	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(elements.getIter() * 0.1, elements.getIter() * 0.1);
	const char * pause = "      ";
	double dzielnik = 1;

	dzielnik = 1 / globalnyDzielnik;
	dzielnikZ = 1 / dzielnikZ;


	/*
	for(int i=0,ileP=points.getIter();i<ileP;++i){

		if(points.getElement(i).getZ()>maxZ){maxZ=points.getElement(i).getZ();}

	}
	*/

	for (int i = 0, e1, e2, e3, e4, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; i++) {

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();


		if (e1 <= -1) {
			face.setElement(p4, p2, p1, 'b');
		}
		if (e2 <= -1) {
			face.setElement(p4, p3, p2, 'b');
		}
		if (e3 <= -1) {
			face.setElement(p4, p1, p3, 'b');
		}
		if (e4 <= -1) {
			face.setElement(p1, p2, p3, 'b');
		}

	}



	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM



	for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
		if (face.getElement(i).getGranica() == 'b') {

			p1 = face.getElement(i).getX();
			p2 = face.getElement(i).getY();
			p3 = face.getElement(i).getZ();

			x1 = points.getElement(p1).getX();
			x2 = points.getElement(p2).getX();
			x3 = points.getElement(p3).getX();

			y1 = points.getElement(p1).getY();
			y2 = points.getElement(p2).getY();
			y3 = points.getElement(p3).getY();

			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();

			if (scianyProstopadl) {

				if (x1 == x2 && x2 == x3) {
					if (x1 == 0) {
						face.getElement(i).setGranica('c');
					}
					else if (x1 == maxX) {
						face.getElement(i).setGranica('d');
					}
				}
				else if (y1 == y2 && y2 == y3) {
					if (y1 == 0) {
						face.getElement(i).setGranica('e');
					}
					else if (y1 == maxY) {
						face.getElement(i).setGranica('f');
					}
				}
				else if (z1 == z2 && z2 == z3) {
					if (z1 == 0) {
						face.getElement(i).setGranica('g');
					}
					else if (z1 == maxZ) {
						face.getElement(i).setGranica('h');
					}
				}
			}
			else {


				/*
				if(y1 < x1*3.5087 + 40.35087 && y2 < x2*3.5087 + 40.35087 && y3 < x3*3.5087 + 40.35087){

				   if(z1>45){
					   face.getElement(i).setGranica('g');
				   }
				   else{
					   face.getElement(i).setGranica('h');
				   }

				}
				*/

				/*
				//z zastawka
				if(z1 >= x1*(-4.25531) + 461.702 && z2 >= x2*(-4.25531) + 461.702 && z3 >= x3*(-4.25531) + 461.702){

					if(y1>45){
						face.getElement(i).setGranica('g');;
					}
					else{
						face.getElement(i).setGranica('h');;
					}

				}
				*/



				if ((z1 + z2 + z3) / 3 > -1.5) { //52 tetra
					face.getElement(i).setGranica('h');
				}
				else {
					face.getElement(i).setGranica('g');
				}


			}
		}
	}


	ofstream zapis(nazwa);
	zapis.precision(12);


	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << (points.getElement(i).getX() - transX)*dzielnik << " " << pause << (points.getElement(i).getY() - transY)*dzielnik << " " << pause << (points.getElement(i).getZ() - transZ)*dzielnik * dzielnikZ << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	zapis << endl << endl;



	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = war[2];
				break;
			case 'd':
				rodzaj = war[3];
				break;
			case 'e':
				rodzaj = war[4];
				break;
			case 'f':
				rodzaj = war[5];
				break;
			case 'g':
				rodzaj = war[0];
				break;
			case 'h':
				rodzaj = war[1];
				break;
			case 'b':
				rodzaj = 9;
				break;
		}

		if (rodzaj) {
			zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;
		}

	}


	zapis.close();
}

void Delanouy3D::zapiszDoPlikuZageszczeniePunktowNaPow(bool zastosujGlobPrzesuniecia, int ileRazy, const char * nazwa) {

	double skal = 1, trX = 0, trY = 0, trZ = 0;

	if (zastosujGlobPrzesuniecia) {

		skal = 1 / globalSkal;
		trX = globalTranX;
		trY = globalTranY;
		trZ = globalTranZ;


	}
	IntList scianyPryzm(100000, 100000);
	IntList scianyPryzm1(100000, 100000);
	IntList scianyPryzm2(10, 10);

	wyszukajScianyNaPryzmyNormalna(scianyPryzm);
	//m1->Lines->Add(scianyPryzm.getIter());
	PunktList dodPunkty(10, 10);


	for (int i = 0, p1, p2, ileE = scianyPryzm.getIter(); i < ileE; i += 2) {

		p1 = scianyPryzm.getElement(i);
		p2 = scianyPryzm.getElement(i + 1);


		switch (p2) {

			case 1:
				scianyPryzm1.setElement(elements.getElement(p1).getP4());
				scianyPryzm1.setElement(elements.getElement(p1).getP2());
				scianyPryzm1.setElement(elements.getElement(p1).getP1());

				break;
			case 2:
				scianyPryzm1.setElement(elements.getElement(p1).getP4());
				scianyPryzm1.setElement(elements.getElement(p1).getP3());
				scianyPryzm1.setElement(elements.getElement(p1).getP2());

				break;
			case 3:
				scianyPryzm1.setElement(elements.getElement(p1).getP4());
				scianyPryzm1.setElement(elements.getElement(p1).getP1());
				scianyPryzm1.setElement(elements.getElement(p1).getP3());

				break;
			case 4:
				scianyPryzm1.setElement(elements.getElement(p1).getP1());
				scianyPryzm1.setElement(elements.getElement(p1).getP2());
				scianyPryzm1.setElement(elements.getElement(p1).getP3());

				break;

		}

	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3, wx, wy, wz;
	int ilePStart = points.getIter();


	int zlicz = 0;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'g') {
			zlicz++;
		}
	}

	dodPunkty.czysc(zlicz + 10, zlicz + 10);

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'g') {
			dodPunkty.setElement(points.getElement(i));
		}
	}



	for (int j = 0; j < ileRazy; ++j) {

		scianyPryzm2.czysc(scianyPryzm1.getIter() * (3 / 2) + 100, scianyPryzm1.getIter() * (3 / 2) + 100);

		for (int i = 0, p1, p2, p3, ileE = scianyPryzm1.getIter(); i < ileE; i += 3) {

			p1 = scianyPryzm1.getElement(i);
			p2 = scianyPryzm1.getElement(i + 1);
			p3 = scianyPryzm1.getElement(i + 2);

			x1 = points.getElement(p1).getX();
			x2 = points.getElement(p2).getX();
			x3 = points.getElement(p3).getX();

			y1 = points.getElement(p1).getY();
			y2 = points.getElement(p2).getY();
			y3 = points.getElement(p3).getY();

			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();

			wx = (x1 + x2 + x3) / 3;
			wy = (y1 + y2 + y3) / 3;
			wz = (z1 + z2 + z3) / 3;

			points.setElement(wx, wy, wz);

			scianyPryzm2.setElement(points.getIter() - 1);
			scianyPryzm2.setElement(p1);
			scianyPryzm2.setElement(p2);

			scianyPryzm2.setElement(points.getIter() - 1);
			scianyPryzm2.setElement(p2);
			scianyPryzm2.setElement(p3);

			scianyPryzm2.setElement(points.getIter() - 1);
			scianyPryzm2.setElement(p1);
			scianyPryzm2.setElement(p3);
		}

		scianyPryzm1.zamienIntListy(scianyPryzm2);

	}


	ofstream zapis(nazwa);
	zapis.precision(10);

	zapis << (points.getIter() - ilePStart) + dodPunkty.getIter() << endl << endl;

	if (!zastosujGlobPrzesuniecia) {

		int j = 0;
		for (int i = 0, ileP = dodPunkty.getIter(); i < ileP; ++i, ++j) {
			zapis << j + 1 << "  " << dodPunkty.getElement(i).getX() << "  " << dodPunkty.getElement(i).getY() << "  " << dodPunkty.getElement(i).getZ() << endl;
		}

		for (int i = ilePStart, ileP = points.getIter(); i < ileP; ++i, ++j) {
			zapis << j + 1 << "  " << points.getElement(i).getX() << "  " << points.getElement(i).getY() << "  " << points.getElement(i).getZ() << endl;
		}
	}
	else {

		int j = 0;
		for (int i = 0, ileP = dodPunkty.getIter(); i < ileP; ++i, ++j) {
			zapis << j + 1 << "  " << (dodPunkty.getElement(i).getX() - trX)*skal << "  " << (dodPunkty.getElement(i).getY() - trY)*skal << "  " << (dodPunkty.getElement(i).getZ() - trZ)*skal << endl;
		}

		for (int i = ilePStart, ileP = points.getIter(); i < ileP; ++i, ++j) {
			zapis << j + 1 << "  " << (points.getElement(i).getX() - trX)*skal << "  " << (points.getElement(i).getY() - trY)*skal << "  " << (points.getElement(i).getZ() - trZ)*skal << endl;
		}

	}

	zapis.close();

//m1->Lines->Add("Punktow na Gr.");

//m1->Lines->Add((points.getIter()-ilePStart)+dodPunkty.getIter());

	points.ustawIter(ilePStart);

}

void Delanouy3D::zapiszPunktyGraniczne() {

	int zlicz = 0;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'g') {
			zlicz++;
		}
	}

	wzorGranicy.czysc(zlicz + 10, zlicz + 10);

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'g') {
			wzorGranicy.setElement(i);
		}
	}
}

void Delanouy3D::wczytajPunktyGraniczne() {

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		points.getElement(i).setGranica('a');
	}

	for (int i = 0, g, ileP = wzorGranicy.getIter(); i < ileP; ++i) {
		g = wzorGranicy.getElement(i);
		points.getElement(g).setGranica('g');
	}

}

double Delanouy3D::powSciany(int p1, int p2, int p3) {

	double x1, x2, x3, y1, y2, y3, z1, z2, z3, a, b, c, p;

	x1 = points.getElement(p1).getX();
	x2 = points.getElement(p2).getX();
	x3 = points.getElement(p3).getX();

	y1 = points.getElement(p1).getY();
	y2 = points.getElement(p2).getY();
	y3 = points.getElement(p3).getY();

	z1 = points.getElement(p1).getZ();
	z2 = points.getElement(p2).getZ();
	z3 = points.getElement(p3).getZ();

	a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
	b = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
	c = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));

	p = 0.5 * (a + b + c);

	return sqrt(p * (p - a) * (p - b) * (p - c));

}


void Delanouy3D::oznaczWszytkieWezlyGraniczneNaBrzegach() {

	double x, y, z;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		x = points.getElement(i).getX();
		y = points.getElement(i).getY();
		z = points.getElement(i).getZ();

		if (( x == 50 || x == 0 ||
						y == 50 || y == 0 ||
						z == 50 || z == 0    )) {

			points.getElement(i).setGranica('g');

		}
	}
}



void Delanouy3D::tworzeniePryzm(double gruboscWarstwyPryzmatycznej, int ileWarstw, int procentowaWielkoscW, double wagaWygladzania, double doklOptyWygladzania) {

	zapiszPunktyGraniczne();

	double grubosc = gruboscWarstwyPryzmatycznej / ileWarstw;
	pryzmy.setGrow(ileScianyNaPryzmy() * (ileWarstw + 1));

	for (int i = 0; i < ileWarstw; ++i) {

		tworzeniePryzmWygladzanie(grubosc, procentowaWielkoscW);
		sasiedniePunkty();

		for (int j = 0; j < 10; ++j) {
			wygladzanieLaplaceWaga(2, wagaWygladzania, 1);
		}
		//m1->Lines->Add(grubosc*(i+1));
		sasiednieElementy();
	}

	sasiednieElementy();
	for (int j = 0; j < 10; ++j) {
		optymalizacjaStosunekRdor(doklOptyWygladzania);
	}

	wczytajPunktyGraniczne();
}

void Delanouy3D::tworzeniePryzmWygladzanie(double grubosc, int procentowaWielkoscW) {

	int ileScian = ileScianyNaPryzmy();
	//int ilePryzm = ileScian;

	int licznikBezPomniejszania = 0;
	int licznikZPomniejszaniem = 0;


	IntList scianyPryzm(100, ileScian + 10);
	IntList zamianaPunktowWTetra(100, ileScian + 10);
	PunktList wektoryNormalne;

	wyszukajScianyNaPryzmyNormalna(scianyPryzm);


	//pryzmy.czysc(100,ilePryzm+10);

	/*
	if(ilePunktowCzworosciany){points.ustawIter(ilePunktowCzworosciany);}
	*/

	//ilePunktowCzworosciany=points.getIter();


	int ilePun = points.getIter();
	int * tabPunktow = new int[ilePun];
	double dlPryzmy = grubosc;



	for (int i = 0; i < ilePun; ++i) {
		tabPunktow[i] = -1;
	}

	wyszukajNormalnePowierzchniDlaPunktowGranicznych();


	double x, y, z, x1, y1, z1, dl;

	sasiednieElementy();



	for (int i = 0, ileS = scianyPryzm.getIter(), e, nrS, p1, p2, p3; i < ileS; i += 2) {

		dl = dlPryzmy;
		e = scianyPryzm.getElement(i);
		nrS = scianyPryzm.getElement(i + 1);

		switch (nrS) {

			case 1:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP1();
				elements.getElement(e).setE1(-2);
				break;

			case 2:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP3();
				p3 = elements.getElement(e).getP2();
				elements.getElement(e).setE2(-2);
				break;

			case 3:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP1();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE3(-2);
				break;

			case 4:
				p1 = elements.getElement(e).getP1();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE4(-2);
				break;

			default:
				break;//m1->Lines->Add("zadna sciana");
		}

		if (tabPunktow[p1] <= -1) {
			x = points.getElement(p1).getX();
			y = points.getElement(p1).getY();
			z = points.getElement(p1).getZ();

			x1 = wektoryNormalne.getElement(p1).getX() * 0.01;
			y1 = wektoryNormalne.getElement(p1).getY() * 0.01;
			z1 = wektoryNormalne.getElement(p1).getZ() * 0.01;

			tabPunktow[p1] = points.getIter();


			points.setElement(x + x1, y + y1, z + z1);



			if (dlPryzmy > 0) {
				//points.zamianaMiejsc(p1,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p1);
				zamianaPunktowWTetra.setElement(tabPunktow[p1]);
				points.getElement(points.getIter() - 1).setGranica('g');
			}
		}
		if (tabPunktow[p2] <= -1) {
			x = points.getElement(p2).getX();
			y = points.getElement(p2).getY();
			z = points.getElement(p2).getZ();

			x1 = wektoryNormalne.getElement(p2).getX() * 0.01;
			y1 = wektoryNormalne.getElement(p2).getY() * 0.01;
			z1 = wektoryNormalne.getElement(p2).getZ() * 0.01;

			tabPunktow[p2] = points.getIter();


			points.setElement(x + x1, y + y1, z + z1);

			if (dlPryzmy > 0) {
				//points.zamianaMiejsc(p2,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p2);
				zamianaPunktowWTetra.setElement(tabPunktow[p2]);
				points.getElement(points.getIter() - 1).setGranica('g');
			}

		}
		if (tabPunktow[p3] <= -1) {

			x = points.getElement(p3).getX();
			y = points.getElement(p3).getY();
			z = points.getElement(p3).getZ();

			x1 = wektoryNormalne.getElement(p3).getX() * 0.01;
			y1 = wektoryNormalne.getElement(p3).getY() * 0.01;
			z1 = wektoryNormalne.getElement(p3).getZ() * 0.01;

			tabPunktow[p3] = points.getIter();


			points.setElement(x + x1, y + y1, z + z1);



			if (dlPryzmy > 0) {
				//points.zamianaMiejsc(p3,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p3);
				zamianaPunktowWTetra.setElement(tabPunktow[p3]);
				points.getElement(points.getIter() - 1).setGranica('g');
			}

		}



		pryzmy.setElement(p1, p3, p2, tabPunktow[p1], tabPunktow[p3], tabPunktow[p2]);


	}

	for (int i = 0, pS, pN, ileP = zamianaPunktowWTetra.getIter(); i < ileP; i += 2) {
		pS = zamianaPunktowWTetra.getElement(i);
		pN = zamianaPunktowWTetra.getElement(i + 1);
		elements.podmienPunktWElementach(pS, pN, laplas[pS]);

	}

	elements.ustawSasiadow(-2, -1);


	IntList doLaplasa(zamianaPunktowWTetra.getIter() * 0.5 + 2, zamianaPunktowWTetra.getIter() * 0.5 + 2);

	for (int i = 0, ileP = zamianaPunktowWTetra.getIter(); i < ileP; i += 2) {
		doLaplasa.setElement(zamianaPunktowWTetra.getElement(i + 1));
	}


	sasiedniePunktyPryzm();

	for (int ii = 0; ii < 6; ++ii) {
		wygladzanieLaplaceWybranePunkty(doLaplasa);
	}
	//for(int ii=0;ii<6;++ii){wygladzanieLaplaceWagaWybranePunkty(2,0.1,doLaplasa,1);}

	for (int i = 0, pS, pN, ileP = zamianaPunktowWTetra.getIter(); i < ileP; i += 2) {

		pS = zamianaPunktowWTetra.getElement(i);
		pN = zamianaPunktowWTetra.getElement(i + 1);

		//wyliczamy wektor
		x = points.getElement(pS).getX();
		y = points.getElement(pS).getY();
		z = points.getElement(pS).getZ();

		x1 = points.getElement(pN).getX() - x;
		y1 = points.getElement(pN).getY() - y;
		z1 = points.getElement(pN).getZ() - z;

		//wektor normalny
		double temp = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
		x1 /= temp;
		y1 /= temp;
		z1 /= temp;

		temp *= procentowaWielkoscW;

		if (temp < dlPryzmy) {
			points.getElement(pN).setPunkt(points.getElement(pS).getX() + x1 * temp, points.getElement(pS).getY() + y1 * temp, points.getElement(pS).getZ() + z1 * temp, 'g');
			++licznikZPomniejszaniem;
		}
		else {

			points.getElement(pN).setPunkt(points.getElement(pS).getX() + x1 * dlPryzmy, points.getElement(pS).getY() + y1 * dlPryzmy, points.getElement(pS).getZ() + z1 * dlPryzmy, 'g');
			++licznikBezPomniejszania;
		}

	}
	/*
	m1->Lines->Add("Pomniejszaj nieP");
	m1->Lines->Add(licznikZPomniejszaniem);
	m1->Lines->Add(licznikBezPomniejszania);
	*/


}


/*
void Delanouy3D::tworzeniePryzmNormalna(double grubosc,int ileWarstw,bool flagaJakosci){

	int ileScian = ileScianyNaPryzmy();
	int ilePryzm = ileScian+ileWarstw;

	IntList scianyPryzm(100,ileScian+10);
	IntList zamianaPunktowWTetra(100,ileScian+10);
	PunktList wektoryNormalne;

	wyszukajScianyNaPryzmyNormalna(scianyPryzm);

	pryzmy.czysc(100,ilePryzm+10);

	if(ilePunktowCzworosciany){points.ustawIter(ilePunktowCzworosciany);}
	ilePunktowCzworosciany=points.getIter();

	int ilePun=points.getIter();
	int *tabPunktow = new int[ilePun];
	double dlPryzmy = grubosc / ileWarstw;

	for(int i=0;i<ilePun;++i){tabPunktow[i]=-1;}

	wyszukajNormalnePowierzchniDlaPunktowGranicznych(wektoryNormalne,flagaJakosci);

	double x,y,z,x1,y1,z1,dl;

	sasiednieElementy();

	for(int i=0,ileS=scianyPryzm.getIter(),e,nrS,p1,p2,p3;i<ileS;i+=2){

		dl = dlPryzmy;
		e = scianyPryzm.getElement(i);
		nrS = scianyPryzm.getElement(i+1);

		switch(nrS){

			case 1:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP1();
				elements.getElement(e).setE1(-2);
				break;

			case 2:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP3();
				p3 = elements.getElement(e).getP2();
				elements.getElement(e).setE2(-2);
				break;

			case 3:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP1();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE3(-2);
				break;

			case 4:
				p1 = elements.getElement(e).getP1();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE4(-2);
				break;

			default: m1->Lines->Add("zadna sciana");
		}

		if(tabPunktow[p1]==-1){
			x = points.getElement(p1).getX();
			y = points.getElement(p1).getY();
			z = points.getElement(p1).getZ();

			x1 = wektoryNormalne.getElement(p1).getX()*dlPryzmy;
			y1 = wektoryNormalne.getElement(p1).getY()*dlPryzmy;
			z1 = wektoryNormalne.getElement(p1).getZ()*dlPryzmy;

			tabPunktow[p1]=points.getIter();

			for(int w=0;w<ileWarstw;++w){
				points.setElement(x+x1*(w+1),y+y1*(w+1),z+z1*(w+1));

			}

			if(dlPryzmy>0){
				//points.zamianaMiejsc(p1,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p1);zamianaPunktowWTetra.setElement(points.getIter()-1);
				points.getElement(points.getIter()-1).setGranica('g');
			}
		}
		if(tabPunktow[p2]==-1){
			x = points.getElement(p2).getX();
			y = points.getElement(p2).getY();
			z = points.getElement(p2).getZ();

			x1 = wektoryNormalne.getElement(p2).getX()*dlPryzmy;
			y1 = wektoryNormalne.getElement(p2).getY()*dlPryzmy;
			z1 = wektoryNormalne.getElement(p2).getZ()*dlPryzmy;

			tabPunktow[p2]=points.getIter();

			for(int w=0;w<ileWarstw;++w){
				points.setElement(x+x1*(w+1),y+y1*(w+1),z+z1*(w+1));
			}
			if(dlPryzmy>0){
				//points.zamianaMiejsc(p2,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p2);zamianaPunktowWTetra.setElement(points.getIter()-1);
				points.getElement(points.getIter()-1).setGranica('g');
			}

		}
		if(tabPunktow[p3]==-1){

			x = points.getElement(p3).getX();
			y = points.getElement(p3).getY();
			z = points.getElement(p3).getZ();

			x1 = wektoryNormalne.getElement(p3).getX()*dlPryzmy;
			y1 = wektoryNormalne.getElement(p3).getY()*dlPryzmy;
			z1 = wektoryNormalne.getElement(p3).getZ()*dlPryzmy;

			tabPunktow[p3]=points.getIter();

			for(int w=0;w<ileWarstw;++w){
				points.setElement(x+x1*(w+1),y+y1*(w+1),z+z1*(w+1));
			}


			if(dlPryzmy>0){
				//points.zamianaMiejsc(p3,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p3);zamianaPunktowWTetra.setElement(points.getIter()-1);
				points.getElement(points.getIter()-1).setGranica('g');
			}

		}

		for(int w=0;w<ileWarstw;++w){

			if(w==0){pryzmy.setElement(p1,p3,p2,tabPunktow[p1],tabPunktow[p3],tabPunktow[p2]);}
			else{pryzmy.setElement(tabPunktow[p1]+w-1,tabPunktow[p3]+w-1,tabPunktow[p2]+w-1,tabPunktow[p1]+w,tabPunktow[p3]+w,tabPunktow[p2]+w);}
			//else{pryzmy.setElement(tabPunktow[p1]+w-1,tabPunktow[p2]+w-1,tabPunktow[p3]+w-1,tabPunktow[p1]+w,tabPunktow[p2]+w,tabPunktow[p3]+w);}

		}

	}

	for(int i=0,pS,pN,ileP=zamianaPunktowWTetra.getIter();i<ileP;i+=2){
		pS=zamianaPunktowWTetra.getElement(i);
		pN=zamianaPunktowWTetra.getElement(i+1);
		elements.podmienPunktWElementach(pS,pN,laplas[pS]);

	}


	IntList doLaplasa(zamianaPunktowWTetra.getIter()*0.5+2,zamianaPunktowWTetra.getIter()*0.5+2);

	for(int i=0,pS,pN,ileP=zamianaPunktowWTetra.getIter();i<ileP;i+=2){
		doLaplasa.setElement(zamianaPunktowWTetra.getElement(i+1));
	}


	sasiedniePunktyPryzm();

	wygladzanieLaplaceWybranePunkty(doLaplasa);
	wygladzanieLaplaceWybranePunkty(doLaplasa);
	wygladzanieLaplaceWybranePunkty(doLaplasa);
	wygladzanieLaplaceWybranePunkty(doLaplasa);
	wygladzanieLaplaceWybranePunkty(doLaplasa);
	wygladzanieLaplaceWybranePunkty(doLaplasa);



	for(int i=0,ileP=doLaplasa.getIter(),p;i<ileP;++i){
		p = doLaplasa.getElement(i);
		points.getElement(p).setGranica('g');

	}





	//wygladzanieLaplaceWOparciuOWybranePunktu(doLaplasa,10);
	//wygladzanieLaplaceWagaWOparciuOWybranePunktu(doLaplasa,2,0.3,1,10);

}



*/


void Delanouy3D::tworzeniePryzmNormalna(double grubosc, int ileWarstw, bool flagaJakosci) {

	int ileScian = ileScianyNaPryzmy();
	int ilePryzm = ileScian + ileWarstw;

	IntList scianyPryzm(100, ileScian + 10);
	IntList zamianaPunktowWTetra(100, ileScian + 10);
	PunktList wektoryNormalne;
	wyszukajScianyNaPryzmyNormalna(scianyPryzm);

	pryzmy.czysc(100, ilePryzm + 10);

	if (ilePunktowCzworosciany) {
		points.ustawIter(ilePunktowCzworosciany);
	}
	ilePunktowCzworosciany = points.getIter();

	int ilePun = points.getIter();
	int * tabPunktow = new int[ilePun];
	double dlPryzmy = grubosc / ileWarstw;

	for (int i = 0; i < ilePun; ++i) {
		tabPunktow[i] = -1;
	}

	wyszukajNormalnePowierzchniDlaPunktowGranicznych();

	double x, y, z, x1, y1, z1, dl;

	sasiednieElementy();

	for (int i = 0, ileS = scianyPryzm.getIter(), e, nrS, p1, p2, p3; i < ileS; i += 2) {

		dl = dlPryzmy;
		e = scianyPryzm.getElement(i);
		nrS = scianyPryzm.getElement(i + 1);

		switch (nrS) {

			case 1:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP1();
				elements.getElement(e).setE1(-2);
				break;

			case 2:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP3();
				p3 = elements.getElement(e).getP2();
				elements.getElement(e).setE2(-2);
				break;

			case 3:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP1();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE3(-2);
				break;

			case 4:
				p1 = elements.getElement(e).getP1();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE4(-2);
				break;

			default:
				break;//m1->Lines->Add("zadna sciana");
		}

		if (tabPunktow[p1] <= -1) {
			x = points.getElement(p1).getX();
			y = points.getElement(p1).getY();
			z = points.getElement(p1).getZ();

			x1 = wektoryNormalne.getElement(p1).getX() * dlPryzmy;
			y1 = wektoryNormalne.getElement(p1).getY() * dlPryzmy;
			z1 = wektoryNormalne.getElement(p1).getZ() * dlPryzmy;

			tabPunktow[p1] = points.getIter();

			for (int w = 0; w < ileWarstw; ++w) {
				points.setElement(x + x1 * (w + 1), y + y1 * (w + 1), z + z1 * (w + 1));

			}

			if (dlPryzmy > 0) {
				//points.zamianaMiejsc(p1,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p1);
				zamianaPunktowWTetra.setElement(points.getIter() - 1);
				points.getElement(points.getIter() - 1).setGranica('g');
			}
		}
		if (tabPunktow[p2] == -1) {
			x = points.getElement(p2).getX();
			y = points.getElement(p2).getY();
			z = points.getElement(p2).getZ();

			x1 = wektoryNormalne.getElement(p2).getX() * dlPryzmy;
			y1 = wektoryNormalne.getElement(p2).getY() * dlPryzmy;
			z1 = wektoryNormalne.getElement(p2).getZ() * dlPryzmy;

			tabPunktow[p2] = points.getIter();

			for (int w = 0; w < ileWarstw; ++w) {
				points.setElement(x + x1 * (w + 1), y + y1 * (w + 1), z + z1 * (w + 1));
			}
			if (dlPryzmy > 0) {
				//points.zamianaMiejsc(p2,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p2);
				zamianaPunktowWTetra.setElement(points.getIter() - 1);
				points.getElement(points.getIter() - 1).setGranica('g');
			}

		}
		if (tabPunktow[p3] <= -1) {

			x = points.getElement(p3).getX();
			y = points.getElement(p3).getY();
			z = points.getElement(p3).getZ();

			x1 = wektoryNormalne.getElement(p3).getX() * dlPryzmy;
			y1 = wektoryNormalne.getElement(p3).getY() * dlPryzmy;
			z1 = wektoryNormalne.getElement(p3).getZ() * dlPryzmy;

			tabPunktow[p3] = points.getIter();

			for (int w = 0; w < ileWarstw; ++w) {
				points.setElement(x + x1 * (w + 1), y + y1 * (w + 1), z + z1 * (w + 1));
			}


			if (dlPryzmy > 0) {
				//points.zamianaMiejsc(p3,points.getIter()-1);
				zamianaPunktowWTetra.setElement(p3);
				zamianaPunktowWTetra.setElement(points.getIter() - 1);
				points.getElement(points.getIter() - 1).setGranica('g');
			}



		}

		for (int w = 0; w < ileWarstw; ++w) {

			if (w == 0) {
				pryzmy.setElement(p1, p3, p2, tabPunktow[p1], tabPunktow[p3], tabPunktow[p2]);
			}
			else {
				pryzmy.setElement(tabPunktow[p1] + w - 1, tabPunktow[p3] + w - 1, tabPunktow[p2] + w - 1, tabPunktow[p1] + w, tabPunktow[p3] + w, tabPunktow[p2] + w);
			}
			//else{pryzmy.setElement(tabPunktow[p1]+w-1,tabPunktow[p2]+w-1,tabPunktow[p3]+w-1,tabPunktow[p1]+w,tabPunktow[p2]+w,tabPunktow[p3]+w);}

		}

	}

	for (int i = 0, pS, pN, ileP = zamianaPunktowWTetra.getIter(); i < ileP; i += 2) {
		pS = zamianaPunktowWTetra.getElement(i);
		pN = zamianaPunktowWTetra.getElement(i + 1);
		elements.podmienPunktWElementach(pS, pN, laplas[pS]);

	}


}

/*
int Delanouy3D::ileScianyNaPryzmy(){

	int ileScian=0;

	for(int i=0,ileE=elements.getIter();i<ileE;++i){

		if(elements.getElement(i).getE1()==-1){


		elements.getElement(i).get
		++ileScian;

		}
		if(elements.getElement(i).getE2()==-1){++ileScian;}
		if(elements.getElement(i).getE3()==-1){++ileScian;}
		if(elements.getElement(i).getE4()==-1){++ileScian;}

	}

	return ileScian;
}
*/


int Delanouy3D::ileScianyNaPryzmy() {

	int ileScian = 0;
	int p1, p2, p3, p4;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		x1 = points.getElement(p1).getX();
		y1 = points.getElement(p1).getY();
		z1 = points.getElement(p1).getZ();

		x2 = points.getElement(p2).getX();
		y2 = points.getElement(p2).getY();
		z2 = points.getElement(p2).getZ();

		x3 = points.getElement(p3).getX();
		y3 = points.getElement(p3).getY();
		z3 = points.getElement(p3).getZ();

		x4 = points.getElement(p4).getX();
		y4 = points.getElement(p4).getY();
		z4 = points.getElement(p4).getZ();

		if (elements.getElement(i).getE1() <= -1) {

			points.getElement(p4).setGranica('g');
			points.getElement(p1).setGranica('g');
			points.getElement(p2).setGranica('g');

			++ileScian;
		}

		if (elements.getElement(i).getE2() <= -1) {



			points.getElement(p4).setGranica('g');
			points.getElement(p3).setGranica('g');
			points.getElement(p2).setGranica('g');

			++ileScian;
		}

		if (elements.getElement(i).getE3() <= -1) {



			points.getElement(p4).setGranica('g');
			points.getElement(p1).setGranica('g');
			points.getElement(p3).setGranica('g');

			++ileScian;
		}
		if (elements.getElement(i).getE4() <= -1) {


			points.getElement(p1).setGranica('g');
			points.getElement(p2).setGranica('g');
			points.getElement(p3).setGranica('g');

			++ileScian;
		}

	}

	return ileScian;
}



void Delanouy3D::wyszukajScianyNaPryzmyNormalna(IntList & scianyPryzm) {

	int ileScian = 0;

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		if (elements.getElement(i).getE1() <= -1) {

			scianyPryzm.setElement(i);
			scianyPryzm.setElement(1);
			++ileScian;

		}
		if (elements.getElement(i).getE2() <= -1) {

			scianyPryzm.setElement(i);
			scianyPryzm.setElement(2);
			++ileScian;

		}
		if (elements.getElement(i).getE3() <= -1) {

			scianyPryzm.setElement(i);
			scianyPryzm.setElement(3);
			++ileScian;

		}
		if (elements.getElement(i).getE4() <= -1) {

			scianyPryzm.setElement(i);
			scianyPryzm.setElement(4);
			++ileScian;

		}

	}


}

//void Delanouy3D::wyszukajNormalnePowierzchniDlaPunktowGranicznych(PunktList &wektoryNormalne,bool flagaJakosci){
void Delanouy3D::wyszukajNormalnePowierzchniDlaPunktowGranicznych() {
	int ilePG = 0;
	PunktList sprawdz(100, 100);
//wektoryNormalne.czysc(10,points.getIter()+10);
	IntList elementySas;
	elementySas.czysc(100, 100);
	DoubleList polaEleSas;
	polaEleSas.czysc(50, 50);

	sasiednieElementy();
	Punkt tempG, tempP;
	tempP.setPunkt(-10000.0, -10000.0, -10000.0, 'a');

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() == 'g') {

			wyszukajWlementyScianGranicznychZPunktem(i, elementySas, polaEleSas);

			wyliczWektorNormalnyDlaPunkty(tempG, i, elementySas, polaEleSas, sprawdz, true);

			points.getElement(i).setVecNorm(tempG.getX(), tempG.getY(), tempG.getZ());
			//wektoryNormalne.setElement(tempG);
			++ilePG;

			points.getElement(i).setGranica('g');
		}
		else {

			//points.getElement(i).setVecNorm(tempP.getX(),tempP.getY(),tempP.getZ());
			points.getElement(i).setGranica('a');
			//wektoryNormalne.setElement(tempP);

		}

	}

}

void Delanouy3D::wyszukajWlementyScianGranicznychZPunktem(int numerWezla, IntList & elementySas, DoubleList & polaEleSas) {

	elementySas.ustawIter(0);


	for ( int i = 0, e, ileE = laplas[numerWezla].getIter() ; i < ileE ; ++i ) {
		e =  laplas[numerWezla].getElement(i);

		if (elements.getElement(e).getE1() <= -1 && elements.getElement(e).getP3() != numerWezla) {
			elementySas.setElement(e);
			elementySas.setElement(1);
		}
		if (elements.getElement(e).getE2() <= -1 && elements.getElement(e).getP1() != numerWezla) {
			elementySas.setElement(e);
			elementySas.setElement(2);
		}
		if (elements.getElement(e).getE3() <= -1 && elements.getElement(e).getP2() != numerWezla) {
			elementySas.setElement(e);
			elementySas.setElement(3);
		}
		if (elements.getElement(e).getE4() <= -1 && elements.getElement(e).getP4() != numerWezla) {
			elementySas.setElement(e);
			elementySas.setElement(4);
		}

	}


	polaEleSas.ustawIter(0);

	for (int i = 0, e, p1 = -1, p2 = -1, p3 = -1, ileE = elementySas.getIter(); i < ileE; i += 2) {
		e = elementySas.getElement(i);

		switch (elementySas.getElement(i + 1)) {

			case 1:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP1();
				break;

			case 2:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP3();
				p3 = elements.getElement(e).getP2();
				break;

			case 3:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP1();
				p3 = elements.getElement(e).getP3();
				break;

			case 4:
				p1 = elements.getElement(e).getP1();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP3();
				break;

		}


		polaEleSas.setElement(powSciany(p1, p2, p3));

	}


}

void Delanouy3D::wyliczWektorNormalnyDlaPunkty(Punkt & szukWektor, int numerWezla, IntList & elementySas, DoubleList & polaEleSas, PunktList & sprawdz, bool flagaJakosci) {


	double sumaPow = polaEleSas.getSuma();
	double sredPow = polaEleSas.getSrednia();
	double ratPow = 1.0;

	double x, x1 = 0, x2 = 0, y, y1 = 0, y2 = 0, z, z1 = 0, z2 = 0, tempx1, tempy1, tempz1, szukX = 0, szukY = 0, szukZ = 0;
	sprawdz.ustawIter(0);

	x = points.getElement(numerWezla).getX();
	y = points.getElement(numerWezla).getY();
	z = points.getElement(numerWezla).getZ();


	for (int i = 0, ii = 0, e, p1 = -1, p2 = -1, p3 = -1, ileE = elementySas.getIter(); i < ileE; i += 2, ++ii) {
		e = elementySas.getElement(i);

		switch (elementySas.getElement(i + 1)) {

			case 1:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP1();
				break;

			case 2:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP3();
				p3 = elements.getElement(e).getP2();
				break;

			case 3:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP1();
				p3 = elements.getElement(e).getP3();
				break;

			case 4:
				p1 = elements.getElement(e).getP1();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP3();
				break;

		}

		if (p1 == numerWezla) {

			x1 = points.getElement(p2).getX() - x;
			y1 = points.getElement(p2).getY() - y;
			z1 = points.getElement(p2).getZ() - z;

			x2 = points.getElement(p3).getX() - x;
			y2 = points.getElement(p3).getY() - y;
			z2 = points.getElement(p3).getZ() - z;
		}
		else if (p2 == numerWezla) {

			x1 = points.getElement(p3).getX() - x;
			y1 = points.getElement(p3).getY() - y;
			z1 = points.getElement(p3).getZ() - z;

			x2 = points.getElement(p1).getX() - x;
			y2 = points.getElement(p1).getY() - y;
			z2 = points.getElement(p1).getZ() - z;
		}
		else if (p3 == numerWezla) {

			x1 = points.getElement(p1).getX() - x;
			y1 = points.getElement(p1).getY() - y;
			z1 = points.getElement(p1).getZ() - z;

			x2 = points.getElement(p2).getX() - x;
			y2 = points.getElement(p2).getY() - y;
			z2 = points.getElement(p2).getZ() - z;
		}
		else {;}
		//wektor normalny

		tempx1 = y1 * z2 - z1 * y2;
		tempy1 = z1 * x2 - x1 * z2;
		tempz1 = x1 * y2 - y1 * x2;

		wyliczWektorNormalny(tempx1, tempy1, tempz1);

		/*
		tempx1=elements.getElement(e).getSx()-x;
		tempy1=elements.getElement(e).getSy()-y;
		tempz1=elements.getElement(e).getSz()-z;

		wyliczWektorNormalny(tempx1,tempy1,tempz1);
		*/


		if (sprawdz.setElementUniqat(tempx1, tempy1, tempz1, 'g')) {


			//ratPow = polaEleSas.getElement(ii)/sumaPow;
			//ratPow = polaEleSas.getElement(ii)/sredPow;
			//m1->Lines->Add(ratPow);
			//ratPow = polaEleSas.getElement(ii);
			ratPow = 1.0;
			szukX = szukX + tempx1 * ratPow;
			szukY = szukY + tempy1 * ratPow;
			szukZ = szukZ + tempz1 * ratPow;

		}




	}
	//m1->Lines->Add("----------");

	wyliczWektorNormalny(szukX, szukY, szukZ);


	//double av=1,bv=2,new_szukZ=(-(szukX+2*szukY))/szukZ;
	//szukX=1;szukY=2;szukZ=new_szukZ;

	//wyliczWektorNormalny(szukX,szukY,szukZ);


	szukWektor.setPunkt(-szukX, -szukY, -szukZ, 'g');

}

void Delanouy3D::wyliczWektorNormalny(double & x1_zwWyn, double & y1_zwWyn, double & z1_zwWyn) {

	double temp = sqrt(x1_zwWyn * x1_zwWyn + y1_zwWyn * y1_zwWyn + z1_zwWyn * z1_zwWyn);
	x1_zwWyn /= temp;
	y1_zwWyn /= temp;
	z1_zwWyn /= temp;
}

void Delanouy3D::z() {

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
		//m1->Lines->Add(elements.getElement(i).getRodzajZiarna());
	}

}


int Delanouy3D::przygotojTabliceDoDelanoya(int ** & tabElementow2, int WETE2, double WEws2) {

	int TE2 = WETE2;

	double ws2 = WEws2;

	int wymiarTE2X = (int)(dx * ws2) + 2 * TE2;
	int wymiarTE2Y = (int)(dy * ws2) + 2 * TE2;
	int wymiarTE2Z = (int)(dz * ws2) + 2 * TE2;

	tabElementow2 = new int * [wymiarTE2Z];
	for (int i = 0; i < wymiarTE2Z; i++) {
		tabElementow2[i] = new int[wymiarTE2X * wymiarTE2Y];
	}
	for (int i = 0; i < wymiarTE2Z; i++) {
		for (int j = 0; j < wymiarTE2X * wymiarTE2Y; j++) {
			tabElementow2[i][j] = -1;
		}
	}

	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, np1x, np1y, np1z;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		x1 = points.getElement(p1).getX();
		y1 = points.getElement(p1).getY();
		z1 = points.getElement(p1).getZ();

		x2 = points.getElement(p2).getX();
		y2 = points.getElement(p2).getY();
		z2 = points.getElement(p2).getZ();

		x3 = points.getElement(p3).getX();
		y3 = points.getElement(p3).getY();
		z3 = points.getElement(p3).getZ();

		x4 = points.getElement(p4).getX();
		y4 = points.getElement(p4).getY();
		z4 = points.getElement(p4).getZ();

		np1x = (x1 + x2 + x3 + x4) * 0.25;
		np1y = (y1 + y2 + y3 + y4) * 0.25;
		np1z = (z1 + z2 + z3 + z4) * 0.25;

		tabElementow2[int(np1z * ws2) + TE2][(int(np1y * ws2) + TE2)*wymiarTE2X + int(np1x * ws2) + TE2] = i;

	}

	return wymiarTE2Z;

}

void Delanouy3D::zageszczenieSiatkiNaGranicy() {

	int ** tabElementow2;
	int wymiarZusunieciaTab = przygotojTabliceDoDelanoya(tabElementow2, 3, 0.3);

	PunktList temp;
	Punkt temp1, temp2, temp3, temp4;
	double nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4, nz1, nz2, nz3, nz4, p1, p2, p3, p4;



	for (int i = 0, e, e1, e2, e3, e4, ileE = elements.getIter(); i < ileE; ++i) {

		temp.czysc(10, 10);
		e = elements.getElement(i).getRodzajZiarna();

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		nx1 = points.getElement(p1).getX();
		ny1 = points.getElement(p1).getY();
		nz1 = points.getElement(p1).getZ();

		nx2 = points.getElement(p2).getX();
		ny2 = points.getElement(p2).getY();
		nz2 = points.getElement(p2).getZ();

		nx3 = points.getElement(p3).getX();
		ny3 = points.getElement(p3).getY();
		nz3 = points.getElement(p3).getZ();

		nx4 = points.getElement(p4).getX();
		ny4 = points.getElement(p4).getY();
		nz4 = points.getElement(p4).getZ();

		if (e != -1) {



			if (e1 != -1 && elements.getElement(e1).getRodzajZiarna() != e ) {

				punktSciana3D(nx4, ny4, nz4, nx2, ny2, nz2, nx1, ny1, nz1, temp1);
				temp.setElement(temp1);

			}

			if (e2 != -1 && elements.getElement(e2).getRodzajZiarna() != e ) {

				punktSciana3D(nx4, ny4, nz4, nx3, ny3, nz3, nx2, ny2, nz2, temp2);
				temp.setElement(temp2);

			}

			if (e3 != -1 && elements.getElement(e3).getRodzajZiarna() != e ) {

				punktSciana3D(nx4, ny4, nz4, nx1, ny1, nz1, nx3, ny3, nz3, temp3);
				temp.setElement(temp3);

			}

			if (e4 != -1 && elements.getElement(e4).getRodzajZiarna() != e ) {

				punktSciana3D(nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, temp4);
				temp.setElement(temp4);

			}


		}

		if (temp.getIter()) {
			//m1->Lines->Add(i);
			delanoyTablicaOPDodajPunktyZTablicaPolozeniaEl(tabElementow2, dx, dy, dz, false, temp, 0.3, 3);
		}


	}


	for (int i = 0; i < wymiarZusunieciaTab; i++) {
		delete []tabElementow2[i];
	}
	delete []tabElementow2;


}

void Delanouy3D::zageszczenieSiatkiPowierzchni() {

	int ** tabElementow2;
	int wymiarZusunieciaTab = przygotojTabliceDoDelanoya(tabElementow2, 3, 0.3);

	PunktList temp;
	Punkt temp1, temp2, temp3, temp4;
	double nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4, nz1, nz2, nz3, nz4, p1, p2, p3, p4;


	for (int i = 0, e1, e2, e3, e4, ileE = elements.getIter(); i < ileE; ++i) {

		temp.czysc(10, 10);

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		nx1 = points.getElement(p1).getX();
		ny1 = points.getElement(p1).getY();
		nz1 = points.getElement(p1).getZ();

		nx2 = points.getElement(p2).getX();
		ny2 = points.getElement(p2).getY();
		nz2 = points.getElement(p2).getZ();

		nx3 = points.getElement(p3).getX();
		ny3 = points.getElement(p3).getY();
		nz3 = points.getElement(p3).getZ();

		nx4 = points.getElement(p4).getX();
		ny4 = points.getElement(p4).getY();
		nz4 = points.getElement(p4).getZ();

		if (e1 <= -1) {
			//punktSciana3D(nx4,ny4,nz4,nx2,ny2,nz2,nx1,ny1,nz1,temp1);temp.setElement(temp1);
			temp.setElement((nx4 + nx2 + nx1) / 3, (ny4 + ny2 + ny1) / 3, (nz4 + nz2 + nz1) / 3, 'g');
		}

		if (e2 <= -1) {
			//punktSciana3D(nx4,ny4,nz4,nx3,ny3,nz3,nx2,ny2,nz2,temp2);temp.setElement(temp2);
			temp.setElement((nx4 + nx2 + nx3) / 3, (ny4 + ny2 + ny3) / 3, (nz4 + nz2 + nz3) / 3, 'g');
		}

		if (e3 <= -1) {
			//punktSciana3D(nx4,ny4,nz4,nx1,ny1,nz1,nx3,ny3,nz3,temp3);temp.setElement(temp3);
			temp.setElement((nx4 + nx3 + nx1) / 3, (ny4 + ny3 + ny1) / 3, (nz4 + nz3 + nz1) / 3, 'g');
		}

		if (e4 <= -1) {
			//punktSciana3D(nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,temp4);temp.setElement(temp4);
			temp.setElement((nx1 + nx2 + nx3) / 3, (ny1 + ny2 + ny3) / 3, (nz1 + nz2 + nz3) / 3, 'g');
		}

		if (temp.getIter()) {
			//m1->Lines->Add(i);
			delanoyTablicaOPDodajPunktyZTablicaPolozeniaEl(tabElementow2, dx, dy, dz, false, temp, 0.3, 3);
		}


	}


	for (int i = 0; i < wymiarZusunieciaTab; i++) {
		delete []tabElementow2[i];
	}
	delete []tabElementow2;


}


int Delanouy3D::wyszukajElement(const double & x, const double & y, const double & z) {


	double X, Y, Z;
//int znalazl=-1;
	float R_od_punktu, Rr;
//int zlicz=0;

	for (int i = 0; i < elements.getIter(); i++) {

		Rr = elements.getElement(i).getR();

		X = elements.getElement(i).getSx();
		Y = elements.getElement(i).getSy();
		Z = elements.getElement(i).getSz();

		R_od_punktu = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);

		if (R_od_punktu <= Rr) {
			return i;
		}
	}

//m1->Lines->Add("nie znalazl");
//text->Lines->Add(zlicz);
	return -1;
}


/*
int Delanouy3D::rysujPryzmy(bool ujemne){

double maxX = dx;
double maxY = dy;
double maxZ = dz;
glPushMatrix();
int zlicz=0;

glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);

glBegin(GL_LINES);
glColor3f (1, 1, 1);

if(!ujemne){
	for(int i=0;i<pryzmy.getIter();i++){



		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());


		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());


		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());


	}
}
else{
	for(int i=0;i<pryzmy.getIter();i++){
		if(!V_objetoscP(i)){
		zlicz++;
		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());


		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());


		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());
		}

	}
}
glEnd(); //koniec rysowania

glPopMatrix();
glPopMatrix();


return zlicz;
}
*/

/*

void Delanouy3D::rysujPryzme(int ktora){
int i= ktora;

double maxX = dx;
double maxY = dy;
double maxZ = dz;
glPushMatrix();


glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);

glPointSize(5);
glBegin(GL_POINTS);
//glPointSize(3);


		glColor3f (1,0,0);

		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());


		glColor3f (1,1,0);

		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());

glEnd(); //koniec rysowania

glBegin(GL_LINES);
glColor3f (1, 1, 1);

		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());


		glVertex3f(points.getElement(pryzmy.getElement(i).getP1()).getX(),points.getElement(pryzmy.getElement(i).getP1()).getY(), points.getElement(pryzmy.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP2()).getX(),points.getElement(pryzmy.getElement(i).getP2()).getY(), points.getElement(pryzmy.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP3()).getX(),points.getElement(pryzmy.getElement(i).getP3()).getY(), points.getElement(pryzmy.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());


		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP5()).getX(),points.getElement(pryzmy.getElement(i).getP5()).getY(), points.getElement(pryzmy.getElement(i).getP5()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());

		glVertex3f(points.getElement(pryzmy.getElement(i).getP6()).getX(),points.getElement(pryzmy.getElement(i).getP6()).getY(), points.getElement(pryzmy.getElement(i).getP6()).getZ());
		glVertex3f(points.getElement(pryzmy.getElement(i).getP4()).getX(),points.getElement(pryzmy.getElement(i).getP4()).getY(), points.getElement(pryzmy.getElement(i).getP4()).getZ());


glEnd(); //koniec rysowania


glPopMatrix();
glPopMatrix();



}
*/

void Delanouy3D::wyszukajScianyNaPryzmy(IntList & scianyPryzm) {

	//wyszukaj sciany

	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		if (elements.getElement(i).getE1() <= -1) {

			p1 =  elements.getElement(i).getP4();
			p2 =  elements.getElement(i).getP2();
			p3 =  elements.getElement(i).getP1();
			p4 =  elements.getElement(i).getP3();

			if (points.getElement(p1).getY() == points.getElement(p2).getY() && points.getElement(p1).getY() == points.getElement(p3).getY() && points.getElement(p3).getY() != 0) {
				scianyPryzm.setElement(i);
				scianyPryzm.setElement(1);
				if (points.getElement(p1).getY() > points.getElement(p4).getY()) {
					scianyPryzm.setElement(1);
				}
				else {
					scianyPryzm.setElement(-1);
				}

			}

		}
		if (elements.getElement(i).getE2() <= -1) {

			p1 =  elements.getElement(i).getP4();
			p2 =  elements.getElement(i).getP3();
			p3 =  elements.getElement(i).getP2();
			p4 =  elements.getElement(i).getP1();

			if (points.getElement(p1).getY() == points.getElement(p2).getY() && points.getElement(p1).getY() == points.getElement(p3).getY() && points.getElement(p3).getY() != 0) {
				scianyPryzm.setElement(i);
				scianyPryzm.setElement(2);
				if (points.getElement(p1).getY() > points.getElement(p4).getY()) {
					scianyPryzm.setElement(1);
				}
				else {
					scianyPryzm.setElement(-1);
				}

			}

		}
		if (elements.getElement(i).getE3() <= -1) {

			p1 =  elements.getElement(i).getP4();
			p2 =  elements.getElement(i).getP1();
			p3 =  elements.getElement(i).getP3();
			p4 =  elements.getElement(i).getP2();

			if (points.getElement(p1).getY() == points.getElement(p2).getY() && points.getElement(p1).getY() == points.getElement(p3).getY() && points.getElement(p3).getY() != 0) {

				scianyPryzm.setElement(i);
				scianyPryzm.setElement(3);
				if (points.getElement(p1).getY() > points.getElement(p4).getY()) {
					scianyPryzm.setElement(1);
				}
				else {
					scianyPryzm.setElement(-1);
				}

			}

		}
		if (elements.getElement(i).getE4() <= -1) {

			p1 =  elements.getElement(i).getP1();
			p2 =  elements.getElement(i).getP2();
			p3 =  elements.getElement(i).getP3();
			p4 =  elements.getElement(i).getP4();

			if (points.getElement(p1).getY() == points.getElement(p2).getY() && points.getElement(p1).getY() == points.getElement(p3).getY() && points.getElement(p3).getY() != 0) {
				scianyPryzm.setElement(i);
				scianyPryzm.setElement(4);
				if (points.getElement(p1).getY() > points.getElement(p4).getY()) {
					scianyPryzm.setElement(1);
				}
				else {
					scianyPryzm.setElement(-1);
				}

			}

		}
	}
	//

}

void Delanouy3D::skalujIPrzesunCzworosciany(double wartoscZmiejszenia) {

	double wspSkali = (dy - wartoscZmiejszenia * 2) / dy, y;
	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		y = points.getElement(i).getY() * wspSkali;

		points.getElement(i).setY(wartoscZmiejszenia + y);

	}

}

void Delanouy3D::skalujCzworosciany(double gX, double gY, double gZ) {



	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		double x = points.getElement(i).getX() * gX;
		double y = points.getElement(i).getY() * gY;
		double z = points.getElement(i).getZ() * gZ;

		points.getElement(i).setXYZ(x, y, z);

	}


}
/*
void Delanouy3D::tworzeniePryzm(double grubosc,int ileWarstw,bool zeSkalowaniem){


	if(zeSkalowaniem){skalujIPrzesunCzworosciany(grubosc);}


	pryzmy.czysc(10000,10000);
	int ilePun=points.getIter();
	int *tabPunktow = new int[ilePun];
	IntList scianyPryzm(30000,30000);

	double dlPryzmy = grubosc / ileWarstw;

	for(int i=0;i<ilePun;++i){tabPunktow[i]=-1;}
	wyszukajScianyNaPryzmy(scianyPryzm);

	double x,y,z,dl;
	for(int i=0,ileS=scianyPryzm.getIter(),e,nrS,p1,p2,p3;i<ileS;i+=3){

		dl = dlPryzmy*scianyPryzm.getElement(i+2);
		e = scianyPryzm.getElement(i);
		nrS = scianyPryzm.getElement(i+1);

		switch(nrS){

			case 1:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP1();
				elements.getElement(e).setE1(-2);
				break;

			case 2:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP3();
				p3 = elements.getElement(e).getP2();
				elements.getElement(e).setE2(-2);
				break;

			case 3:
				p1 = elements.getElement(e).getP4();
				p2 = elements.getElement(e).getP1();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE3(-2);
				break;

			case 4:
				p1 = elements.getElement(e).getP1();
				p2 = elements.getElement(e).getP2();
				p3 = elements.getElement(e).getP3();
				elements.getElement(e).setE4(-2);
				break;

			default: m1->Lines->Add("zadna sciana");
		}

		if(tabPunktow[p1]==-1){
			x = points.getElement(p1).getX();
			y = points.getElement(p1).getY();
			z = points.getElement(p1).getZ();
			tabPunktow[p1]=points.getIter();

			for(int w=0;w<ileWarstw;++w){
				points.setElement(x,y+dl*(w+1),z);
			}

		}
		if(tabPunktow[p2]==-1){
			x = points.getElement(p2).getX();
			y = points.getElement(p2).getY();
			z = points.getElement(p2).getZ();
			tabPunktow[p2]=points.getIter();

			for(int w=0;w<ileWarstw;++w){
				points.setElement(x,y+dl*(w+1),z);
			}

		}
		if(tabPunktow[p3]==-1){
			x = points.getElement(p3).getX();
			y = points.getElement(p3).getY();
			z = points.getElement(p3).getZ();
			tabPunktow[p3]=points.getIter();

			for(int w=0;w<ileWarstw;++w){
				points.setElement(x,y+dl*(w+1),z);
			}

		}

		for(int w=0;w<ileWarstw;++w){

			if(w==0){pryzmy.setElement(p1,p2,p3,tabPunktow[p1],tabPunktow[p2],tabPunktow[p3]);}
			else{pryzmy.setElement(tabPunktow[p1]+w-1,tabPunktow[p2]+w-1,tabPunktow[p3]+w-1,tabPunktow[p1]+w,tabPunktow[p2]+w,tabPunktow[p3]+w);}

		}

	}


}

*/

void Delanouy3D::dopasujElementDoZiarnaNaPodstawieFunkcji(double wartoscKata, double rZasat, double rTube) {

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
		elements.getElement(i).setRodzajZiarna(0);
	}

	double ax, ay, az;
	double sX = 20, sY = 20, sZ = 35, doZastZmin = 34, doZastZmax = 36;
	double promien = rTube;
	double promienZas = rZasat;
	double alfa = wartoscKata * (M_PI / 180);
	double pX, pZ;
	//double x0,y0,z0,pY;

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		ax = elements.getElement(i).getAx();
		ay = elements.getElement(i).getAy();
		az = elements.getElement(i).getAz();

		pX = sX + (ax - sX) * cos(alfa) + (az - sZ) * sin(alfa);
		pZ = sZ + (az - sZ) * cos(alfa) - (ax - sX) * sin(alfa);

		if (sqrt((pX - sX) * (pX - sX) + (ay - sY) * (ay - sY)) <= promienZas && pZ >= doZastZmin && pZ <= doZastZmax) {
			elements.getElement(i).setRodzajZiarna(2);
		}
		//if(pZ>23 && pZ<27){elements.getElement(i).setRodzajZiarna(2);}

		else {

			if (sqrt((ax - sX) * (ax - sX) + (ay - sY) * (ay - sY)) <= promien) {
				elements.getElement(i).setRodzajZiarna(3);
			}
			else {
				elements.getElement(i).setRodzajZiarna(1);
			}

		}
	}


}

void Delanouy3D::pokazMinMaxV() {

	double V;
	double min = 99999999, max = 0;

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		V = V_objetoscT(i);

		if (min > V) {
			min = V;
		}
		if (max < V) {
			max = V;
		}

	}
	//double mnoz = 1/min;
	//m1->Lines->Add(" Min Max stos min:max");
	//m1->Lines->Add(min);
	//m1->Lines->Add(max);

	//AnsiString stosunek;
	//stosunek=mnoz*min;
	//stosunek+=" : ";
	//stosunek+=mnoz*max;
	//m1->Lines->Add(stosunek);


}

bool Delanouy3D::sprawdzCzyNalezyDoE(int nrE, const double & x, const double & y, const double & z) {

	double wynik, wynik1, wynik2, wynik3, wynik4;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
//bool flaga=true;
	int p1, p2, p3, p4;

	p1 = elements.getElement(nrE).getP1();
	p2 = elements.getElement(nrE).getP2();
	p3 = elements.getElement(nrE).getP3();
	p4 = elements.getElement(nrE).getP4();

	x1 = points.getElement(p1).getX();
	y1 = points.getElement(p1).getY();
	z1 = points.getElement(p1).getZ();

	x2 = points.getElement(p2).getX();
	y2 = points.getElement(p2).getY();
	z2 = points.getElement(p2).getZ();

	x3 = points.getElement(p3).getX();
	y3 = points.getElement(p3).getY();
	z3 = points.getElement(p3).getZ();

	x4 = points.getElement(p4).getX();
	y4 = points.getElement(p4).getY();
	z4 = points.getElement(p4).getZ();

	wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
			(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

	wynik1 = (x2 - x1) * (y3 - y1) * (z - z1) + (x3 - x1) * (y - y1) * (z2 - z1) + (x - x1) * (y2 - y1) * (z3 - z1) -
			(z2 - z1) * (y3 - y1) * (x - x1) - (z3 - z1) * (y - y1) * (x2 - x1) - (z - z1) * (y2 - y1) * (x3 - x1) ;

	wynik2 = (x2 - x1) * (y - y1) * (z4 - z1) + (x - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z - z1) -
			(z2 - z1) * (y - y1) * (x4 - x1) - (z - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x - x1) ;

	wynik3 = (x - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z - z1) + (x4 - x1) * (y - y1) * (z3 - z1) -
			(z - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x - x1) - (z4 - z1) * (y - y1) * (x3 - x1) ;

	wynik4 = (x2 - x ) * (y3 - y ) * (z4 - z ) + (x3 - x ) * (y4 - y ) * (z2 - z ) + (x4 - x ) * (y2 - y ) * (z3 - z ) -
			(z2 - z ) * (y3 - y ) * (x4 - x ) - (z3 - z ) * (y4 - y ) * (x2 - x ) - (z4 - z ) * (y2 - y ) * (x3 - x ) ;


	if (wynik1 < 0) {
		wynik1 *= -1;
	}
	if (wynik2 < 0) {
		wynik2 *= -1;
	}
	if (wynik3 < 0) {
		wynik3 *= -1;
	}
	if (wynik4 < 0) {
		wynik4 *= -1;
	}

	if (wynik + 0.00000001 > wynik1 + wynik2 + wynik3 + wynik4) {

		return true;

	}

	return false;
}

void Delanouy3D::uzupelnienieElementow() {

//bool flaga=false;
//double p,a,b,c,d,e,f,V,r,R,pol;
	double a2, b2, c2, l2, a3, c3, b3, l3, a4, b4, c4, l4;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double X, Y, Z, W, wynik, R, p;


	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		//R=sqrt(elements.getElement(i).getR());


		if (p4 < p1) {
			p = p4;
			p4 = p1;
			p1 = p;
		}
		if (p4 < p2) {
			p = p4;
			p4 = p2;
			p2 = p;
		}
		if (p4 < p3) {
			p = p4;
			p4 = p3;
			p3 = p;
		}

		if (p1 < p2) {
			p = p1;
			p1 = p2;
			p2 = p;
		}
		if (p1 < p3) {
			p = p1;
			p1 = p3;
			p3 = p;
		}

		x1 = points.getElement(p1).getX();
		y1 = points.getElement(p1).getY();
		z1 = points.getElement(p1).getZ();

		x2 = points.getElement(p2).getX();
		y2 = points.getElement(p2).getY();
		z2 = points.getElement(p2).getZ();

		x3 = points.getElement(p3).getX();
		y3 = points.getElement(p3).getY();
		z3 = points.getElement(p3).getZ();

		x4 = points.getElement(p4).getX();
		y4 = points.getElement(p4).getY();
		z4 = points.getElement(p4).getZ();

		wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
				(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

		if (wynik < 0) {
			p = p2;
			p2 = p3;
			p3 = p;
		}

		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);

		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =   (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) ;


		elements.getElement(i).setPSO(X, Y, Z, R);
		elements.getElement(i).setAx((x1 + x2 + x3 + x4) * 0.25);
		elements.getElement(i).setAy((y1 + y2 + y3 + y4) * 0.25);
		elements.getElement(i).setAz((z1 + z2 + z3 + z4) * 0.25);

		elements.getElement(i).setP1(p1);
		elements.getElement(i).setP2(p2);
		elements.getElement(i).setP3(p3);
		elements.getElement(i).setP4(p4);


	}


}

void Delanouy3D::uzupelnienieElementowWybranych(IntList & wybraneEl) {

//bool flaga=false;
//double p,a,b,c,d,e,f,V,r,R,pol;
	double a2, b2, c2, l2, a3, c3, b3, l3, a4, b4, c4, l4;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double X, Y, Z, W, wynik, R, p;


	for (int i = 0, p1, p2, p3, p4, ktory, ileE = wybraneEl.getIter(); i < ileE; ++i) {

		ktory = wybraneEl.getElement(i);

		p1 = elements.getElement(ktory).getP1();
		p2 = elements.getElement(ktory).getP2();
		p3 = elements.getElement(ktory).getP3();
		p4 = elements.getElement(ktory).getP4();
		//R=sqrt(elements.getElement(i).getR());


		if (p4 < p1) {
			p = p4;
			p4 = p1;
			p1 = p;
		}
		if (p4 < p2) {
			p = p4;
			p4 = p2;
			p2 = p;
		}
		if (p4 < p3) {
			p = p4;
			p4 = p3;
			p3 = p;
		}

		if (p1 < p2) {
			p = p1;
			p1 = p2;
			p2 = p;
		}
		if (p1 < p3) {
			p = p1;
			p1 = p3;
			p3 = p;
		}

		x1 = points.getElement(p1).getX();
		y1 = points.getElement(p1).getY();
		z1 = points.getElement(p1).getZ();

		x2 = points.getElement(p2).getX();
		y2 = points.getElement(p2).getY();
		z2 = points.getElement(p2).getZ();

		x3 = points.getElement(p3).getX();
		y3 = points.getElement(p3).getY();
		z3 = points.getElement(p3).getZ();

		x4 = points.getElement(p4).getX();
		y4 = points.getElement(p4).getY();
		z4 = points.getElement(p4).getZ();

		wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
				(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

		if (wynik < 0) {
			p = p2;
			p2 = p3;
			p3 = p;
		}

		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);

		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =   (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) ;


		elements.getElement(ktory).setPSO(X, Y, Z, R);
		elements.getElement(ktory).setAx((x1 + x2 + x3 + x4) * 0.25);
		elements.getElement(ktory).setAy((y1 + y2 + y3 + y4) * 0.25);
		elements.getElement(ktory).setAz((z1 + z2 + z3 + z4) * 0.25);

		elements.getElement(ktory).setP1(p1);
		elements.getElement(ktory).setP2(p2);
		elements.getElement(ktory).setP3(p3);
		elements.getElement(ktory).setP4(p4);


	}


}

void Delanouy3D::uzupelnienieElementuWybranego(int ktory) {

//bool flaga=false;
//double p,a,b,c,d,e,f,V,r,R,pol;
	double a2, b2, c2, l2, a3, c3, b3, l3, a4, b4, c4, l4;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double X, Y, Z, W, wynik, R, p;


	int p1 = elements.getElement(ktory).getP1();
	int p2 = elements.getElement(ktory).getP2();
	int p3 = elements.getElement(ktory).getP3();
	int p4 = elements.getElement(ktory).getP4();
	//R=sqrt(elements.getElement(i).getR());


	if (p4 < p1) {
		p = p4;
		p4 = p1;
		p1 = p;
	}
	if (p4 < p2) {
		p = p4;
		p4 = p2;
		p2 = p;
	}
	if (p4 < p3) {
		p = p4;
		p4 = p3;
		p3 = p;
	}

	if (p1 < p2) {
		p = p1;
		p1 = p2;
		p2 = p;
	}
	if (p1 < p3) {
		p = p1;
		p1 = p3;
		p3 = p;
	}

	x1 = points.getElement(p1).getX();
	y1 = points.getElement(p1).getY();
	z1 = points.getElement(p1).getZ();

	x2 = points.getElement(p2).getX();
	y2 = points.getElement(p2).getY();
	z2 = points.getElement(p2).getZ();

	x3 = points.getElement(p3).getX();
	y3 = points.getElement(p3).getY();
	z3 = points.getElement(p3).getZ();

	x4 = points.getElement(p4).getX();
	y4 = points.getElement(p4).getY();
	z4 = points.getElement(p4).getZ();

	wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
			(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

	if (wynik < 0) {
		p = p2;
		p2 = p3;
		p3 = p;
	}

	a2 = (-x3 + x1) * 2;
	b2 = (-y3 + y1) * 2;
	c2 = (-z3 + z1) * 2;
	l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

	a3 = (-x4 + x1) * 2;
	b3 = (-y4 + y1) * 2;
	c3 = (-z4 + z1) * 2;
	l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);

	a4 = (-x4 + x2) * 2;
	b4 = (-y4 + y2) * 2;
	c4 = (-z4 + z2) * 2;
	l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

	W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
	W = 1 / W;
	l2 *= W;
	l3 *= W;
	l4 *= W;

	X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
	Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
	Y = Y * -1;
	Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

	R =   (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) ;


	elements.getElement(ktory).setPSO(X, Y, Z, R);
	elements.getElement(ktory).setAx((x1 + x2 + x3 + x4) * 0.25);
	elements.getElement(ktory).setAy((y1 + y2 + y3 + y4) * 0.25);
	elements.getElement(ktory).setAz((z1 + z2 + z3 + z4) * 0.25);

	elements.getElement(ktory).setP1(p1);
	elements.getElement(ktory).setP2(p2);
	elements.getElement(ktory).setP3(p3);
	elements.getElement(ktory).setP4(p4);

}


void Delanouy3D::sprawdzWaznoscPunktow() {

	double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, wynik;

	for (int i = 0, p1, p2, p3, p4, p, ileE = elements.getIter(); i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		if (p4 < p1) {
			p = p4;
			p4 = p1;
			p1 = p;
		}
		if (p4 < p2) {
			p = p4;
			p4 = p2;
			p2 = p;
		}
		if (p4 < p3) {
			p = p4;
			p4 = p3;
			p3 = p;
		}

		if (p1 < p2) {
			p = p1;
			p1 = p2;
			p2 = p;
		}
		if (p1 < p3) {
			p = p1;
			p1 = p3;
			p3 = p;
		}

		x1 = points.getElement(p1).getX();
		y1 = points.getElement(p1).getY();
		z1 = points.getElement(p1).getZ();

		x2 = points.getElement(p2).getX();
		y2 = points.getElement(p2).getY();
		z2 = points.getElement(p2).getZ();

		x3 = points.getElement(p3).getX();
		y3 = points.getElement(p3).getY();
		z3 = points.getElement(p3).getZ();

		x4 = points.getElement(p4).getX();
		y4 = points.getElement(p4).getY();
		z4 = points.getElement(p4).getZ();

		wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
				(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

		if (wynik < 0) {
			p = p2;
			p2 = p3;
			p3 = p;
		}

		elements.getElement(i).setP1(p1);
		elements.getElement(i).setP2(p2);
		elements.getElement(i).setP3(p3);
		elements.getElement(i).setP4(p4);

	}


}

void Delanouy3D::dodajNowePunktyNaKrawedziach() {

	PunktList nowyP(elements.getIter(), elements.getIter() + 10);

	//double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
	//double np1x,np1y,np1z;
	//int p1,p2,p3,p4;

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
		/*
		p1=elements.getElement(i).getP1();
		p2=elements.getElement(i).getP2();
		p3=elements.getElement(i).getP3();
		p4=elements.getElement(i).getP4();

		x1=points.getElement(p1).getX();
		y1=points.getElement(p1).getY();
		z1=points.getElement(p1).getZ();

		x2=points.getElement(p2).getX();
		y2=points.getElement(p2).getY();
		z2=points.getElement(p2).getZ();

		x3=points.getElement(p3).getX();
		y3=points.getElement(p3).getY();
		z3=points.getElement(p3).getZ();

		x4=points.getElement(p4).getX();
		y4=points.getElement(p4).getY();
		z4=points.getElement(p4).getZ();

		np1x=x1+x2+x3+x4;
		np1y=y1+y2+y3+y4;
		np1z=z1+z2+z3+z4;
		*/


		nowyP.setElement(elements.getElement(i).getSx(), elements.getElement(i).getSy(), elements.getElement(i).getSz());

	}

	delanoyTablicaOPDodajPunkty(dx, dy, dz, false, nowyP, 0.3, 3);

}

void Delanouy3D::delanoyTablicaOPDodajPunkty(double xdl, double ydl, double zdl, bool reMes, PunktList & Tab, double ws, int te) {


	int TE2 = te;

	double ws2 = ws;

	int wymiarTE2X = (int)(xdl * ws2) + 2 * TE2;
	int wymiarTE2Y = (int)(ydl * ws2) + 2 * TE2;
	int wymiarTE2Z = (int)(zdl * ws2) + 2 * TE2;


//sprawdzNumeracje();

	IntList delElement(2000, 1000);
	IntList zlyElement(200, 100);
	IntList pomoc(2000, 1000);
	IntList sprawdz(50, 20);

	int ** tabElementow2;
	tabElementow2 = new int * [wymiarTE2Z];
	for (int i = 0; i < wymiarTE2Z; i++) {
		tabElementow2[i] = new int[wymiarTE2X * wymiarTE2Y];
	}
	for (int i = 0; i < wymiarTE2Z; i++) {
		for (int j = 0; j < wymiarTE2X * wymiarTE2Y; j++) {
			tabElementow2[i][j] = -1;
		}
	}



//int naliczanie1=0;
//int naliczanie2=0;
	int naliczanie4 = 0;
	char granica;
	double np1x, np1y, np1z;
	double x0, x1, x2, x3, x4, y0, y1, y2, y3, y4, z0, z1, z2, z3, z4, a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, WW, W, X, Y, Z;
	int wszystkiePunkty = Tab.getIter();
	int wszystkiePunktyPomoc = wszystkiePunkty;
	int ilePowtorzenPunktow = 10;

	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		x1 = points.getElement(p1).getX();
		y1 = points.getElement(p1).getY();
		z1 = points.getElement(p1).getZ();

		x2 = points.getElement(p2).getX();
		y2 = points.getElement(p2).getY();
		z2 = points.getElement(p2).getZ();

		x3 = points.getElement(p3).getX();
		y3 = points.getElement(p3).getY();
		z3 = points.getElement(p3).getZ();

		x4 = points.getElement(p4).getX();
		y4 = points.getElement(p4).getY();
		z4 = points.getElement(p4).getZ();

		np1x = (x1 + x2 + x3 + x4) * 0.25;
		np1y = (y1 + y2 + y3 + y4) * 0.25;
		np1z = (z1 + z2 + z3 + z4) * 0.25;

		tabElementow2[int(np1z * ws2) + TE2][(int(np1y * ws2) + TE2)*wymiarTE2X + int(np1x * ws2) + TE2] = i;

	}

	int punktReMes = 0;
	if (reMes) {
		punktReMes = 8;
	}

	for (int pun = punktReMes; pun < wszystkiePunkty; pun++) {


		granica = Tab.getElement(pun).getGranica();
		x0 = Tab.getElement(pun).getX();
		y0 = Tab.getElement(pun).getY();
		z0 = Tab.getElement(pun).getZ();

		delElement.ustawIter(0);

		int nrNowegoEl = elements.getIter();
		int iloscEStart = elements.getIter();
		int ileElUsunieto = 0;
		int znaleziony = -1;

		int wx2 = ((int)(x0 * ws2)) + TE2;
		int wy2 = ((int)(y0 * ws2)) + TE2;
		int wz2 = ((int)(z0 * ws2)) + TE2;


		if (tabElementow2[wz2][wy2 * wymiarTE2X + wx2] != -1) {
			znaleziony = tabElementow2[wz2][wy2 * wymiarTE2X + wx2];
			//naliczanie1++;
		}

		else {

			for (int t2z = -TE2, max = 100000000; t2z < TE2; t2z++) {

				for (int t2y = -TE2; t2y < TE2; t2y++) {

					for (int t2x = -TE2; t2x < TE2; t2x++) {

						if (tabElementow2[wz2 + t2z][(wy2 + t2y)*wymiarTE2X + (wx2 + t2x)] != -1) {
							if (max > t2x * t2x + t2y * t2y + t2z * t2z) {
								znaleziony = tabElementow2[wz2 + t2z][(wy2 + t2y) * wymiarTE2X + (wx2 + t2x)];
								max = t2x * t2x + t2y * t2y + t2z * t2z;
							}
						}

					}
				}

			}

		}




		if (znaleziony == -1) {;}


		znaleziony = wyszukajElementOP(znaleziony, x0, y0, z0, zlyElement, pomoc);
		//elementyDoDelnoyaZeSprObj



		if (elementyDoDelnoyaZeSprObj(znaleziony, x0, y0, z0, sprawdz)) {

			points.setElement(x0, y0, z0, granica);
			ileElUsunieto = sprawdz.getIter();



			for (int i = 0, wybrany, iter; i < ileElUsunieto; i++) {

				wybrany = sprawdz.getElement(i);


// 0 - nowa pozycja powstalego elementu
// 1 - pozycja elementu niszczonego
// 2 - pierwszy wezel nowego el
// 3 - drugi wezel nowego elementu
// 4 - trzeci wezel nowego elementu
// 5 - kolejnosc sciany w elemencie niszczonym
// 6 - 1 sasiedni element
// 7 - 2 sasiedni element
// 8 - 3 sasiedni element
// 9 - 4 sasiedni element

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(2);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(3);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(4);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				iter = delElement.getIter() - 40;

				for (int a = 0; a < iter; a += 10) {

					if (delElement.getElement(iter + 2) == delElement.getElement(a + 2) && delElement.getElement(iter + 3) == delElement.getElement(a + 4) && delElement.getElement(iter + 4) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 12) == delElement.getElement(a + 2) && delElement.getElement(iter + 13) == delElement.getElement(a + 4) && delElement.getElement(iter + 14) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 10, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 22) == delElement.getElement(a + 2) && delElement.getElement(iter + 23) == delElement.getElement(a + 4) && delElement.getElement(iter + 24) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 20, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 32) == delElement.getElement(a + 2) && delElement.getElement(iter + 33) == delElement.getElement(a + 4) && delElement.getElement(iter + 34) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 30, -1);
						delElement.podmienElement(a, -1);
					}

				}
			}




			int iter = delElement.getIter();

			for (int i = 0; i < iter; i += 10) {


				if (delElement.getElement(i) != -1) {

					double Rr;
					int e1 = delElement.getElement(i + 2), e2 = delElement.getElement(i + 3), e3 = delElement.getElement(i + 4);

					x1 = points.getElement(e1).getX();
					x2 = points.getElement(e2).getX();
					x3 = points.getElement(e3).getX();

					y1 = points.getElement(e1).getY();
					y2 = points.getElement(e2).getY();
					y3 = points.getElement(e3).getY();

					z1 = points.getElement(e1).getZ();
					z2 = points.getElement(e2).getZ();
					z3 = points.getElement(e3).getZ();


					WW = (x0 - x1) * (y2 - y1) * (z3 - z1) + (x2 - x1) * (y3 - y1) * (z0 - z1) + (x3 - x1) * (y0 - y1) * (z2 - z1)
							- (z0 - z1) * (y2 - y1) * (x3 - x1) - (z2 - z1) * (y3 - y1) * (x0 - x1) - (z3 - z1) * (y0 - y1) * (x2 - x1);


					//WW  =x0*y2*z3+x2*y3*z0+x3*y0*z2-z0*y2*x3-z2*y3*x0-z3*y0*x2;
					if (WW) {

						a2 = (-x3 + x1) * 2;
						b2 = (-y3 + y1) * 2;
						c2 = (-z3 + z1) * 2;
						l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

						a3 = (-x0 + x1) * 2;
						b3 = (-y0 + y1) * 2;
						c3 = (-z0 + z1) * 2;
						l3 = -(x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1);

						a4 = (-x0 + x2) * 2;
						b4 = (-y0 + y2) * 2;
						c4 = (-z0 + z2) * 2;
						l4 = -(x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 - z2 * z2);

						W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
						W = 1 / W;
						l2 *= W;
						l3 *= W;
						l4 *= W;
						if (W) {
							X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
							Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
							Y = Y * -1;
							Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

							Rr =  (X - x0) * (X - x0) + (Y - y0) * (Y - y0) + (Z - z0) * (Z - z0);
							//Rr = X*X-2*X*x0+x0*x0 + Y*Y-2*Y*y0+y0*y0 + Z*Z-2*Z*z0+z0*z0;

							elements.setElement(e1, e2, e3, points.getIter() - 1, X, Y, Z, Rr, (x1 + x2 + x3 + x0) * 0.25, (y1 + y2 + y3 + y0) * 0.25, (z1 + z2 + z3 + z0) * 0.25);
							//elements.setElement(e1,e2,e3,points.getIter()-1,X,Y,Z,Rr,(x1+x2+x3+x0)*0.25,(y1+y2+y3+y0)*0.25,(z1+z2+z3+z0)*0.25);
						}
						else {
							delElement.podmienElement(i, -1);
						}
					}
					else {
						delElement.podmienElement(i, -1);
					}
				}
			}






			for (int a = 0; a < iter; a += 10) {

				if (delElement.getElement(a) != -1) {
					delElement.podmienElement(a, nrNowegoEl++);
				}

			}




			int iter2 = iter;

			for ( int a = iter - 10, i = 1; i < iter; a -= 10) {

				if (a >= 0) {
					if (delElement.getElement(a) != -1) {

						delElement.podmienElement(a, delElement.getElement(i));
						i += 40;

					}
				}
				else {
					iter2 += int(((elements.getIter() - iloscEStart) - iter * 0.025) * 40);

					break;
				}


			}

			nrNowegoEl -= elements.getIter();

//////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {

					int sasiad = -1;

					switch (delElement.getElement(i + 5)) {

						case 1:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE1();

							break;
						case 2:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE2();

							break;
						case 3:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE3();

							break;
						case 4:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE4();

							break;
					}

					//i+9 ostatnia sciana 123
					delElement.podmienElement(i + 9, sasiad);

					//podejscie=0;
					//if(sasiad!=-1 && podejscie==0){
					if (sasiad != -1) {
						//if(0){


						if (elements.getElement(sasiad).getP4() == delElement.getElement(i + 2)) {

							if (elements.getElement(sasiad).getP1() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE1(delElement.getElement(i));
							}
							else if (elements.getElement(sasiad).getP2() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE2(delElement.getElement(i));
							}
							else {
								elements.getElement(sasiad).setE3(delElement.getElement(i));
							}

						}

						else {
							elements.getElement(sasiad).setE4(delElement.getElement(i));
						}
					}


					for (int j = i + 10; j < iter; j += 10) {

						if (delElement.getElement(j) != -1) {

							if (delElement.getElement(i + 2) == delElement.getElement(j + 2)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 3)) {
								if (delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+3 = j+2");
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 4)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+4 = j+2");

								}
							}
							else {

								if (delElement.getElement(i + 3) == delElement.getElement(j + 2) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 3) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 4) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}

							}

						}
					}

				}

			}

			/*
			if(podejscie){
				optymalizacjaStosunekRdorDODelanoya(0.1,iloscEStart-1);
				Tab.getElement(pun).setPunkt(points.getLastElement());
				points.zmienIterO(-1);pun--;
				elements.ustawIter(iloscEStart);
				podejscie=0;
			}
			else{
				podejscie=1;
			*/
			//usuwanie elementow
			//usuwanie wpisu elementu w tablicy

			for (int i = 1; i < iter2; i += 40) {
				elements.deleteElement(delElement.getElement(i));
			}

			for (int i = 0, ww; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {
					ww = delElement.getElement(i);
					elements.getElement(ww).setElementySasiad(delElement.getElement(i + 6), delElement.getElement(i + 7), delElement.getElement(i + 8), delElement.getElement(i + 9));

					tabElementow2[int(elements.getElement(ww).getAz()*ws2) + TE2][(int(elements.getElement(ww).getAy()*ws2) + TE2)*wymiarTE2X + int(elements.getElement(ww).getAx()*ws2) + TE2] = ww;

				}

			}

			for (int i = iter2 + 1; i < iter; i += 40) {

				elements.deleteElementPodmienS(delElement.getElement(i));

			}
			//}

		}

		else {

			--wszystkiePunktyPomoc;

			if (wszystkiePunktyPomoc > pun) {
				Tab.zamianaMiejsc(pun, wszystkiePunktyPomoc);
				--pun;
				++naliczanie4;
			}
			else if (ilePowtorzenPunktow > 0 && wszystkiePunktyPomoc < pun) {
				wszystkiePunktyPomoc = wszystkiePunkty - 1;
				--pun;
				--ilePowtorzenPunktow;

			}

		}
	}
//wymiarTE2++;
	for (int i = 0; i < wymiarTE2Z; i++) {
		delete []tabElementow2[i];
	}
	delete []tabElementow2;


//m1->Lines->Add("+++++++++");
	if (naliczanie4) {;}
//m1->Lines->Add(naliczanie2);
//m1->Lines->Add(naliczanie3-(naliczanie2+naliczanie1));

}



void Delanouy3D::delanoyTablicaOPDodajPunktyZTablicaPolozeniaEl(int ** tabElementow2, double xdl, double ydl, double zdl, bool reMes, PunktList & Tab, double ws, int te) {


	int TE2 = te;

	double ws2 = ws;

	int wymiarTE2X = (int)(xdl * ws2) + 2 * TE2;
//int wymiarTE2Y = (int)(ydl*ws2)+2*TE2;
//int wymiarTE2Z = (int)(zdl*ws2)+2*TE2;


//sprawdzNumeracje();

	IntList delElement(2000, 1000);
	IntList zlyElement(200, 100);
	IntList pomoc(2000, 1000);
	IntList sprawdz(50, 20);



//int naliczanie1=0;
//int naliczanie2=0;
	int naliczanie4 = 0;
	char granica;
//double np1x,np1y,np1z;
	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3, a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, WW, W, X, Y, Z;
	int wszystkiePunkty = Tab.getIter();
	int wszystkiePunktyPomoc = wszystkiePunkty;
	int ilePowtorzenPunktow = 10;



	int punktReMes = 0;
	if (reMes) {
		punktReMes = 8;
	}

	for (int pun = punktReMes; pun < wszystkiePunkty; pun++) {


		granica = Tab.getElement(pun).getGranica();
		x0 = Tab.getElement(pun).getX();
		y0 = Tab.getElement(pun).getY();
		z0 = Tab.getElement(pun).getZ();

		delElement.ustawIter(0);

		int nrNowegoEl = elements.getIter();
		int iloscEStart = elements.getIter();
		int ileElUsunieto = 0;
		int znaleziony = -1;

		int wx2 = ((int)(x0 * ws2)) + TE2;
		int wy2 = ((int)(y0 * ws2)) + TE2;
		int wz2 = ((int)(z0 * ws2)) + TE2;


		if (tabElementow2[wz2][wy2 * wymiarTE2X + wx2] != -1) {
			znaleziony = tabElementow2[wz2][wy2 * wymiarTE2X + wx2];
			//naliczanie1++;
		}

		else {

			for (int t2z = -TE2, max = 100000000; t2z < TE2; t2z++) {

				for (int t2y = -TE2; t2y < TE2; t2y++) {

					for (int t2x = -TE2; t2x < TE2; t2x++) {

						if (tabElementow2[wz2 + t2z][(wy2 + t2y)*wymiarTE2X + (wx2 + t2x)] != -1) {
							if (max > t2x * t2x + t2y * t2y + t2z * t2z) {
								znaleziony = tabElementow2[wz2 + t2z][(wy2 + t2y) * wymiarTE2X + (wx2 + t2x)];
								max = t2x * t2x + t2y * t2y + t2z * t2z;
							}
						}

					}
				}

			}

		}




		if (znaleziony <= -1) {;}


		znaleziony = wyszukajElementOP(znaleziony, x0, y0, z0, zlyElement, pomoc);
		//elementyDoDelnoyaZeSprObj



		if (elementyDoDelnoyaZeSprObj(znaleziony, x0, y0, z0, sprawdz)) {

			points.setElement(x0, y0, z0, granica);
			ileElUsunieto = sprawdz.getIter();



			for (int i = 0, wybrany, iter; i < ileElUsunieto; i++) {

				wybrany = sprawdz.getElement(i);


// 0 - nowa pozycja powstalego elementu
// 1 - pozycja elementu niszczonego
// 2 - pierwszy wezel nowego el
// 3 - drugi wezel nowego elementu
// 4 - trzeci wezel nowego elementu
// 5 - kolejnosc sciany w elemencie niszczonym
// 6 - 1 sasiedni element
// 7 - 2 sasiedni element
// 8 - 3 sasiedni element
// 9 - 4 sasiedni element

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(2);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(3);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(4);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				iter = delElement.getIter() - 40;

				for (int a = 0; a < iter; a += 10) {

					if (delElement.getElement(iter + 2) == delElement.getElement(a + 2) && delElement.getElement(iter + 3) == delElement.getElement(a + 4) && delElement.getElement(iter + 4) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 12) == delElement.getElement(a + 2) && delElement.getElement(iter + 13) == delElement.getElement(a + 4) && delElement.getElement(iter + 14) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 10, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 22) == delElement.getElement(a + 2) && delElement.getElement(iter + 23) == delElement.getElement(a + 4) && delElement.getElement(iter + 24) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 20, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 32) == delElement.getElement(a + 2) && delElement.getElement(iter + 33) == delElement.getElement(a + 4) && delElement.getElement(iter + 34) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 30, -1);
						delElement.podmienElement(a, -1);
					}

				}
			}




			int iter = delElement.getIter();

			for (int i = 0; i < iter; i += 10) {


				if (delElement.getElement(i) != -1) {

					double Rr;
					int e1 = delElement.getElement(i + 2), e2 = delElement.getElement(i + 3), e3 = delElement.getElement(i + 4);

					x1 = points.getElement(e1).getX();
					x2 = points.getElement(e2).getX();
					x3 = points.getElement(e3).getX();

					y1 = points.getElement(e1).getY();
					y2 = points.getElement(e2).getY();
					y3 = points.getElement(e3).getY();

					z1 = points.getElement(e1).getZ();
					z2 = points.getElement(e2).getZ();
					z3 = points.getElement(e3).getZ();


					WW = (x0 - x1) * (y2 - y1) * (z3 - z1) + (x2 - x1) * (y3 - y1) * (z0 - z1) + (x3 - x1) * (y0 - y1) * (z2 - z1)
							- (z0 - z1) * (y2 - y1) * (x3 - x1) - (z2 - z1) * (y3 - y1) * (x0 - x1) - (z3 - z1) * (y0 - y1) * (x2 - x1);


					//WW  =x0*y2*z3+x2*y3*z0+x3*y0*z2-z0*y2*x3-z2*y3*x0-z3*y0*x2;
					if (WW) {

						a2 = (-x3 + x1) * 2;
						b2 = (-y3 + y1) * 2;
						c2 = (-z3 + z1) * 2;
						l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

						a3 = (-x0 + x1) * 2;
						b3 = (-y0 + y1) * 2;
						c3 = (-z0 + z1) * 2;
						l3 = -(x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1);

						a4 = (-x0 + x2) * 2;
						b4 = (-y0 + y2) * 2;
						c4 = (-z0 + z2) * 2;
						l4 = -(x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 - z2 * z2);

						W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;

						if (W > 0.000001) {
							W = 1 / W;
							l2 *= W;
							l3 *= W;
							l4 *= W;
							if (W) {
								X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
								Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
								Y = Y * -1;
								Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

								Rr =  (X - x0) * (X - x0) + (Y - y0) * (Y - y0) + (Z - z0) * (Z - z0);
								//Rr = X*X-2*X*x0+x0*x0 + Y*Y-2*Y*y0+y0*y0 + Z*Z-2*Z*z0+z0*z0;

								elements.setElement(e1, e2, e3, points.getIter() - 1, X, Y, Z, Rr, (x1 + x2 + x3 + x0) * 0.25, (y1 + y2 + y3 + y0) * 0.25, (z1 + z2 + z3 + z0) * 0.25);
								//elements.setElement(e1,e2,e3,points.getIter()-1,X,Y,Z,Rr,(x1+x2+x3+x0)*0.25,(y1+y2+y3+y0)*0.25,(z1+z2+z3+z0)*0.25);
							}

						}
						else {
							delElement.podmienElement(i, -1);
						}
					}
					else {
						delElement.podmienElement(i, -1);
					}
				}
			}






			for (int a = 0; a < iter; a += 10) {

				if (delElement.getElement(a) != -1) {
					delElement.podmienElement(a, nrNowegoEl++);
				}

			}




			int iter2 = iter;

			for ( int a = iter - 10, i = 1; i < iter; a -= 10) {

				if (a >= 0) {
					if (delElement.getElement(a) != -1) {

						delElement.podmienElement(a, delElement.getElement(i));
						i += 40;

					}
				}
				else {
					iter2 += int(((elements.getIter() - iloscEStart) - iter * 0.025) * 40);

					break;
				}


			}

			nrNowegoEl -= elements.getIter();

//////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {

					int sasiad = -1;

					switch (delElement.getElement(i + 5)) {

						case 1:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE1();

							break;
						case 2:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE2();

							break;
						case 3:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE3();

							break;
						case 4:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE4();

							break;
					}

					//i+9 ostatnia sciana 123
					delElement.podmienElement(i + 9, sasiad);

					//podejscie=0;
					//if(sasiad!=-1 && podejscie==0){
					if (sasiad != -1) {
						//if(0){


						if (elements.getElement(sasiad).getP4() == delElement.getElement(i + 2)) {

							if (elements.getElement(sasiad).getP1() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE1(delElement.getElement(i));
							}
							else if (elements.getElement(sasiad).getP2() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE2(delElement.getElement(i));
							}
							else {
								elements.getElement(sasiad).setE3(delElement.getElement(i));
							}

						}

						else {
							elements.getElement(sasiad).setE4(delElement.getElement(i));
						}
					}


					for (int j = i + 10; j < iter; j += 10) {

						if (delElement.getElement(j) != -1) {

							if (delElement.getElement(i + 2) == delElement.getElement(j + 2)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 3)) {
								if (delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+3 = j+2");
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 4)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+4 = j+2");

								}
							}
							else {

								if (delElement.getElement(i + 3) == delElement.getElement(j + 2) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 3) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 4) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}

							}

						}
					}

				}

			}

			/*
			if(podejscie){
				optymalizacjaStosunekRdorDODelanoya(0.1,iloscEStart-1);
				Tab.getElement(pun).setPunkt(points.getLastElement());
				points.zmienIterO(-1);pun--;
				elements.ustawIter(iloscEStart);
				podejscie=0;
			}
			else{
				podejscie=1;
			*/
			//usuwanie elementow
			//usuwanie wpisu elementu w tablicy




			for (int i = 1; i < iter2; i += 40) {
				elements.deleteElement(delElement.getElement(i));
			}

			for (int i = 0, ww; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {
					ww = delElement.getElement(i);
					elements.getElement(ww).setElementySasiad(delElement.getElement(i + 6), delElement.getElement(i + 7), delElement.getElement(i + 8), delElement.getElement(i + 9));

					tabElementow2[int(elements.getElement(ww).getAz()*ws2) + TE2][(int(elements.getElement(ww).getAy()*ws2) + TE2)*wymiarTE2X + int(elements.getElement(ww).getAx()*ws2) + TE2] = ww;

				}

			}

			for (int i = iter2 + 1; i < iter; i += 40) {

				elements.deleteElementPodmienS(delElement.getElement(i));

			}
			//}

		}

		else {

			--wszystkiePunktyPomoc;

			if (wszystkiePunktyPomoc > pun) {
				Tab.zamianaMiejsc(pun, wszystkiePunktyPomoc);
				--pun;
				++naliczanie4;
			}
			else if (ilePowtorzenPunktow > 0 && wszystkiePunktyPomoc < pun) {
				wszystkiePunktyPomoc = wszystkiePunkty - 1;
				--pun;
				--ilePowtorzenPunktow;
				//m1->Lines->Add(ilePowtorzenPunktow);
			}

		}
	}
//wymiarTE2++;



//m1->Lines->Add("+++++++++");
	if (naliczanie4) {;}
//m1->Lines->Add(naliczanie2);
//m1->Lines->Add(naliczanie3-(naliczanie2+naliczanie1));



}

void Delanouy3D::delanoyTablicaOPDodajPunktyZTablicaPolozeniaElUstalJakoscEl(int ** tabElementow2, double xdl, double ydl, double zdl, bool reMes, PunktList & Tab, double ws, int te) {


//int iterOstatniegoElementu = elements.getIter();

	int TE2 = te;

	double ws2 = ws;

	int wymiarTE2X = (int)(xdl * ws2) + 2 * TE2;
//int wymiarTE2Y = (int)(ydl*ws2)+2*TE2;
//int wymiarTE2Z = (int)(zdl*ws2)+2*TE2;


//sprawdzNumeracje();

	IntList delElement(2000, 1000);
	IntList zlyElement(200, 100);
	IntList pomoc(2000, 1000);
	IntList sprawdz(50, 20);



//int naliczanie1=0;
//int naliczanie2=0;
	int naliczanie4 = 0;
	char granica;
//double np1x,np1y,np1z;
	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3, a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, WW, W, X, Y, Z;
	double r_r, a_r, b_r, c_r, d_r, e_r, f_r, pol_r, p_r, V_r;

	int wszystkiePunkty = Tab.getIter();
	int wszystkiePunktyPomoc = wszystkiePunkty;
	int ilePowtorzenPunktow = 10;



	int punktReMes = 0;
	if (reMes) {
		punktReMes = 8;
	}

	for (int pun = punktReMes; pun < wszystkiePunkty; pun++) {


		granica = Tab.getElement(pun).getGranica();
		x0 = Tab.getElement(pun).getX();
		y0 = Tab.getElement(pun).getY();
		z0 = Tab.getElement(pun).getZ();

		delElement.ustawIter(0);

		int nrNowegoEl = elements.getIter();
		int iloscEStart = elements.getIter();
		int ileElUsunieto = 0;
		int znaleziony = -1;

		int wx2 = ((int)(x0 * ws2)) + TE2;
		int wy2 = ((int)(y0 * ws2)) + TE2;
		int wz2 = ((int)(z0 * ws2)) + TE2;


		if (tabElementow2[wz2][wy2 * wymiarTE2X + wx2] != -1) {
			znaleziony = tabElementow2[wz2][wy2 * wymiarTE2X + wx2];
			//naliczanie1++;
		}

		else {

			for (int t2z = -TE2, max = 100000000; t2z < TE2; t2z++) {

				for (int t2y = -TE2; t2y < TE2; t2y++) {

					for (int t2x = -TE2; t2x < TE2; t2x++) {

						if (tabElementow2[wz2 + t2z][(wy2 + t2y)*wymiarTE2X + (wx2 + t2x)] != -1) {
							if (max > t2x * t2x + t2y * t2y + t2z * t2z) {
								znaleziony = tabElementow2[wz2 + t2z][(wy2 + t2y) * wymiarTE2X + (wx2 + t2x)];
								max = t2x * t2x + t2y * t2y + t2z * t2z;
							}
						}

					}
				}

			}

		}




		if (znaleziony == -1) {;}


		znaleziony = wyszukajElementOP(znaleziony, x0, y0, z0, zlyElement, pomoc);
		//elementyDoDelnoyaZeSprObj



		if (elementyDoDelnoyaZeSprObj(znaleziony, x0, y0, z0, sprawdz)) {

			points.setElement(x0, y0, z0, granica);
			ileElUsunieto = sprawdz.getIter();



			for (int i = 0, wybrany, iter; i < ileElUsunieto; i++) {

				wybrany = sprawdz.getElement(i);


// 0 - nowa pozycja powstalego elementu
// 1 - pozycja elementu niszczonego
// 2 - pierwszy wezel nowego el
// 3 - drugi wezel nowego elementu
// 4 - trzeci wezel nowego elementu
// 5 - kolejnosc sciany w elemencie niszczonym
// 6 - 1 sasiedni element
// 7 - 2 sasiedni element
// 8 - 3 sasiedni element
// 9 - 4 sasiedni element

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(2);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(3);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(4);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				iter = delElement.getIter() - 40;

				for (int a = 0; a < iter; a += 10) {

					if (delElement.getElement(iter + 2) == delElement.getElement(a + 2) && delElement.getElement(iter + 3) == delElement.getElement(a + 4) && delElement.getElement(iter + 4) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 12) == delElement.getElement(a + 2) && delElement.getElement(iter + 13) == delElement.getElement(a + 4) && delElement.getElement(iter + 14) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 10, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 22) == delElement.getElement(a + 2) && delElement.getElement(iter + 23) == delElement.getElement(a + 4) && delElement.getElement(iter + 24) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 20, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 32) == delElement.getElement(a + 2) && delElement.getElement(iter + 33) == delElement.getElement(a + 4) && delElement.getElement(iter + 34) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 30, -1);
						delElement.podmienElement(a, -1);
					}

				}
			}




			int iter = delElement.getIter();

			for (int i = 0; i < iter; i += 10) {


				if (delElement.getElement(i) != -1) {

					double Rr;
					int e1 = delElement.getElement(i + 2), e2 = delElement.getElement(i + 3), e3 = delElement.getElement(i + 4);

					x1 = points.getElement(e1).getX();
					x2 = points.getElement(e2).getX();
					x3 = points.getElement(e3).getX();

					y1 = points.getElement(e1).getY();
					y2 = points.getElement(e2).getY();
					y3 = points.getElement(e3).getY();

					z1 = points.getElement(e1).getZ();
					z2 = points.getElement(e2).getZ();
					z3 = points.getElement(e3).getZ();


					WW = (x0 - x1) * (y2 - y1) * (z3 - z1) + (x2 - x1) * (y3 - y1) * (z0 - z1) + (x3 - x1) * (y0 - y1) * (z2 - z1)
							- (z0 - z1) * (y2 - y1) * (x3 - x1) - (z2 - z1) * (y3 - y1) * (x0 - x1) - (z3 - z1) * (y0 - y1) * (x2 - x1);

					//WW  =x0*y2*z3+x2*y3*z0+x3*y0*z2-z0*y2*x3-z2*y3*x0-z3*y0*x2;
					if (WW) {

						a2 = (-x3 + x1) * 2;
						b2 = (-y3 + y1) * 2;
						c2 = (-z3 + z1) * 2;
						l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

						a3 = (-x0 + x1) * 2;
						b3 = (-y0 + y1) * 2;
						c3 = (-z0 + z1) * 2;
						l3 = -(x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1);

						a4 = (-x0 + x2) * 2;
						b4 = (-y0 + y2) * 2;
						c4 = (-z0 + z2) * 2;
						l4 = -(x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 - z2 * z2);

						W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
						if (W) {

							W = 1 / W;
							l2 *= W;
							l3 *= W;
							l4 *= W;

							X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
							Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
							Y = Y * -1;
							Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

							Rr =  (X - x0) * (X - x0) + (Y - y0) * (Y - y0) + (Z - z0) * (Z - z0);

							a_r = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
							b_r = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
							c_r = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
							d_r = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0));
							e_r = sqrt((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2) + (z0 - z2) * (z0 - z2));
							f_r = sqrt((x0 - x3) * (x0 - x3) + (y0 - y3) * (y0 - y3) + (z0 - z3) * (z0 - z3));

							pol_r = 0;
							p_r = 0.5 * (a_r + d_r + e_r);

							pol_r += sqrt(p_r * (p_r - a_r) * (p_r - d_r) * (p_r - e_r));
							p_r = 0.5 * (b_r + f_r + e_r);
							pol_r += sqrt(p_r * (p_r - b_r) * (p_r - f_r) * (p_r - e_r));
							p_r = 0.5 * (c_r + d_r + f_r);
							pol_r += sqrt(p_r * (p_r - c_r) * (p_r - d_r) * (p_r - f_r));
							p_r = 0.5 * (a_r + b_r + c_r);
							pol_r += sqrt(p_r * (p_r - a_r) * (p_r - b_r) * (p_r - c_r));

							V_r = ((x2 - x1) * (y3 - y1) * (z0 - z1) + (x3 - x1) * (y0 - y1) * (z2 - z1) + (x0 - x1) * (y2 - y1) * (z3 - z1) -
											(z2 - z1) * (y3 - y1) * (x0 - x1) - (z3 - z1) * (y0 - y1) * (x2 - x1) - (z0 - z1) * (y2 - y1) * (x3 - x1)) / 6;

							r_r = (V_r * 3) / pol_r;
							r_r = sqrt(Rr) / r_r;
							//Rr = X*X-2*X*x0+x0*x0 + Y*Y-2*Y*y0+y0*y0 + Z*Z-2*Z*z0+z0*z0;

							elements.setElement(e1, e2, e3, points.getIter() - 1, X, Y, Z, Rr, (x1 + x2 + x3 + x0) * 0.25, (y1 + y2 + y3 + y0) * 0.25, (z1 + z2 + z3 + z0) * 0.25, r_r); //
							//elements.setElement(e1,e2,e3,points.getIter()-1,X,Y,Z,Rr,(x1+x2+x3+x0)*0.25,(y1+y2+y3+y0)*0.25,(z1+z2+z3+z0)*0.25);
						}
						else {
							delElement.podmienElement(i, -1);
						}
					}
					else {
						delElement.podmienElement(i, -1);
					}
				}
			}






			for (int a = 0; a < iter; a += 10) {

				if (delElement.getElement(a) != -1) {
					delElement.podmienElement(a, nrNowegoEl++);
				}

			}




			int iter2 = iter;

			for ( int a = iter - 10, i = 1; i < iter; a -= 10) {

				if (a >= 0) {
					if (delElement.getElement(a) != -1) {

						delElement.podmienElement(a, delElement.getElement(i));
						i += 40;

					}
				}
				else {
					iter2 += int(((elements.getIter() - iloscEStart) - iter * 0.025) * 40);

					break;
				}


			}

			nrNowegoEl -= elements.getIter();

//////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {

					int sasiad = -1;

					switch (delElement.getElement(i + 5)) {

						case 1:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE1();

							break;
						case 2:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE2();

							break;
						case 3:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE3();

							break;
						case 4:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE4();

							break;
					}

					//i+9 ostatnia sciana 123
					delElement.podmienElement(i + 9, sasiad);

					//podejscie=0;
					//if(sasiad!=-1 && podejscie==0){
					if (sasiad != -1) {
						//if(0){


						if (elements.getElement(sasiad).getP4() == delElement.getElement(i + 2)) {

							if (elements.getElement(sasiad).getP1() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE1(delElement.getElement(i));
							}
							else if (elements.getElement(sasiad).getP2() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE2(delElement.getElement(i));
							}
							else {
								elements.getElement(sasiad).setE3(delElement.getElement(i));
							}

						}

						else {
							elements.getElement(sasiad).setE4(delElement.getElement(i));
						}
					}


					for (int j = i + 10; j < iter; j += 10) {

						if (delElement.getElement(j) != -1) {

							if (delElement.getElement(i + 2) == delElement.getElement(j + 2)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 3)) {
								if (delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+3 = j+2");
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 4)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+4 = j+2");

								}
							}
							else {

								if (delElement.getElement(i + 3) == delElement.getElement(j + 2) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 3) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 4) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}

							}

						}
					}

				}

			}

			/*
			if(podejscie){
				optymalizacjaStosunekRdorDODelanoya(0.1,iloscEStart-1);
				Tab.getElement(pun).setPunkt(points.getLastElement());
				points.zmienIterO(-1);pun--;
				elements.ustawIter(iloscEStart);
				podejscie=0;
			}
			else{
				podejscie=1;
			*/
			//usuwanie elementow
			//usuwanie wpisu elementu w tablicy




			for (int i = 1; i < iter2; i += 40) {
				elements.deleteElement(delElement.getElement(i));
			}

			for (int i = 0, ww; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {
					ww = delElement.getElement(i);
					elements.getElement(ww).setElementySasiad(delElement.getElement(i + 6), delElement.getElement(i + 7), delElement.getElement(i + 8), delElement.getElement(i + 9));

					tabElementow2[int(elements.getElement(ww).getAz()*ws2) + TE2][(int(elements.getElement(ww).getAy()*ws2) + TE2)*wymiarTE2X + int(elements.getElement(ww).getAx()*ws2) + TE2] = ww;

				}

			}

			for (int i = iter2 + 1; i < iter; i += 40) {

				elements.deleteElementPodmienS(delElement.getElement(i));

			}
			//}

		}

		else {

			--wszystkiePunktyPomoc;

			if (wszystkiePunktyPomoc > pun) {
				Tab.zamianaMiejsc(pun, wszystkiePunktyPomoc);
				--pun;
				++naliczanie4;
			}
			else if (ilePowtorzenPunktow > 0 && wszystkiePunktyPomoc < pun) {
				wszystkiePunktyPomoc = wszystkiePunkty - 1;
				--pun;
				--ilePowtorzenPunktow;
				//m1->Lines->Add(ilePowtorzenPunktow);
			}

		}
	}
//wymiarTE2++;



//m1->Lines->Add("+++++++++");
	if (naliczanie4) {;}
//m1->Lines->Add(naliczanie2);
//m1->Lines->Add(naliczanie3-(naliczanie2+naliczanie1));








}

void Delanouy3D::delanoyTablicaDodajPunktyWOparciuOTabliceElementow(PunktList & Tab, IntList & zElementu) {

//sprawdzNumeracje();
	IntList delElement(2000, 1000);
	IntList zlyElement(200, 100);
	IntList pomoc(2000, 1000);
	IntList sprawdz(50, 20);


//int naliczanie1=0;
//int naliczanie2=0;
	int naliczanie4 = 0;
	char granica;

	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3, a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, WW, W, X, Y, Z;
	int wszystkiePunkty = Tab.getIter();
	int wszystkiePunktyPomoc = wszystkiePunkty;
	int ilePowtorzenPunktow = 10;
	int punktReMes = 0;

	for (int pun = punktReMes; pun < wszystkiePunkty; pun++) {


		x0 =  Tab.getElement(pun).getX();
		y0 =  Tab.getElement(pun).getY();
		z0 =  Tab.getElement(pun).getZ();
		granica = Tab.getElement(pun).getGranica();
		delElement.ustawIter(0);

		int nrNowegoEl = elements.getIter();
		int iloscEStart = elements.getIter();
		int ileElUsunieto = 0;
		int znaleziony = -1;

		znaleziony = zElementu.getElement(pun);

		//if(znaleziony==-1){wyszukajElementOP(znaleziony,x0,y0,z0,zlyElement,pomoc);}
		//znaleziony = wyszukajElement(x0,y0,z0,sprawdz);

		if (elementyDoDelnoya(znaleziony, x0, y0, z0, sprawdz)) {

			points.setElement(x0, y0, z0, granica);
			ileElUsunieto = sprawdz.getIter();

			for (int i = 0, wybrany, iter; i < ileElUsunieto; i++) {

				wybrany = sprawdz.getElement(i);

// 0 - nowa pozycja powstalego elementu
// 1 - pozycja elementu niszczonego
// 2 - pierwszy wezel nowego el
// 3 - drugi wezel nowego elementu
// 4 - trzeci wezel nowego elementu
// 5 - kolejnosc sciany w elemencie niszczonym
// 6 - 1 sasiedni element
// 7 - 2 sasiedni element
// 8 - 3 sasiedni element
// 9 - 4 sasiedni element

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(2);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(3);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(4);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				iter = delElement.getIter() - 40;

				for (int a = 0; a < iter; a += 10) {

					if (delElement.getElement(iter + 2) == delElement.getElement(a + 2) && delElement.getElement(iter + 3) == delElement.getElement(a + 4) && delElement.getElement(iter + 4) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 12) == delElement.getElement(a + 2) && delElement.getElement(iter + 13) == delElement.getElement(a + 4) && delElement.getElement(iter + 14) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 10, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 22) == delElement.getElement(a + 2) && delElement.getElement(iter + 23) == delElement.getElement(a + 4) && delElement.getElement(iter + 24) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 20, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 32) == delElement.getElement(a + 2) && delElement.getElement(iter + 33) == delElement.getElement(a + 4) && delElement.getElement(iter + 34) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 30, -1);
						delElement.podmienElement(a, -1);
					}

				}
			}




			int iter = delElement.getIter();

			for (int i = 0; i < iter; i += 10) {


				if (delElement.getElement(i) != -1) {

					double Rr;
					int e1 = delElement.getElement(i + 2), e2 = delElement.getElement(i + 3), e3 = delElement.getElement(i + 4);

					x1 = points.getElement(e1).getX();
					x2 = points.getElement(e2).getX();
					x3 = points.getElement(e3).getX();

					y1 = points.getElement(e1).getY();
					y2 = points.getElement(e2).getY();
					y3 = points.getElement(e3).getY();

					z1 = points.getElement(e1).getZ();
					z2 = points.getElement(e2).getZ();
					z3 = points.getElement(e3).getZ();


					WW = (x0 - x1) * (y2 - y1) * (z3 - z1) + (x2 - x1) * (y3 - y1) * (z0 - z1) + (x3 - x1) * (y0 - y1) * (z2 - z1)
							- (z0 - z1) * (y2 - y1) * (x3 - x1) - (z2 - z1) * (y3 - y1) * (x0 - x1) - (z3 - z1) * (y0 - y1) * (x2 - x1);


					//WW  =x0*y2*z3+x2*y3*z0+x3*y0*z2-z0*y2*x3-z2*y3*x0-z3*y0*x2;
					if (WW) {

						a2 = (-x3 + x1) * 2;
						b2 = (-y3 + y1) * 2;
						c2 = (-z3 + z1) * 2;
						l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

						a3 = (-x0 + x1) * 2;
						b3 = (-y0 + y1) * 2;
						c3 = (-z0 + z1) * 2;
						l3 = -(x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1);

						a4 = (-x0 + x2) * 2;
						b4 = (-y0 + y2) * 2;
						c4 = (-z0 + z2) * 2;
						l4 = -(x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 - z2 * z2);

						W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
						W = 1 / W;
						l2 *= W;
						l3 *= W;
						l4 *= W;
						if (W) {
							X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
							Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
							Y = Y * -1;
							Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

							Rr =  (X - x0) * (X - x0) + (Y - y0) * (Y - y0) + (Z - z0) * (Z - z0);
							//Rr = X*X-2*X*x0+x0*x0 + Y*Y-2*Y*y0+y0*y0 + Z*Z-2*Z*z0+z0*z0;

							elements.setElement(e1, e2, e3, points.getIter() - 1, X, Y, Z, Rr, (x1 + x2 + x3 + x0) * 0.25, (y1 + y2 + y3 + y0) * 0.25, (z1 + z2 + z3 + z0) * 0.25);
							//elements.setElement(e1,e2,e3,points.getIter()-1,X,Y,Z,Rr,(x1+x2+x3+x0)*0.25,(y1+y2+y3+y0)*0.25,(z1+z2+z3+z0)*0.25);
						}
						else {
							delElement.podmienElement(i, -1);
						}
					}
					else {
						delElement.podmienElement(i, -1);
					}
				}
			}






			for (int a = 0; a < iter; a += 10) {

				if (delElement.getElement(a) != -1) {
					delElement.podmienElement(a, nrNowegoEl++);
				}

			}




			int iter2 = iter;

			for ( int a = iter - 10, i = 1; i < iter; a -= 10) {

				if (a >= 0) {
					if (delElement.getElement(a) != -1) {

						delElement.podmienElement(a, delElement.getElement(i));
						i += 40;

					}
				}
				else {
					iter2 += int(((elements.getIter() - iloscEStart) - iter * 0.025) * 40);

					break;
				}


			}

			nrNowegoEl -= elements.getIter();

//////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {

					int sasiad = -1;

					switch (delElement.getElement(i + 5)) {

						case 1:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE1();

							break;
						case 2:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE2();

							break;
						case 3:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE3();

							break;
						case 4:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE4();

							break;
					}

					//i+9 ostatnia sciana 123
					delElement.podmienElement(i + 9, sasiad);

					if (sasiad != -1) {
						//if(0){


						if (elements.getElement(sasiad).getP4() == delElement.getElement(i + 2)) {

							if (elements.getElement(sasiad).getP1() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE1(delElement.getElement(i));
							}
							else if (elements.getElement(sasiad).getP2() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE2(delElement.getElement(i));
							}
							else {
								elements.getElement(sasiad).setE3(delElement.getElement(i));
							}

						}

						else {
							elements.getElement(sasiad).setE4(delElement.getElement(i));
						}
					}


					for (int j = i + 10; j < iter; j += 10) {

						if (delElement.getElement(j) != -1) {

							if (delElement.getElement(i + 2) == delElement.getElement(j + 2)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 3)) {
								if (delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+3 = j+2");
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 4)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+4 = j+2");

								}
							}
							else {

								if (delElement.getElement(i + 3) == delElement.getElement(j + 2) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 3) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 4) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}

							}

						}
					}

				}

			}


			//usuwanie elementow
			//usuwanie wpisu elementu w tablicy
			for (int i = 1; i < iter2; i += 40) {
				//ww = delElement.getElement(i);
				//tabElementow2[(int(elementy.getElement(ww).getAx()*ws2))+TE2][(int(elementy.getElement(ww).getAy()*ws2))+TE2] = -1;
				//tabElementow2[(int(elements.getElement(ww).getAz()*ws2))+TE2][((int(elements.getElement(ww).getAy()*ws2))+TE2)*wymiarTE2 + (int(elements.getElement(ww).getAx()*ws2))+TE2] = -1;
				elements.deleteElement(delElement.getElement(i));
			}

			for (int i = 0, ww; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {
					ww = delElement.getElement(i);
					elements.getElement(ww).setElementySasiad(delElement.getElement(i + 6), delElement.getElement(i + 7), delElement.getElement(i + 8), delElement.getElement(i + 9));

					//tabElementow2[int(elements.getElement(ww).getAz()*ws2)+TE2][(int(elements.getElement(ww).getAy()*ws2)+TE2)*wymiarTE2 + int(elements.getElement(ww).getAx()*ws2)+TE2] = ww;
					//tabElementow1[int(elementy.getElement(delElement.getElement(i)).getAx()*ws1)][int(elementy.getElement(delElement.getElement(i)).getAy()*ws1)]++;

				}

			}

			for (int i = iter2 + 1; i < iter; i += 40) {

				elements.deleteElementPodmienS(delElement.getElement(i));


			}

		}
		else {

			--wszystkiePunktyPomoc;

			if (wszystkiePunktyPomoc > pun) {
				Tab.zamianaMiejsc(pun, wszystkiePunktyPomoc);
				--pun;
				++naliczanie4;
			}
			else if (ilePowtorzenPunktow > 0 && wszystkiePunktyPomoc < pun) {
				wszystkiePunktyPomoc = wszystkiePunkty - 1;
				--pun;
				--ilePowtorzenPunktow;
				//m1->Lines->Add(ilePowtorzenPunktow);
			}

		}
	}
//wymiarTE2++;
//for(int i=0;i<wymiarTE2;i++){delete []tabElementow2[i];}
//delete []tabElementow2;


//m1->Lines->Add("+++++++++");
	if (naliczanie4) {;}
//m1->Lines->Add(naliczanie2);
//m1->Lines->Add(naliczanie3-(naliczanie2+naliczanie1));


}

void Delanouy3D::poprawGranice() {
	sasiedniePunkty();

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double a2, b2, c2, l2, a3, b3, c3, l3, a4, b4, c4, l4, W;
	double  X, Y, Z;

	for (int i = 0, p, p1, p2, flaga, ilePG = wzorGranicy.getIter(); i < ilePG; i += 3) {

		p = wzorGranicy.getElement(i);
		p1 = wzorGranicy.getElement(i + 1);
		p2 = wzorGranicy.getElement(i + 2);
		flaga = 0;

		for (int j = 0, ileP = laplas[p].getIter(); j < ileP; ++j) {

			if (laplas[p].getElement(j) == p1) {
				++flaga;
			}
			if (laplas[p].getElement(j) == p2) {
				++flaga;
			}

		}

		if (flaga != 2) {

			x1 = points.getElement(p).getX();
			y1 = points.getElement(p).getY();
			z1 = points.getElement(p).getZ();

			x2 = points.getElement(p1).getX();
			y2 = points.getElement(p1).getY();
			z2 = points.getElement(p1).getZ();

			x3 = points.getElement(p2).getX();
			y3 = points.getElement(p2).getY();
			z3 = points.getElement(p2).getZ();

			a2 = (-x3 + x1) * 2;
			b2 = (-y3 + y1) * 2;
			c2 = (-z3 + z1) * 2;
			l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

			a3 = (-x2 + x1) * 2;
			b3 = (-y2 + y1) * 2;
			c3 = (-z2 + z1) * 2;
			l3 = -(x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + z2 * z2 - z1 * z1);

			a4 = y1 * z2 + y2 * z3 + y3 * z1 - z2 * y3 - z3 * y1 - z1 * y2;
			b4 = -(x1 * z2 + x2 * z3 + x3 * z1 - z2 * x3 - z3 * x1 - z1 * x2);
			c4 = x1 * y2 + x2 * y3 + x3 * y1 - y2 * x3 - y3 * x1 - y1 * x2;
			l4 = (x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - z1 * y2 * x3 - z2 * y3 * x1 - z3 * y1 * x2);

			W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
			W = 1 / W;
			l2 *= W;
			l3 *= W;
			l4 *= W;

			if (W) {

				X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
				Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
				Y = Y * -1;
				Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

				//points.setElement((x1+x2+x3)*0.33333,(y1+y2+y3)*0.33333,(z1+z2+z3)*0.33333,'g');

				if (X <= maxX && X >= 0 && Y <= maxY && Y >= 0 && Z <= maxZ && Z >= 0 ) {
					points.setElement(X, Y, Z, 'g');
				}
				else {

					//m1->Lines->Add("poprawa poza obszar");
					points.setElement((x1 + x2 + x3) * 0.33333, (y1 + y2 + y3) * 0.33333, (z1 + z2 + z3) * 0.33333, 'g');

				}

			}



		}

	}

}


//void Delanouy3D::usunieciePunktowKoloG(){



//}

void Delanouy3D::poprawGraniceUsuniecieP() {
	sasiedniePunkty();

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double a2, b2, c2, l2, a3, b3, c3, l3, a4, b4, c4, l4, W;
	double  X, Y, Z;
	IntList scianyDoPoprawy(10000, 10000);

	for (int i = 0, p, p1, p2, flaga, ilePG = wzorGranicy.getIter(); i < ilePG; i += 3) {

		p = wzorGranicy.getElement(i);
		p1 = wzorGranicy.getElement(i + 1);
		p2 = wzorGranicy.getElement(i + 2);
		flaga = 0;

		for (int j = 0, ileP = laplas[p].getIter(); j < ileP; ++j) {

			if (laplas[p].getElement(j) == p1) {
				++flaga;
			}
			if (laplas[p].getElement(j) == p2) {
				++flaga;
			}

		}

		if (flaga != 2) {

			scianyDoPoprawy.setElement(i);

			x1 = points.getElement(p).getX();
			y1 = points.getElement(p).getY();
			z1 = points.getElement(p).getZ();

			x2 = points.getElement(p1).getX();
			y2 = points.getElement(p1).getY();
			z2 = points.getElement(p1).getZ();

			x3 = points.getElement(p2).getX();
			y3 = points.getElement(p2).getY();
			z3 = points.getElement(p2).getZ();

			a2 = (-x3 + x1) * 2;
			b2 = (-y3 + y1) * 2;
			c2 = (-z3 + z1) * 2;
			l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

			a3 = (-x2 + x1) * 2;
			b3 = (-y2 + y1) * 2;
			c3 = (-z2 + z1) * 2;
			l3 = -(x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + z2 * z2 - z1 * z1);

			a4 = y1 * z2 + y2 * z3 + y3 * z1 - z2 * y3 - z3 * y1 - z1 * y2;
			b4 = -(x1 * z2 + x2 * z3 + x3 * z1 - z2 * x3 - z3 * x1 - z1 * x2);
			c4 = x1 * y2 + x2 * y3 + x3 * y1 - y2 * x3 - y3 * x1 - y1 * x2;
			l4 = (x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - z1 * y2 * x3 - z2 * y3 * x1 - z3 * y1 * x2);

			W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
			W = 1 / W;
			l2 *= W;
			l3 *= W;
			l4 *= W;

			if (W) {

				X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
				Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
				Y = Y * -1;
				Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

				//points.setElement((x1+x2+x3)*0.33333,(y1+y2+y3)*0.33333,(z1+z2+z3)*0.33333,'g');

				if (X <= maxX && X >= 0 && Y <= maxY && Y >= 0 && Z <= maxZ && Z >= 0 ) {
					points.setElement(X, Y, Z, 'g');
				}
				else {

					//m1->Lines->Add("poprawa poza obszar");
					points.setElement((x1 + x2 + x3) * 0.33333, (y1 + y2 + y3) * 0.33333, (z1 + z2 + z3) * 0.33333, 'g');

				}

			}



		}

	}
	/*
	for(int i=0,in=0,p,p1,p2,flaga,ilePG=scianyDoPoprawy.getIter();i<ilePG;++i){

		in = scianyDoPoprawy.getElement(i);

		p=wzorGranicy.getElement(in);
		p1=wzorGranicy.getElement(in+1);
		p2=wzorGranicy.getElement(in+2);

	}
	*/

}

void Delanouy3D::poprawGranice(PunktList & noweP) {

	noweP.ustawIter(0);
	sasiedniePunkty();


	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	for (int i = 0, p, p1, p2, flaga, ilePG = wzorGranicy.getIter(); i < ilePG; i += 3) {

		p = wzorGranicy.getElement(i);
		p1 = wzorGranicy.getElement(i + 1);
		p2 = wzorGranicy.getElement(i + 2);
		flaga = 0;

		for (int j = 0, ileP = laplas[p].getIter(); j < ileP; ++j) {

			if (laplas[p].getElement(j) == p1) {
				++flaga;
			}
			if (laplas[p].getElement(j) == p2) {
				++flaga;
			}

		}

		if (flaga != 2) {

			x1 = points.getElement(p).getX();
			y1 = points.getElement(p).getY();
			z1 = points.getElement(p).getZ();

			x2 = points.getElement(p1).getX();
			y2 = points.getElement(p1).getY();
			z2 = points.getElement(p1).getZ();

			x3 = points.getElement(p2).getX();
			y3 = points.getElement(p2).getY();
			z3 = points.getElement(p2).getZ();

			noweP.setElement((x1 + x2 + x3) * 0.33333, (y1 + y2 + y3) * 0.33333, (z1 + z2 + z3) * 0.33333, 'g');

		}

	}

}

void Delanouy3D::poprawGranice(IntList & wzorGranicy, PunktList & noweP) {

	noweP.ustawIter(0);
	sasiedniePunkty();


	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	for (int i = 0, p, p1, p2, flaga, ilePG = wzorGranicy.getIter(); i < ilePG; i += 3) {

		p = wzorGranicy.getElement(i);
		p1 = wzorGranicy.getElement(i + 1);
		p2 = wzorGranicy.getElement(i + 2);
		flaga = 0;

		for (int j = 0, ileP = laplas[p].getIter(); j < ileP; ++j) {

			if (laplas[p].getElement(j) == p1) {
				++flaga;
			}
			if (laplas[p].getElement(j) == p2) {
				++flaga;
			}

		}

		if (flaga != 2) {

			x1 = points.getElement(p).getX();
			y1 = points.getElement(p).getY();
			z1 = points.getElement(p).getZ();

			x2 = points.getElement(p1).getX();
			y2 = points.getElement(p1).getY();
			z2 = points.getElement(p1).getZ();

			x3 = points.getElement(p2).getX();
			y3 = points.getElement(p2).getY();
			z3 = points.getElement(p2).getZ();

			noweP.setElement((x1 + x2 + x3) * 0.33333, (y1 + y2 + y3) * 0.33333, (z1 + z2 + z3) * 0.33333, 'g');

		}

	}

}

void Delanouy3D::zapiszWzorGranicy() {

	wzorGranicy.ustawIter(0);

	int ileE = elements.getIter();
	bool * flagiE = new bool[ileE];
	for (int i = 0; i < ileE; ++i) {
		flagiE[i] = true;
	}


	for (int i = 0, ziarno, e1, e2, e3, e4, p1, p2, p3; i < ileE; ++i) {

		if (flagiE[i]) {

			ziarno = elements.getElement(i).getRodzajZiarna();

			e1 = elements.getElement(i).getE1();
			e2 = elements.getElement(i).getE2();
			e3 = elements.getElement(i).getE3();
			e4 = elements.getElement(i).getE4();

			if (e1 != -1 && elements.getElement(e1).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP4();
				p2 = elements.getElement(i).getP2();
				p3 = elements.getElement(i).getP1();
				flagiE[e1] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}

			if (e2 != -1 && elements.getElement(e2).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP4();
				p2 = elements.getElement(i).getP3();
				p3 = elements.getElement(i).getP2();
				flagiE[e2] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}
			if (e3 != -1 && elements.getElement(e3).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP4();
				p2 = elements.getElement(i).getP1();
				p3 = elements.getElement(i).getP3();
				flagiE[e3] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}
			if (e4 != -1 && elements.getElement(e4).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP1();
				p2 = elements.getElement(i).getP2();
				p3 = elements.getElement(i).getP3();
				flagiE[e4] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}

			//if
		}


	}





	delete[]flagiE;
}

void Delanouy3D::zapiszWzorGranicy(IntList & wzorGranicy) {

	wzorGranicy.ustawIter(0);

	int ileE = elements.getIter();
	bool * flagiE = new bool[ileE];
	for (int i = 0; i < ileE; ++i) {
		flagiE[i] = true;
	}


	for (int i = 0, ziarno, e1, e2, e3, e4, p1, p2, p3; i < ileE; ++i) {

		if (flagiE[i]) {

			ziarno = elements.getElement(i).getRodzajZiarna();

			e1 = elements.getElement(i).getE1();
			e2 = elements.getElement(i).getE2();
			e3 = elements.getElement(i).getE3();
			e4 = elements.getElement(i).getE4();

			if (e1 != -1 && elements.getElement(e1).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP4();
				p2 = elements.getElement(i).getP2();
				p3 = elements.getElement(i).getP1();
				flagiE[e1] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}

			if (e2 != -1 && elements.getElement(e2).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP4();
				p2 = elements.getElement(i).getP3();
				p3 = elements.getElement(i).getP2();
				flagiE[e2] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}
			if (e3 != -1 && elements.getElement(e3).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP4();
				p2 = elements.getElement(i).getP1();
				p3 = elements.getElement(i).getP3();
				flagiE[e3] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}
			if (e4 != -1 && elements.getElement(e4).getRodzajZiarna() != ziarno) {

				p1 = elements.getElement(i).getP1();
				p2 = elements.getElement(i).getP2();
				p3 = elements.getElement(i).getP3();
				flagiE[e4] = false;
				wzorGranicy.setElement(p1);
				wzorGranicy.setElement(p2);
				wzorGranicy.setElement(p3);

			}

			//if
		}


	}


	delete[]flagiE;
}



void Delanouy3D::optymalizacjaStosunekRdorNajlepszayWynikCalosc(double dlKroku, bool bezGranicy) {

	char znakG = 'g';
	if (!bezGranicy) {
		znakG = 'w';
	}
	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R = 0;
	double maxR = 0;
	int ktory;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() != znakG) {
			bool flaga = true;
			int razy = 1000;

			while (flaga && razy > 0) {

				tX = points.getElement(i).getX();
				tY = points.getElement(i).getY();
				tZ = points.getElement(i).getZ();

				//if(tX==0 || tX==maxX || tY==0 || tY==maxY || tZ==0 || tZ==maxZ){break;}

				--razy;
				maxR = stosunekRdorNajlepszayWynikDlaPunktu(i);
				ktory = 0;


				if (0 < tX - dlKroku && tX != maxX) {

					points.getElement(i).setPunktG(tX - dlKroku, tY, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 1;
					}
				}

				if (maxX > tX + dlKroku && tX != 0) {

					points.getElement(i).setPunktG(tX + dlKroku, tY, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 2;
					}
				}

				if (0 < tY - dlKroku && tY != maxY) {

					points.getElement(i).setPunktG(tX, tY - dlKroku, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 3;
					}
				}

				if (maxY > tY + dlKroku && tY != 0) {

					points.getElement(i).setPunktG(tX, tY + dlKroku, tZ);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 4;
					}
				}


				if (0 < tZ - dlKroku && tZ != maxZ) {

					points.getElement(i).setPunktG(tX, tY, tZ - dlKroku);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 5;
					}
				}

				if (maxZ > tZ + dlKroku && tZ != 0) {

					points.getElement(i).setPunktG(tX, tY, tZ + dlKroku);
					R = stosunekRdorNajlepszayWynikDlaPunktu(i);
					if (maxR > R) {
						ktory = 6;
					}
				}

				switch (ktory) {

					case 0:
						points.getElement(i).setPunkt(tX, tY, tZ);
						flaga = false;
						break;
					case 1:
						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						break;
					case 2:
						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						break;
					case 3:
						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						break;
					case 4:
						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						break;
					case 5:
						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						break;
					case 6:
						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						break;

				}

				//while
			}
		}
	}

}


void Delanouy3D::optymalizacjaStosunekRdorCalosc_zmienna_dl_kroku(double mnoznik_test, bool bezGranicy) {



	char znakG = 'g';
	if (!bezGranicy) {
		znakG = 'w';
	}
	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R = 0;
	double maxR;
	double dlKroku;
	int ktory;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() != znakG) {

			dlKroku = srednieRamiePunktu(i, mnoznik_test);
			//m1->Lines->Add(dlKroku);
			//if(1){
			bool flaga = true;
			int razy = 10000;

			while (flaga && razy > 0) {

				tX = points.getElement(i).getX();
				tY = points.getElement(i).getY();
				tZ = points.getElement(i).getZ();

				//if(tX==0 || tX==maxX || tY==0 || tY==maxY || tZ==0 || tZ==maxZ){break;}

				--razy;
				maxR = stosunekRdorNajwiekszyDlaPunktu(i);
				ktory = 0;

				////////

				//if(0<tX-dlKroku && tX!=maxX){

				points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
				R = stosunekRdorNajwiekszyDlaPunktu(i);
				if (maxR > R) {
					maxR = R;
					ktory = 1;
				}
				//}

				//if(maxX>tX+dlKroku && tX!=0){

				points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
				R = stosunekRdorNajwiekszyDlaPunktu(i);
				if (maxR > R) {
					maxR = R;
					ktory = 2;
				}
				//}

				//if(0<tY-dlKroku && tY!=maxY){

				points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
				R = stosunekRdorNajwiekszyDlaPunktu(i);
				if (maxR > R) {
					maxR = R;
					ktory = 3;
				}
				//}

				//if(maxY>tY+dlKroku && tY!=0){

				points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
				R = stosunekRdorNajwiekszyDlaPunktu(i);
				if (maxR > R) {
					maxR = R;
					ktory = 4;
				}
				//}

				///////////

				//if(0<tZ-dlKroku && tZ!=maxZ){

				points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
				R = stosunekRdorNajwiekszyDlaPunktu(i);
				if (maxR > R) {
					maxR = R;
					ktory = 5;
				}
				//}

				//if(maxZ>tZ+dlKroku && tZ!=0){

				points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
				R = stosunekRdorNajwiekszyDlaPunktu(i);
				if (maxR > R) {
					ktory = 6;
				}
				//}




				switch (ktory) {

					case 0:
						points.getElement(i).setPunkt(tX, tY, tZ);
						flaga = false;
						break;
					case 1:
						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						break;
					case 2:
						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						break;
					case 3:
						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						break;
					case 4:
						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						break;
					case 5:
						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						break;
					case 6:
						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						break;

				}

				//while
			}
		}
	}


}


double Delanouy3D::dlugosc_do_sciany(int punktID, int elID) {
	double x[5], y[5], z[5], wynik = 0;


	for (int i = 0, ii = 1; i < 4; i++) {
		if (punktID == elements.getElement(elID).getP(i)) {
			x[0] = points.getElement(elements.getElement(elID).getP(i)).getX();
			y[0] = points.getElement(elements.getElement(elID).getP(i)).getY();
			z[0] = points.getElement(elements.getElement(elID).getP(i)).getZ();

		}
		else {
			x[ii] = points.getElement(elements.getElement(elID).getP(i)).getX();
			y[ii] = points.getElement(elements.getElement(elID).getP(i)).getY();
			z[ii] = points.getElement(elements.getElement(elID).getP(i)).getZ();
			ii++;
		}

	}

	x[4] = (x[1] + x[2] + x[3]) / 3;
	y[4] = (y[1] + y[2] + y[3]) / 3;
	z[4] = (z[1] + z[2] + z[3]) / 3;



	return wynik = sqrt((x[0] - x[4]) * (x[0] - x[4]) + (y[0] - y[4]) * (y[0] - y[4]) + (z[0] - z[4]) * (z[0] - z[4]));

}

double Delanouy3D::srednieRamiePunktu(int punktID, double mnoznik) {

	double l_Patch = 10000000.0;
	//double l_Patch=0.0;
	double dlugosc;

	for (int i = 0, nrE, ile = laplas[punktID].getIter(); i < ile; ++i) {

		nrE = laplas[punktID].getElement(i);
		dlugosc = dlugosc_do_sciany(punktID, nrE);

		//l_Patch += dlugosc;

		if (l_Patch > dlugosc) {
			l_Patch = dlugosc;
		}


	}

	//return (l_Patch/laplas[punktID].getIter())*mnoznik;
	return l_Patch * mnoznik;

}

void Delanouy3D::optymalizacjaStosunekRdorCalosc(double dlKroku, bool bezGranicy) {



	char znakG = 'g';
	if (!bezGranicy) {
		znakG = 'w';
	}
	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R = 0;
	double maxR;
	int ktory;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() != znakG) {
			//if(1){
			bool flaga = true;
			int razy = 10000;

			while (flaga && razy > 0) {

				tX = points.getElement(i).getX();
				tY = points.getElement(i).getY();
				tZ = points.getElement(i).getZ();

				//if(tX==0 || tX==maxX || tY==0 || tY==maxY || tZ==0 || tZ==maxZ){break;}

				--razy;
				maxR = stosunekRdorNajwiekszyDlaPunktu(i);
				ktory = 0;

				////////

				if (0 < tX - dlKroku && tX != maxX) {

					points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 1;
					}
				}

				if (maxX > tX + dlKroku && tX != 0) {

					points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 2;
					}
				}

				if (0 < tY - dlKroku && tY != maxY) {

					points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 3;
					}
				}

				if (maxY > tY + dlKroku && tY != 0) {

					points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 4;
					}
				}

				///////////

				if (0 < tZ - dlKroku && tZ != maxZ) {

					points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 5;
					}
				}

				if (maxZ > tZ + dlKroku && tZ != 0) {

					points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						ktory = 6;
					}
				}




				switch (ktory) {

					case 0:
						points.getElement(i).setPunkt(tX, tY, tZ);
						flaga = false;
						break;
					case 1:
						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						break;
					case 2:
						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						break;
					case 3:
						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						break;
					case 4:
						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						break;
					case 5:
						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						break;
					case 6:
						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						break;

				}

				//while
			}
		}
	}


}


double Delanouy3D::stosunekRdorNajlepszayWynikDlaPunktu(int ktoryP) {

	double V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, nx, ny, nz;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;
//double maxX = dx;
	double maxR = 0;
	double maxV = 0;

	for (int i = 0, ktoryE, ileE = laplas[ktoryP].getIter(), p1, p2, p3, p4; i < ileE; ++i) {

		ktoryE = laplas[ktoryP].getElement(i);

		p1 = elements.getElement(ktoryE).getP1();
		p2 = elements.getElement(ktoryE).getP2();
		p3 = elements.getElement(ktoryE).getP3();
		p4 = elements.getElement(ktoryE).getP4();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);
		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);
		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );


		pol = 0;
		Ax = x2 - x1;
		Bx = x3 - x1;
		Cx = x4 - x1;
		Ay = y2 - y1;
		By = y3 - y1;
		Cy = y4 - y1;
		Az = z2 - z1;
		Bz = z3 - z1;
		Cz = z4 - z1;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = Ay * Cz - Cy * Az;
		ny = Cx * Az - Ax * Cz;
		nz = Ax * Cy - Cx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = By * Cz - Cy * Bz;
		ny = Cx * Bz - Bx * Cz;
		nz = Bx * Cy - Cx * By;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		V = (nx * Ax + ny * Ay + nz * Az) * 0.5;

		Ax = x2 - x4;
		Bx = x3 - x4;
		Ay = y2 - y4;
		By = y3 - y4;
		Az = z2 - z4;
		Bz = z3 - z4;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;


		r = (V) / pol;

		if (V <= 0) {
			V = (V * -1) + 100000;
			if (maxV < V) {
				maxV = V;
			}
		}

		R = R / r;

		//if(maxR<R){maxR=R;}
		maxR += R;
	}

	if (maxV) {
		return maxV;
	}
	return maxR;

}

/*
double Delanouy3D::stosunekRdorNajlepszayWynikDlaPunktu(int ktoryP){

double V,r,R,pol;
double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,nx,ny,nz;
double a2,a3,a4,b2,b3,b4,c2,c3,c4,l2,l3,l4,X,Y,Z,W;
//double maxX = dx;
double maxR=0;
double maxV=0;

	for(int i=0,ktoryE,ileE=laplas[ktoryP].getIter(),p1,p2,p3,p4;i<ileE;++i){



		ktoryE=laplas[ktoryP].getElement(i);

		p1=elements.getElement(ktoryE).getP1();
		p2=elements.getElement(ktoryE).getP2();
		p3=elements.getElement(ktoryE).getP3();
		p4=elements.getElement(ktoryE).getP4();

		x1=points.getElement(p1).getX();
		x2=points.getElement(p2).getX();
		x3=points.getElement(p3).getX();
		x4=points.getElement(p4).getX();
		y1=points.getElement(p1).getY();
		y2=points.getElement(p2).getY();
		y3=points.getElement(p3).getY();
		y4=points.getElement(p4).getY();
		z1=points.getElement(p1).getZ();
		z2=points.getElement(p2).getZ();
		z3=points.getElement(p3).getZ();
		z4=points.getElement(p4).getZ();


		a2 = (-x3+x1)*2;b2 = (-y3+y1)*2;c2 = (-z3+z1)*2;l2 = -(x3*x3-x1*x1 + y3*y3-y1*y1 + z3*z3-z1*z1);
		a3 = (-x4+x1)*2;b3 = (-y4+y1)*2;c3 = (-z4+z1)*2;l3 = -(x4*x4-x1*x1 + y4*y4-y1*y1 + z4*z4-z1*z1);
		a4 = (-x4+x2)*2;b4 = (-y4+y2)*2;c4 = (-z4+z2)*2;l4 = -(x4*x4-x2*x2 + y4*y4-y2*y2 + z4*z4-z2*z2);

		W =  a2*b3*c4+a3*b4*c2+a4*b2*c3-c2*b3*a4-c3*b4*a2-c4*b2*a3;
		W=1/W;
		l2*=W;l3*=W;l4*=W;

		X = (b2*c3*l4+b3*c4*l2+b4*c2*l3-l2*c3*b4-l3*c4*b2-l4*c2*b3);
		Y = (a2*c3*l4+a3*c4*l2+a4*c2*l3-l2*c3*a4-l3*c4*a2-l4*c2*a3);Y=Y*-1;
		Z = (a2*b3*l4+a3*b4*l2+a4*b2*l3-l2*b3*a4-l3*b4*a2-l4*b2*a3);

		R =  sqrt( (X-x4)*(X-x4) + (Y-y4)*(Y-y4) + (Z-z4)*(Z-z4) );


		pol=0;
		Ax=x2-x1; Bx=x3-x1; Cx=x4-x1;
		Ay=y2-y1; By=y3-y1; Cy=y4-y1;
		Az=z2-z1; Bz=z3-z1; Cz=z4-z1;

		nx=Ay*Bz-By*Az; ny=Bx*Az-Ax*Bz; nz=Ax*By-Bx*Ay;
		pol += sqrt(nx*nx+ny*ny+nz*nz)*0.5;

		nx=Ay*Cz-Cy*Az; ny=Cx*Az-Ax*Cz; nz=Ax*Cy-Cx*Ay;
		pol += sqrt(nx*nx+ny*ny+nz*nz)*0.5;

		nx=By*Cz-Cy*Bz; ny=Cx*Bz-Bx*Cz; nz=Bx*Cy-Cx*By;
		pol += sqrt(nx*nx+ny*ny+nz*nz)*0.5;

		V = (nx*Ax+ny*Ay+nz*Az)*0.5;

		Ax=x2-x4; Bx=x3-x4;
		Ay=y2-y4; By=y3-y4;
		Az=z2-z4; Bz=z3-z4;

		nx=Ay*Bz-By*Az; ny=Bx*Az-Ax*Bz; nz=Ax*By-Bx*Ay;
		pol += sqrt(nx*nx+ny*ny+nz*nz)*0.5;


		r = (V)/pol;

		if(V<=0){V = (V*-1)+100000; if(maxV<V){maxV=V;}}

		R= R/r;

		//if(maxR<R){maxR=R;}
		maxR+=R;
	}

if(maxV){return maxV;}



return maxR;

}
*/
double Delanouy3D::stosunekRdorNajlepszayWynikDlaPunktu1(int ktoryP) {

	double V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, nx, ny, nz;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;
//double maxX = dx;
	double maxR = 0;
	double maxV = 0;

	for (int i = 0, ktoryE = 0, ileE = laplas[ktoryP].getIter(), p1, p2, p3, p4; i < ileE; ++i) {

		if (i == 10000) {
			cout << "ble1";
		}



		ktoryE = laplas[ktoryP].getElement(i);

		p1 = elements.getElement(ktoryE).getP1();
		p2 = elements.getElement(ktoryE).getP2();
		p3 = elements.getElement(ktoryE).getP3();
		p4 = elements.getElement(ktoryE).getP4();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);
		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);
		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );


		pol = 0;
		Ax = x2 - x1;
		Bx = x3 - x1;
		Cx = x4 - x1;
		Ay = y2 - y1;
		By = y3 - y1;
		Cy = y4 - y1;
		Az = z2 - z1;
		Bz = z3 - z1;
		Cz = z4 - z1;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = Ay * Cz - Cy * Az;
		ny = Cx * Az - Ax * Cz;
		nz = Ax * Cy - Cx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = By * Cz - Cy * Bz;
		ny = Cx * Bz - Bx * Cz;
		nz = Bx * Cy - Cx * By;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		V = (nx * Ax + ny * Ay + nz * Az) * 0.5;

		Ax = x2 - x4;
		Bx = x3 - x4;
		Ay = y2 - y4;
		By = y3 - y4;
		Az = z2 - z4;
		Bz = z3 - z4;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;


		r = (V) / pol;

		if (V <= 0) {
			V = (V * -1) + 100000;
			if (maxV < V) {
				maxV = V;
			}
		}

		R = R / r;

		//if(maxR<R){maxR=R;}

		maxR += R;


	}

	if (maxV) {
		return maxV;
	}


	return maxR;

}




void Delanouy3D::optymalizacjaStosunekRdor(double dlKroku) {


	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R = 0;
	double maxR;
	int ktory;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() != 'g') {
			bool flaga = true;
			int razy = 10000;

			while (flaga && razy > 0) {

				tX = points.getElement(i).getX();
				tY = points.getElement(i).getY();
				tZ = points.getElement(i).getZ();

				if (tX == 0 || tX == maxX || tY == 0 || tY == maxY || tZ == 0 || tZ == maxZ) {
					break;
				}

				--razy;


				maxR = stosunekRdorNajwiekszyDlaPunktu(i);
				ktory = 0;

				if (0 < tX - dlKroku) {

					points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 1;
					}
				}

				if (maxX > tX + dlKroku) {

					points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 2;
					}
				}

				if (0 < tY - dlKroku) {

					points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 3;
					}
				}

				if (maxY > tY + dlKroku) {

					points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 4;
					}
				}

				if (0 < tZ - dlKroku) {

					points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 5;
					}
				}

				if (maxZ > tZ + dlKroku) {

					points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
					R = stosunekRdorNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						ktory = 6;
					}
				}

				switch (ktory) {

					case 0:
						points.getElement(i).setPunkt(tX, tY, tZ);
						flaga = false;
						break;
					case 1:
						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						break;
					case 2:
						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						break;
					case 3:
						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						break;
					case 4:
						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						break;
					case 5:
						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						break;
					case 6:
						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						break;

				}

				//while
			}
		}
	}


}



void Delanouy3D::optymalizacjaStosunekGMSH(double dlKroku) {

	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R;
	double maxR;
	int ktory;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() != 'g') {
			bool flaga = true;
			int razy = 10000;



			while (flaga && razy > 0) {

				tX = points.getElement(i).getX();
				tY = points.getElement(i).getY();
				tZ = points.getElement(i).getZ();

				if (tX == 0 || tX == maxX || tY == 0 || tY == maxY || tZ == 0 || tZ == maxZ) {
					break;
				}

				--razy;

				maxR = stosunekGMSHNajwiekszyDlaPunktu(i);
				ktory = 0;

				if (0 < tX - dlKroku) {

					points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
					R = stosunekGMSHNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 1;
					}
				}

				if (maxX > tX + dlKroku) {

					points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
					R = stosunekGMSHNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 2;
					}
				}

				if (0 < tY - dlKroku) {

					points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
					R = stosunekGMSHNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 3;
					}
				}

				if (maxY > tY + dlKroku) {

					points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
					R = stosunekGMSHNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 4;
					}
				}

				if (0 < tZ - dlKroku) {

					points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
					R = stosunekGMSHNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						maxR = R;
						ktory = 5;
					}
				}

				if (maxZ > tZ + dlKroku) {

					points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
					R = stosunekGMSHNajwiekszyDlaPunktu(i);
					if (maxR > R) {
						ktory = 6;
					}
				}

				switch (ktory) {

					case 0:
						points.getElement(i).setPunkt(tX, tY, tZ);
						flaga = false;
						break;
					case 1:
						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						break;
					case 2:
						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						break;
					case 3:
						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						break;
					case 4:
						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						break;
					case 5:
						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						break;
					case 6:
						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						break;

				}

				//while
			}
		}
	}


}

void Delanouy3D::optymalizacjaV(double dlKroku) {


	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double V;
	double stosunekV;
	int ktory;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() != 'g') {
			bool flaga = true;
			int razy = 10000;



			while (flaga && razy > 0) {

				tX = points.getElement(i).getX();
				tY = points.getElement(i).getY();
				tZ = points.getElement(i).getZ();

				if (tX == 0 || tX == maxX || tY == 0 || tY == maxY || tZ == 0 || tZ == maxZ) {
					break;
				}

				--razy;

				stosunekV = stosunekVDlaPunktu(i);
				ktory = 0;

				if (0 < tX - dlKroku) {

					points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
					V = stosunekVDlaPunktu(i);
					if (stosunekV > V) {
						stosunekV = V;
						ktory = 1;
					}
				}

				if (maxX > tX + dlKroku) {

					points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
					V = stosunekVDlaPunktu(i);
					if (stosunekV > V) {
						stosunekV = V;
						ktory = 2;
					}
				}

				if (0 < tY - dlKroku) {

					points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
					V = stosunekVDlaPunktu(i);
					if (stosunekV > V) {
						stosunekV = V;
						ktory = 3;
					}
				}

				if (maxY > tY + dlKroku) {

					points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
					V = stosunekVDlaPunktu(i);
					if (stosunekV > V) {
						stosunekV = V;
						ktory = 4;
					}
				}

				if (0 < tZ - dlKroku) {

					points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
					V = stosunekVDlaPunktu(i);
					if (stosunekV > V) {
						stosunekV = V;
						ktory = 5;
					}
				}

				if (maxZ > tZ + dlKroku) {

					points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
					V = stosunekVDlaPunktu(i);
					if (stosunekV > V) {
						ktory = 6;
					}
				}

				switch (ktory) {

					case 0:
						points.getElement(i).setPunkt(tX, tY, tZ);
						flaga = false;
						break;
					case 1:
						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						break;
					case 2:
						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						break;
					case 3:
						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						break;
					case 4:
						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						break;
					case 5:
						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						break;
					case 6:
						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						break;

				}

				//while
			}
		}
	}


}

void Delanouy3D::optymalizacjaStosunekRdorDODelanoya(double dlKroku, int odElementu) {

//sasiednieElementy();
	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R;
	double maxR;
	int ktory;


	int i = points.getIter() - 1;

	if (points.getElement(i).getGranica() != 'g') {
		bool flaga = true;
		int razy = 10000;

		while (flaga && razy > 0) {

			tX = points.getElement(i).getX();
			tY = points.getElement(i).getY();
			tZ = points.getElement(i).getZ();

			if (tX == 0 || tX == maxX || tY == 0 || tY == maxY || tZ == 0 || tZ == maxZ) {
				break;
			}

			--razy;

			maxR = stosunekRdorNajwiekszyDlaPunktuDODelanoya(i, odElementu);
			ktory = 0;

			if (0 < tX - dlKroku) {

				points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
				R = stosunekRdorNajwiekszyDlaPunktuDODelanoya(i, odElementu);
				if (maxR > R) {
					maxR = R;
					ktory = 1;
				}
			}

			if (maxX > tX + dlKroku) {

				points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
				R = stosunekRdorNajwiekszyDlaPunktuDODelanoya(i, odElementu);
				if (maxR > R) {
					maxR = R;
					ktory = 2;
				}
			}

			if (0 < tY - dlKroku) {

				points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
				R = stosunekRdorNajwiekszyDlaPunktuDODelanoya(i, odElementu);
				if (maxR > R) {
					maxR = R;
					ktory = 3;
				}
			}

			if (maxY > tY + dlKroku) {

				points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
				R = stosunekRdorNajwiekszyDlaPunktuDODelanoya(i, odElementu);
				if (maxR > R) {
					maxR = R;
					ktory = 4;
				}
			}

			if (0 < tZ - dlKroku) {

				points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
				R = stosunekRdorNajwiekszyDlaPunktuDODelanoya(i, odElementu);
				if (maxR > R) {
					maxR = R;
					ktory = 5;
				}
			}

			if (maxZ > tZ + dlKroku) {

				points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
				R = stosunekRdorNajwiekszyDlaPunktuDODelanoya(i, odElementu);
				if (maxR > R) {
					ktory = 6;
				}
			}

			switch (ktory) {

				case 0:
					points.getElement(i).setPunkt(tX, tY, tZ);
					flaga = false;
					break;
				case 1:
					points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
					break;
				case 2:
					points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
					break;
				case 3:
					points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
					break;
				case 4:
					points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
					break;
				case 5:
					points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
					break;
				case 6:
					points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
					break;

			}

			//while
		}
	}

}

void Delanouy3D::optymalizacjaStosunekRdorWybranePunkty(double dlKroku, IntList & wybraneP, int ileRazy) {

	sasiednieElementy();
	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double R;
	double maxR;
	int ktory;

	for (int ile = 0; ile < ileRazy; ++ile) {

		for (int ip = 0, ileP = wybraneP.getIter(), i; ip < ileP; ++ip) {
			i = wybraneP.getElement(ip);
			if (points.getElement(i).getGranica() != 'g') {
				bool flaga = true;
				int razy = 10000;

				while (flaga && razy > 0) {

					tX = points.getElement(i).getX();
					tY = points.getElement(i).getY();
					tZ = points.getElement(i).getZ();

					if (tX == 0 || tX == maxX || tY == 0 || tY == maxY || tZ == 0 || tZ == maxZ) {
						break;
					}

					--razy;

					maxR = stosunekRdorNajwiekszyDlaPunktu(i);
					ktory = 0;

					if (0 < tX - dlKroku) {

						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						R = stosunekRdorNajwiekszyDlaPunktu(i);
						if (maxR > R) {
							maxR = R;
							ktory = 1;
						}
					}

					if (maxX > tX + dlKroku) {

						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						R = stosunekRdorNajwiekszyDlaPunktu(i);
						if (maxR > R) {
							maxR = R;
							ktory = 2;
						}
					}

					if (0 < tY - dlKroku) {

						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						R = stosunekRdorNajwiekszyDlaPunktu(i);
						if (maxR > R) {
							maxR = R;
							ktory = 3;
						}
					}

					if (maxY > tY + dlKroku) {

						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						R = stosunekRdorNajwiekszyDlaPunktu(i);
						if (maxR > R) {
							maxR = R;
							ktory = 4;
						}
					}

					if (0 < tZ - dlKroku) {

						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						R = stosunekRdorNajwiekszyDlaPunktu(i);
						if (maxR > R) {
							maxR = R;
							ktory = 5;
						}
					}

					if (maxZ > tZ + dlKroku) {

						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						R = stosunekRdorNajwiekszyDlaPunktu(i);
						if (maxR > R) {
							ktory = 6;
						}
					}

					switch (ktory) {

						case 0:
							points.getElement(i).setPunkt(tX, tY, tZ);
							flaga = false;
							break;
						case 1:
							points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
							break;
						case 2:
							points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
							break;
						case 3:
							points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
							break;
						case 4:
							points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
							break;
						case 5:
							points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
							break;
						case 6:
							points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
							break;

					}

					//while
				}
			}
		}
	}

}



void Delanouy3D::optymalizacjaVWybranePunkty(double dlKroku, IntList & wybraneP, int ileRazy) {

	sasiednieElementy();
	double tX, tY, tZ;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double V;
	double stosunekV;
	int ktory;

	for (int ile = 0; ile < ileRazy; ++ile) {

		for (int ip = 0, ileP = wybraneP.getIter(), i; ip < ileP; ++ip) {
			i = wybraneP.getElement(ip);
			if (points.getElement(i).getGranica() != 'g') {
				bool flaga = true;
				int razy = 10000;

				while (flaga && razy > 0) {

					tX = points.getElement(i).getX();
					tY = points.getElement(i).getY();
					tZ = points.getElement(i).getZ();

					if (tX == 0 || tX == maxX || tY == 0 || tY == maxY || tZ == 0 || tZ == maxZ) {
						break;
					}

					--razy;

					stosunekV = stosunekVDlaPunktu(i);
					ktory = 0;

					if (0 < tX - dlKroku) {

						points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
						V = stosunekVDlaPunktu(i);
						if (stosunekV > V) {
							stosunekV = V;
							ktory = 1;
						}
					}

					if (maxX > tX + dlKroku) {

						points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
						V = stosunekVDlaPunktu(i);
						if (stosunekV > V) {
							stosunekV = V;
							ktory = 2;
						}
					}

					if (0 < tY - dlKroku) {

						points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
						V = stosunekVDlaPunktu(i);
						if (stosunekV > V) {
							stosunekV = V;
							ktory = 3;
						}
					}

					if (maxY > tY + dlKroku) {

						points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
						V = stosunekVDlaPunktu(i);
						if (stosunekV > V) {
							stosunekV = V;
							ktory = 4;
						}
					}

					if (0 < tZ - dlKroku) {

						points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
						V = stosunekVDlaPunktu(i);
						if (stosunekV > V) {
							stosunekV = V;
							ktory = 5;
						}
					}

					if (maxZ > tZ + dlKroku) {

						points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
						V = stosunekVDlaPunktu(i);
						if (stosunekV > V) {
							ktory = 6;
						}
					}

					switch (ktory) {

						case 0:
							points.getElement(i).setPunkt(tX, tY, tZ);
							flaga = false;
							break;
						case 1:
							points.getElement(i).setPunkt(tX - dlKroku, tY, tZ);
							break;
						case 2:
							points.getElement(i).setPunkt(tX + dlKroku, tY, tZ);
							break;
						case 3:
							points.getElement(i).setPunkt(tX, tY - dlKroku, tZ);
							break;
						case 4:
							points.getElement(i).setPunkt(tX, tY + dlKroku, tZ);
							break;
						case 5:
							points.getElement(i).setPunkt(tX, tY, tZ - dlKroku);
							break;
						case 6:
							points.getElement(i).setPunkt(tX, tY, tZ + dlKroku);
							break;

					}

					//while
				}
			}
		}
	}

}

void Delanouy3D::wyszukajPunktuDoWygladzenia(IntList & punk, int numerZiarna, int ilePok) {
	punk.ustawIter(0);
	sasiedniePunkty();

	for (int i = 0, ileE = elements.getIter(), p1, p2, p3, p4; i < ileE; ++i) {
		if (elements.getElement(i).getRodzajZiarna() == numerZiarna) {

			p1 = elements.getElement(i).getP1();
			p2 = elements.getElement(i).getP2();
			p3 = elements.getElement(i).getP3();
			p4 = elements.getElement(i).getP4();

			if (points.getElement(p1).getGranica() == 'g') {
				points.getElement(p1).setGranica('r');
			}
			else if (points.getElement(p1).getGranica() == 'a') {
				points.getElement(p1).setGranica('m');
			}
			if (points.getElement(p2).getGranica() == 'g') {
				points.getElement(p2).setGranica('r');
			}
			else if (points.getElement(p2).getGranica() == 'a') {
				points.getElement(p2).setGranica('m');
			}
			if (points.getElement(p3).getGranica() == 'g') {
				points.getElement(p3).setGranica('r');
			}
			else if (points.getElement(p3).getGranica() == 'a') {
				points.getElement(p3).setGranica('m');
			}
			if (points.getElement(p4).getGranica() == 'g') {
				points.getElement(p4).setGranica('r');
			}
			else if (points.getElement(p4).getGranica() == 'a') {
				points.getElement(p4).setGranica('m');
			}

		}

	}


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'r') {

			for (int j = 0, ilePP = laplas[i].getIter(), p; j < ilePP; ++j) {

				p = laplas[i].getElement(j);
				if (points.getElement(p).getGranica() == 'a') {
					points.getElement(p).setGranica('m');
					punk.setElement(p);
				}

			}

		}
	}

	int i = 0, ii; //ii stara pozycja i uzywana pozycja z poczatku 0
	if (ilePok < 0) {
		ilePok = 0;
	}
	while (ilePok--) {
		ii = punk.getIter();

		for (int p0, ileP = punk.getIter(); i < ileP; ++i) {

			p0 = punk.getElement(i);
			for (int j = 0, ilePP = laplas[p0].getIter(), p1; j < ilePP; ++j) {

				p1 = laplas[p0].getElement(j);
				if (points.getElement(p1).getGranica() == 'a') {
					points.getElement(p1).setGranica('m');
					punk.setElement(p1);
				}

			}


		}

		i = ii;
	}


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() == 'm') {
			points.getElement(i).setGranica('a');
		}
		if (points.getElement(i).getGranica() == 'r') {
			points.getElement(i).setGranica('g');
		}
	}

}

void Delanouy3D::uwtorzWartstwePrzysciennaDlaZiarna(int nrZ, double dlOdsuniecia) {

	//sciany.czysc(1000,1000);

	//znajdzScianyGraniczneDlaZiarna(sciany,nrZ,false);
	//obliczGradiety(sciany,dlOdsuniecia);

	/*
	zmianaWartGrad(sciany);
	zmianaWartGrad(sciany);
	zmianaWartGrad(sciany);
	zmianaWartGrad(sciany);
	zmianaWartGrad(sciany);
	*/

	//wyznaczOdleglosc(sciany,dlOdsuniecia);
	//tworzenieWarstwyPrzysciennej(sciany);
}

/*
void Delanouy3D::znajdzScianyGraniczneDlaZiarna(IntList &sciany,int nrZiarna,bool czyKrawedz){

	for(int i=0,e1,e2,e3,e4,ileE=elementy.getIter();i<ileE;++i){


		if(elementy.getElement(i).getRodzajZiarna()==nrZiarna){

			e1=elementy.getElement(i).getE1();
			e2=elementy.getElement(i).getE2();
			e3=elementy.getElement(i).getE3();
			e4=elementy.getElement(i).getE4();

			if(e1!=-1){
				if(elementy.getElement(e1).getRodzajZiarna()!=nrZiarna){
					sciany.setElement(i);
					sciany.setElement(1);
				}
			}
			else{
				if(czyKrawedz){
					sciany.setElement(i);
					sciany.setElement(1);
				}
			}

			if(e2!=-1){
				if(elementy.getElement(e2).getRodzajZiarna()!=nrZiarna){
					sciany.setElement(i);
					sciany.setElement(2);
				}
			}
			else{
				if(czyKrawedz){
					sciany.setElement(i);
					sciany.setElement(2);
				}
			}

			if(e3!=-1){
				if(elementy.getElement(e3).getRodzajZiarna()!=nrZiarna){
					sciany.setElement(i);
					sciany.setElement(3);
				}
			}
			else{
				if(czyKrawedz){
					sciany.setElement(i);
					sciany.setElement(3);
				}
			}

		}


	}



}


void Delanouy3D::rozpoznawaniePlaszczyzn(){

wybraneZleEl.ustawIter(0);
int bardzozle=0;

//ustalenie elementu zerowego
int elZero=0;

for(int w=0,ileE=elements.getIter();w<ileE;w++){
	if(elements.getElement(w).getP2()==0 || elements.getElement(w).getP3()==0){elZero=w;}
	elements.getElement(w).setRodzajZiarna(-1);
}

int numerZ=0;

IntList runda(1450,500);
IntList pomoc(1300,500);
elZero=0;
for(int w=0,e1,e2,e3,e4,ileE=elements.getIter();w<ileE;w++){


	if(elements.getElement(w).getRodzajZiarna()==-1){
	runda.ustawIter(0);

	e1 = elements.getElement(w).getE1();
	e2 = elements.getElement(w).getE2();
	e3 = elements.getElement(w).getE3();
	e4 = elements.getElement(w).getE4();

   if( points.getElement(p1).getGranica()!='a' && points.getElement(p2).getGranica()!='a' && points.getElement(p3).getGranica()!='a' && points.getElement(p4).getGranica()!='a'){
   bardzozle++;wybraneZleEl.setElement(w);
   }
   else{ runda.setElement(w);numerZ++; }

	while(runda.getIter()){

		for(int i=0,el;i<runda.getIter();i++){
			el = runda.getElement(i);

			if(elements.getElement(el).getRodzajZiarna()==-1){
				elements.getElement(el).setRodzajZiarna(numerZ);

				int n1=0,n2=0,n3=0,n4=0;

				p1 = elements.getElement(el).getP1();
				p2 = elements.getElement(el).getP2();
				p3 = elements.getElement(el).getP3();
				p4 = elements.getElement(el).getP4();

				if(points.getElement(p1).getGranica()!='a'){n1=1;}
				if(points.getElement(p2).getGranica()!='a'){n2=1;}
				if(points.getElement(p3).getGranica()!='a'){n3=1;}
				if(points.getElement(p4).getGranica()!='a'){n4=1;}


				if(n4+n2+n1<3){
					if(elements.getElement(el).getE1()!=-1){pomoc.setElement(elements.getElement(el).getE1());}
				}
				if(n4+n3+n2<3){
					if(elements.getElement(el).getE2()!=-1){pomoc.setElement(elements.getElement(el).getE2());}
				}
				if(n4+n1+n3<3){
					if(elements.getElement(el).getE3()!=-1){pomoc.setElement(elements.getElement(el).getE3());}
				}
				if(n1+n2+n3<3){
					if(elements.getElement(el).getE4()!=-1){pomoc.setElement(elements.getElement(el).getE4());}
				}
			}

		}

		runda.ustawIter(0);

		for(int i=0;i<pomoc.getIter();i++){runda.setElement(pomoc.getElement(i));}

		pomoc.ustawIter(0);

	//koniec while
	}

	}

elZero=w;
//koniec for
}


IntList zleE(10000,5000);
int zlychE=0;
//m1->Lines->Add(numerZ);

///

 for(int i=0;i<elements.getIter();i++){
	if(elements.getElement(i).getRodzajZiarna()==-1){zleE.setElement(i);}
 }
///

//poprawaDopasowaniaZiaren2Sas(wybraneZleEl);
//poprawaDopasowaniaZiarenNajPow();
//poprawaDopasowaniaZiaren2Sas(wybraneZleEl);
//poprawaDopasowaniaZiarenNajPow();
//poprawaDopasowaniaZiaren2Sas(wybraneZleEl);



 for(int i=0;i<elements.getIter();i++){
	if(elements.getElement(i).getRodzajZiarna()==-1){++zlychE;}
 }

m1->Lines->Add(bardzozle);
m1->Lines->Add(zlychE);

return numerZ;

}
*/

double Delanouy3D::ileWezlowG() {
	int zlicz = 0;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() == 'g') {
			++zlicz;
		}

	}

	return zlicz;
}


void Delanouy3D::oznaczWezlyGraniczneNaPodstawieScianG() {

	for (int i = 0, p1, p2, p3, ileE = elements.getIter(); i < ileE; ++i) {

		if (elements.getElement(i).getE1() <= -1) {
			p1 = elements.getElement(i).getP4();
			p2 = elements.getElement(i).getP2();
			p3 = elements.getElement(i).getP1();

			points.getElement(p1).setGranica('r');
			points.getElement(p2).setGranica('r');
			points.getElement(p3).setGranica('r');
		}

		if (elements.getElement(i).getE2() <= -1) {
			p1 = elements.getElement(i).getP4();
			p2 = elements.getElement(i).getP3();
			p3 = elements.getElement(i).getP2();

			points.getElement(p1).setGranica('r');
			points.getElement(p2).setGranica('r');
			points.getElement(p3).setGranica('r');
		}

		if (elements.getElement(i).getE3() <= -1) {
			p1 = elements.getElement(i).getP4();
			p2 = elements.getElement(i).getP1();
			p3 = elements.getElement(i).getP3();

			points.getElement(p1).setGranica('r');
			points.getElement(p2).setGranica('r');
			points.getElement(p3).setGranica('r');
		}

		if (elements.getElement(i).getE4() <= -1) {
			p1 = elements.getElement(i).getP1();
			p2 = elements.getElement(i).getP2();
			p3 = elements.getElement(i).getP3();

			points.getElement(p1).setGranica('r');
			points.getElement(p2).setGranica('r');
			points.getElement(p3).setGranica('r');
		}

	}

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() == 'r') {
			points.getElement(i).setGranica('g');
		}
		else {
			points.getElement(i).setGranica('a');
		}

	}

}

void Delanouy3D::wyszukanieSasidniegoElementu(int nrEl, int numerSciany) {


	IntList wspEl(1000, 1000);

	int p1 = -1, p2 = -1, p3 = -1;
	switch (numerSciany) {

		case 1:

			p1 = elements.getElement(nrEl).getP4();
			p2 = elements.getElement(nrEl).getP2();
			p3 = elements.getElement(nrEl).getP1();
			elements.getElement(nrEl).setE1(-1);
			break;
		case 2:

			p1 = elements.getElement(nrEl).getP4();
			p2 = elements.getElement(nrEl).getP3();
			p3 = elements.getElement(nrEl).getP2();
			elements.getElement(nrEl).setE2(-1);
			break;
		case 3:

			p1 = elements.getElement(nrEl).getP4();
			p2 = elements.getElement(nrEl).getP1();
			p3 = elements.getElement(nrEl).getP3();
			elements.getElement(nrEl).setE3(-1);
			break;
		case 4:

			p1 = elements.getElement(nrEl).getP1();
			p2 = elements.getElement(nrEl).getP2();
			p3 = elements.getElement(nrEl).getP3();
			elements.getElement(nrEl).setE4(-1);
			break;

	}

	for (int j = 0, p, ileSas = laplas[p1].getIter(); j < ileSas; ++j) {

		p =  laplas[p1].getElement(j);

		for (int jj = 0, ileSas2 = laplas[p2].getIter(); jj < ileSas2; ++j) {
			if (p == laplas[p2].getElement(jj)) {
				wspEl.setElement(p);
				break;
			}
		}
	}

	for (int j = 0, p, ileSas = wspEl.getIter(); j < ileSas; ++j) {

		p =  wspEl.getElement(j);

		for (int jj = 0, ileSas2 = laplas[p3].getIter(); jj < ileSas2; ++j) {
			if (p == laplas[p3].getElement(jj) && p != nrEl) {


				switch (numerSciany) {
					case 1:
						elements.getElement(nrEl).setE1(p);
						break;
					case 2:
						elements.getElement(nrEl).setE2(p);
						break;
					case 3:
						elements.getElement(nrEl).setE3(p);
						break;
					case 4:
						elements.getElement(nrEl).setE4(p);
						break;
				}

				j = ileSas;
				break;
			}
		}
	}

}



void Delanouy3D::wyszukanieSasidniegoElementu(int nrEl) {

	IntList wspEl1(1000, 1000);
	IntList wspEl2(1000, 1000);

	elements.getElement(nrEl).setE1(-1);
	elements.getElement(nrEl).setE2(-1);
	elements.getElement(nrEl).setE3(-1);
	elements.getElement(nrEl).setE4(-1);

	int p1, p2, p3, p4;

	p1 = elements.getElement(nrEl).getP1();
	p2 = elements.getElement(nrEl).getP2();
	p3 = elements.getElement(nrEl).getP3();
	p4 = elements.getElement(nrEl).getP4();


	for (int j = 0, p, ileSas = laplas[p4].getIter(); j < ileSas; ++j) {

		p =  laplas[p4].getElement(j);

		for (int jj = 0, ileSas2 = laplas[p2].getIter(); jj < ileSas2; ++j) {
			if (p == laplas[p2].getElement(jj)) {
				wspEl1.setElement(p);
				break;
			}
		}

	}

	for (int j = 0, p, ileSas = laplas[p1].getIter(); j < ileSas; ++j) {

		p =  laplas[p1].getElement(j);

		for (int jj = 0, ileSas1 = laplas[p3].getIter(); jj < ileSas1; ++j) {
			if (p == laplas[p3].getElement(jj)) {
				wspEl2.setElement(p);
				break;
			}
		}
	}

	//421
	for (int j = 0, p, ileSas = wspEl1.getIter(); j < ileSas; ++j) {

		p =  wspEl1.getElement(j);

		for (int jj = 0, ileSas2 = laplas[p1].getIter(); jj < ileSas2; ++j) {
			if (p == laplas[p1].getElement(jj) && p != nrEl) {

				elements.getElement(nrEl).setE1(p);
				j = ileSas;
				break;
			}
		}
	}

	//432
	for (int j = 0, p, ileSas = wspEl1.getIter(); j < ileSas; ++j) {

		p =  wspEl1.getElement(j);

		for (int jj = 0, ileSas2 = laplas[p3].getIter(); jj < ileSas2; ++j) {
			if (p == laplas[p3].getElement(jj) && p != nrEl) {

				elements.getElement(nrEl).setE2(p);
				j = ileSas;
				break;
			}
		}
	}

	//413
	for (int j = 0, p, ileSas = wspEl2.getIter(); j < ileSas; ++j) {

		p =  wspEl2.getElement(j);

		for (int jj = 0, ileSas2 = laplas[p4].getIter(); jj < ileSas2; ++j) {
			if (p == laplas[p4].getElement(jj) && p != nrEl) {

				elements.getElement(nrEl).setE3(p);
				j = ileSas;
				break;
			}
		}
	}

	//123
	for (int j = 0, p, ileSas = wspEl2.getIter(); j < ileSas; ++j) {

		p =  wspEl2.getElement(j);

		for (int jj = 0, ileSas2 = laplas[p2].getIter(); jj < ileSas2; ++j) {
			if (p == laplas[p2].getElement(jj) && p != nrEl) {

				elements.getElement(nrEl).setE4(p);
				j = ileSas;
				break;
			}
		}
	}


}

void Delanouy3D::wyszukanieSasidnichElementowE() {

	sasiednieElementy();
	IntList wspEl1(1000, 1000);
	IntList wspEl2(1000, 1000);

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		elements.getElement(i).setE1(-1);
		elements.getElement(i).setE2(-1);
		elements.getElement(i).setE3(-1);
		elements.getElement(i).setE4(-1);

	}

	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		wspEl1.ustawIter(0);
		wspEl2.ustawIter(0);

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		for (int j = 0, p, ileSas = laplas[p4].getIter(); j < ileSas; ++j) {

			p =  laplas[p4].getElement(j);

			for (int jj = 0, ileSas2 = laplas[p2].getIter(); jj < ileSas2; ++jj) {
				if (p == laplas[p2].getElement(jj)) {
					wspEl1.setElement(p);
					break;
				}
			}
		}

		for (int j = 0, p, ileSas = laplas[p1].getIter(); j < ileSas; ++j) {

			p =  laplas[p1].getElement(j);

			for (int jj = 0, ileSas1 = laplas[p3].getIter(); jj < ileSas1; ++jj) {
				if (p == laplas[p3].getElement(jj)) {
					wspEl2.setElement(p);
					break;
				}
			}
		}

		//421
		for (int j = 0, p, ileSas = wspEl1.getIter(); j < ileSas; ++j) {

			p =  wspEl1.getElement(j);

			for (int jj = 0, ileSas2 = laplas[p1].getIter(); jj < ileSas2; ++jj) {
				if (p == laplas[p1].getElement(jj) && p != i) {

					elements.getElement(i).setE1(p);
					j = ileSas;
					break;
				}
			}
		}

		//432
		for (int j = 0, p, ileSas = wspEl1.getIter(); j < ileSas; ++j) {

			p =  wspEl1.getElement(j);

			for (int jj = 0, ileSas2 = laplas[p3].getIter(); jj < ileSas2; ++jj) {
				if (p == laplas[p3].getElement(jj) && p != i) {

					elements.getElement(i).setE2(p);
					j = ileSas;
					break;
				}
			}
		}

		//413
		for (int j = 0, p, ileSas = wspEl2.getIter(); j < ileSas; ++j) {

			p =  wspEl2.getElement(j);

			for (int jj = 0, ileSas2 = laplas[p4].getIter(); jj < ileSas2; ++jj) {
				if (p == laplas[p4].getElement(jj) && p != i) {

					elements.getElement(i).setE3(p);
					j = ileSas;
					break;
				}
			}
		}

		//123
		for (int j = 0, p, ileSas = wspEl2.getIter(); j < ileSas; ++j) {

			p =  wspEl2.getElement(j);

			for (int jj = 0, ileSas2 = laplas[p2].getIter(); jj < ileSas2; ++jj) {
				if (p == laplas[p2].getElement(jj) && p != i) {

					elements.getElement(i).setE4(p);
					j = ileSas;
					break;
				}
			}
		}



	}

}


void Delanouy3D::oznaczWezlyGraniczneNaPodstSasP() {
	int zlicz = 0;
	for (int i = 0, flaga, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getGranica() == 'g') {
			flaga = 0;

			for (int j = 0, ilePP = laplas[i].getIter(); j < ilePP; ++j) {

				if (points.getElement(laplas[i].getElement(j)).getGranica() == 'a') {
					flaga = 1;
					break;
				}

			}

			if (!flaga) {
				++zlicz;
				points.getElement(i).setGranica('a');
			}
		}

	}
//m1->Lines->Add(zlicz);
}

void Delanouy3D::zaznaczElementyZPPocz() {


	for (int i = 0, flaga, ileE = elements.getIter(); i < ileE; ++i) {
		flaga = 0;

		if (elements.getElement(i).getP1() < 8) {
			flaga = 1;
		}
		else if (elements.getElement(i).getP2() < 8) {
			flaga = 1;
		}
		else if (elements.getElement(i).getP3() < 8) {
			flaga = 1;
		}
		else if (elements.getElement(i).getP4() < 8) {
			flaga = 1;
		}

		if (flaga) {
			elements.getElement(i).setRodzajZiarna(2);
		}

	}

}

void Delanouy3D::usunWybraneElementy(bool zaznaczone, bool usunluboznacz_niepotrzebneWezly, IntList & nrZ) {
//usunluboznacz_niepotrzebneWezly==true usun

	//nrZazZiaren=nrZ;

	if (zaznaczone) {
		usunZaznaczoneEl(nrZ);
	}
	else {
		usunNieZaznaczoneEl(nrZ);
	}

	sasiednieElementy();


	if (usunluboznacz_niepotrzebneWezly) {
		usunPunktyNiepotrzebne();
		sasiednieElementy();
	}
	else {
		oznaczPunktyNiepotrzebne();
	}



	////sasiedniePunkty();

}


void Delanouy3D::oznaczPunktyNiepotrzebne() {

	for (int i = 0, j; i < points.getIter(); ++i) {
		if (laplas[i].getIter() == 0) {

			points.getElement(i).setGranica('z');

		}
	}

}

void Delanouy3D::usunPunktyNiepotrzebne() {


	for (int i = 0, j; i < points.getIter(); ++i) {
		if (laplas[i].getIter() == 0) {

			j = points.getIter() - 1;
			points.usunElement(i);
			laplas[i].zamienIntListy(laplas[j]);

			for (int z = 0, ileE = laplas[i].getIter(); z < ileE; ++z) {
				elements.getElement(laplas[i].getElement(z)).podmienPunkt(j, i);

			}

			--i;

		}
	}

}










void Delanouy3D::ZapiszDoPlikuSiatkaPow(const char * nazwa) {

	IntList elementy(3000, 3000);
	bool b1, b2, b3, b4;
	bool * zazEl = new bool[elements.getIter()];
	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
		zazEl[i] = true;
	}


	for (int i = 0, ileE = elements.getIter(), p1, p2, p3, p4; i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		if (points.getElement(p1).getGranica() == 'g') {
			b1 = true;
		}
		else {
			b1 = false;
		}
		if (points.getElement(p2).getGranica() == 'g') {
			b2 = true;
		}
		else {
			b2 = false;
		}
		if (points.getElement(p3).getGranica() == 'g') {
			b3 = true;
		}
		else {
			b3 = false;
		}
		if (points.getElement(p4).getGranica() == 'g') {
			b4 = true;
		}
		else {
			b4 = false;
		}

		if (b4 && b2 && b1 && zazEl[elements.getElement(i).getE1()]) {
			elementy.setElement(p4);
			elementy.setElement(p2);
			elementy.setElement(p1);
			zazEl[i] = false;
		}
		if (b4 && b3 && b2 && zazEl[elements.getElement(i).getE2()]) {
			elementy.setElement(p4);
			elementy.setElement(p3);
			elementy.setElement(p2);
			zazEl[i] = false;
		}
		if (b4 && b1 && b3 && zazEl[elements.getElement(i).getE3()]) {
			elementy.setElement(p4);
			elementy.setElement(p1);
			elementy.setElement(p3);
			zazEl[i] = false;
		}
		if (b1 && b2 && b3 && zazEl[elements.getElement(i).getE4()]) {
			elementy.setElement(p1);
			elementy.setElement(p2);
			elementy.setElement(p3);
			zazEl[i] = false;
		}

	}

	delete []zazEl;

	ofstream zapis(nazwa);
	zapis.precision(12);


	int * pozPunktow = new int[points.getIter()];

	int iP = 0;
	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'g') {
			pozPunktow[i] = iP;
			++iP;
		}
	}

//zapisywanie punktow
	zapis << "Points" << endl;
	zapis << iP << endl << endl;

	for (int i = 0, j = 0, ileP = points.getIter(); i < ileP; ++i) {
		if (points.getElement(i).getGranica() == 'g') {
			zapis << j << "       " << points.getElement(i).getX() << "   " << points.getElement(i).getY() << "   " << points.getElement(i).getZ() << endl;
			++j;
		}
	}

//-------------------
//zapisywanie elementow
	zapis << endl << "Elements" << endl;
	zapis << elementy.getIter() / 3 << endl << endl;

	for (int i = 0, j = 0, e1, e2, e3, ileE = elementy.getIter(); i < ileE; i = i + 3) {

		e1 = elementy.getElement(i);
		e2 = elementy.getElement(i + 1);
		e3 = elementy.getElement(i + 2);

		zapis << j << "       " << pozPunktow[e1] << "   " << pozPunktow[e2] << "   " << pozPunktow[e3] << endl;
		++j;
	}


//-------------------
	zapis.close();

	delete []pozPunktow;

}

void Delanouy3D::oznaczWezlyNaPow_b() {

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double x, y, z;

	for (int i = 0, ileP = points.getIter(); ileP < i; ++i) {

		x = points.getElement(i).getX();
		y = points.getElement(i).getY();
		z = points.getElement(i).getZ();

		if (0 == x || maxX == x || 0 == y || maxY == y || 0 == z || maxZ == z) {
			if (points.getElement(i).getGranica() != 'g') {
				points.getElement(i).setGranica('b');
			}
		}

	}

}

void Delanouy3D::odznaczWezlyNaPow_b() {


	for (int i = 0, ileP = points.getIter(); ileP < i; ++i) {

		if (points.getElement(i).getGranica() == 'b') {
			points.getElement(i).setGranica('a');
		}

	}
}


void Delanouy3D::dodajPunktWSrodkuElementow() {


	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
		if (elements.getElement(i).getRodzajZiarna() == -1) {
			points.setElement(elements.getElement(i).getAx(), elements.getElement(i).getAy(), elements.getElement(i).getAz(), 'a');
		}
	}

}

void Delanouy3D::dodajPunktPopElement() {

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {
		if (elements.getElement(i).getRodzajZiarna() == -1) {
			points.setElement(elements.getElement(i).getSx(), elements.getElement(i).getSy(), elements.getElement(i).getSz(), 'a');
		}
	}

}


bool Delanouy3D::sprawdzElement(int jaki, IntList & lista) {

	if (jaki != -1) {
		for (int i = 0; i < lista.getIter(); i++) {
			if (jaki == lista.getElement(i)) {
				return false;
			}
		}
		return true;
	}

	return false;
}

int Delanouy3D::sprawdzenieElementow(const double & x, const double & y, const double & z, IntList & sprawEle, const int & od) {

	double X, Y, Z;
	float R, temp;
	float max = 10000000;
	int wybrany = 0;
	for (int i = od, ktory, ileE = sprawEle.getIter(); i < ileE; ++i) {

		ktory = sprawEle.getElement(i);

		X = elements.getElement(ktory).getSx();
		Y = elements.getElement(ktory).getSy();
		Z = elements.getElement(ktory).getSz();
		R = elements.getElement(ktory).getR();
		temp  = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);

		if (temp <= R) {
			return ktory;
		}

		X = elements.getElement(ktory).getAx();
		Y = elements.getElement(ktory).getAy();
		Z = elements.getElement(ktory).getAz();
		temp  = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);

		if (max > temp) {
			max = temp;
			wybrany = ktory;
		}

	}

	return wybrany;
}





bool Delanouy3D::elementyDoDelnoya(int wybrany, double x, double y, double z, IntList & delElementy) {

	if (wybrany == -1) {
		return false;
	}

	int e1, e2, e3, e4;

	double X, Y, Z;
	float Rr;

	delElementy.ustawIter(0);

	delElementy.setElement(wybrany);

	int start = 0, koniec = 1;

	while (start != koniec) {

		for (int i = start; i < koniec; i++) {

			wybrany = delElementy.getElement(i);

			e1 = elements.getElement(wybrany).getE1();
			e2 = elements.getElement(wybrany).getE2();
			e3 = elements.getElement(wybrany).getE3();
			e4 = elements.getElement(wybrany).getE4();

			if (sprawdzElement(e1, delElementy)) {



				X = elements.getElement(e1).getSx();
				Y = elements.getElement(e1).getSy();
				Z = elements.getElement(e1).getSz();
				//R = elements.getElement(e1).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;

				if (Rr <= elements.getElement(e1).getR()) {

					delElementy.setElement(e1);
				}
				/*
						else{
						if(sprawdzCzyNalezyDoE(e1,x,y)){
						m1->Lines->Add("dodano nowy e1");
						delElementy.setElement(e1);
						}
						}
				*/

			}
			if (sprawdzElement(e2, delElementy)) {

				X = elements.getElement(e2).getSx();
				Y = elements.getElement(e2).getSy();
				Z = elements.getElement(e2).getSz();
				//R = elements.getElement(e2).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;
				if (Rr <= elements.getElement(e2).getR()) {

					delElementy.setElement(e2);
				}

				/*
						else{
						if(sprawdzCzyNalezyDoE(e2,x,y)){
						m1->Lines->Add("dodano nowy e2");
						delElementy.setElement(e2);
						}
						}
				*/

			}
			if (sprawdzElement(e3, delElementy)) {

				X = elements.getElement(e3).getSx();
				Y = elements.getElement(e3).getSy();
				Z = elements.getElement(e3).getSz();
				//R = elements.getElement(e3).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;
				if (Rr <= elements.getElement(e3).getR()) {

					delElementy.setElement(e3);
				}
				/*
						else{
						if(sprawdzCzyNalezyDoE(e3,x,y)){
						m1->Lines->Add("dodano nowy e3");
						delElementy.setElement(e3);
						}
						}
				*/

			}
			if (sprawdzElement(e4, delElementy)) {

				X = elements.getElement(e4).getSx();
				Y = elements.getElement(e4).getSy();
				Z = elements.getElement(e4).getSz();
				//R = elements.getElement(e4).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;
				if (Rr <= elements.getElement(e4).getR()) {

					delElementy.setElement(e4);
				}
				/*
						else{
						if(sprawdzCzyNalezyDoE(e4,x,y)){
						m1->Lines->Add("dodano nowy e4");
						delElementy.setElement(e4);
						}
						}
				*/

			}

//koniec for
		}

		start = koniec;
		koniec = delElementy.getIter();

//koniec while
	}

	bool flag = true;
	/*
	for(int i=0;i<delElementy.getIter();i++){

	if(sprawdzCzyNalezyDoE(delElementy.getElement(i),x,y)){flag = true;break;}


	}

	if(!flag){text->Lines->Add(" nie nalezy do e ");}
	*/


	return flag;

}


bool Delanouy3D::elementyDoDelnoyaZeSprObj(int wybrany, double x, double y, double z, IntList & delElementy) {

	if (wybrany == -1) {
		return false;
	}

	int e1, e2, e3, e4;

	double X, Y, Z;
	float Rr;

	delElementy.ustawIter(0);

	delElementy.setElement(wybrany);

	int start = 0, koniec = 1;

	while (start != koniec) {

		for (int i = start; i < koniec; i++) {

			wybrany = delElementy.getElement(i);

			e1 = elements.getElement(wybrany).getE1();
			e2 = elements.getElement(wybrany).getE2();
			e3 = elements.getElement(wybrany).getE3();
			e4 = elements.getElement(wybrany).getE4();

			if (sprawdzElement(e1, delElementy)) {



				X = elements.getElement(e1).getSx();
				Y = elements.getElement(e1).getSy();
				Z = elements.getElement(e1).getSz();
				//R = elements.getElement(e1).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;

				if (Rr <= elements.getElement(e1).getR()) {

					delElementy.setElement(e1);
				}
				/*
						else{
						if(sprawdzCzyNalezyDoE(e1,x,y)){
						m1->Lines->Add("dodano nowy e1");
						delElementy.setElement(e1);
						}
						}
				*/

			}
			if (sprawdzElement(e2, delElementy)) {

				X = elements.getElement(e2).getSx();
				Y = elements.getElement(e2).getSy();
				Z = elements.getElement(e2).getSz();
				//R = elements.getElement(e2).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;
				if (Rr <= elements.getElement(e2).getR()) {

					delElementy.setElement(e2);
				}

				/*
						else{
						if(sprawdzCzyNalezyDoE(e2,x,y)){
						m1->Lines->Add("dodano nowy e2");
						delElementy.setElement(e2);
						}
						}
				*/

			}
			if (sprawdzElement(e3, delElementy)) {

				X = elements.getElement(e3).getSx();
				Y = elements.getElement(e3).getSy();
				Z = elements.getElement(e3).getSz();
				//R = elements.getElement(e3).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;
				if (Rr <= elements.getElement(e3).getR()) {

					delElementy.setElement(e3);
				}
				/*
						else{
						if(sprawdzCzyNalezyDoE(e3,x,y)){
						m1->Lines->Add("dodano nowy e3");
						delElementy.setElement(e3);
						}
						}
				*/

			}
			if (sprawdzElement(e4, delElementy)) {

				X = elements.getElement(e4).getSx();
				Y = elements.getElement(e4).getSy();
				Z = elements.getElement(e4).getSz();
				//R = elements.getElement(e4).getR();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
				//Rr = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;
				if (Rr <= elements.getElement(e4).getR()) {

					delElementy.setElement(e4);
				}
				/*
						else{
						if(sprawdzCzyNalezyDoE(e4,x,y,z)){
						m1->Lines->Add("dodano nowy e4");
						delElementy.setElement(e4);
						}
						}
				*/

			}

//koniec for
		}

		start = koniec;
		koniec = delElementy.getIter();

//koniec while
	}

	bool flag = false;

	for (int i = 0; i < delElementy.getIter(); i++) {
		if (sprawdzCzyNalezyDoE(delElementy.getElement(i), x, y, z)) {
			flag = true;
			break;
		}
	}

	if (!flag) {
		//m1->Lines->Add(" nie nalezy do e ");
	}

	return flag;

}

bool Delanouy3D::elementyDoDelnoyaPoprawaE(int wybrany, IntList & delElementy) {
	delElementy.ustawIter(0);

	if (wybrany == -1) {
		return false;
	}

	const double x = elements.getElement(wybrany).getSx();
	const double y = elements.getElement(wybrany).getSy();
	const double z = elements.getElement(wybrany).getSz();


	int e1, e2, e3, e4;

	double X, Y, Z;
	float Rr;

	delElementy.setElement(wybrany);

	int start = 0, koniec = 1;

	while (start != koniec) {

		for (int i = start; i < koniec; ++i) {

			wybrany = delElementy.getElement(i);

			e1 = elements.getElement(wybrany).getE1();
			e2 = elements.getElement(wybrany).getE2();
			e3 = elements.getElement(wybrany).getE3();
			e4 = elements.getElement(wybrany).getE4();

			if (sprawdzElement(e1, delElementy)) {

				X = elements.getElement(e1).getSx();
				Y = elements.getElement(e1).getSy();
				Z = elements.getElement(e1).getSz();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);

				if (Rr <= elements.getElement(e1).getR()) {

					delElementy.setElement(e1);

				}

			}
			if (sprawdzElement(e2, delElementy)) {

				X = elements.getElement(e2).getSx();
				Y = elements.getElement(e2).getSy();
				Z = elements.getElement(e2).getSz();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);

				if (Rr <= elements.getElement(e2).getR()) {

					delElementy.setElement(e2);
				}

			}
			if (sprawdzElement(e3, delElementy)) {

				X = elements.getElement(e3).getSx();
				Y = elements.getElement(e3).getSy();
				Z = elements.getElement(e3).getSz();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);

				if (Rr <= elements.getElement(e3).getR()) {

					delElementy.setElement(e3);

				}


			}
			if (sprawdzElement(e4, delElementy)) {

				X = elements.getElement(e4).getSx();
				Y = elements.getElement(e4).getSy();
				Z = elements.getElement(e4).getSz();
				Rr = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);

				if (Rr <= elements.getElement(e4).getR()) {

					delElementy.setElement(e4);

				}


			}

			return true;
//koniec for
		}

		start = koniec;
		koniec = delElementy.getIter();

//koniec while
	}

	return true;
}


void Delanouy3D::Start(double xdl, double ydl, double zdl) {



	dx = xdl;
	dy = ydl;
	dz = zdl;

	wzorGranicy.czysc(1000, 1000);
	points.czysc(10000, 5000);
	elements.czysc(70000, 10000);

	Punkt start1(0, 0, 0), start2(xdl, 0, 0), start3(0, ydl, 0), start4(xdl, ydl, 0),
		  start5(0, 0, zdl), start6(xdl, 0, zdl), start7(0, ydl, zdl), start8(xdl, ydl, zdl);

//Element el1(6,3,7,4),el2(2,6,0,3),el3(4,0,6,3),el4(5,7,1,4),el5(1,3,0,4),el6(7,4,3,1);
//Element el1(3,7,6,4),el2(0,2,6,3),el3(0,6,4,3),el4(1,5,7,4),el5(0,1,3,4),el6(3,7,4,1);

//Element el1(7,6,3,4),el2(6,0,2,3),el3(6,4,0,3),el4(7,1,5,4),el5(3,0,1,4),el6(7,3,1,4);

	Element el1(6, 4, 3, 7), el2(3, 2, 0, 6), el3(4, 3, 0, 6), el4(5, 1, 4, 7), el5(3, 0, 1, 4), el6(4, 1, 3, 7);

//Element el1,el2,el3,el4,el5,el6;

	/*

	gora

		2  3
		0  1

	dol
		6  7
		4  5

	*/

	points.setElement(start1);
	points.setElement(start2);
	points.setElement(start3);
	points.setElement(start4);
	points.setElement(start5);
	points.setElement(start6);
	points.setElement(start7);
	points.setElement(start8);


	elements.setElement(el1);
	elements.setElement(el2);
	elements.setElement(el3);
	elements.setElement(el4);
	elements.setElement(el5);
	elements.setElement(el6);

	ilePunktowCzworosciany = 0;
}

void Delanouy3D::Start(double xdl, double ydl, double zdl, int ileP) {

	dx = xdl;
	dy = ydl;
	dz = zdl;
	ilePunktowCzworosciany = 0;
	points.czysc(ileP * 2, ileP * 0.5);
	elements.czysc(ileP * 7, ileP);

	Punkt start1(0, 0, 0), start2(xdl, 0, 0), start3(0, ydl, 0), start4(xdl, ydl, 0),
		  start5(0, 0, zdl), start6(xdl, 0, zdl), start7(0, ydl, zdl), start8(xdl, ydl, zdl);

	double sx = xdl * 0.5, sy = ydl * 0.5, sz = zdl * 0.5, r = (xdl - sx) * (xdl - sx) + (ydl - sy) * (ydl - sy) + (zdl - sz) * (zdl - sz);



	points.setElement(start1);
	points.setElement(start2);
	points.setElement(start3);
	points.setElement(start4);
	points.setElement(start5);
	points.setElement(start6);
	points.setElement(start7);
	points.setElement(start8);

//Element el1(6,4,3,7),el2(3,2,0,6),el3(4,3,0,6),el4(5,1,4,7),el5(3,0,1,4),el6(4,1,3,7);

	elements.setElement(6, 4, 3, 7, sx, sy, sz, r, (xdl + xdl) * 0.25, (ydl + ydl + ydl) * 0.25, (zdl + zdl + zdl) * 0.25);
	elements.setElement(3, 2, 0, 6, sx, sy, sz, r, xdl * 0.25, (ydl + ydl + ydl) * 0.25, zdl * 0.25);
	elements.setElement(4, 3, 0, 6, sx, sy, sz, r, xdl * 0.25, (ydl + ydl) * 0.25, (zdl + zdl) * 0.25);
	elements.setElement(5, 1, 4, 7, sx, sy, sz, r, (xdl + xdl + xdl) * 0.25, ydl * 0.25, (zdl + zdl + zdl) * 0.25);
	elements.setElement(3, 0, 1, 4, sx, sy, sz, r, (xdl + xdl) * 0.25, ydl * 0.25, zdl * 0.25);
	elements.setElement(4, 1, 3, 7, sx, sy, sz, r, (xdl + xdl + xdl) * 0.25, (ydl + ydl) * 0.25, (zdl + zdl) * 0.25);

	elements.getElement(0).setElementySasiad(-1, 5, -1, 2);
	elements.getElement(1).setElementySasiad(-1, -1, 2, -1);
	elements.getElement(2).setElementySasiad(0, 1, -1, 4);
	elements.getElement(3).setElementySasiad(-1, 5, -1, -1);
	elements.getElement(4).setElementySasiad(2, -1, 5, -1);
	elements.getElement(5).setElementySasiad(3, -1, 0, 4);


}




void Delanouy3D::delanoyTablicaOP(double xdl, double ydl, double zdl, bool reMes, PunktList & Tab, double ws, int te, bool dokladneWyszukanie) {


	Start(xdl, ydl, zdl, Tab.getIter());

	int TE2 = te;

	double ws2 = ws;
	int wymiarTE2X = (int)(xdl * ws2) + 2 * TE2;
	int wymiarTE2Y = (int)(ydl * ws2) + 2 * TE2;
	int wymiarTE2Z = (int)(zdl * ws2) + 2 * TE2;


//sprawdzNumeracje();
	IntList delElement(2000, 1000);
	IntList zlyElement(200, 100);
	IntList pomoc(2000, 1000);
	IntList sprawdz(50, 20);

	int ** tabElementow2;
	tabElementow2 = new int * [wymiarTE2Z];
	for (int i = 0; i < wymiarTE2Z; i++) {
		tabElementow2[i] = new int[wymiarTE2X * wymiarTE2Y];
	}
	for (int i = 0; i < wymiarTE2Z; i++) {
		for (int j = 0; j < wymiarTE2X * wymiarTE2Y; j++) {
			tabElementow2[i][j] = -1;
		}
	}


//int naliczanie1=0;
//int naliczanie2=0;
	int naliczanie4 = 0;
	char granica;

	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3, a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, WW, W, X, Y, Z;
	int wszystkiePunkty = Tab.getIter();
	int wszystkiePunktyPomoc = wszystkiePunkty;
	int ilePowtorzenPunktow = 10;


	int punktReMes = 0;
	if (reMes) {
		punktReMes = 8;
	}

	for (int pun = punktReMes; pun < wszystkiePunkty; pun++) {

		granica = Tab.getElement(pun).getGranica();
		x0 = Tab.getElement(pun).getX();
		y0 = Tab.getElement(pun).getY();
		z0 = Tab.getElement(pun).getZ();

		delElement.ustawIter(0);

		int nrNowegoEl = elements.getIter();
		int iloscEStart = elements.getIter();
		int ileElUsunieto = 0;
		int znaleziony = -1;

		/*
		if(pun==0){
		 m1->Lines->Add("start");
		 m1->Lines->Add(wymiarTE2X);
		 m1->Lines->Add(TE2);
		 m1->Lines->Add("start");

			for(int w1=0;w1<zdl;w1++){
				for(int w2=0;w2<ydl;w2++){
					for(int w3=0;w3<xdl;w3++){

						int wz2r = (int((w1+0.3)*ws2))+TE2;
						int wy2r = (int((w2+0.3)*ws2))+TE2;
						int wx2r = (int((w3+0.3)*ws2))+TE2;

						m1->Lines->Add(wz2r);
						m1->Lines->Add(wy2r*(wymiarTE2X)+wx2r);

					}
				}
			}
		}
		*/

		int wx2 = ((int)(x0 * ws2)) + TE2;
		int wy2 = ((int)(y0 * ws2)) + TE2;
		int wz2 = ((int)(z0 * ws2)) + TE2;


		if (tabElementow2[wz2][wy2 * wymiarTE2X + wx2] != -1) {

			znaleziony = tabElementow2[wz2][wy2 * wymiarTE2X + wx2];
			//naliczanie1++;

		}

		else {

			for (int t2z = -TE2, max = 100000000; t2z < TE2; t2z++) {

				for (int t2y = -TE2; t2y < TE2; t2y++) {

					for (int t2x = -TE2; t2x < TE2; t2x++) {

						if (tabElementow2[wz2 + t2z][(wy2 + t2y)*wymiarTE2X + (wx2 + t2x)] != -1) {
							if (max > t2x * t2x + t2y * t2y + t2z * t2z) {
								znaleziony = tabElementow2[wz2 + t2z][(wy2 + t2y) * wymiarTE2X + (wx2 + t2x)];
								max = t2x * t2x + t2y * t2y + t2z * t2z;
							}
						}

					}
				}

			}

		}




		//if(znaleziony==-1){naliczanie2++;}
		//naliczanie3++;

		znaleziony = wyszukajElementOP(znaleziony, x0, y0, z0, zlyElement, pomoc);

		if (znaleziony == -1 && dokladneWyszukanie) {
			znaleziony = wyszukajElement(x0, y0, z0);;
		}
		//znaleziony = wyszukajElement(x0,y0,z0,sprawdz);

		if (elementyDoDelnoya(znaleziony, x0, y0, z0, sprawdz)) {

			points.setElement(x0, y0, z0, granica);
			ileElUsunieto = sprawdz.getIter();

			for (int i = 0, wybrany, iter; i < ileElUsunieto; i++) {

				wybrany = sprawdz.getElement(i);

// 0 - nowa pozycja powstalego elementu
// 1 - pozycja elementu niszczonego
// 2 - pierwszy wezel nowego el
// 3 - drugi wezel nowego elementu
// 4 - trzeci wezel nowego elementu
// 5 - kolejnosc sciany w elemencie niszczonym
// 6 - 1 sasiedni element
// 7 - 2 sasiedni element
// 8 - 3 sasiedni element
// 9 - 4 sasiedni element

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(2);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP4());
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(3);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				delElement.setElement(1);
				delElement.setElement(wybrany);
				delElement.setElement(elements.getElement(wybrany).getP1());
				delElement.setElement(elements.getElement(wybrany).getP2());
				delElement.setElement(elements.getElement(wybrany).getP3());
				delElement.setElement(4);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);
				delElement.setElement(-1);

				iter = delElement.getIter() - 40;

				for (int a = 0; a < iter; a += 10) {

					if (delElement.getElement(iter + 2) == delElement.getElement(a + 2) && delElement.getElement(iter + 3) == delElement.getElement(a + 4) && delElement.getElement(iter + 4) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 12) == delElement.getElement(a + 2) && delElement.getElement(iter + 13) == delElement.getElement(a + 4) && delElement.getElement(iter + 14) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 10, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 22) == delElement.getElement(a + 2) && delElement.getElement(iter + 23) == delElement.getElement(a + 4) && delElement.getElement(iter + 24) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 20, -1);
						delElement.podmienElement(a, -1);
					}
					if (delElement.getElement(iter + 32) == delElement.getElement(a + 2) && delElement.getElement(iter + 33) == delElement.getElement(a + 4) && delElement.getElement(iter + 34) == delElement.getElement(a + 3)) {
						delElement.podmienElement(iter + 30, -1);
						delElement.podmienElement(a, -1);
					}

				}
			}




			int iter = delElement.getIter();

			for (int i = 0; i < iter; i += 10) {


				if (delElement.getElement(i) != -1) {

					double Rr;
					int e1 = delElement.getElement(i + 2), e2 = delElement.getElement(i + 3), e3 = delElement.getElement(i + 4);

					x1 = points.getElement(e1).getX();
					y1 = points.getElement(e1).getY();
					z1 = points.getElement(e1).getZ();

					x2 = points.getElement(e2).getX();
					y2 = points.getElement(e2).getY();
					z2 = points.getElement(e2).getZ();

					x3 = points.getElement(e3).getX();
					y3 = points.getElement(e3).getY();
					z3 = points.getElement(e3).getZ();


					WW = (x0 - x1) * (y2 - y1) * (z3 - z1) + (x2 - x1) * (y3 - y1) * (z0 - z1) + (x3 - x1) * (y0 - y1) * (z2 - z1)
							- (z0 - z1) * (y2 - y1) * (x3 - x1) - (z2 - z1) * (y3 - y1) * (x0 - x1) - (z3 - z1) * (y0 - y1) * (x2 - x1);

					//WW  =x0*y2*z3+x2*y3*z0+x3*y0*z2-z0*y2*x3-z2*y3*x0-z3*y0*x2;
					if (WW) {

						a2 = (-x3 + x1) * 2;
						b2 = (-y3 + y1) * 2;
						c2 = (-z3 + z1) * 2;
						l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

						a3 = (-x0 + x1) * 2;
						b3 = (-y0 + y1) * 2;
						c3 = (-z0 + z1) * 2;
						l3 = -(x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1);

						a4 = (-x0 + x2) * 2;
						b4 = (-y0 + y2) * 2;
						c4 = (-z0 + z2) * 2;
						l4 = -(x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 - z2 * z2);

						W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
						W = 1 / W;
						l2 *= W;
						l3 *= W;
						l4 *= W;
						if (W) {
							X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
							Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
							Y = Y * -1;
							Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

							Rr =  (X - x0) * (X - x0) + (Y - y0) * (Y - y0) + (Z - z0) * (Z - z0);
							//Rr = X*X-2*X*x0+x0*x0 + Y*Y-2*Y*y0+y0*y0 + Z*Z-2*Z*z0+z0*z0;

							elements.setElement(e1, e2, e3, points.getIter() - 1, X, Y, Z, Rr, (x1 + x2 + x3 + x0) * 0.25, (y1 + y2 + y3 + y0) * 0.25, (z1 + z2 + z3 + z0) * 0.25);
							//elements.setElement(e1,e2,e3,points.getIter()-1,X,Y,Z,Rr,(x1+x2+x3+x0)*0.25,(y1+y2+y3+y0)*0.25,(z1+z2+z3+z0)*0.25);
						}
						else {
							delElement.podmienElement(i, -1);
						}
					}
					else {
						delElement.podmienElement(i, -1);
					}
				}
			}




			for (int a = 0; a < iter; a += 10) {

				if (delElement.getElement(a) != -1) {
					delElement.podmienElement(a, nrNowegoEl++);
				}

			}




			int iter2 = iter;

			for ( int a = iter - 10, i = 1; i < iter; a -= 10) {

				if (a >= 0) {
					if (delElement.getElement(a) != -1) {

						delElement.podmienElement(a, delElement.getElement(i));
						i += 40;

					}
				}
				else {

					iter2 += int(((elements.getIter() - iloscEStart) - iter * 0.025) * 40);
					break;
				}


			}

			nrNowegoEl -= elements.getIter();

//////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {

					int sasiad = -1;

					switch (delElement.getElement(i + 5)) {

						case 1:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE1();

							break;
						case 2:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE2();

							break;
						case 3:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE3();

							break;
						case 4:
							sasiad = elements.getElement(delElement.getElement(i + 1)).getE4();

							break;
					}

					//i+9 ostatnia sciana 123
					delElement.podmienElement(i + 9, sasiad);

					if (sasiad != -1) {
						//if(0){


						if (elements.getElement(sasiad).getP4() == delElement.getElement(i + 2)) {

							if (elements.getElement(sasiad).getP1() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE1(delElement.getElement(i));
							}
							else if (elements.getElement(sasiad).getP2() == delElement.getElement(i + 3)) {
								elements.getElement(sasiad).setE2(delElement.getElement(i));
							}
							else {
								elements.getElement(sasiad).setE3(delElement.getElement(i));
							}

						}

						else {
							elements.getElement(sasiad).setE4(delElement.getElement(i));
						}
					}


					for (int j = i + 10; j < iter; j += 10) {

						if (delElement.getElement(j) != -1) {

							if (delElement.getElement(i + 2) == delElement.getElement(j + 2)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 3)) {
								if (delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 8, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+3 = j+2");
								}
							}
							else if (delElement.getElement(i + 2) == delElement.getElement(j + 4)) {
								if (delElement.getElement(i + 3) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 6, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									//m1->Lines->Add("i+4 = j+2");

								}
							}
							else {

								if (delElement.getElement(i + 3) == delElement.getElement(j + 2) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 4)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 8, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 3) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 2)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 6, delElement.getElement(i));
								}
								else if (delElement.getElement(i + 3) == delElement.getElement(j + 4) &&
										delElement.getElement(i + 4) == delElement.getElement(j + 3)) {

									delElement.podmienElement(i + 7, delElement.getElement(j));
									delElement.podmienElement(j + 7, delElement.getElement(i));
								}

							}

						}
					}

				}

			}


			//usuwanie elementow
			//usuwanie wpisu elementu w tablicy
			for (int i = 1; i < iter2; i += 40) {
				//ww = delElement.getElement(i);
				//tabElementow2[(int(elementy.getElement(ww).getAx()*ws2))+TE2][(int(elementy.getElement(ww).getAy()*ws2))+TE2] = -1;
				//tabElementow2[(int(elements.getElement(ww).getAz()*ws2))+TE2][((int(elements.getElement(ww).getAy()*ws2))+TE2)*wymiarTE2 + (int(elements.getElement(ww).getAx()*ws2))+TE2] = -1;
				elements.deleteElement(delElement.getElement(i));

			}

			for (int i = 0, ww; i < iter; i += 10) {

				if (delElement.getElement(i) != -1) {
					ww = delElement.getElement(i);
					elements.getElement(ww).setElementySasiad(delElement.getElement(i + 6), delElement.getElement(i + 7), delElement.getElement(i + 8), delElement.getElement(i + 9));

					tabElementow2[int(elements.getElement(ww).getAz()*ws2) + TE2][(int(elements.getElement(ww).getAy()*ws2) + TE2)*wymiarTE2X + int(elements.getElement(ww).getAx()*ws2) + TE2] = ww;
					//tabElementow1[int(elementy.getElement(delElement.getElement(i)).getAx()*ws1)][int(elementy.getElement(delElement.getElement(i)).getAy()*ws1)]++;

				}

			}

			for (int i = iter2 + 1; i < iter; i += 40) {

				elements.deleteElementPodmienS(delElement.getElement(i));

			}

		}
		else {

			--wszystkiePunktyPomoc;

			if (wszystkiePunktyPomoc > pun) {
				Tab.zamianaMiejsc(pun, wszystkiePunktyPomoc);
				--pun;
				++naliczanie4;
			}
			else if (ilePowtorzenPunktow > 0 && wszystkiePunktyPomoc < pun) {
				wszystkiePunktyPomoc = wszystkiePunkty - 1;
				--pun;
				--ilePowtorzenPunktow;
				//m1->Lines->Add(ilePowtorzenPunktow);
			}

		}
	}

//wymiarTE2++;
	for (int i = 0; i < wymiarTE2Z; i++) {
		delete []tabElementow2[i];
	}
	delete []tabElementow2;


//m1->Lines->Add("+++++++++");
	if (naliczanie4) { //m1->Lines->Add(naliczanie4);
	}
//m1->Lines->Add(naliczanie2);
//m1->Lines->Add(naliczanie3-(naliczanie2+naliczanie1));



}

int Delanouy3D::wyszukajElementOP(const int & wybranyTE, const double & x, const double & y, const double & z, IntList & zlyElement, IntList & pomoc) {
//int Siatka2D::wyszukajElement(double x,double y){

//int coWchodzi= wybranyTE;
	double xx, yy, zz;
//-----------
	zlyElement.ustawIter(0);

	int ileE = elements.getIter(), ktory, wybrany;
	double sX, sY, sZ;
	float max = 99999999, maxStart = 99999999;
	double aX, aY, aZ;
	float tempR;

	if (wybranyTE != -1) {
		wybrany = wybranyTE;
		sX = elements.getElement(wybrany).getSx();
		sY = elements.getElement(wybrany).getSy();
		sZ = elements.getElement(wybrany).getSz();
		tempR = (sX - x) * (sX - x) + (sY - y) * (sY - y) + (sZ - z) * (sZ - z);
		xx = sX;
		yy = sY;
		zz = sZ;
		//tempR = sX*sX-2*sX*x+x*x + sY*sY-2*sY*y+y*y + sZ*sZ-2*sZ*z+z*z;
		float R = elements.getElement(wybrany).getR();

		if (tempR <= R) {
			return wybrany;
		}





	}
	else {

		for (int i = 0; i < 25; i++) {


			ktory = (double)rand() / (RAND_MAX) * ileE;

			sX = elements.getElement(ktory).getAx();
			sY = elements.getElement(ktory).getAy();
			sZ = elements.getElement(ktory).getAz();
			tempR = (sX - x) * (sX - x) + (sY - y) * (sY - y) + (sZ - z) * (sZ - z);

			//tempR = sX*sX-2*sX*x+x*x + sY*sY-2*sY*y+y*y + sZ*sZ-2*sZ*z+z*z;
			if (max > tempR) {
				max = tempR;
				wybrany = ktory;
			}
		}


	}

	int war = 1;
	int ww = 0;
	bool drugiZly = false;
	int zly = -1;
	max = 99999999;



	while (war && ww < 200) {






		ww++;

		while (war) {


			war = 0;
			int e1 = elements.getElement(wybrany).getE1();
			int e2 = elements.getElement(wybrany).getE2();
			int e3 = elements.getElement(wybrany).getE3();
			int e4 = elements.getElement(wybrany).getE4();

			;
			//m1->Lines->Add(e1);
			//m1->Lines->Add(e2);
			//m1->Lines->Add(e3);
			//m1->Lines->Add(e4);

			float temp;
			//int punk;

			if (zlyElement.sprawdzElement(e1)) {


				aX = elements.getElement(e1).getAx();
				aY = elements.getElement(e1).getAy();
				aZ = elements.getElement(e1).getAz();
				/*
				m1->Lines->Add(e1);
				m1->Lines->Add(x);
				m1->Lines->Add(aX);
				m1->Lines->Add(y);
				m1->Lines->Add(aY);
				m1->Lines->Add(z);
				m1->Lines->Add(aZ);

				m1->Lines->Add(sqrt((x-aX)*(x-aX)+(y-aY)*(y-aY)+(z-aZ)*(z-aZ)));
				*/
				/*
				punk = elements.getElement(e1).getP4();
				aX = points.getElement(punk).getX();
				aY = points.getElement(punk).getY();
				aZ = points.getElement(punk).getZ();
				*/

				temp = (aX - x) * (aX - x) + (aY - y) * (aY - y) + (aZ - z) * (aZ - z);
				//temp = aX*aX-2*aX*x+x*x + aY*aY-2*aY*y+y*y + aZ*aZ-2*aZ*z+z*z;

				if (max > temp) {
					max = temp;
					wybrany = e1;
					war = 1;
				}

			}
			if (zlyElement.sprawdzElement(e2)) {

				aX = elements.getElement(e2).getAx();
				aY = elements.getElement(e2).getAy();
				aZ = elements.getElement(e2).getAz();

				/*
				punk = elements.getElement(e2).getP4();
				aX = points.getElement(punk).getX();
				aY = points.getElement(punk).getY();
				aZ = points.getElement(punk).getZ();
				*/
				temp = (aX - x) * (aX - x) + (aY - y) * (aY - y) + (aZ - z) * (aZ - z);
				//temp = aX*aX-2*aX*x+x*x + aY*aY-2*aY*y+y*y + aZ*aZ-2*aZ*z+z*z;
				if (max > temp) {
					max = temp;
					wybrany = e2;
					war = 1;
				}


			}
			if (zlyElement.sprawdzElement(e3)) {

				aX = elements.getElement(e3).getAx();
				aY = elements.getElement(e3).getAy();
				aZ = elements.getElement(e3).getAz();

				/*
				punk = elements.getElement(e3).getP4();
				aX = points.getElement(punk).getX();
				aY = points.getElement(punk).getY();
				aZ = points.getElement(punk).getZ();
				*/
				temp = (aX - x) * (aX - x) + (aY - y) * (aY - y) + (aZ - z) * (aZ - z);
				//temp = aX*aX-2*aX*x+x*x + aY*aY-2*aY*y+y*y + aZ*aZ-2*aZ*z+z*z;
				if (max > temp) {
					max = temp;
					wybrany = e3;
					war = 1;
				}

			}
			if (zlyElement.sprawdzElement(e4)) {

				aX = elements.getElement(e4).getAx();
				aY = elements.getElement(e4).getAy();
				aZ = elements.getElement(e4).getAz();

				/*
				punk = elements.getElement(e4).getP4();
				aX = points.getElement(punk).getX();
				aY = points.getElement(punk).getY();
				aZ = points.getElement(punk).getZ();
				*/
				temp = (aX - x) * (aX - x) + (aY - y) * (aY - y) + (aZ - z) * (aZ - z);
				//temp = aX*aX-2*aX*x+x*x + aY*aY-2*aY*y+y*y + aZ*aZ-2*aZ*z+z*z;
				if (max > temp) {
					max = temp;
					wybrany = e4;
					war = 1;
				}


			}

			if (drugiZly) {
				zlyElement.setElement(wybrany);
				drugiZly = false;
			}
		}

		double X = elements.getElement(wybrany).getSx();
		double Y = elements.getElement(wybrany).getSy();
		double Z = elements.getElement(wybrany).getSz();
		float R = elements.getElement(wybrany).getR();
		float temp = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
//double temp = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;
		if (temp <= R) {
			return wybrany;
		}
		else {
			//pomoc.ustawIter(0);
			//pomoc.setElement(wybrany);
			//wybrany = wyszukajElementGwiazda(x,y,z,pomoc,3);
			//X = elements.getElement(wybrany).getSx();
			//Y = elements.getElement(wybrany).getSy();
			//Z = elements.getElement(wybrany).getSz();
			//R = elements.getElement(wybrany).getR();
			//if(porowRO(temp,R)){return wybrany;}


			//if(porowL(temp,R)){

			zlyElement.setElement(wybrany);
			max = maxStart + 10;
			war = 1;
			drugiZly = true;

			if (zly != -1) {
				wybrany = zlyElement.getElement(zly);
			}
			zly = zlyElement.getIter() - 1;

			/*
			   pomoc.ustawIter(0);
			   pomoc.setElement(wybrany);
			 wybrany = wyszukajElementGwiazda(x,y,z,pomoc,5);
			  zlyElement.setElement(wybrany);
			  max=maxStart+10;war=1;
			 */
		}



	}


	double X = elements.getElement(wybrany).getSx();
	double Y = elements.getElement(wybrany).getSy();
	double Z = elements.getElement(wybrany).getSz();
	float R = elements.getElement(wybrany).getR();
	float temp = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
//double temp = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;

	if (temp <= R) {
		return wybrany;
	}



	zlyElement.ustawIter(0);
	ww = 0;
	drugiZly = false;
	zly = -1;
	int e1, e2, e3, e4;

	max = 99999999;
	while (war && ww < 200) {

		ww++;
		while (war) {

			war = 0;
			e1 = elements.getElement(wybrany).getE1();
			e2 = elements.getElement(wybrany).getE2();
			e3 = elements.getElement(wybrany).getE3();
			e4 = elements.getElement(wybrany).getE4();


			if (zlyElement.sprawdzElement(e1)) {

				sX = elements.getElement(e1).getSx();
				sY = elements.getElement(e1).getSy();
				sZ = elements.getElement(e1).getSz();
				temp = (sX - x) * (sX - x) + (sY - y) * (sY - y) + (sZ - z) * (sZ - z);
				//temp=sX*sX-2*sX*x+x*x + sY*sY-2*sY*y+y*y + sZ*sZ-2*sZ*z+z*z;
				if (temp <= elements.getElement(e1).getR()) {
					return e1;
				}
				if (max > temp) {
					max = temp;
					wybrany = e1;
					war = 1;
				}

			}
			if (zlyElement.sprawdzElement(e2)) {

				sX = elements.getElement(e2).getSx();
				sY = elements.getElement(e2).getSy();
				sZ = elements.getElement(e2).getSz();
				temp = (sX - x) * (sX - x) + (sY - y) * (sY - y) + (sZ - z) * (sZ - z);
				//temp=sX*sX-2*sX*x+x*x + sY*sY-2*sY*y+y*y + sZ*sZ-2*sZ*z+z*z;
				if (temp <= elements.getElement(e2).getR()) {
					return e2;
				}
				if (max > temp) {
					max = temp;
					wybrany = e2;
					war = 1;
				}
			}
			if (zlyElement.sprawdzElement(e3)) {

				sX = elements.getElement(e3).getSx();
				sY = elements.getElement(e3).getSy();
				sZ = elements.getElement(e3).getSz();
				temp = (sX - x) * (sX - x) + (sY - y) * (sY - y) + (sZ - z) * (sZ - z);
				//temp=sX*sX-2*sX*x+x*x + sY*sY-2*sY*y+y*y + sZ*sZ-2*sZ*z+z*z;
				if (temp <= elements.getElement(e3).getR()) {
					return e3;
				}
				if (max > temp) {
					max = temp;
					wybrany = e3;
					war = 1;
				}
			}
			if (zlyElement.sprawdzElement(e4)) {

				sX = elements.getElement(e4).getSx();
				sY = elements.getElement(e4).getSy();
				sZ = elements.getElement(e4).getSz();
				temp = (sX - x) * (sX - x) + (sY - y) * (sY - y) + (sZ - z) * (sZ - z);
				//temp =sX*sX-2*sX*x+x*x + sY*sY-2*sY*y+y*y + sZ*sZ-2*sZ*z+z*z;
				if (temp <= elements.getElement(e4).getR()) {
					return e4;
				}
				if (max > temp) {
					max = temp;
					wybrany = e4;
					war = 1;
				}
			}
			if (drugiZly) {
				zlyElement.setElement(wybrany);
				drugiZly = false;
			}

		}

		double X = elements.getElement(wybrany).getSx();
		double Y = elements.getElement(wybrany).getSy();
		double Z = elements.getElement(wybrany).getSz();
		float R = elements.getElement(wybrany).getR();
		float temp  = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
//double temp = X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;

		if (temp <= R) {
			return wybrany;
		}
		else {
			zlyElement.setElement(wybrany);
			max = maxStart + 10;
			war = 1;
			drugiZly = true;

			if (zly != -1) {
				wybrany = zlyElement.getElement(zly);
			}
			zly = zlyElement.getIter() - 1;
		}

	}


	X = elements.getElement(wybrany).getSx();
	Y = elements.getElement(wybrany).getSy();
	Z = elements.getElement(wybrany).getSz();
	R = elements.getElement(wybrany).getR();
//int e;


	tempR = (X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z);
//tempR=X*X-2*X*x+x*x + Y*Y-2*Y*y+y*y + Z*Z-2*Z*z+z*z;

	if (tempR <= R) {
		return wybrany;
	}

	else {
		/*
			int nowyW=wybrany;

			pomoc.ustawIter(0);
			pomoc.setElement(nowyW);
			nowyW=wyszukajElementGwiazda(x,y,z,pomoc,20);

			X = elements.getElement(nowyW).getSx();
			Y = elements.getElement(nowyW).getSy();
			Z = elements.getElement(nowyW).getSz();
			R = elements.getElement(nowyW).getR();
			tempR = (X-x) * (X-x) + (Y-y) * (Y-y) + (Z-z) * (Z-z);

			zlyElement.ustawIter(0);

			if(porowRO(tempR,R)){
			//text->Lines->Add(" Pomogl ");
			return nowyW;}
			else{ */
		//m1->Lines->Add(" !!! ERROR  ");
		//m1->Lines->Add(coWchodzi);

		//historia.setElement(x,y,z);
		//historia.setElement(X,Y,Z);

		return -1;
		//}
	}
}


int Delanouy3D::wyszukajElementPoprawaOP(const double & x, const double & y, const double & z, IntList & zlyElement) {
	/*
		for(int i=0;i<100;i++){



		}
	*/
	return 0;
}

int Delanouy3D::wyszukajElementGwiazda(const double & x, const double & y, const double & z, IntList & sprawEle, const int & wielkoscG) {


	int el;
	for (int i = 0, ileE, ileEStare = 0; i < wielkoscG; ++i) {
		ileE = sprawEle.getIter();

		for (int j = ileEStare; j < ileE; ++j) {

			el = sprawEle.getElement(j);
			if (el != -1) {
				sprawEle.setElementUniqat(elements.getElement(el).getE1());
				sprawEle.setElementUniqat(elements.getElement(el).getE2());
				sprawEle.setElementUniqat(elements.getElement(el).getE3());
				sprawEle.setElementUniqat(elements.getElement(el).getE4());
			}
		}

		ileEStare = ileE;
	}

	for (int i = 0, ileElem = elements.getIter(), ileE = sprawEle.getIter(); i < ileE; ++i) {
		if (sprawEle.getElement(i) == -1) {
			sprawEle.usunElement(i);
		}
		if (sprawEle.getElement(i) > ileElem) {
			sprawEle.usunElement(i);
		}
	}

	return  sprawdzenieElementow(x, y, z, sprawEle, 0);
}

void Delanouy3D::reMes() {

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	PunktList temp;
	temp.copyPunktList(points);
//for(int i=8,ileP=points.getIter();i<ileP;i++){temp.setElement(points.getElement(i));}

	delanoyTablicaOP(maxX, maxY, maxZ, true, temp, 0.3, 3, false);

}

void Delanouy3D::reMes(double wX, double wY, double wZ) {


	PunktList temp;
	temp.copyPunktList(points);
//for(int i=8,ileP=points.getIter();i<ileP;i++){temp.setElement(points.getElement(i));}

	delanoyTablicaOP(wX, wY, wZ, true, temp, 0.3, 3, true);

}

void Delanouy3D::wygladzanieLaplace() {

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	int ileP = points.getIter();

	for (int i = 8; i < ileP; i++) {
		if (points.getElement(i).getGranica() != 'g') {

			if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

				double y = 0, z = 0;
				int ile = laplas[i].getIter();

				if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

					for (int j = 0; j < ile; j++) {
						z += points.getElement(laplas[i].getElement(j)).getZ();
					}
					if (ile) {
						points.getElement(i).setZ(z / ile);
					}

				}
				else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

					for (int j = 0; j < ile; j++) {
						y += points.getElement(laplas[i].getElement(j)).getY();
					}
					if (ile) {
						points.getElement(i).setY(y / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {
						y += points.getElement(laplas[i].getElement(j)).getY();
						z += points.getElement(laplas[i].getElement(j)).getZ();
					}
					if (ile) {
						points.getElement(i).setY(y / ile);
						points.getElement(i).setZ(z / ile);
					}
				}

			}
			else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY) {

				double x = 0, z = 0;
				int ile = laplas[i].getIter();

				if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

					for (int j = 0; j < ile; j++) {
						z += points.getElement(laplas[i].getElement(j)).getZ();
					}
					if (ile) {
						points.getElement(i).setZ(z / ile);
					}

				}
				else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
						z += points.getElement(laplas[i].getElement(j)).getZ();
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
						points.getElement(i).setZ(z / ile);
					}

				}

			}
			else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ) {

				double x = 0, y = 0;
				int ile = laplas[i].getIter();

				if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

					for (int j = 0; j < ile; j++) {
						y += points.getElement(laplas[i].getElement(j)).getY();
					}
					if (ile) {
						points.getElement(i).setY(y / ile);
					}

				}
				else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
						y += points.getElement(laplas[i].getElement(j)).getY();
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
						points.getElement(i).setY(y / ile);
					}

				}

			}

			else {

				//if(points.getElement(i).getGranica()=='a'){

				double x = 0, y = 0, z = 0;
				int ile = laplas[i].getIter();

				for (int j = 0; j < ile; j++) {

					x += points.getElement(laplas[i].getElement(j)).getX();
					y += points.getElement(laplas[i].getElement(j)).getY();
					z += points.getElement(laplas[i].getElement(j)).getZ();

				}

				if (ile) {
					points.getElement(i).setPunkt(x / ile, y / ile, z / ile);
				}
				//}
			}

		}//
	}

}



void Delanouy3D::wygladzanieLaplaceWybranePunkty(IntList & punk) {

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	for (int ip = 0, i, ileP = punk.getIter(); ip < ileP; ip++) {

		i = punk.getElement(ip);


		if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

			double y = 0, z = 0;
			int ile = laplas[i].getIter();

			if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

				for (int j = 0; j < ile; j++) {
					z += points.getElement(laplas[i].getElement(j)).getZ();
				}
				if (ile) {
					points.getElement(i).setZ(z / ile);
				}

			}
			else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

				for (int j = 0; j < ile; j++) {
					y += points.getElement(laplas[i].getElement(j)).getY();
				}
				if (ile) {
					points.getElement(i).setY(y / ile);
				}

			}
			else {

				for (int j = 0; j < ile; j++) {
					y += points.getElement(laplas[i].getElement(j)).getY();
					z += points.getElement(laplas[i].getElement(j)).getZ();
				}
				if (ile) {
					points.getElement(i).setY(y / ile);
					points.getElement(i).setZ(z / ile);
				}
			}

		}
		else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY) {

			double x = 0, z = 0;
			int ile = laplas[i].getIter();

			if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

				for (int j = 0; j < ile; j++) {
					z += points.getElement(laplas[i].getElement(j)).getZ();
				}
				if (ile) {
					points.getElement(i).setZ(z / ile);
				}

			}
			else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

				for (int j = 0; j < ile; j++) {
					x += points.getElement(laplas[i].getElement(j)).getX();
				}
				if (ile) {
					points.getElement(i).setX(x / ile);
				}

			}
			else {

				for (int j = 0; j < ile; j++) {
					x += points.getElement(laplas[i].getElement(j)).getX();
					z += points.getElement(laplas[i].getElement(j)).getZ();
				}
				if (ile) {
					points.getElement(i).setX(x / ile);
					points.getElement(i).setZ(z / ile);
				}

			}

		}
		else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ) {

			double x = 0, y = 0;
			int ile = laplas[i].getIter();

			if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

				for (int j = 0; j < ile; j++) {
					y += points.getElement(laplas[i].getElement(j)).getY();
				}
				if (ile) {
					points.getElement(i).setY(y / ile);
				}

			}
			else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

				for (int j = 0; j < ile; j++) {
					x += points.getElement(laplas[i].getElement(j)).getX();
				}
				if (ile) {
					points.getElement(i).setX(x / ile);
				}

			}
			else {

				for (int j = 0; j < ile; j++) {
					x += points.getElement(laplas[i].getElement(j)).getX();
					y += points.getElement(laplas[i].getElement(j)).getY();
				}
				if (ile) {
					points.getElement(i).setX(x / ile);
					points.getElement(i).setY(y / ile);
				}

			}

		}

		else {

			//if(points.getElement(i).getGranica()=='a'){

			double x = 0, y = 0, z = 0;
			int ile = laplas[i].getIter();

			for (int j = 0; j < ile; j++) {

				x += points.getElement(laplas[i].getElement(j)).getX();
				y += points.getElement(laplas[i].getElement(j)).getY();
				z += points.getElement(laplas[i].getElement(j)).getZ();

			}

			if (ile) {
				points.getElement(i).setPunkt(x / ile, y / ile, z / ile);
			}
			//}

		}


	}

}

void Delanouy3D::wygladzanieLaplaceWOparciuOWybranePunktu(IntList & punk, int ilePow) {

	sasiedniePunkty();
	/*
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	int ileP = points.getIter();
	*/
	//zaznaczam punkty nalezoce do wyciagnietej warstwy
	for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {

		p = punk.getElement(i);
		points.getElement(p).setGranica('w');

	}

	//z zaznaczonych punktow usuwam nie zaznaczonych sasiadow
	for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {

		p = punk.getElement(i);
		for (int j = 0, ilel = laplas[p].getIter(), pp; j < ilel; ++j) {

			pp = laplas[p].getElement(j);

			if (points.getElement(pp).getGranica() != 'w') {

				laplas[p].usunElement(j);
				--j;
				--ilel;

			}

		}

	}

	for (int s = 0; s < ilePow; ++s) {

		//wygladzanie w oparciu o wybrane punkty
		for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {
			p = punk.getElement(i);


			double x = 0, y = 0, z = 0;
			int ileL = laplas[p].getIter();

			for (int j = 0; j < ileL; j++) {

				x += points.getElement(laplas[p].getElement(j)).getX();
				y += points.getElement(laplas[p].getElement(j)).getY();
				z += points.getElement(laplas[p].getElement(j)).getZ();

			}

			if (ileL) {
				points.getElement(p).setPunkt(x / ileL, y / ileL, z / ileL);
			}


		}
	}


	//set granica 'a'
	for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {

		p = punk.getElement(i);
		points.getElement(p).setGranica('g');

	}

	sasiedniePunkty();

}

void Delanouy3D::wygladzanieLaplaceWaga(int sposobZbieraniaWag, double waga, int ileWag) {

	DoubleList wagi(100, 100);
	IntList sas0(100, 100), sas1(100, 100);
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	int ileP = points.getIter();

	for (int i = 8; i < ileP; i++) {
		if (points.getElement(i).getGranica() != 'g') {

			if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

				double y = 0, z = 0, yP, zP;;
				int ile = laplas[i].getIter();

				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

					for (int j = 0; j < ile; j++) {
						z += points.getElement(laplas[i].getElement(j)).getZ();
						//    z+=(zP*waga+zP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setZ(z / ile);
					}

				}
				else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

					for (int j = 0; j < ile; j++) {
						y += points.getElement(laplas[i].getElement(j)).getY();
						//   y+=(yP*waga+yP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setY(y / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {
						yP = points.getElement(laplas[i].getElement(j)).getY();
						zP = points.getElement(laplas[i].getElement(j)).getZ();

						y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);
						z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);

					}
					if (ile) {
						points.getElement(i).setY(y / ile);
						points.getElement(i).setZ(z / ile);
					}
				}

			}
			else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY) {

				double x = 0, z = 0, xP, zP;
				int ile = laplas[i].getIter();

				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

					for (int j = 0; j < ile; j++) {
						z += points.getElement(laplas[i].getElement(j)).getZ();
						//   z+=(zP*waga+zP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setZ(z / ile);
					}

				}
				else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
						// x+=(xP*waga+xP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {
						xP = points.getElement(laplas[i].getElement(j)).getX();
						zP = points.getElement(laplas[i].getElement(j)).getZ();

						x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
						z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
						points.getElement(i).setZ(z / ile);
					}

				}

			}
			else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ) {

				double x = 0, y = 0, xP, yP;
				int ile = laplas[i].getIter();

				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

					for (int j = 0; j < ile; j++) {
						y += points.getElement(laplas[i].getElement(j)).getY();
						//y+=(yP*waga+yP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setY(y / ile);
					}

				}
				else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
						//x+=(xP*waga+xP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {

						xP = points.getElement(laplas[i].getElement(j)).getX();
						yP = points.getElement(laplas[i].getElement(j)).getY();

						x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
						y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);

					}
					if (ile) {
						points.getElement(i).setX(x / ile);
						points.getElement(i).setY(y / ile);
					}

				}

			}


			else {

				//if(points.getElement(i).getGranica()!='g'){

				double x = 0, y = 0, z = 0, xP, yP, zP;
				int ile = laplas[i].getIter();


				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				for (int j = 0; j < ile; j++) {

					xP = points.getElement(laplas[i].getElement(j)).getX();
					yP = points.getElement(laplas[i].getElement(j)).getY();
					zP = points.getElement(laplas[i].getElement(j)).getZ();

					x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
					y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);
					z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);
					//x+=xP*wagi.getElement(j);
					//y+=yP*wagi.getElement(j);

				}
				if (ile) {
					points.getElement(i).setPunkt(x / ile, y / ile, z / ile);
				}

				//}
			}

		}//
	}

}

void Delanouy3D::wygladzanieLaplaceWagaWybraneP(int sposobZbieraniaWag, double waga, int ODwybrane, int ileWag) {

	DoubleList wagi(100, 100);
	IntList sas0(100, 100), sas1(100, 100);
	double maxX =  dx;
	double maxY =  dy;
	double maxZ =  dz;

	int ileP = points.getIter();

	for (int i = ODwybrane; i < ileP; i++) {
		if (points.getElement(i).getGranica() != 'g') {

			if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

				double y = 0, z = 0, yP, zP;;
				int ile = laplas[i].getIter();

				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

					for (int j = 0; j < ile; j++) {
						z += points.getElement(laplas[i].getElement(j)).getZ();
						//    z+=(zP*waga+zP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setZ(z / ile);
					}

				}
				else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

					for (int j = 0; j < ile; j++) {
						y += points.getElement(laplas[i].getElement(j)).getY();
						//   y+=(yP*waga+yP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setY(y / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {
						yP = points.getElement(laplas[i].getElement(j)).getY();
						zP = points.getElement(laplas[i].getElement(j)).getZ();

						y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);
						z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);

					}
					if (ile) {
						points.getElement(i).setY(y / ile);
						points.getElement(i).setZ(z / ile);
					}
				}

			}
			else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY) {

				double x = 0, z = 0, xP, zP;
				int ile = laplas[i].getIter();

				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

					for (int j = 0; j < ile; j++) {
						z += points.getElement(laplas[i].getElement(j)).getZ();
						//   z+=(zP*waga+zP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setZ(z / ile);
					}

				}
				else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
						// x+=(xP*waga+xP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {
						xP = points.getElement(laplas[i].getElement(j)).getX();
						zP = points.getElement(laplas[i].getElement(j)).getZ();

						x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
						z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
						points.getElement(i).setZ(z / ile);
					}

				}

			}
			else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ) {

				double x = 0, y = 0, xP, yP;
				int ile = laplas[i].getIter();

				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

					for (int j = 0; j < ile; j++) {
						y += points.getElement(laplas[i].getElement(j)).getY();
						//y+=(yP*waga+yP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setY(y / ile);
					}

				}
				else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

					for (int j = 0; j < ile; j++) {
						x += points.getElement(laplas[i].getElement(j)).getX();
						//x+=(xP*waga+xP*wagi.getElement(j))/(waga+1);
					}
					if (ile) {
						points.getElement(i).setX(x / ile);
					}

				}
				else {

					for (int j = 0; j < ile; j++) {

						xP = points.getElement(laplas[i].getElement(j)).getX();
						yP = points.getElement(laplas[i].getElement(j)).getY();

						x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
						y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);

					}
					if (ile) {
						points.getElement(i).setX(x / ile);
						points.getElement(i).setY(y / ile);
					}

				}

			}


			else {

				//if(points.getElement(i).getGranica()!='g'){

				double x = 0, y = 0, z = 0, xP, yP, zP;
				int ile = laplas[i].getIter();

				zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

				for (int j = 0; j < ile; j++) {

					xP = points.getElement(laplas[i].getElement(j)).getX();
					yP = points.getElement(laplas[i].getElement(j)).getY();
					zP = points.getElement(laplas[i].getElement(j)).getZ();

					x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
					y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);
					z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);
					//x+=xP*wagi.getElement(j);
					//y+=yP*wagi.getElement(j);

				}
				if (ile) {
					points.getElement(i).setPunkt(x / ile, y / ile, z / ile);
				}

				//}
			}

		}//
	}

}

void Delanouy3D::wygladzanieLaplaceWagaWybranePunkty(int sposobZbieraniaWag, double waga, IntList & punk, int ileWag) {

	DoubleList wagi(100, 100);
	IntList sas0(100, 100), sas1(100, 100);
	double maxX =  dx;
	double maxY =  dy;
	double maxZ =  dz;

	//dla wygody dlalem i ktore przyjmuje numer z &punk
	for (int ip = 0, i, ileP = punk.getIter(); ip < ileP; ip++) {

		i = punk.getElement(ip);

		if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

			double y = 0, z = 0, yP, zP;;
			int ile = laplas[i].getIter();

			zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

			if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

				for (int j = 0; j < ile; j++) {
					z += points.getElement(laplas[i].getElement(j)).getZ();
					//    z+=(zP*waga+zP*wagi.getElement(j))/(waga+1);
				}
				if (ile) {
					points.getElement(i).setZ(z / ile);
				}

			}
			else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

				for (int j = 0; j < ile; j++) {
					y += points.getElement(laplas[i].getElement(j)).getY();
					//   y+=(yP*waga+yP*wagi.getElement(j))/(waga+1);
				}
				if (ile) {
					points.getElement(i).setY(y / ile);
				}

			}
			else {

				for (int j = 0; j < ile; j++) {
					yP = points.getElement(laplas[i].getElement(j)).getY();
					zP = points.getElement(laplas[i].getElement(j)).getZ();

					y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);
					z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);

				}
				if (ile) {
					points.getElement(i).setY(y / ile);
					points.getElement(i).setZ(z / ile);
				}
			}

		}
		else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY) {

			double x = 0, z = 0, xP, zP;
			int ile = laplas[i].getIter();

			zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

			if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

				for (int j = 0; j < ile; j++) {
					z += points.getElement(laplas[i].getElement(j)).getZ();
					//   z+=(zP*waga+zP*wagi.getElement(j))/(waga+1);
				}
				if (ile) {
					points.getElement(i).setZ(z / ile);
				}

			}
			else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ ) {

				for (int j = 0; j < ile; j++) {
					x += points.getElement(laplas[i].getElement(j)).getX();
					// x+=(xP*waga+xP*wagi.getElement(j))/(waga+1);
				}
				if (ile) {
					points.getElement(i).setX(x / ile);
				}

			}
			else {

				for (int j = 0; j < ile; j++) {
					xP = points.getElement(laplas[i].getElement(j)).getX();
					zP = points.getElement(laplas[i].getElement(j)).getZ();

					x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
					z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);
				}
				if (ile) {
					points.getElement(i).setX(x / ile);
					points.getElement(i).setZ(z / ile);
				}

			}

		}
		else if (points.getElement(i).getZ() == 0 || points.getElement(i).getZ() == maxZ) {

			double x = 0, y = 0, xP, yP;
			int ile = laplas[i].getIter();

			zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

			if (points.getElement(i).getX() == 0 || points.getElement(i).getX() == maxX ) {

				for (int j = 0; j < ile; j++) {
					y += points.getElement(laplas[i].getElement(j)).getY();
					//y+=(yP*waga+yP*wagi.getElement(j))/(waga+1);
				}
				if (ile) {
					points.getElement(i).setY(y / ile);
				}

			}
			else if (points.getElement(i).getY() == 0 || points.getElement(i).getY() == maxY ) {

				for (int j = 0; j < ile; j++) {
					x += points.getElement(laplas[i].getElement(j)).getX();
					//x+=(xP*waga+xP*wagi.getElement(j))/(waga+1);
				}
				if (ile) {
					points.getElement(i).setX(x / ile);
				}

			}
			else {

				for (int j = 0; j < ile; j++) {

					xP = points.getElement(laplas[i].getElement(j)).getX();
					yP = points.getElement(laplas[i].getElement(j)).getY();

					x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
					y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);

				}
				if (ile) {
					points.getElement(i).setX(x / ile);
					points.getElement(i).setY(y / ile);
				}

			}

		}


		else {

			//if(points.getElement(i).getGranica()!='g'){

			double x = 0, y = 0, z = 0, xP, yP, zP;
			int ile = laplas[i].getIter();


			zbierzWagi(sposobZbieraniaWag, i, wagi, ileWag);

			for (int j = 0; j < ile; j++) {

				xP = points.getElement(laplas[i].getElement(j)).getX();
				yP = points.getElement(laplas[i].getElement(j)).getY();
				zP = points.getElement(laplas[i].getElement(j)).getZ();

				x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
				y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);
				z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);
				//x+=xP*wagi.getElement(j);
				//y+=yP*wagi.getElement(j);

			}
			if (ile) {
				points.getElement(i).setPunkt(x / ile, y / ile, z / ile);
			}

			//}
		}
	}

}

void Delanouy3D::wygladzanieLaplaceWagaWOparciuOWybranePunktu(IntList & punk, int sposobZbieraniaWag, double waga, int ileWag, int ilePow) {

	sasiedniePunkty();

	DoubleList wagi(100, 100);
	IntList sas0(100, 100), sas1(100, 100);
	/*
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	int ileP = points.getIter();
	*/
	//zaznaczam punkty nalezoce do wyciagnietej warstwy
	for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {

		p = punk.getElement(i);
		points.getElement(p).setGranica('w');

	}

	//z zaznaczonych punktow usuwam nie zaznaczonych sasiadow
	for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {

		p = punk.getElement(i);
		for (int j = 0, ilel = laplas[p].getIter(), pp; j < ilel; ++j) {

			pp = laplas[p].getElement(j);

			if (points.getElement(pp).getGranica() != 'w') {

				laplas[p].usunElement(j);
				--j;
				--ilel;

			}

		}

	}


	for (int s = 0; s < ilePow; ++s) {

		//wygladzanie w oparciu o wybrane punkty
		for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {
			p = punk.getElement(i);

			double x = 0, y = 0, z = 0, xP, yP, zP;
			int ileL = laplas[p].getIter();

			zbierzWagi(sposobZbieraniaWag, p, wagi, ileWag);

			for (int j = 0; j < ileL; j++) {

				xP = points.getElement(laplas[p].getElement(j)).getX();
				yP = points.getElement(laplas[p].getElement(j)).getY();
				zP = points.getElement(laplas[p].getElement(j)).getZ();

				x += (xP * waga + xP * wagi.getElement(j)) / (waga + 1);
				y += (yP * waga + yP * wagi.getElement(j)) / (waga + 1);
				z += (zP * waga + zP * wagi.getElement(j)) / (waga + 1);


			}

			if (ileL) {
				points.getElement(p).setPunkt(x / ileL, y / ileL, z / ileL);
			}

		}

	}

	for (int i = 0, ile = punk.getIter(), p; i < ile; ++i) {

		p = punk.getElement(i);
		points.getElement(p).setGranica('g');

	}

	sasiedniePunkty();

}


void Delanouy3D::zbierzWagiOdleglosciKrawedzi(int ktoryP, DoubleList & wagi) {

	wagi.ustawIter(0);
	int p1 = -1, p2 = -1, p3 = -1, p = 0;
	double waga = 0;

	for (int i = 0, ileP = laplas[ktoryP].getIter(), pk; i < ileP; ++i) {

		pk = laplas[ktoryP].getElement(i);
		p1 = -1;
		p2 = -1;
		p = 0;

		for (int j = 0, ilePP = laplas[pk].getIter(); j < ilePP; ++j) {

			if (!laplas[ktoryP].sprCzyJestEleWList(laplas[pk].getElement(j))) {

				if (p1 == -1) {
					p1 = laplas[pk].getElement(j);
				}
				else if (p2 == -1) {
					p2 = laplas[pk].getElement(j);
				}
				else {
					p3 = laplas[pk].getElement(j);
					waga = points.getElement(pk).odlegloscPunktu(points.getElement(p1)) + points.getElement(pk).odlegloscPunktu(points.getElement(p2)) + points.getElement(pk).odlegloscPunktu(points.getElement(p3));
					waga = waga * waga * waga * waga * waga;
					wagi.setElement(waga);
					break;
				}


			}

		}

	}

	waga = wagi.getSrednia();
	wagi.podzielPrzezL(waga);


}

void Delanouy3D::zbierzWagiOdleglosciKrawedziPoprawne(int ktoryP, DoubleList & wagi, int ileWag) {

	wagi.ustawIter(0);
	double waga = 0;

	for (int i = 0, ileP = laplas[ktoryP].getIter(), pk; i < ileP; ++i) {
		pk = laplas[ktoryP].getElement(i);
		waga = 0;

		for (int j = 0, ilePP = laplas[pk].getIter(), p; j < ilePP; ++j) {
			p = laplas[pk].getElement(j);
			if (laplas[ktoryP].sprawdzCzyJest(p)) {

				waga += points.getElement(ktoryP).odlegloscPunktuSqrt(points.getElement(p));

			}
		}

		//waga += points.getElement(ktoryP).odlegloscPunktuSqrt(points.getElement(pk));

		waga = pow(waga, ileWag);
		wagi.setElement(waga);

	}

	waga = wagi.getSrednia();
	wagi.podzielPrzezL(waga);

}

void Delanouy3D::zbierzWagiOdleglosciOdPunktowWagi(int ktoryP, DoubleList & wagi, int ileWag) {

	//IntList krawedzie(50,50);
	wagi.ustawIter(0);
	double waga;

	for (int i = 0, ileP = laplas[ktoryP].getIter(), pk; i < ileP; ++i) {
		pk = laplas[ktoryP].getElement(i);
		waga = points.getElement(ktoryP).odlegloscPunktuSqrt(points.getElement(pk));

		waga = pow(waga, ileWag);
		wagi.setElement(waga);

	}

	waga = wagi.getSrednia();
	wagi.podzielPrzezL(waga);


}



void Delanouy3D::sasiedniePunkty() {
	flagaLaplasa = 2;
	delete []laplas;
	int ileP = points.getIter();
	int ileE = elements.getIter();

	laplas = new IntList[ileP];

	for (int i = 0; i < ileP; i++) {
		laplas[i].czysc(15, 15);
	}

	for (int i = 0, p1, p2, p3, p4; i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		laplas[p1].setElementUniqat(p4);
		laplas[p1].setElementUniqat(p3);
		laplas[p1].setElementUniqat(p2);

		laplas[p2].setElementUniqat(p4);
		laplas[p2].setElementUniqat(p1);
		laplas[p2].setElementUniqat(p3);

		laplas[p3].setElementUniqat(p4);
		laplas[p3].setElementUniqat(p2);
		laplas[p3].setElementUniqat(p1);

		laplas[p4].setElementUniqat(p1);
		laplas[p4].setElementUniqat(p2);
		laplas[p4].setElementUniqat(p3);

	}

}

void Delanouy3D::sasiedniePunktyPryzm() {

	delete []laplas;
	int ileP = points.getIter();
	int ileE = elements.getIter();
	int ilePryzm = pryzmy.getIter();

	laplas = new IntList[ileP];

	for (int i = 0; i < ileP; i++) {
		laplas[i].czysc(100, 100);
	}

	for (int i = 0, p1, p2, p3, p4; i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		laplas[p1].setElementUniqat(p4);
		laplas[p1].setElementUniqat(p3);
		laplas[p1].setElementUniqat(p2);

		laplas[p2].setElementUniqat(p4);
		laplas[p2].setElementUniqat(p1);
		laplas[p2].setElementUniqat(p3);

		laplas[p3].setElementUniqat(p4);
		laplas[p3].setElementUniqat(p2);
		laplas[p3].setElementUniqat(p1);

		laplas[p4].setElementUniqat(p1);
		laplas[p4].setElementUniqat(p2);
		laplas[p4].setElementUniqat(p3);

	}

	for (int i = 0, p1, p2, p3, p4, p5, p6; i < ilePryzm; ++i) {

		p1 = pryzmy.getElement(i).getP1();
		p2 = pryzmy.getElement(i).getP2();
		p3 = pryzmy.getElement(i).getP3();
		p4 = pryzmy.getElement(i).getP4();
		p5 = pryzmy.getElement(i).getP5();
		p6 = pryzmy.getElement(i).getP6();

		laplas[p1].setElementUniqat(p4);
		laplas[p1].setElementUniqat(p2);
		laplas[p1].setElementUniqat(p3);

		laplas[p2].setElementUniqat(p5);
		laplas[p2].setElementUniqat(p1);
		laplas[p2].setElementUniqat(p3);

		laplas[p3].setElementUniqat(p6);
		laplas[p3].setElementUniqat(p1);
		laplas[p3].setElementUniqat(p2);

		laplas[p4].setElementUniqat(p1);
		laplas[p4].setElementUniqat(p5);
		laplas[p4].setElementUniqat(p6);

		laplas[p5].setElementUniqat(p2);
		laplas[p5].setElementUniqat(p4);
		laplas[p5].setElementUniqat(p6);

		laplas[p6].setElementUniqat(p3);
		laplas[p6].setElementUniqat(p4);
		laplas[p6].setElementUniqat(p5);



	}

}



void Delanouy3D::sasiednieElementy() {
	flagaLaplasa = 1;
	delete []laplas;

	int ileP = points.getIter();
	int ileE = elements.getIter();
	laplas = new IntList[ileP];

//m1->Lines->Add(ileP);

	for (int i = 0; i < ileP; i++) {
		laplas[i].czysc(100, 100);
	}

	for (int i = 0, p1, p2, p3, p4; i < ileE; i++) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		//if(p1>ileP || p2>ileP || p3>ileP|| p4>ileP){m1->Lines->Add(p1);}

		laplas[p1].setElement(i);
		laplas[p2].setElement(i);
		laplas[p3].setElement(i);
		laplas[p4].setElement(i);

		//m1->Lines->Add(p1);

	}
}

/*
int Delanouy3D::rysuj(bool ujemne){


double maxX=  dx;
double maxY=  dy;
double maxZ=  dz;
glPushMatrix();


glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);

int zlicz=0;

glBegin(GL_LINES);
glColor3f (1, 1, 1);

if(!ujemne){
	for(int i=0;i<elements.getIter();i++){


		zlicz++;
		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());


		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

	}
}
else{

	for(int i=0;i<elements.getIter();i++){

		if(V_objetoscT(i)<0){
			zlicz++;

			glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
			glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());

			glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
			glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());

			glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
			glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());


			glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
			glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

			glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
			glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

			glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
			glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());
		}
	}
}

glEnd(); //koniec rysowania

glPopMatrix();
glPopMatrix();

return zlicz;

}
*/

/*
int Delanouy3D::rysujWybraneZiarna(IntList &niePokZiarn){

double maxX=  dx;
double maxY=  dy;
double maxZ=  dz;
glPushMatrix();


glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);

int zlicz=0;
glLineWidth(1.5);

//glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
//glEnable(GL_BLEND);

//glEnable(GL_LINE_SMOOTH);
//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
//glEnable(GL_POLYGON_SMOOTH);
//glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);


glBegin(GL_LINES);

glColor3f (1, 1, 1);

for(int i=0;i<elements.getIter();i++){

	if(niePokZiarn.sprCzyJestEleWList(elements.getElement(i).getRodzajZiarna())){

		zlicz++;
		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());


		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

	}

}

glEnd(); //koniec rysowania

glPopMatrix();
glPopMatrix();
glLineWidth(1);


//glDisable(GL_LINE_SMOOTH);
//glDisable(GL_BLEND);
//glDisable(GL_POLYGON_SMOOTH);

return zlicz;


}
*/
/*
int Delanouy3D::rysujWybraneZiarnaPrzekroj(IntList &niePokZiarn,int Px,int Py,int Pz,bool rysuj){

double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
double maxX=  dx;
double maxY=  dy;
double maxZ=  dz;
glPushMatrix();


glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);
glLineWidth(2);
int zlicz=0;

glBegin(GL_LINES);

glColor3f (1, 1, 1);

for(int i=0;i<elements.getIter();i++){

	if(niePokZiarn.sprCzyJestEleWList(elements.getElement(i).getRodzajZiarna())){

		x1 = points.getElement(elements.getElement(i).getP1()).getX();
		x2 = points.getElement(elements.getElement(i).getP2()).getX();
		x3 = points.getElement(elements.getElement(i).getP3()).getX();
		x4 = points.getElement(elements.getElement(i).getP4()).getX();

		y1 = points.getElement(elements.getElement(i).getP1()).getY();
		y2 = points.getElement(elements.getElement(i).getP2()).getY();
		y3 = points.getElement(elements.getElement(i).getP3()).getY();
		y4 = points.getElement(elements.getElement(i).getP4()).getY();

		z1 = points.getElement(elements.getElement(i).getP1()).getZ();
		z2 = points.getElement(elements.getElement(i).getP2()).getZ();
		z3 = points.getElement(elements.getElement(i).getP3()).getZ();
		z4 = points.getElement(elements.getElement(i).getP4()).getZ();


		if(elements.getElement(i).getRodzajZiarna()==2 && rysuj){
			//glLineWidth(1);
			glColor3f (1, 0.5, 0.5);
			glVertex3f(x1,y1,z1);glVertex3f(x2,y2,z2);glVertex3f(x2,y2,z2);glVertex3f(x3,y3,z3);
			glVertex3f(x3,y3,z3);glVertex3f(x1,y1,z1);glVertex3f(x1,y1,z1);glVertex3f(x4,y4,z4);
			glVertex3f(x2,y2,z2);glVertex3f(x4,y4,z4);glVertex3f(x3,y3,z3);glVertex3f(x4,y4,z4);
			glColor3f (1, 1, 1);
			//glLineWidth(2);

		}
		else if(x1<Px && x2<Px && x3<Px && x4<Px && y1<Py && y2<Py && y3<Py && y4<Py && z1<Pz && z2<Pz && z3<Pz && z4<Pz){


			glVertex3f(x1,y1,z1);glVertex3f(x2,y2,z2);glVertex3f(x2,y2,z2);glVertex3f(x3,y3,z3);
			glVertex3f(x3,y3,z3);glVertex3f(x1,y1,z1);glVertex3f(x1,y1,z1);glVertex3f(x4,y4,z4);
			glVertex3f(x2,y2,z2);glVertex3f(x4,y4,z4);glVertex3f(x3,y3,z3);glVertex3f(x4,y4,z4);

		zlicz++;
		}
	}

}

glEnd(); //koniec rysowania

glPopMatrix();
glPopMatrix();

glLineWidth(1);

return zlicz;


}
*/

double Delanouy3D::V_objetosc(FILE * out) {

	return V_objetoscT(out) + V_objetoscP();

}

double Delanouy3D::V_objetoscT(FILE * out) {
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, wynik, wynikAll = 0;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		x1 = points.getElement(elements.getElement(i).getP1()).getX();
		x2 = points.getElement(elements.getElement(i).getP2()).getX();
		x3 = points.getElement(elements.getElement(i).getP3()).getX();
		x4 = points.getElement(elements.getElement(i).getP4()).getX();

		y1 = points.getElement(elements.getElement(i).getP1()).getY();
		y2 = points.getElement(elements.getElement(i).getP2()).getY();
		y3 = points.getElement(elements.getElement(i).getP3()).getY();
		y4 = points.getElement(elements.getElement(i).getP4()).getY();

		z1 = points.getElement(elements.getElement(i).getP1()).getZ();
		z2 = points.getElement(elements.getElement(i).getP2()).getZ();
		z3 = points.getElement(elements.getElement(i).getP3()).getZ();
		z4 = points.getElement(elements.getElement(i).getP4()).getZ();



		wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
				(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

		if (wynik < 0) {
			fprintf(out, "error volume < 0 id Element :%d \n", i + 1);
			wynik *= -1;
		}

		wynikAll += wynik / 6;

	}

	return wynikAll;
}

double Delanouy3D::V_objetoscTPopraw() {

	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, wynik, wynikAll = 0;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		x1 = points.getElement(elements.getElement(i).getP1()).getX();
		x2 = points.getElement(elements.getElement(i).getP2()).getX();
		x3 = points.getElement(elements.getElement(i).getP3()).getX();
		x4 = points.getElement(elements.getElement(i).getP4()).getX();

		y1 = points.getElement(elements.getElement(i).getP1()).getY();
		y2 = points.getElement(elements.getElement(i).getP2()).getY();
		y3 = points.getElement(elements.getElement(i).getP3()).getY();
		y4 = points.getElement(elements.getElement(i).getP4()).getY();

		z1 = points.getElement(elements.getElement(i).getP1()).getZ();
		z2 = points.getElement(elements.getElement(i).getP2()).getZ();
		z3 = points.getElement(elements.getElement(i).getP3()).getZ();
		z4 = points.getElement(elements.getElement(i).getP4()).getZ();



		wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
				(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

		if (wynik < 0) {

			//m1->Lines->Add("zleTetra");

			x3 = elements.getElement(i).getP2();
			elements.getElement(i).setP2(elements.getElement(i).getP3());
			elements.getElement(i).setP3(x3);

			wynik *= -1;
		}

		wynikAll += wynik / 6;

	}

	return wynikAll;
}

double Delanouy3D::V_objetoscP() {
	double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6, z1, z2, z3, z4, z5, z6, wynik1, wynik2, wynik3, wynikAll = 0;
	double xa, ya, za, wynik4, wynik5;

	for (int i = 0, ileE = pryzmy.getIter(); i < ileE; i++) {

		x1 = points.getElement(pryzmy.getElement(i).getP1()).getX();
		x2 = points.getElement(pryzmy.getElement(i).getP3()).getX();
		x3 = points.getElement(pryzmy.getElement(i).getP2()).getX();
		x4 = points.getElement(pryzmy.getElement(i).getP4()).getX();
		x5 = points.getElement(pryzmy.getElement(i).getP6()).getX();
		x6 = points.getElement(pryzmy.getElement(i).getP5()).getX();

		y1 = points.getElement(pryzmy.getElement(i).getP1()).getY();
		y2 = points.getElement(pryzmy.getElement(i).getP3()).getY();
		y3 = points.getElement(pryzmy.getElement(i).getP2()).getY();
		y4 = points.getElement(pryzmy.getElement(i).getP4()).getY();
		y5 = points.getElement(pryzmy.getElement(i).getP6()).getY();
		y6 = points.getElement(pryzmy.getElement(i).getP5()).getY();

		z1 = points.getElement(pryzmy.getElement(i).getP1()).getZ();
		z2 = points.getElement(pryzmy.getElement(i).getP3()).getZ();
		z3 = points.getElement(pryzmy.getElement(i).getP2()).getZ();
		z4 = points.getElement(pryzmy.getElement(i).getP4()).getZ();
		z5 = points.getElement(pryzmy.getElement(i).getP6()).getZ();
		z6 = points.getElement(pryzmy.getElement(i).getP5()).getZ();

		xa = (x1 + x4 + x6 + x3) * 0.25;
		ya = (y1 + y4 + y6 + y3) * 0.25;
		za = (z1 + z4 + z6 + z3) * 0.25;
		/*
			3 2 1 - 5
			0 7 8 - 9

			wynik = (x7-x0)*(y8-y0)*(z9-z0)+(x8-x0)*(y9-y0)*(z7-z0)+(x9-x0)*(y7-y0)*(z8-z0)-
			(z7-z0)*(y8-y0)*(x9-x0)-(z8-z0)*(y9-y0)*(x7-x0)-(z9-z0)*(y7-y0)*(x8-x0);

			tetra
			wynik = (x2-x1)*(y3-y1)*(z4-z1)+(x3-x1)*(y4-y1)*(z2-z1)+(x4-x1)*(y2-y1)*(z3-z1)-
			(z2-z1)*(y3-y1)*(x4-x1)-(z3-z1)*(y4-y1)*(x2-x1)-(z4-z1)*(y2-y1)*(x3-x1) ;

		*/

		wynik1 = (x2 - x3) * (y1 - y3) * (z5 - z3) + (x1 - x3) * (y5 - y3) * (z2 - z3) + (x5 - x3) * (y2 - y3) * (z1 - z3) -
				(z2 - z3) * (y1 - y3) * (x5 - x3) - (z1 - z3) * (y5 - y3) * (x2 - x3) - (z5 - z3) * (y2 - y3) * (x1 - x3);

		/*
			wynik2 = (x5-x4)*(y6-y4)*(z1-z4)+(x6-x4)*(y1-y4)*(z5-z4)+(x1-x4)*(y5-y4)*(z6-z4)-
			(z5-z4)*(y6-y4)*(x1-x4)-(z6-z4)*(y1-y4)*(x5-x4)-(z1-z4)*(y5-y4)*(x6-x4);

			wynik3 = (x5-x3)*(y1-y3)*(z6-z3)+(x1-x3)*(y6-y3)*(z5-z3)+(x6-x3)*(y5-y3)*(z1-z3)-
			(z5-z3)*(y1-y3)*(x6-x3)-(z1-z3)*(y6-y3)*(x5-x3)-(z6-z3)*(y5-y3)*(x1-x3);
		*/

		wynik2 = (x4 - x1) * (ya - y1) * (z5 - z1) + (xa - x1) * (y5 - y1) * (z4 - z1) + (x5 - x1) * (y4 - y1) * (za - z1) -
				(z4 - z1) * (ya - y1) * (x5 - x1) - (za - z1) * (y5 - y1) * (x4 - x1) - (z5 - z1) * (y4 - y1) * (xa - x1);

		wynik3 = (x6 - x4) * (ya - y4) * (z5 - z4) + (xa - x4) * (y5 - y4) * (z6 - z4) + (x5 - x4) * (y6 - y4) * (za - z4) -
				(z6 - z4) * (ya - y4) * (x5 - x4) - (za - z4) * (y5 - y4) * (x6 - x4) - (z5 - z4) * (y6 - y4) * (xa - x4);

		wynik4 = (x6 - xa) * (y3 - ya) * (z5 - za) + (x3 - xa) * (y5 - ya) * (z6 - za) + (x5 - xa) * (y6 - ya) * (z3 - za) -
				(z6 - za) * (y3 - ya) * (x5 - xa) - (z3 - za) * (y5 - ya) * (x6 - xa) - (z5 - za) * (y6 - ya) * (x3 - xa);

		wynik5 = (xa - x1) * (y3 - y1) * (z5 - z1) + (x3 - x1) * (y5 - y1) * (za - z1) + (x5 - x1) * (ya - y1) * (z3 - z1) -
				(za - z1) * (y3 - y1) * (x5 - x1) - (z3 - z1) * (y5 - y1) * (xa - x1) - (z5 - z1) * (ya - y1) * (x3 - x1);


		if (wynik1 < 0) {
			//m1->Lines->Add("zle1");
			wynik1 *= -1;
		}
		else {;}

		if (wynik2 < 0) {
			//m1->Lines->Add("zle2");
			wynik2 *= -1;
		}
		else {;}

		if (wynik3 < 0) {
			//m1->Lines->Add("zle3");
			wynik3 *= -1;
		}
		else {;}

		if (wynik4 < 0) {
			//m1->Lines->Add("zle4");
			wynik4 *= -1;
		}
		else {;}

		if (wynik5 < 0) {
			//m1->Lines->Add("zle5");
			wynik5 *= -1;
		}
		else {;}

		wynikAll += (wynik1 / 6 + wynik2 / 6 + wynik3 / 6 + wynik4 / 6 + wynik5 / 6);

	}

	return wynikAll;
}

double Delanouy3D::V_objetoscT(int i) {

	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, wynik;

	x1 = points.getElement(elements.getElement(i).getP1()).getX();
	x2 = points.getElement(elements.getElement(i).getP2()).getX();
	x3 = points.getElement(elements.getElement(i).getP3()).getX();
	x4 = points.getElement(elements.getElement(i).getP4()).getX();

	y1 = points.getElement(elements.getElement(i).getP1()).getY();
	y2 = points.getElement(elements.getElement(i).getP2()).getY();
	y3 = points.getElement(elements.getElement(i).getP3()).getY();
	y4 = points.getElement(elements.getElement(i).getP4()).getY();

	z1 = points.getElement(elements.getElement(i).getP1()).getZ();
	z2 = points.getElement(elements.getElement(i).getP2()).getZ();
	z3 = points.getElement(elements.getElement(i).getP3()).getZ();
	z4 = points.getElement(elements.getElement(i).getP4()).getZ();


	wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
			(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

	return wynik / 6;

}

bool Delanouy3D::V_objetoscP(int i) {
	bool poprawny = true;
	double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6, z1, z2, z3, z4, z5, z6, wynik1, wynik2, wynik3, wynikAll = 0;
	double xa, ya, za, wynik4, wynik5;

	x1 = points.getElement(pryzmy.getElement(i).getP1()).getX();
	x2 = points.getElement(pryzmy.getElement(i).getP3()).getX();
	x3 = points.getElement(pryzmy.getElement(i).getP2()).getX();
	x4 = points.getElement(pryzmy.getElement(i).getP4()).getX();
	x5 = points.getElement(pryzmy.getElement(i).getP6()).getX();
	x6 = points.getElement(pryzmy.getElement(i).getP5()).getX();

	y1 = points.getElement(pryzmy.getElement(i).getP1()).getY();
	y2 = points.getElement(pryzmy.getElement(i).getP3()).getY();
	y3 = points.getElement(pryzmy.getElement(i).getP2()).getY();
	y4 = points.getElement(pryzmy.getElement(i).getP4()).getY();
	y5 = points.getElement(pryzmy.getElement(i).getP6()).getY();
	y6 = points.getElement(pryzmy.getElement(i).getP5()).getY();

	z1 = points.getElement(pryzmy.getElement(i).getP1()).getZ();
	z2 = points.getElement(pryzmy.getElement(i).getP3()).getZ();
	z3 = points.getElement(pryzmy.getElement(i).getP2()).getZ();
	z4 = points.getElement(pryzmy.getElement(i).getP4()).getZ();
	z5 = points.getElement(pryzmy.getElement(i).getP6()).getZ();
	z6 = points.getElement(pryzmy.getElement(i).getP5()).getZ();

	xa = (x1 + x4 + x6 + x3) * 0.25;
	ya = (y1 + y4 + y6 + y3) * 0.25;
	za = (z1 + z4 + z6 + z3) * 0.25;
	/*
		3 2 1 - 5
		0 7 8 - 9

		wynik = (x7-x0)*(y8-y0)*(z9-z0)+(x8-x0)*(y9-y0)*(z7-z0)+(x9-x0)*(y7-y0)*(z8-z0)-
		(z7-z0)*(y8-y0)*(x9-x0)-(z8-z0)*(y9-y0)*(x7-x0)-(z9-z0)*(y7-y0)*(x8-x0);

		tetra
		wynik = (x2-x1)*(y3-y1)*(z4-z1)+(x3-x1)*(y4-y1)*(z2-z1)+(x4-x1)*(y2-y1)*(z3-z1)-
		(z2-z1)*(y3-y1)*(x4-x1)-(z3-z1)*(y4-y1)*(x2-x1)-(z4-z1)*(y2-y1)*(x3-x1) ;

	*/

	wynik1 = (x2 - x3) * (y1 - y3) * (z5 - z3) + (x1 - x3) * (y5 - y3) * (z2 - z3) + (x5 - x3) * (y2 - y3) * (z1 - z3) -
			(z2 - z3) * (y1 - y3) * (x5 - x3) - (z1 - z3) * (y5 - y3) * (x2 - x3) - (z5 - z3) * (y2 - y3) * (x1 - x3);

	/*
		wynik2 = (x5-x4)*(y6-y4)*(z1-z4)+(x6-x4)*(y1-y4)*(z5-z4)+(x1-x4)*(y5-y4)*(z6-z4)-
		(z5-z4)*(y6-y4)*(x1-x4)-(z6-z4)*(y1-y4)*(x5-x4)-(z1-z4)*(y5-y4)*(x6-x4);

		wynik3 = (x5-x3)*(y1-y3)*(z6-z3)+(x1-x3)*(y6-y3)*(z5-z3)+(x6-x3)*(y5-y3)*(z1-z3)-
		(z5-z3)*(y1-y3)*(x6-x3)-(z1-z3)*(y6-y3)*(x5-x3)-(z6-z3)*(y5-y3)*(x1-x3);
	*/

	wynik2 = (x4 - x1) * (ya - y1) * (z5 - z1) + (xa - x1) * (y5 - y1) * (z4 - z1) + (x5 - x1) * (y4 - y1) * (za - z1) -
			(z4 - z1) * (ya - y1) * (x5 - x1) - (za - z1) * (y5 - y1) * (x4 - x1) - (z5 - z1) * (y4 - y1) * (xa - x1);

	wynik3 = (x6 - x4) * (ya - y4) * (z5 - z4) + (xa - x4) * (y5 - y4) * (z6 - z4) + (x5 - x4) * (y6 - y4) * (za - z4) -
			(z6 - z4) * (ya - y4) * (x5 - x4) - (za - z4) * (y5 - y4) * (x6 - x4) - (z5 - z4) * (y6 - y4) * (xa - x4);

	wynik4 = (x6 - xa) * (y3 - ya) * (z5 - za) + (x3 - xa) * (y5 - ya) * (z6 - za) + (x5 - xa) * (y6 - ya) * (z3 - za) -
			(z6 - za) * (y3 - ya) * (x5 - xa) - (z3 - za) * (y5 - ya) * (x6 - xa) - (z5 - za) * (y6 - ya) * (x3 - xa);

	wynik5 = (xa - x1) * (y3 - y1) * (z5 - z1) + (x3 - x1) * (y5 - y1) * (za - z1) + (x5 - x1) * (ya - y1) * (z3 - z1) -
			(za - z1) * (y3 - y1) * (x5 - x1) - (z3 - z1) * (y5 - y1) * (xa - x1) - (z5 - z1) * (ya - y1) * (x3 - x1);


	if (wynik1 < 0) {
		//m1->Lines->Add("zle1");
		wynik1 *= -1;
	}
	else {
		poprawny = false;
	}

	if (wynik2 < 0) {
		//m1->Lines->Add("zle2");
		wynik2 *= -1;
	}
	else {
		poprawny = false;
	}

	if (wynik3 < 0) {
		//m1->Lines->Add("zle3");
		wynik3 *= -1;
	}
	else {
		poprawny = false;
	}

	if (wynik4 < 0) {
		//m1->Lines->Add("zle4");
		wynik4 *= -1;
	}
	else {
		poprawny = false;
	}

	if (wynik5 < 0) {
		//m1->Lines->Add("zle5");
		wynik5 *= -1;
	}
	else {
		poprawny = false;
	}

	wynikAll += (wynik1 / 6 + wynik2 / 6 + wynik3 / 6 + wynik4 / 6 + wynik5 / 6);



	return poprawny;

}


void Delanouy3D::ustalSSE() {

	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		x1 = points.getElement(elements.getElement(i).getP1()).getX();
		x2 = points.getElement(elements.getElement(i).getP2()).getX();
		x3 = points.getElement(elements.getElement(i).getP3()).getX();
		x4 = points.getElement(elements.getElement(i).getP4()).getX();

		y1 = points.getElement(elements.getElement(i).getP1()).getY();
		y2 = points.getElement(elements.getElement(i).getP2()).getY();
		y3 = points.getElement(elements.getElement(i).getP3()).getY();
		y4 = points.getElement(elements.getElement(i).getP4()).getY();

		z1 = points.getElement(elements.getElement(i).getP1()).getZ();
		z2 = points.getElement(elements.getElement(i).getP2()).getZ();
		z3 = points.getElement(elements.getElement(i).getP3()).getZ();
		z4 = points.getElement(elements.getElement(i).getP4()).getZ();

		elements.getElement(i).setSSE((x1 + x2 + x3 + x4) * 0.25, (y1 + y2 + y3 + y4) * 0.25, (z1 + z2 + z3 + z4) * 0.25);
	}

}

int Delanouy3D::dopasujElementDoZiarna(bool pokazWynik) {

	wybraneZleEl.ustawIter(0);
	int bardzozle = 0;

//ustalenie elementu zerowego
	int elZero = 0;

	for (int w = 0, ileE = elements.getIter(); w < ileE; w++) {
		if (elements.getElement(w).getP2() == 0 || elements.getElement(w).getP3() == 0) {
			elZero = w;
		}
		elements.getElement(w).setRodzajZiarna(-1);
	}

	int numerZ = 0;

	IntList runda(1450, 500);
	IntList pomoc(1300, 500);
	elZero = 0;
	for (int w = 0, p1, p2, p3, p4, ileE = elements.getIter(); w < ileE; w++) {


		if (elements.getElement(w).getRodzajZiarna() == -1) {
			runda.ustawIter(0);

			p1 = elements.getElement(w).getP1();
			p2 = elements.getElement(w).getP2();
			p3 = elements.getElement(w).getP3();
			p4 = elements.getElement(w).getP4();

			if ( points.getElement(p1).getGranica() != 'a' && points.getElement(p2).getGranica() != 'a' && points.getElement(p3).getGranica() != 'a' && points.getElement(p4).getGranica() != 'a') {
				bardzozle++;
				wybraneZleEl.setElement(w);
			}
			else {
				runda.setElement(w);
				numerZ++;
			}

			while (runda.getIter()) {

				for (int i = 0, el; i < runda.getIter(); i++) {
					el = runda.getElement(i);

					if (elements.getElement(el).getRodzajZiarna() == -1) {
						elements.getElement(el).setRodzajZiarna(numerZ);

						int n1 = 0, n2 = 0, n3 = 0, n4 = 0;

						p1 = elements.getElement(el).getP1();
						p2 = elements.getElement(el).getP2();
						p3 = elements.getElement(el).getP3();
						p4 = elements.getElement(el).getP4();

						if (points.getElement(p1).getGranica() != 'a') {
							n1 = 1;
						}
						if (points.getElement(p2).getGranica() != 'a') {
							n2 = 1;
						}
						if (points.getElement(p3).getGranica() != 'a') {
							n3 = 1;
						}
						if (points.getElement(p4).getGranica() != 'a') {
							n4 = 1;
						}


						if (n4 + n2 + n1 < 3) {
							if (elements.getElement(el).getE1() != -1) {
								pomoc.setElement(elements.getElement(el).getE1());
							}
						}
						if (n4 + n3 + n2 < 3) {
							if (elements.getElement(el).getE2() != -1) {
								pomoc.setElement(elements.getElement(el).getE2());
							}
						}
						if (n4 + n1 + n3 < 3) {
							if (elements.getElement(el).getE3() != -1) {
								pomoc.setElement(elements.getElement(el).getE3());
							}
						}
						if (n1 + n2 + n3 < 3) {
							if (elements.getElement(el).getE4() != -1) {
								pomoc.setElement(elements.getElement(el).getE4());
							}
						}
					}

				}

				runda.ustawIter(0);

				for (int i = 0; i < pomoc.getIter(); i++) {
					runda.setElement(pomoc.getElement(i));
				}

				pomoc.ustawIter(0);

				//koniec while
			}

		}

		elZero = w;
//koniec for
	}


	IntList zleE(10000, 5000);
	int zlychE = 0;
	if (pokazWynik) { ;}

	/*

	 for(int i=0;i<elements.getIter();i++){
		if(elements.getElement(i).getRodzajZiarna()==-1){zleE.setElement(i);}
	 }
	*/
//poprawaDopasowaniaZiaren2Sas(wybraneZleEl);
//poprawaDopasowaniaZiarenNajPow();
//poprawaDopasowaniaZiaren2Sas(wybraneZleEl);
//poprawaDopasowaniaZiarenNajPow();
//poprawaDopasowaniaZiaren2Sas(wybraneZleEl);



	for (int i = 0; i < elements.getIter(); i++) {
		if (elements.getElement(i).getRodzajZiarna() == -1) {
			++zlychE;
		}
	}

	if (pokazWynik) { ;}
	if (pokazWynik) { ;}

	return numerZ;

}

int Delanouy3D::dopasujElementDoZiarnaPunktyStartowe(PunktList & startP) {

	wybraneZleEl.ustawIter(0);
	int bardzozle = 0;

//ustalenie elementu zerowego
	int elZero = 0;

	for (int w = 0, ileE = elements.getIter(); w < ileE; w++) {
		if (elements.getElement(w).getP2() == 0 || elements.getElement(w).getP3() == 0) {
			elZero = w;
		}
		elements.getElement(w).setRodzajZiarna(-1);
	}

	int numerZ = 0;

	IntList runda(1450, 500);
	IntList pomoc(1300, 500);
//dziwne
	IntList zlyElement(1300, 500);

	elZero = 0;

	for (int i = 0, w, p1, p2, p3, p4, ileE = startP.getIter(); i < ileE; i++) {

		//element
		w = wyszukajElementOP(-1, startP.getElement(i).getX(), startP.getElement(i).getY(), startP.getElement(i).getZ(), zlyElement, pomoc);

		if (elements.getElement(w).getRodzajZiarna() == -1) {
			runda.ustawIter(0);

			p1 = elements.getElement(w).getP1();
			p2 = elements.getElement(w).getP2();
			p3 = elements.getElement(w).getP3();
			p4 = elements.getElement(w).getP4();

			if ( points.getElement(p1).getGranica() != 'a' && points.getElement(p2).getGranica() != 'a' && points.getElement(p3).getGranica() != 'a' && points.getElement(p4).getGranica() != 'a') {
				bardzozle++;
				wybraneZleEl.setElement(w);
			}
			else {
				runda.setElement(w);
				numerZ++;
			}

			while (runda.getIter()) {

				for (int i = 0, el; i < runda.getIter(); i++) {
					el = runda.getElement(i);

					if (elements.getElement(el).getRodzajZiarna() == -1) {
						elements.getElement(el).setRodzajZiarna(numerZ);

						int n1 = 0, n2 = 0, n3 = 0, n4 = 0;

						p1 = elements.getElement(el).getP1();
						p2 = elements.getElement(el).getP2();
						p3 = elements.getElement(el).getP3();
						p4 = elements.getElement(el).getP4();

						if (points.getElement(p1).getGranica() != 'a') {
							n1 = 1;
						}
						if (points.getElement(p2).getGranica() != 'a') {
							n2 = 1;
						}
						if (points.getElement(p3).getGranica() != 'a') {
							n3 = 1;
						}
						if (points.getElement(p4).getGranica() != 'a') {
							n4 = 1;
						}


						if (n4 + n2 + n1 < 3) {
							if (elements.getElement(el).getE1() != -1) {
								pomoc.setElement(elements.getElement(el).getE1());
							}
						}
						if (n4 + n3 + n2 < 3) {
							if (elements.getElement(el).getE2() != -1) {
								pomoc.setElement(elements.getElement(el).getE2());
							}
						}
						if (n4 + n1 + n3 < 3) {
							if (elements.getElement(el).getE3() != -1) {
								pomoc.setElement(elements.getElement(el).getE3());
							}
						}
						if (n1 + n2 + n3 < 3) {
							if (elements.getElement(el).getE4() != -1) {
								pomoc.setElement(elements.getElement(el).getE4());
							}
						}
					}

				}

				runda.ustawIter(0);

				for (int i = 0; i < pomoc.getIter(); i++) {
					runda.setElement(pomoc.getElement(i));
				}

				pomoc.ustawIter(0);

				//koniec while
			}

		}

		elZero = w;
//koniec for
	}

	int zlychE = 0;

	for (int i = 0, ileZ = startP.getIter() + 1; i < elements.getIter(); i++) {

		if (elements.getElement(i).getRodzajZiarna() == -1) {
			++zlychE;
			elements.getElement(i).setRodzajZiarna(ileZ);
		}

	}



//m1->Lines->Add(bardzozle);
//m1->Lines->Add(zlychE);

	return numerZ;

}

void Delanouy3D::poprawaDopasowaniaZiaren2Sas(IntList & zleE) {

	for (int i = 0, a, el, z1, z2, z3, z4, ileE = zleE.getIter(); i < ileE; ++i) {

		el =  zleE.getElement(i);

		z1 = elements.getElement(elements.getElement(el).getE1()).getRodzajZiarna();
		z2 = elements.getElement(elements.getElement(el).getE2()).getRodzajZiarna();
		z3 = elements.getElement(elements.getElement(el).getE3()).getRodzajZiarna();
		z4 = elements.getElement(elements.getElement(el).getE4()).getRodzajZiarna();

		a = 1;
		while (a--) {
			if (z1 == z2) {
				if (z1 != -1) {
					elements.getElement(el).setRodzajZiarna(z1);
					break;
				}
			}
			if (z1 == z3) {
				if (z1 != -1) {
					elements.getElement(el).setRodzajZiarna(z1);
					break;
				}
			}
			if (z1 == z4) {
				if (z1 != -1) {
					elements.getElement(el).setRodzajZiarna(z1);
					break;
				}
			}
			if (z2 == z3) {
				if (z2 != -1) {
					elements.getElement(el).setRodzajZiarna(z2);
					break;
				}
			}
			if (z2 == z4) {
				if (z2 != -1) {
					elements.getElement(el).setRodzajZiarna(z2);
					break;
				}
			}
			if (z3 == z4) {
				if (z3 != -1) {
					elements.getElement(el).setRodzajZiarna(z3);
					break;
				}
			}
		}
	}


	for (int i = 0, w; i < zleE.getIter(); ++i) {
		w = zleE.getElement(i);
		if (elements.getElement(w).getRodzajZiarna() == -1) {
			zleE.usunElement(i);
			--i;
		}
	}

}


void Delanouy3D::objetoscSrodAryt(double obj[4], const int & p1, const int & p2, const int & p3, const int & p4) {


	double x1, x2, x3, x4, x5, y1, y2, y3, y4, y5, z1, z2, z3, z4, z5, wynik;

	x1 = points.getElement(p1).getX();
	x2 = points.getElement(p2).getX();
	x3 = points.getElement(p3).getX();
	x4 = points.getElement(p4).getX();
	x5 = (x1 + x2 + x3 + x4) * 0.25;

	y1 = points.getElement(p1).getY();
	y2 = points.getElement(p2).getY();
	y3 = points.getElement(p3).getY();
	y4 = points.getElement(p4).getY();
	y5 = (y1 + y2 + y3 + y4) * 0.25;

	z1 = points.getElement(p1).getZ();
	z2 = points.getElement(p2).getZ();
	z3 = points.getElement(p3).getZ();
	z4 = points.getElement(p4).getZ();
	z5 = (z1 + z2 + z3 + z4) * 0.25;



	//4215
	obj[0] = ((x2 - x4) * (y1 - y4) * (z5 - z4) + (x1 - x4) * (y5 - y4) * (z2 - z4) + (x5 - x4) * (y2 - y4) * (z1 - z4) -
					(z2 - z4) * (y1 - y4) * (x5 - x4) - (z1 - z4) * (y5 - y4) * (x2 - x4) - (z5 - z4) * (y2 - y4) * (x1 - x4)) / 6 ;

	//4325
	obj[1] = ((x3 - x4) * (y2 - y4) * (z5 - z4) + (x2 - x4) * (y5 - y4) * (z3 - z4) + (x5 - x4) * (y3 - y4) * (z2 - z4) -
					(z3 - z4) * (y2 - y4) * (x5 - x4) - (z2 - z4) * (y5 - y4) * (x3 - x4) - (z5 - z4) * (y3 - y4) * (x2 - x4)) / 6 ;

	//4135
	obj[2] = ((x1 - x4) * (y3 - y4) * (z5 - z4) + (x3 - x4) * (y5 - y4) * (z1 - z4) + (x5 - x4) * (y1 - y4) * (z3 - z4) -
					(z1 - z4) * (y3 - y4) * (x5 - x4) - (z3 - z4) * (y5 - y4) * (x1 - x4) - (z5 - z4) * (y1 - y4) * (x3 - x4)) / 6 ;

	//1235
	obj[3] = ((x2 - x1) * (y3 - y1) * (z5 - z1) + (x3 - x1) * (y5 - y1) * (z2 - z1) + (x5 - x1) * (y2 - y1) * (z3 - z1) -
					(z2 - z1) * (y3 - y1) * (x5 - x1) - (z3 - z1) * (y5 - y1) * (x2 - x1) - (z5 - z1) * (y2 - y1) * (x3 - x1)) / 6 ;

	wynik = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
					(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6 ;

	if (wynik == obj[0] + obj[1] + obj[2] + obj[3]) {}
	else {;}

}

void Delanouy3D::krawedzieSrodAryt(double kra[4], const int & p1, const int & p2, const int & p3, const int & p4) {

	//1-2 a 2-3 b 3-1 c 1-4 d 2-4 e 3-4 f
	//double pol421,pol432,pol413,pol123;

	double x1, x2, x3, x4, x5, y1, y2, y3, y4, y5, z1, z2, z3, z4, z5;
	double dkp1, dkp2, dkp3, dkp4;

	x1 = points.getElement(p1).getX();
	x2 = points.getElement(p2).getX();
	x3 = points.getElement(p3).getX();
	x4 = points.getElement(p4).getX();
	x5 = (x1 + x2 + x3 + x4) * 0.25;

	y1 = points.getElement(p1).getY();
	y2 = points.getElement(p2).getY();
	y3 = points.getElement(p3).getY();
	y4 = points.getElement(p4).getY();
	y5 = (y1 + y2 + y3 + y4) * 0.25;

	z1 = points.getElement(p1).getZ();
	z2 = points.getElement(p2).getZ();
	z3 = points.getElement(p3).getZ();
	z4 = points.getElement(p4).getZ();
	z5 = (z1 + z2 + z3 + z4) * 0.25;

	dkp1 = sqrt((x5 - x1) * (x5 - x1) + (y5 - y1) * (y5 - y1) + (z5 - z1) * (z5 - z1));
	dkp2 = sqrt((x5 - x2) * (x5 - x2) + (y5 - y2) * (y5 - y2) + (z5 - z2) * (z5 - z2));
	dkp3 = sqrt((x5 - x3) * (x5 - x3) + (y5 - y3) * (y5 - y3) + (z5 - z3) * (z5 - z3));
	dkp4 = sqrt((x5 - x4) * (x5 - x4) + (y5 - y4) * (y5 - y4) + (z5 - z4) * (z5 - z4));

	kra[0] = dkp4 + dkp2 + dkp1;
	kra[1] = dkp4 + dkp3 + dkp2;
	kra[2] = dkp4 + dkp1 + dkp3;
	kra[3] = dkp1 + dkp2 + dkp3;

}

void Delanouy3D::powierzchnieScian(double pol[4], const int & p1, const int & p2, const int & p3, const int & p4) {

	//1-2 a 2-3 b 3-1 c 1-4 d 2-4 e 3-4 f
	//double pol421,pol432,pol413,pol123;

	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a, b, c, d, e, f, p;

	x1 = points.getElement(p1).getX();
	x2 = points.getElement(p2).getX();
	x3 = points.getElement(p3).getX();
	x4 = points.getElement(p4).getX();
	y1 = points.getElement(p1).getY();
	y2 = points.getElement(p2).getY();
	y3 = points.getElement(p3).getY();
	y4 = points.getElement(p4).getY();
	z1 = points.getElement(p1).getZ();
	z2 = points.getElement(p2).getZ();
	z3 = points.getElement(p3).getZ();
	z4 = points.getElement(p4).getZ();

	a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
	b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
	c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
	d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
	e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
	f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

	p = 0.5 * (a + d + e);
	pol[0] = p * (p - a) * (p - d) * (p - e);
	p = 0.5 * (b + f + e);
	pol[1] = p * (p - b) * (p - f) * (p - e);
	p = 0.5 * (c + d + f);
	pol[2] = p * (p - c) * (p - d) * (p - f);
	p = 0.5 * (a + b + c);
	pol[3] = p * (p - a) * (p - b) * (p - c);

}


void Delanouy3D::poprawaDopasowaniaZiaren(const int & dokladDop, const int & rodzaj) {


	double pol[5];
	pol[4] = 0;
	int kolej[5];//kolejnosc kolej[4] = max kolej[0] = 0;
	kolej[0] = 4;

	int licz = 0;
	for (int i = 0, p1, p2, p3, p4, eS, ileE = elements.getIter(); i < ileE; ++i) {

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();

		if ((points.getElement(p1).getGranica() != 'a' && points.getElement(p2).getGranica() != 'a' && points.getElement(p3).getGranica() != 'a' && points.getElement(p4).getGranica() != 'a') || elements.getElement(i).getRodzajZiarna() == -1) {
			if (elements.getElement(i).getRodzajZiarna() == -1) {
				++licz;
			}
			//++licz;

			switch (rodzaj) {

				case 0:
					powierzchnieScian(pol, p1, p2, p3, p4);
					break;
				case 1:
					objetoscSrodAryt(pol, p1, p2, p3, p4);
					break;
				case 2:
					krawedzieSrodAryt(pol, p1, p2, p3, p4);
					break;

			}


			//powierzchnieScian(pol,p1,p2,p3,p4);
			//objetoscSrodAryt(pol,p1,p2,p3,p4);
			//krawedzieSrodAryt(pol,p1,p2,p3,p4);

			//  pol[kolej[3]]> pol[kolej[1]] > pol[kolej[2]] > pol[kolej[0]]
			segreguj4(pol, kolej);


			int zliczE = 4;
			int dokladnoscDop = dokladDop;
			while (zliczE && dokladnoscDop) {

				--dokladnoscDop;


				switch (kolej[zliczE]) {

					case 0:
						eS = elements.getElement(i).getE1();
						if (eS != -1 && pol[kolej[zliczE]] >= pol[kolej[zliczE - 1]]) {
							if (elements.getElement(eS).getRodzajZiarna() != -1) {
								elements.getElement(i).setRodzajZiarna(elements.getElement(eS).getRodzajZiarna());
								zliczE = 1;
							}
						}

						break;

					case 1:
						eS = elements.getElement(i).getE2();
						if (eS != -1 && pol[kolej[zliczE]] >= pol[kolej[zliczE - 1]]) {
							if (elements.getElement(eS).getRodzajZiarna() != -1) {
								elements.getElement(i).setRodzajZiarna(elements.getElement(eS).getRodzajZiarna());
								zliczE = 1;
							}
						}

						break;

					case 2:
						eS = elements.getElement(i).getE3();
						if (eS != -1 && pol[kolej[zliczE]] >= pol[kolej[zliczE - 1]]) {
							if (elements.getElement(eS).getRodzajZiarna() != -1) {
								elements.getElement(i).setRodzajZiarna(elements.getElement(eS).getRodzajZiarna());
								zliczE = 1;
							}
						}

						break;

					case 3:
						eS = elements.getElement(i).getE4();
						if (eS != -1 && pol[kolej[zliczE]] >= pol[kolej[zliczE - 1]]) {
							if (elements.getElement(eS).getRodzajZiarna() != -1) {
								elements.getElement(i).setRodzajZiarna(elements.getElement(eS).getRodzajZiarna());
								zliczE = 1;
							}
						}

						break;
				}

				--zliczE;
				//while
			}

		}



		//for
	}


//m1->Lines->Add(licz);
}

void Delanouy3D::dodajPunktWSrodkuNajwiekszejPow(IntList & zleE) {

//1-2 a 2-3 b 3-1 c 1-4 d 2-4 e 3-4 f
	double pol421, pol432, pol413, pol123;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a, b, c, d, e, f, p;

	for (int i = 0, el, p1, p2, p3, p4, ileE = zleE.getIter(); i < ileE; ++i) {

		el = zleE.getElement(i);
		p1 = elements.getElement(el).getP1();
		p2 = elements.getElement(el).getP2();
		p3 = elements.getElement(el).getP3();
		p4 = elements.getElement(el).getP4();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();

		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		p = 0.5 * (a + d + e);
		pol421 = p * (p - a) * (p - d) * (p - e);
		p = 0.5 * (b + f + e);
		pol432 = p * (p - b) * (p - f) * (p - e);
		p = 0.5 * (c + d + f);
		pol413 = p * (p - c) * (p - d) * (p - f);
		p = 0.5 * (a + b + c);
		pol123 = p * (p - a) * (p - b) * (p - c);


		switch (maximum4(pol421, pol432, pol413, pol123)) {

			case 1:
				points.setElement( (x4 + x2 + x1) / 3, (y4 + y2 + y1) / 3, (z4 + z2 + z1) / 3, 'a' );
				break;

			case 2:
				points.setElement( (x4 + x3 + x2) / 3, (y4 + y3 + y2) / 3, (z4 + z3 + z2) / 3, 'a' );
				break;

			case 3:
				points.setElement( (x4 + x1 + x3) / 3, (y4 + y1 + y3) / 3, (z4 + z1 + z3) / 3, 'a' );
				break;

			case 4:
				points.setElement( (x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3, (z1 + z2 + z3) / 3, 'a' );
				break;
		}

	}

}

int Delanouy3D::maximum4(const double & m1, const double & m2, const double & m3, const double & m4) {

	if (m1 > m2) {
		if (m1 > m3) {
			if (m1 > m4) {
				return 1;
			}
		}
		else {
			if (m3 > m4) {
				return 3;
			}
		}
	}
	else {
		if (m2 > m3) {
			if (m2 > m4) {
				return 2;
			}
		}
		else {
			if (m3 > m4) {
				return 3;
			}
		}
	}

	return 4;
}

int Delanouy3D::segreguj4(double * m, int * kol) {

	kol[1] = 0;
	kol[2] = 1;
	kol[3] = 2;
	kol[4] = 3;
	int k;


	if (m[0] > m[1]) {
		kol[1] = 1;
		kol[2] = 0;
	}
	if (m[2] > m[3]) {
		kol[3] = 3;
		kol[4] = 2;
	}
	if (m[kol[1]] > m[kol[3]]) {
		k = kol[3];
		kol[3] = kol[1];
		kol[1] = k;
	}

	if (m[kol[2]] < m[kol[1]]) {
		k = kol[1];
		kol[1] = kol[3];
		kol[3] = k;
		return 0;
	}
	if (m[kol[1]] > m[kol[3]]) {
		k = kol[3];
		kol[3] = kol[2];
		kol[2] = kol[1];
		kol[1] = k;
		return 0;
	}
	if (m[kol[2]] > m[kol[3]]) {
		k = kol[3];
		kol[3] = kol[2];
		kol[2] = k;
	}

	return 1;
}

void Delanouy3D::sprawdzSasiadow() {

	int zlychSasiadow = 0;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		double p1 = elements.getElement(i).getP1();
		double p2 = elements.getElement(i).getP2();
		double p3 = elements.getElement(i).getP3();
		double p4 = elements.getElement(i).getP4();

		int e1 = elements.getElement(i).getE1();
		int e2 = elements.getElement(i).getE2();
		int e3 = elements.getElement(i).getE3();
		int e4 = elements.getElement(i).getE4();
		int zlicz = 0;

		if (e1 != -1) {
			if (p1 == elements.getElement(e1).getP1() || p1 == elements.getElement(e1).getP2() || p1 == elements.getElement(e1).getP3() || p1 == elements.getElement(e1).getP4()) {
				++zlicz;
			}
			if (p2 == elements.getElement(e1).getP1() || p2 == elements.getElement(e1).getP2() || p2 == elements.getElement(e1).getP3() || p2 == elements.getElement(e1).getP4()) {
				++zlicz;
			}
			if (p3 == elements.getElement(e1).getP1() || p3 == elements.getElement(e1).getP2() || p3 == elements.getElement(e1).getP3() || p3 == elements.getElement(e1).getP4()) {
				++zlicz;
			}
			if (p4 == elements.getElement(e1).getP1() || p4 == elements.getElement(e1).getP2() || p4 == elements.getElement(e1).getP3() || p4 == elements.getElement(e1).getP4()) {
				++zlicz;
			}
			if (zlicz < 3) {
				++zlychSasiadow;
			}
			zlicz = 0;
		}
		if (e2 != -1) {
			if (p1 == elements.getElement(e2).getP1() || p1 == elements.getElement(e2).getP2() || p1 == elements.getElement(e2).getP3() || p1 == elements.getElement(e2).getP4()) {
				++zlicz;
			}
			if (p2 == elements.getElement(e2).getP1() || p2 == elements.getElement(e2).getP2() || p2 == elements.getElement(e2).getP3() || p2 == elements.getElement(e2).getP4()) {
				++zlicz;
			}
			if (p3 == elements.getElement(e2).getP1() || p3 == elements.getElement(e2).getP2() || p3 == elements.getElement(e2).getP3() || p3 == elements.getElement(e2).getP4()) {
				++zlicz;
			}
			if (p4 == elements.getElement(e2).getP1() || p4 == elements.getElement(e2).getP2() || p4 == elements.getElement(e2).getP3() || p4 == elements.getElement(e2).getP4()) {
				++zlicz;
			}
			if (zlicz < 3) {
				++zlychSasiadow;
			}
			zlicz = 0;
		}
		if (e3 != -1) {
			if (p1 == elements.getElement(e3).getP1() || p1 == elements.getElement(e3).getP2() || p1 == elements.getElement(e3).getP3() || p1 == elements.getElement(e3).getP4()) {
				++zlicz;
			}
			if (p2 == elements.getElement(e3).getP1() || p2 == elements.getElement(e3).getP2() || p2 == elements.getElement(e3).getP3() || p2 == elements.getElement(e3).getP4()) {
				++zlicz;
			}
			if (p3 == elements.getElement(e3).getP1() || p3 == elements.getElement(e3).getP2() || p3 == elements.getElement(e3).getP3() || p3 == elements.getElement(e3).getP4()) {
				++zlicz;
			}
			if (p4 == elements.getElement(e3).getP1() || p4 == elements.getElement(e3).getP2() || p4 == elements.getElement(e3).getP3() || p4 == elements.getElement(e3).getP4()) {
				++zlicz;
			}
			if (zlicz < 3) {
				++zlychSasiadow;
			}
			zlicz = 0;
		}
		if (e4 != -1) {
			if (p1 == elements.getElement(e4).getP1() || p1 == elements.getElement(e4).getP2() || p1 == elements.getElement(e4).getP3() || p1 == elements.getElement(e4).getP4()) {
				++zlicz;
			}
			if (p2 == elements.getElement(e4).getP1() || p2 == elements.getElement(e4).getP2() || p2 == elements.getElement(e4).getP3() || p2 == elements.getElement(e4).getP4()) {
				++zlicz;
			}
			if (p3 == elements.getElement(e4).getP1() || p3 == elements.getElement(e4).getP2() || p3 == elements.getElement(e4).getP3() || p3 == elements.getElement(e4).getP4()) {
				++zlicz;
			}
			if (p4 == elements.getElement(e4).getP1() || p4 == elements.getElement(e4).getP2() || p4 == elements.getElement(e4).getP3() || p4 == elements.getElement(e4).getP4()) {
				++zlicz;
			}
			if (zlicz < 3) {
				++zlychSasiadow;
			}
			zlicz = 0;
		}

//for
	}
	if (zlychSasiadow) {
		//m1->Lines->Add("Z£YCH SASIADOW");
		//m1->Lines->Add(zlychSasiadow);
	}
}



void Delanouy3D::ZapiszDoPlikuAbaqus(const char * nazwa) {

	IntList * MX;
	IntList eDoZ;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {
		eDoZ.setElementUniqat(elements.getElement(i).getRodzajZiarna());
	}

	MX = new IntList[eDoZ.getIter()];

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		for (int j = 0; j < eDoZ.getIter(); j++) {
			if (elements.getElement(i).getRodzajZiarna() == eDoZ.getElement(j)) {
				MX[j].setElement(i + 1);
			}
		}

	}


	textDoPliku text;
	const char * pause = "      ";
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;


	IntList nset1, nset2, nset3, nset4, nset5, nset6;
	bool n1, n2, n3, n4;
	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {
		n1 = true;
		n2 = true;
		n3 = true;
		n4 = true;

		if (points.getElement(i).getX() == 0) {
			nset1.setElement(i + 1);
			n1 = false;

		}
		if (points.getElement(i).getX() == maxX) {
			nset2.setElement(i + 1);
			n2 = false;

		}
		if (points.getElement(i).getY() == 0) {
			if (n1 && n2) {
				nset3.setElement(i + 1);
				n3 = false;
			}

		}
		if (points.getElement(i).getY() == maxY) {
			if (n1 && n2) {
				nset4.setElement(i + 1);
				n4 = false;
			}

		}
		if (points.getElement(i).getZ() == 0) {
			if (n1 && n2 && n3 && n4) {
				nset5.setElement(i + 1);
			}

		}
		if (points.getElement(i).getZ() == maxZ) {
			if (n1 && n2 && n3 && n4) {
				nset6.setElement(i + 1);
			}

		}

	}


	ofstream zapis(nazwa);
	zapis.precision(12);


	zapis << text.poczatekPliku << endl;

	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << pause << i + 1 << "," << pause << points.getElement(i).getX() << "," << pause << points.getElement(i).getY() << "," << pause << points.getElement(i).getZ() << endl;

	}

	zapis << "*Element, type=C3D4" << endl;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		zapis << i + 1 << "," << pause << elements.getElement(i).getP1() + 1 << "," << pause << elements.getElement(i).getP2() + 1 << "," << pause << elements.getElement(i).getP3() + 1 << "," << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	for (int i = 0; i < eDoZ.getIter(); i++) {
		zapis << "*Elset, elset=M" << i + 1 << endl;

		for (int j = 0, l = 1, ileP = MX[i].getIter(); j < ileP; j++, l++) {

			zapis << MX[i].getElement(j);
			if (j < ileP - 1) {
				zapis << ", ";
			}
			if (l > 9) {
				zapis << endl;
				l = 0;
			}

		}

		zapis << endl;

	}

	zapis << "*Nset, nset=N1" << endl;

	for (int i = 0, j = 1, ileP = nset1.getIter(); i < ileP; i++, j++) {

		zapis << nset1.getElement(i);
		if (i < ileP - 1) {
			zapis << ", ";
		}
		if (j > 9) {
			zapis << endl;
			j = 0;
		}

	}

	zapis << endl;
	zapis << "*Nset, nset=N2" << endl;

	for (int i = 0, j = 1, ileP = nset2.getIter(); i < ileP; i++, j++) {

		zapis << nset2.getElement(i);
		if (i < ileP - 1) {
			zapis << ", ";
		}
		if (j > 9) {
			zapis << endl;
			j = 0;
		}

	}

	zapis << endl;
	zapis << "*Nset, nset=N3" << endl;

	for (int i = 0, j = 1, ileP = nset3.getIter(); i < ileP; i++, j++) {

		zapis << nset3.getElement(i);
		if (i < ileP - 1) {
			zapis << ", ";
		}
		if (j > 9) {
			zapis << endl;
			j = 0;
		}

	}

	zapis << endl;
	zapis << "*Nset, nset=N4" << endl;

	for (int i = 0, j = 1, ileP = nset4.getIter(); i < ileP; i++, j++) {

		zapis << nset4.getElement(i);
		if (i < ileP - 1) {
			zapis << ", ";
		}
		if (j > 9) {
			zapis << endl;
			j = 0;
		}

	}

	zapis << endl;
	zapis << "*Nset, nset=N5" << endl;

	for (int i = 0, j = 1, ileP = nset5.getIter(); i < ileP; i++, j++) {

		zapis << nset5.getElement(i);
		if (i < ileP - 1) {
			zapis << ", ";
		}
		if (j > 9) {
			zapis << endl;
			j = 0;
		}

	}

	zapis << endl;
	zapis << "*Nset, nset=N6" << endl;

	for (int i = 0, j = 1, ileP = nset6.getIter(); i < ileP; i++, j++) {

		zapis << nset6.getElement(i);
		if (i < ileP - 1) {
			zapis << ", ";
		}
		if (j > 9) {
			zapis << endl;
			j = 0;
		}

	}


	zapis << endl;

	for (int i = 0; i < eDoZ.getIter(); i++) {
		zapis << "** Section: Section-" << i + 1 << "-M" << i + 1 << endl;
		zapis << "*Solid Section, elset=M" << i + 1 << ", material=A" << endl;
		zapis << "1.," << endl;
	}

	zapis << text.koncowkaPliku << endl;

	zapis.close();

	delete[] MX;

}

void Delanouy3D::ZapiszDoPlikuNas(const char * nazwa, double dzielnikZ, bool wartoscieDo1, bool scianyProstopadl) {


	double maxX =  dx;
	double maxY =  dy;
	double maxZ =  dz;

	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(elements.getIter() * 0.1, elements.getIter() * 0.1);
	const char * pause = "      ";
	double dzielnik = 1;

	if (wartoscieDo1) {
		dzielnik = 1 / maxXX;
	}
	dzielnikZ = 1 / dzielnikZ;

	/*
	for(int i=0,ileP=points.getIter();i<ileP;++i){

		if(points.getElement(i).getZ()>maxZ){maxZ=points.getElement(i).getZ();}

	}
	*/

	for (int i = 0, e1, e2, e3, e4, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; i++) {

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();


		if (e1 != -2) {
			if (e1 != -1) {
				face.setElement(p4, p2, p1);

				if (elements.getElement(e1).getE1() == i) {
					elements.getElement(e1).setE1(-2);
				}
				else if (elements.getElement(e1).getE2() == i) {
					elements.getElement(e1).setE2(-2);
				}
				else if (elements.getElement(e1).getE3() == i) {
					elements.getElement(e1).setE3(-2);
				}
				else {
					elements.getElement(e1).setE4(-2);
				}

			}
			else {
				face.setElement(p4, p2, p1, 'b');

			}
		}

		if (e2 != -2) {
			if (e2 != -1) {

				face.setElement(p4, p3, p2);
				if (elements.getElement(e2).getE1() == i) {
					elements.getElement(e2).setE1(-2);
				}
				else if (elements.getElement(e2).getE2() == i) {
					elements.getElement(e2).setE2(-2);
				}
				else if (elements.getElement(e2).getE3() == i) {
					elements.getElement(e2).setE3(-2);
				}
				else {
					elements.getElement(e2).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p3, p2, 'b');

			}
		}

		if (e3 != -2) {
			if (e3 != -1) {

				face.setElement(p4, p1, p3);
				if (elements.getElement(e3).getE1() == i) {
					elements.getElement(e3).setE1(-2);
				}
				else if (elements.getElement(e3).getE2() == i) {
					elements.getElement(e3).setE2(-2);
				}
				else if (elements.getElement(e3).getE3() == i) {
					elements.getElement(e3).setE3(-2);
				}
				else {
					elements.getElement(e3).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p1, p3, 'b');

			}
		}

		if (e4 != -2) {
			if (e4 != -1) {

				face.setElement(p1, p2, p3);
				if (elements.getElement(e4).getE1() == i) {
					elements.getElement(e4).setE1(-2);
				}
				else if (elements.getElement(e4).getE2() == i) {
					elements.getElement(e4).setE2(-2);
				}
				else if (elements.getElement(e4).getE3() == i) {
					elements.getElement(e4).setE3(-2);
				}
				else {
					elements.getElement(e4).setE4(-2);
				}
			}
			else {
				face.setElement(p1, p2, p3, 'b');

			}
		}

	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM

	if (scianyProstopadl) {
		for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
			if (face.getElement(i).getGranica() == 'b') {

				p1 = face.getElement(i).getX();
				p2 = face.getElement(i).getY();
				p3 = face.getElement(i).getZ();

				x1 = points.getElement(p1).getX();
				x2 = points.getElement(p2).getX();
				x3 = points.getElement(p3).getX();

				y1 = points.getElement(p1).getY();
				y2 = points.getElement(p2).getY();
				y3 = points.getElement(p3).getY();

				z1 = points.getElement(p1).getZ();
				z2 = points.getElement(p2).getZ();
				z3 = points.getElement(p3).getZ();


				if (x1 == x2 && x2 == x3) {
					if (x1 == 0) {
						face.getElement(i).setGranica('c');
					}
					else if (x1 == maxX) {
						face.getElement(i).setGranica('d');
					}
				}
				else if (y1 == y2 && y2 == y3) {
					if (y1 == 0) {
						face.getElement(i).setGranica('e');
					}
					else if (y1 == maxY) {
						face.getElement(i).setGranica('f');
					}
				}
				else if (z1 == z2 && z2 == z3) {
					if (z1 == 0) {
						face.getElement(i).setGranica('g');
					}
					else if (z1 == maxZ) {
						face.getElement(i).setGranica('h');
					}
				}

				/*
				if(y1 < x1*3.5087 + 40.35087 && y2 < x2*3.5087 + 40.35087 && y3 < x3*3.5087 + 40.35087){

					if(z1>45){
						face.getElement(i).setGranica('g');
					}
					else{
						face.getElement(i).setGranica('h');
					}

				}

				*/

			}
		}
	}

	ofstream zapis(nazwa);
	zapis.precision(12);


	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << points.getElement(i).getX()*dzielnik << " " << pause << points.getElement(i).getY()*dzielnik << " " << pause << points.getElement(i).getZ()*dzielnik * dzielnikZ << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	zapis << endl << endl;


	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = 3;
				break;
			case 'd':
				rodzaj = 4;
				break;
			case 'e':
				rodzaj = 5;
				break;
			case 'f':
				rodzaj = 6;
				break;
			case 'g':
				rodzaj = 1;
				break;
			case 'h':
				rodzaj = 2;
				break;
			case 'b':
				rodzaj = 9;
				break;
		}


		zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;

	}


	zapis.close();
}


void Delanouy3D::ZapiszDoPlikuNasPryzm(const char * nazwa, double dzielnikZ, bool wartoscieDo1, bool scianyProstopadl) {

//double dlaZetowej = 50;
//double dlaZetowej = 2.5;

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(elements.getIter() * 0.1, elements.getIter() * 0.1);
	IntList face2(50000, 50000);
	const char * pause = "      ";
	double dzielnik = 1;

	if (wartoscieDo1) {
		dzielnik = 1 / maxXX;
	}
	dzielnikZ = 1 / dzielnikZ;


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getZ() > maxZ) {
			maxZ = points.getElement(i).getZ();
		}

	}

	for (int i = 0, e1, e2, e3, e4, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; i++) {

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();


		if (e1 != -2) {
			if (e1 != -1) {
				face.setElement(p4, p2, p1);

				if (elements.getElement(e1).getE1() == i) {
					elements.getElement(e1).setE1(-2);
				}
				else if (elements.getElement(e1).getE2() == i) {
					elements.getElement(e1).setE2(-2);
				}
				else if (elements.getElement(e1).getE3() == i) {
					elements.getElement(e1).setE3(-2);
				}
				else {
					elements.getElement(e1).setE4(-2);
				}

			}
			else {
				face.setElement(p4, p2, p1, 'b');

			}
		}

		if (e2 != -2) {
			if (e2 != -1) {

				face.setElement(p4, p3, p2);
				if (elements.getElement(e2).getE1() == i) {
					elements.getElement(e2).setE1(-2);
				}
				else if (elements.getElement(e2).getE2() == i) {
					elements.getElement(e2).setE2(-2);
				}
				else if (elements.getElement(e2).getE3() == i) {
					elements.getElement(e2).setE3(-2);
				}
				else {
					elements.getElement(e2).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p3, p2, 'b');

			}
		}

		if (e3 != -2) {
			if (e3 != -1) {

				face.setElement(p4, p1, p3);
				if (elements.getElement(e3).getE1() == i) {
					elements.getElement(e3).setE1(-2);
				}
				else if (elements.getElement(e3).getE2() == i) {
					elements.getElement(e3).setE2(-2);
				}
				else if (elements.getElement(e3).getE3() == i) {
					elements.getElement(e3).setE3(-2);
				}
				else {
					elements.getElement(e3).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p1, p3, 'b');

			}
		}

		if (e4 != -2) {
			if (e4 != -1) {

				face.setElement(p1, p2, p3);
				if (elements.getElement(e4).getE1() == i) {
					elements.getElement(e4).setE1(-2);
				}
				else if (elements.getElement(e4).getE2() == i) {
					elements.getElement(e4).setE2(-2);
				}
				else if (elements.getElement(e4).getE3() == i) {
					elements.getElement(e4).setE3(-2);
				}
				else {
					elements.getElement(e4).setE4(-2);
				}
			}
			else {
				face.setElement(p1, p2, p3, 'b');

			}
		}

	}

	double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6, z1, z2, z3, z4, z5, z6;

	for (int i = 0, p1, p2, p3, p4, p5, p6, ileE = pryzmy.getIter(); i < ileE; i++) {

		p1 = pryzmy.getElement(i).getP1();
		p2 = pryzmy.getElement(i).getP2();
		p3 = pryzmy.getElement(i).getP3();
		p4 = pryzmy.getElement(i).getP4();
		p5 = pryzmy.getElement(i).getP5();
		p6 = pryzmy.getElement(i).getP6();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		x5 = points.getElement(p5).getX();
		x6 = points.getElement(p6).getX();

		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		y5 = points.getElement(p5).getY();
		y6 = points.getElement(p6).getY();

		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();
		z5 = points.getElement(p5).getZ();
		z6 = points.getElement(p6).getZ();

		if ( y1 == 0 || y1 == maxY ) {
			//m1->Lines->Add("dol");
			if (y1 == 0) {
				face.setElement(p1, p3, p2, 'e');
			}
			else {
				face.setElement(p1, p3, p2, 'f');
			}
		}
		else if ( y4 == 0 || y4 == maxY ) {
			if (y4 == 0) {
				face.setElement(p6, p4, p5, 'e');
			}
			else {
				face.setElement(p6, p4, p5, 'f');
			}
		}

		if ( (x1 == 0 || x1 == maxX) && x1 == x2) {
			if (x1 == 0) {
				face2.setElement(p5);
				face2.setElement(p4);
				face2.setElement(p1);
				face2.setElement(p2);
				face2.setElement(3);
			}
			else {
				face2.setElement(p5);
				face2.setElement(p4);
				face2.setElement(p1);
				face2.setElement(p2);
				face2.setElement(4);
			}
		}
		else if ( (x1 == 0 || x3 == maxX) && x1 == x3 ) {
			if (x1 == 0) {
				face2.setElement(p6);
				face2.setElement(p3);
				face2.setElement(p1);
				face2.setElement(p4);
				face2.setElement(3);
			}
			else {
				face2.setElement(p6);
				face2.setElement(p3);
				face2.setElement(p1);
				face2.setElement(p4);
				face2.setElement(4);
			}
		}
		else if ( (x2 == 0 || x3 == maxX) && x2 == x3 ) {
			if (x2 == 0) {
				face2.setElement(p6);
				face2.setElement(p5);
				face2.setElement(p2);
				face2.setElement(p3);
				face2.setElement(3);
			}
			else {
				face2.setElement(p6);
				face2.setElement(p5);
				face2.setElement(p2);
				face2.setElement(p3);
				face2.setElement(4);
			}

		}

		if ( (z1 == 0 || z1 == maxZ) && z1 == z2) {
			if (z1 == 0) {
				face2.setElement(p5);
				face2.setElement(p4);
				face2.setElement(p1);
				face2.setElement(p2);
				face2.setElement(7);
			}
			else {
				face2.setElement(p5);
				face2.setElement(p4);
				face2.setElement(p1);
				face2.setElement(p2);
				face2.setElement(8);
			}
		}
		else if ( (z1 == 0 || z3 == maxZ) && z1 == z3 ) {
			if (z1 == 0) {
				face2.setElement(p6);
				face2.setElement(p3);
				face2.setElement(p1);
				face2.setElement(p4);
				face2.setElement(7);
			}
			else {
				face2.setElement(p6);
				face2.setElement(p3);
				face2.setElement(p1);
				face2.setElement(p4);
				face2.setElement(8);
			}
		}
		else if ( (z2 == 0 || z3 == maxZ) && z2 == z3 ) {
			if (z2 == 0) {
				face2.setElement(p6);
				face2.setElement(p5);
				face2.setElement(p2);
				face2.setElement(p3);
				face2.setElement(7);
			}
			else {
				face2.setElement(p6);
				face2.setElement(p5);
				face2.setElement(p2);
				face2.setElement(p3);
				face2.setElement(8);
			}

		}


	}


// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM

	if (scianyProstopadl) {
		for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
			if (face.getElement(i).getGranica() == 'b') {

				p1 = face.getElement(i).getX();
				p2 = face.getElement(i).getY();
				p3 = face.getElement(i).getZ();

				x1 = points.getElement(p1).getX();
				x2 = points.getElement(p2).getX();
				x3 = points.getElement(p3).getX();

				y1 = points.getElement(p1).getY();
				y2 = points.getElement(p2).getY();
				y3 = points.getElement(p3).getY();

				z1 = points.getElement(p1).getZ();
				z2 = points.getElement(p2).getZ();
				z3 = points.getElement(p3).getZ();


				if (x1 == x2 && x2 == x3) {
					if (x1 == 0) {
						face.getElement(i).setGranica('c');
					}
					else if (x1 == maxX) {
						face.getElement(i).setGranica('d');
					}
				}
				else if (y1 == y2 && y2 == y3) {
					if (y1 == 0) {
						face.getElement(i).setGranica('e');
					}
					else if (y1 == maxY) {
						face.getElement(i).setGranica('f');
					}
				}
				else if (z1 == z2 && z2 == z3) {
					if (z1 == 0) {
						face.getElement(i).setGranica('g');
					}
					else if (z1 == maxZ) {
						face.getElement(i).setGranica('h');
					}
				}

			}
		}
	}

	ofstream zapis(nazwa);
	zapis.precision(12);


	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << points.getElement(i).getX()*dzielnik << " " << pause << points.getElement(i).getY()*dzielnik << " " << pause << points.getElement(i).getZ()*dzielnik * dzielnikZ << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	zapis << endl << endl;

	for (int ileE = pryzmy.getIter(), j = 0; j < ileE; ++i, ++j) {


		zapis << "CPRIZM" << pause << i + 1 << " " << pause << 1 << "   " << pryzmy.getElement(j).getP1() + 1 << " " << pause << pryzmy.getElement(j).getP2() + 1 << " " << pause << pryzmy.getElement(j).getP3() + 1 << " " << pause << pryzmy.getElement(j).getP4() + 1 << " " << pause << pryzmy.getElement(j).getP5() + 1 << " " << pause << pryzmy.getElement(j).getP6() + 1 << endl;

	}

	zapis << endl << endl;

	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = 3;
				break;
			case 'd':
				rodzaj = 4;
				break;
			case 'e':
				rodzaj = 5;
				break;
			case 'f':
				rodzaj = 6;
				break;
			case 'g':
				rodzaj = 1;
				break;
			case 'h':
				rodzaj = 2;
				break;
			case 'b':
				rodzaj = 9;
				break;
		}


		zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;

	}

	zapis << endl << endl;

	for (int j = 0, ileF = face2.getIter(), rodzaj = -1; j < ileF; ++i, j += 5) {

		switch (face2.getElement(j + 4)) {

			case 3:
				rodzaj = 3;
				break;
			case 4:
				rodzaj = 4;
				break;
			case 7:
				rodzaj = 1;
				break;
			case 8:
				rodzaj = 2;
				break;
		}

		zapis << "CQUADX" << pause << i + 1 << " " << pause << rodzaj << "   " << face2.getElement(j) + 1 << " " << pause << face2.getElement(j + 1) + 1 << " " << pause << face2.getElement(j + 2) + 1 << " " << pause << face2.getElement(j + 3) + 1 << endl;

	}



	zapis.close();
}

void Delanouy3D::ZapiszDoPlikuNasWOparciuOPowTetra(const char * nazwa, double dzielnikZ, bool wartoscieDo1, bool scianyProstopadl) {

	int zliczIleFace = 0;
	int zliczIleNa1 = 0;
	int zliczIleNa2 = 0;

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(elements.getIter() * 0.1, elements.getIter() * 0.1);
	const char * pause = "      ";
	double dzielnik;

	if (wartoscieDo1) {
		dzielnik = 1 / maxXX;
	}
	dzielnikZ = 1 / dzielnikZ;


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getZ() > maxZ) {
			maxZ = points.getElement(i).getZ();
		}

	}

	for (int i = 0, e1, e2, e3, e4, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; i++) {

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();


		if (e1 != -2) {
			if (e1 != -1) {
				face.setElement(p4, p2, p1);

				if (elements.getElement(e1).getE1() == i) {
					elements.getElement(e1).setE1(-2);
				}
				else if (elements.getElement(e1).getE2() == i) {
					elements.getElement(e1).setE2(-2);
				}
				else if (elements.getElement(e1).getE3() == i) {
					elements.getElement(e1).setE3(-2);
				}
				else {
					elements.getElement(e1).setE4(-2);
				}

			}
			else {
				face.setElement(p4, p2, p1, 'b');
				zliczIleFace++;

			}
		}

		if (e2 != -2) {
			if (e2 != -1) {

				face.setElement(p4, p3, p2);
				if (elements.getElement(e2).getE1() == i) {
					elements.getElement(e2).setE1(-2);
				}
				else if (elements.getElement(e2).getE2() == i) {
					elements.getElement(e2).setE2(-2);
				}
				else if (elements.getElement(e2).getE3() == i) {
					elements.getElement(e2).setE3(-2);
				}
				else {
					elements.getElement(e2).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p3, p2, 'b');
				zliczIleFace++;

			}
		}

		if (e3 != -2) {
			if (e3 != -1) {

				face.setElement(p4, p1, p3);
				if (elements.getElement(e3).getE1() == i) {
					elements.getElement(e3).setE1(-2);
				}
				else if (elements.getElement(e3).getE2() == i) {
					elements.getElement(e3).setE2(-2);
				}
				else if (elements.getElement(e3).getE3() == i) {
					elements.getElement(e3).setE3(-2);
				}
				else {
					elements.getElement(e3).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p1, p3, 'b');
				zliczIleFace++;

			}
		}

		if (e4 != -2) {
			if (e4 != -1) {

				face.setElement(p1, p2, p3);
				if (elements.getElement(e4).getE1() == i) {
					elements.getElement(e4).setE1(-2);
				}
				else if (elements.getElement(e4).getE2() == i) {
					elements.getElement(e4).setE2(-2);
				}
				else if (elements.getElement(e4).getE3() == i) {
					elements.getElement(e4).setE3(-2);
				}
				else {
					elements.getElement(e4).setE4(-2);
				}
			}
			else {
				face.setElement(p1, p2, p3, 'b');
				zliczIleFace++;

			}
		}

	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM

	if (scianyProstopadl) {
		for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
			if (face.getElement(i).getGranica() == 'b') {

				p1 = face.getElement(i).getX();
				p2 = face.getElement(i).getY();
				p3 = face.getElement(i).getZ();

				x1 = points.getElement(p1).getX();
				x2 = points.getElement(p2).getX();
				x3 = points.getElement(p3).getX();

				y1 = points.getElement(p1).getY();
				y2 = points.getElement(p2).getY();
				y3 = points.getElement(p3).getY();

				z1 = points.getElement(p1).getZ();
				z2 = points.getElement(p2).getZ();
				z3 = points.getElement(p3).getZ();

				/*
				if(x1 == x2 && x2 == x3){
					if(x1==0){face.getElement(i).setGranica('c');}
					else if(x1==max){face.getElement(i).setGranica('d');}
				}
				else if(y1 == y2 && y2 == y3){
					if(y1==0){face.getElement(i).setGranica('e');}
					else if(y1==max){face.getElement(i).setGranica('f');}
				}
				else if(z1 == z2 && z2 == z3){
					if(z1==0){face.getElement(i).setGranica('g');}
					else if(z1==maxZ){face.getElement(i).setGranica('h');}
				}
				*/
				/*
				//bez zastawki
				if(y1 >= x1*3.5087 + 40.35087 && y2 >= x2*3.5087 + 40.35087 && y3 >= x3*3.5087 + 40.35087){

					if(z1>45){
						face.getElement(i).setGranica('g');zliczIleNa1++;
					}
					else{
						face.getElement(i).setGranica('h');zliczIleNa2++;
					}

				}
				*/

				//z zastawka
				if (z1 >= x1 * (-4.25531) + 461.702 && z2 >= x2 * (-4.25531) + 461.702 && z3 >= x3 * (-4.25531) + 461.702) {

					if (y1 > 45) {
						face.getElement(i).setGranica('g');
						zliczIleNa1++;
					}
					else {
						face.getElement(i).setGranica('h');
						zliczIleNa2++;
					}

				}



			}
		}
	}

	ofstream zapis(nazwa);
	zapis.precision(12);
//globalSkal =1;
	dzielnik = 1 / globalSkal;

	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << (points.getElement(i).getX() - globalTranX)*dzielnik << " " << pause << (points.getElement(i).getY() - globalTranY)*dzielnik << " " << pause << (points.getElement(i).getZ() - globalTranZ)*dzielnik << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}


	zapis << endl << endl;


	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = 3;
				break;
			case 'd':
				rodzaj = 4;
				break;
			case 'e':
				rodzaj = 5;
				break;
			case 'f':
				rodzaj = 6;
				break;
			case 'g':
				rodzaj = 1;
				break;
			case 'h':
				rodzaj = 2;
				break;
			case 'b':
				rodzaj = 9;
				break;
		}


		if (rodzaj) {
			zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;
		}
		else {
			--i;
		}

		/*
		zapis<<"CTRIAX"<<pause<<i+1<<" "<<pause<<rodzaj<<"   "<<face.getElement(j).getX()+1<<" "<<pause<<face.getElement(j).getY()+1<<" "<<pause<<face.getElement(j).getZ()+1<<endl;
		*/
	}

	elements.ustawSasiadow(-2, -1);

	zapis.close();
	/*
	m1->Lines->Add("Ile F = 1,2,9");
	m1->Lines->Add(zliczIleFace);
	m1->Lines->Add(zliczIleNa1);
	m1->Lines->Add(zliczIleNa2);
	m1->Lines->Add(zliczIleFace-zliczIleNa1-zliczIleNa2);
	*/


}



void Delanouy3D::ZapiszDoPlikuNasWOparciuOPowHybrydPSS(const char * nazwa, double dzielnikZ, bool wartoscieDo1, bool scianyProstopadl) {

	int zliczIleFace = 0;
	int zliczIleNa1 = 0;
	int zliczIleNa2 = 0;


	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(pryzmy.getIter() + 100, pryzmy.getIter() + 100);
	const char * pause = "      ";
	double dzielnik;

	if (wartoscieDo1) {
		dzielnik = 1 / maxXX;
	}
	dzielnikZ = 1 / dzielnikZ;


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getZ() > maxZ) {
			maxZ = points.getElement(i).getZ();
		}

	}

	for (int i = 0, p1, p2, p3, ileE = pryzmy.getIter(); i < ileE; i++) {

		p1 = pryzmy.getElement(i).getP1();
		p2 = pryzmy.getElement(i).getP2();
		p3 = pryzmy.getElement(i).getP3();
		/*

		if( (points.getElement(p1).getGranica()=='g' || points.getElement(p1).getGranica()=='r') &&
			(points.getElement(p2).getGranica()=='g' || points.getElement(p2).getGranica()=='r') &&
			(points.getElement(p3).getGranica()=='g' || points.getElement(p3).getGranica()=='r') ){

				points.getElement(p1).setGranica('r');
				points.getElement(p2).setGranica('r');
				points.getElement(p3).setGranica('r');

		*/

		if ( points.getElement(p1).getGranica() == 'g' &&
				points.getElement(p2).getGranica() == 'g' &&
				points.getElement(p3).getGranica() == 'g' ) {

			face.setElement(p1, p2, p3, 'b');
			zliczIleFace++;
		}

	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM

	if (scianyProstopadl) {
		for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
			if (face.getElement(i).getGranica() == 'b') {

				p1 = face.getElement(i).getX();
				p2 = face.getElement(i).getY();
				p3 = face.getElement(i).getZ();

				x1 = points.getElement(p1).getX();
				x2 = points.getElement(p2).getX();
				x3 = points.getElement(p3).getX();

				y1 = points.getElement(p1).getY();
				y2 = points.getElement(p2).getY();
				y3 = points.getElement(p3).getY();

				z1 = points.getElement(p1).getZ();
				z2 = points.getElement(p2).getZ();
				z3 = points.getElement(p3).getZ();



				if (y1 >= x1 * 3.5087 + 40.35087 && y2 >= x2 * 3.5087 + 40.35087 && y3 >= x3 * 3.5087 + 40.35087) {

					if (z1 > 45) {
						face.getElement(i).setGranica('g');
						zliczIleNa1++;
					}
					else {
						face.getElement(i).setGranica('h');
						zliczIleNa2++;
					}

				}

			}
		}
	}

	ofstream zapis(nazwa);
	zapis.precision(12);
//globalSkal=1;
	dzielnik = 1 / globalSkal;

	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << (points.getElement(i).getX() - globalTranX)*dzielnik << " " << pause << (points.getElement(i).getY() - globalTranY)*dzielnik << " " << pause << (points.getElement(i).getZ() - globalTranZ)*dzielnik << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	zapis << endl << endl;

	for (int ileE = pryzmy.getIter(), j = 0; j < ileE; ++i, ++j) {


		zapis << "CPRIZM" << pause << i + 1 << " " << pause << 1 << "   " << pryzmy.getElement(j).getP1() + 1 << " " << pause << pryzmy.getElement(j).getP2() + 1 << " " << pause << pryzmy.getElement(j).getP3() + 1 << " " << pause << pryzmy.getElement(j).getP4() + 1 << " " << pause << pryzmy.getElement(j).getP5() + 1 << " " << pause << pryzmy.getElement(j).getP6() + 1 << endl;

	}
	zapis << endl << endl;


	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = 3;
				break;
			case 'd':
				rodzaj = 4;
				break;
			case 'e':
				rodzaj = 5;
				break;
			case 'f':
				rodzaj = 6;
				break;
			case 'g':
				rodzaj = 1;
				break;
			case 'h':
				rodzaj = 2;
				break;
			case 'b':
				rodzaj = 9;
				break;
		}


		if (rodzaj) {
			zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;
		}
		else {
			--i;
		}

	}



	zapis.close();
	/*
	m1->Lines->Add("Ile F = 1,2,9");
	m1->Lines->Add(zliczIleFace);
	m1->Lines->Add(zliczIleNa1);
	m1->Lines->Add(zliczIleNa2);
	m1->Lines->Add(zliczIleFace-zliczIleNa1-zliczIleNa2);
	*/
}


void Delanouy3D::ZapiszDoPlikuNasWOparciuOPowTetraKanal(const char * nazwa, int wartoscieDo1) {

	int zliczIleFace = 0;
	int zliczIleNa1 = 0;
	int zliczIleNa2 = 0;

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(elements.getIter() * 0.1, elements.getIter() * 0.1);
	const char * pause = "      ";
	double dzielnik = 1;

	if (wartoscieDo1) {
		dzielnik = wartoscieDo1 / maxXX;
	}



	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getZ() > maxZ) {
			maxZ = points.getElement(i).getZ();
		}

	}

	for (int i = 0, e1, e2, e3, e4, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; i++) {

		e1 = elements.getElement(i).getE1();
		e2 = elements.getElement(i).getE2();
		e3 = elements.getElement(i).getE3();
		e4 = elements.getElement(i).getE4();

		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();


		if (e1 != -2) {
			if (e1 != -1) {
				face.setElement(p4, p2, p1);

				if (elements.getElement(e1).getE1() == i) {
					elements.getElement(e1).setE1(-2);
				}
				else if (elements.getElement(e1).getE2() == i) {
					elements.getElement(e1).setE2(-2);
				}
				else if (elements.getElement(e1).getE3() == i) {
					elements.getElement(e1).setE3(-2);
				}
				else {
					elements.getElement(e1).setE4(-2);
				}

			}
			else {
				face.setElement(p4, p2, p1, 'b');
				zliczIleFace++;

			}
		}

		if (e2 != -2) {
			if (e2 != -1) {

				face.setElement(p4, p3, p2);
				if (elements.getElement(e2).getE1() == i) {
					elements.getElement(e2).setE1(-2);
				}
				else if (elements.getElement(e2).getE2() == i) {
					elements.getElement(e2).setE2(-2);
				}
				else if (elements.getElement(e2).getE3() == i) {
					elements.getElement(e2).setE3(-2);
				}
				else {
					elements.getElement(e2).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p3, p2, 'b');
				zliczIleFace++;

			}
		}

		if (e3 != -2) {
			if (e3 != -1) {

				face.setElement(p4, p1, p3);
				if (elements.getElement(e3).getE1() == i) {
					elements.getElement(e3).setE1(-2);
				}
				else if (elements.getElement(e3).getE2() == i) {
					elements.getElement(e3).setE2(-2);
				}
				else if (elements.getElement(e3).getE3() == i) {
					elements.getElement(e3).setE3(-2);
				}
				else {
					elements.getElement(e3).setE4(-2);
				}
			}
			else {
				face.setElement(p4, p1, p3, 'b');
				zliczIleFace++;

			}
		}

		if (e4 != -2) {
			if (e4 != -1) {

				face.setElement(p1, p2, p3);
				if (elements.getElement(e4).getE1() == i) {
					elements.getElement(e4).setE1(-2);
				}
				else if (elements.getElement(e4).getE2() == i) {
					elements.getElement(e4).setE2(-2);
				}
				else if (elements.getElement(e4).getE3() == i) {
					elements.getElement(e4).setE3(-2);
				}
				else {
					elements.getElement(e4).setE4(-2);
				}
			}
			else {
				face.setElement(p1, p2, p3, 'b');
				zliczIleFace++;

			}
		}

	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM


	for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
		if (face.getElement(i).getGranica() == 'b') {

			p1 = face.getElement(i).getX();
			p2 = face.getElement(i).getY();
			p3 = face.getElement(i).getZ();

			x1 = points.getElement(p1).getX();
			x2 = points.getElement(p2).getX();
			x3 = points.getElement(p3).getX();

			y1 = points.getElement(p1).getY();
			y2 = points.getElement(p2).getY();
			y3 = points.getElement(p3).getY();

			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();


			if (z1 == z2 && z2 == z3 ) {

				if (z1 == 0) {
					face.getElement(i).setGranica('g');
					zliczIleNa1++;
				}
				else if (z1 == maxZ) {
					face.getElement(i).setGranica('h');
					zliczIleNa2++;
				}

			}

		}
	}


	ofstream zapis(nazwa);
	zapis.precision(12);


	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << (points.getElement(i).getX())*dzielnik << " " << pause << (points.getElement(i).getY())*dzielnik << " " << pause << (points.getElement(i).getZ())*dzielnik << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}


	zapis << endl << endl;


	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = 3;
				break;
			case 'd':
				rodzaj = 4;
				break;
			case 'e':
				rodzaj = 5;
				break;
			case 'f':
				rodzaj = 6;
				break;
			case 'g':
				rodzaj = 1;
				break;
			case 'h':
				rodzaj = 2;
				break;
			case 'b':
				rodzaj = 9;
				break;
		}


		if (rodzaj) {
			zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;
		}
		else {
			--i;
		}

	}



	zapis.close();


	elements.ustawSasiadow(-2, -1);

}




void Delanouy3D::ZapiszDoPlikuNasWOparciuOPowHybrydKanal(const char * nazwa, int wartoscMaxDo) {

	int zliczIleFace = 0;
	int zliczIleNa1 = 0;
	int zliczIleNa2 = 0;

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(pryzmy.getIter() + 100, pryzmy.getIter() + 100);
	const char * pause = "      ";
	double dzielnik = 1;

	if (wartoscMaxDo) {
		dzielnik = wartoscMaxDo / maxXX;
	}



	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getZ() > maxZ) {
			maxZ = points.getElement(i).getZ();
		}

	}

	for (int i = 0, p1, p2, p3, ileE = pryzmy.getIter(); i < ileE; i++) {

		p1 = pryzmy.getElement(i).getP1();
		p2 = pryzmy.getElement(i).getP2();
		p3 = pryzmy.getElement(i).getP3();


		if ( points.getElement(p1).getGranica() == 'g' &&
				points.getElement(p2).getGranica() == 'g' &&
				points.getElement(p3).getGranica() == 'g' ) {

			face.setElement(p1, p2, p3, 'b');
			zliczIleFace++;
		}

	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM


	for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
		if (face.getElement(i).getGranica() == 'b') {

			p1 = face.getElement(i).getX();
			p2 = face.getElement(i).getY();
			p3 = face.getElement(i).getZ();

			x1 = points.getElement(p1).getX();
			x2 = points.getElement(p2).getX();
			x3 = points.getElement(p3).getX();

			y1 = points.getElement(p1).getY();
			y2 = points.getElement(p2).getY();
			y3 = points.getElement(p3).getY();

			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();

			//wlasna funkcja

			if (z1 == z2 && z2 == z3 ) {

				if (z1 == 0) {
					face.getElement(i).setGranica('g');
					zliczIleNa1++;
				}
				else if (z1 == maxZ) {
					face.getElement(i).setGranica('h');
					zliczIleNa2++;
				}

			}

		}
	}


	ofstream zapis(nazwa);
	zapis.precision(12);


	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << (points.getElement(i).getX())*dzielnik << " " << pause << (points.getElement(i).getY())*dzielnik << " " << pause << (points.getElement(i).getZ())*dzielnik << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	zapis << endl << endl;

	for (int ileE = pryzmy.getIter(), j = 0; j < ileE; ++i, ++j) {


		zapis << "CPRIZM" << pause << i + 1 << " " << pause << 1 << "   " << pryzmy.getElement(j).getP1() + 1 << " " << pause << pryzmy.getElement(j).getP2() + 1 << " " << pause << pryzmy.getElement(j).getP3() + 1 << " " << pause << pryzmy.getElement(j).getP4() + 1 << " " << pause << pryzmy.getElement(j).getP5() + 1 << " " << pause << pryzmy.getElement(j).getP6() + 1 << endl;

	}
	zapis << endl << endl;


	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = 3;
				break;
			case 'd':
				rodzaj = 4;
				break;
			case 'e':
				rodzaj = 5;
				break;
			case 'f':
				rodzaj = 6;
				break;
			case 'g':
				rodzaj = 1;
				break;
			case 'h':
				rodzaj = 2;
				break;
			case 'b':
				rodzaj = 9;
				break;
		}


		if (rodzaj) {
			zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;
		}
		else {
			--i;
		}

	}



	zapis.close();
	/*
	m1->Lines->Add("Ile F = 1,2,9");
	m1->Lines->Add(zliczIleFace);
	m1->Lines->Add(zliczIleNa1);
	m1->Lines->Add(zliczIleNa2);

	m1->Lines->Add(zliczIleFace-zliczIleNa1-zliczIleNa2);
	*/
}

void Delanouy3D::ZapiszDoPlikuNasWOparciuOPowHybrydDowolny(const char * nazwa, double dzielnikZ, bool wartoscieDo1, int warunki) {


	char cc = 'c', dd = 'd', ee = 'b', ff = 'b', gg = 'g', hh = 'h';
	/*
				   ;if(warunki>=100000){cc='c';}
	warunki-=100000;if(warunki>=10000){dd='d';}else{warunki+=100000;}
	warunki-=10000;if(warunki>=1000){ee='e';}else{warunki+=10000;}
	warunki-=1000;if(warunki>=100){ff='f';}else{warunki+=1000;}
	warunki-=100;if(warunki>=10){gg='g';}else{warunki+=100;}
	warunki-=10;if(warunki>=1){hh='h';}

	m1->Lines->Add(cc);
	m1->Lines->Add(dd);
	m1->Lines->Add(ee);
	m1->Lines->Add(ff);
	m1->Lines->Add(gg);
	m1->Lines->Add(hh);
	*/

	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;
	double maxXX;

	if (maxX > maxY) {
		if (maxX > maxZ) {
			maxXX = maxX;
		}
		else {
			maxXX = maxZ;
		}
	}
	else {
		if (maxY > maxZ) {
			maxXX = maxY;
		}
		else {
			maxXX = maxZ;
		}
	}

	PunktList face(pryzmy.getIter() + 100, pryzmy.getIter() + 100);
	const char * pause = "      ";
	double dzielnik;

	if (wartoscieDo1) {
		dzielnik = 1 / maxXX;
	}
	dzielnikZ = 1 / dzielnikZ;


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		if (points.getElement(i).getZ() > maxZ) {
			maxZ = points.getElement(i).getZ();
		}

	}

	for (int i = 0, p1, p2, p3, ileE = pryzmy.getIter(); i < ileE; i++) {

		p1 = pryzmy.getElement(i).getP1();
		p2 = pryzmy.getElement(i).getP2();
		p3 = pryzmy.getElement(i).getP3();


		if ( points.getElement(p1).getGranica() == 'g' &&
				points.getElement(p2).getGranica() == 'g' &&
				points.getElement(p3).getGranica() == 'g' ) {

			face.setElement(p1, p2, p3, 'b');
		}

	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
// b- graniczny(10) c- x0 d-xM e- y0 f-yM g-z0 h-zM


	for (int i = 0, p1, p2, p3, ileF = face.getIter(); i < ileF; ++i) {
		if (face.getElement(i).getGranica() == 'b') {

			p1 = face.getElement(i).getX();
			p2 = face.getElement(i).getY();
			p3 = face.getElement(i).getZ();

			x1 = points.getElement(p1).getX();
			x2 = points.getElement(p2).getX();
			x3 = points.getElement(p3).getX();

			y1 = points.getElement(p1).getY();
			y2 = points.getElement(p2).getY();
			y3 = points.getElement(p3).getY();

			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();

			if (x1 == x2 && x1 == x3) {
				if (x1 == 0) {
					face.getElement(i).setGranica(cc);
				}
				else if (x1 == maxX) {
					face.getElement(i).setGranica(dd);
				}

			}
			else if (y1 == y2 && y1 == y3) {
				if (y1 == 0) {
					face.getElement(i).setGranica(ee);
				}
				else if (y1 == maxY) {
					face.getElement(i).setGranica(ff);
				}

			}
			else if (z1 == z2 && z1 == z3) {
				if (z1 == 0) {
					face.getElement(i).setGranica(gg);
				}
				else if (z1 == maxZ) {
					face.getElement(i).setGranica(hh);
				}

			}


		}
	}


	ofstream zapis(nazwa);
	zapis.precision(12);
	globalSkal = 1;
	if (wartoscieDo1) {
		dzielnik = 1 / maxXX;
	}
	else {
		dzielnik = 1;
	}


	for (int i = 0, ileP = points.getIter(); i < ileP; i++) {

		zapis << "GRID" << pause << i + 1 << " " << pause << points.getElement(i).getX()*dzielnik << " " << pause << points.getElement(i).getY()*dzielnik << " " << pause << points.getElement(i).getZ()*dzielnik << "  " << 0 << endl;

	}

	zapis << endl << endl;

	int i = 0;
	for (int ileE = elements.getIter(); i < ileE; ++i) {

		zapis << "CTETRA" << pause << i + 1 << " " << pause << 1 << "   " << elements.getElement(i).getP1() + 1 << " " << pause << elements.getElement(i).getP2() + 1 << " " << pause << elements.getElement(i).getP3() + 1 << " " << pause << elements.getElement(i).getP4() + 1 << endl;

	}

	zapis << endl << endl;

	for (int ileE = pryzmy.getIter(), j = 0; j < ileE; ++i, ++j) {


		zapis << "CPRIZM" << pause << i + 1 << " " << pause << 1 << "   " << pryzmy.getElement(j).getP1() + 1 << " " << pause << pryzmy.getElement(j).getP2() + 1 << " " << pause << pryzmy.getElement(j).getP3() + 1 << " " << pause << pryzmy.getElement(j).getP4() + 1 << " " << pause << pryzmy.getElement(j).getP5() + 1 << " " << pause << pryzmy.getElement(j).getP6() + 1 << endl;

	}
	zapis << endl << endl;


	for (int j = 0, ileF = face.getIter(), rodzaj = -1; j < ileF; ++i, ++j) {

		switch (face.getElement(j).getGranica()) {
			case 'a':
				rodzaj = 0;
				break;
			case 'c':
				rodzaj = 3;
				break;
			case 'd':
				rodzaj = 4;
				break;
			case 'e':
				rodzaj = 5;
				break;
			case 'f':
				rodzaj = 6;
				break;
			case 'g':
				rodzaj = 1;
				break;
			case 'h':
				rodzaj = 2;
				break;
			case 'b':
				rodzaj = 9;
				break;
		}


		if (rodzaj) {
			zapis << "CTRIAX" << pause << i + 1 << " " << pause << rodzaj << "   " << face.getElement(j).getX() + 1 << " " << pause << face.getElement(j).getY() + 1 << " " << pause << face.getElement(j).getZ() + 1 << endl;
		}
		else {
			--i;
		}

	}



	zapis.close();
}


void Delanouy3D::stosunekRdorCout() {
//bool flaga=false;
	double p, a, b, c, d, e, f, V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;
//double maxX = dx;
	int z5 = 0, z6 = 0, z7 = 0, z8 = 0, z9 = 0, z10 = 0, z11 = 0, z12 = 0, z13 = 0, z15 = 0, z17 = 0, z19 = 0, z21 = 0, z23 = 0, zZ = 0;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {


		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		//R=sqrt(elements.getElement(i).getR());

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);

		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );




		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;
		p = 0.5 * (a + d + e);

		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		r = (V * 3) / pol;

		R = R / r;

		if (R < 5) {
			++z5;
		}
		else if (R < 6) {
			++z6;
		}
		else if (R < 7) {
			++z7;
		}
		else if (R < 8) {
			++z8;
		}
		else if (R < 9) {
			++z9;
		}
		else if (R < 10) {
			++z10;
		}
		else if (R < 11) {
			++z11;
		}
		else if (R < 12) {
			++z12;
		}
		else if (R < 13) {
			++z13;
		}
		else if (R < 15) {
			++z15;
		}
		else if (R < 17) {
			++z17;
		}
		else if (R < 19) {
			++z19;
		}
		else if (R < 21) {
			++z21;
		}
		else if (R < 23) {
			++z23;
		}
		else {
			++zZ;


			if (R < 0) {
				++z5;
			}
			else if (R < 0.2) {
				++z6;
			}
			else if (R < 0.3) {
				++z7;
			}
			else if (R < 0.4) {
				++z8;
			}
			else if (R < 0.5) {
				++z9;
			}
			else if (R < 4.5) {
				++z10;
			}
			else if (R < 4.6) {
				++z11;
			}
			else if (R < 4.7) {
				++z12;
			}
			else if (R < 4.8) {
				++z13;
			}
			else if (R < 4.9) {
				++z15;
			}
			else if (R < 5) {
				++z17;
			}
			else if (R < 19) {
				++z19;
			}
			else if (R < 21) {
				++z21;
			}
			else if (R < 23) {
				++z23;
			}
			else {
				++zZ;
			}

			/*
			x1=elements.getElement(i).getSx();
			y1=elements.getElement(i).getSy();
			z1=elements.getElement(i).getSz();

				if(0<x1 && x1<max && 0<y1 && y1<max && 0<z1 && z1<max){
				points.setElement(x1,y1,z1);
				}
			*/

		}
	}

	cout << " R<5 : " << z5 << endl;
	cout << " R<6 : " << z6 << endl;
	cout << " R<7 : " << z7 << endl;
	cout << " R<8 : " << z8 << endl;
	cout << " R<9 : " << z9 << endl;
	cout << " R<10 : " << z10 << endl;
	cout << " R<11 : " << z11 << endl;
	cout << " R<12 : " << z12 << endl;
	cout << " R<13 : " << z13 << endl;
	cout << " R<15 : " << z15 << endl;
	cout << " R<17 : " << z17 << endl;
	cout << " R<19 : " << z19 << endl;
	cout << " R<21 : " << z21 << endl;
	cout << " R<23 : " << z23 << endl;
	cout << " R>23 : " << zZ << endl;


}


void Delanouy3D::stosunekRdorZElementow() {

//bool flaga=false;
	double R;
//double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,p,a,b,c,d,e,f,V,r,pol;
//double a2,a3,a4,b2,b3,b4,c2,c3,c4,l2,l3,l4,X,Y,Z,W;
//double maxX = dx;
	int z5 = 0, z6 = 0, z7 = 0, z8 = 0, z9 = 0, z10 = 0, z11 = 0, z12 = 0, z13 = 0, z15 = 0, z17 = 0, z19 = 0, z21 = 0, z23 = 0, zZ = 0;
	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {



		R = elements.getElement(i).getStosunekR_r();

		if (R < 5) {
			++z5;
		}
		else if (R < 6) {
			++z6;
		}
		else if (R < 7) {
			++z7;
		}
		else if (R < 8) {
			++z8;
		}
		else if (R < 9) {
			++z9;
		}
		else if (R < 10) {
			++z10;
		}
		else if (R < 11) {
			++z11;
		}
		else if (R < 12) {
			++z12;
		}
		else if (R < 13) {
			++z13;
		}
		else if (R < 15) {
			++z15;
		}
		else if (R < 17) {
			++z17;
		}
		else if (R < 19) {
			++z19;
		}
		else if (R < 21) {
			++z21;
		}
		else if (R < 23) {
			++z23;
		}
		else {
			++zZ;

			/*
			if(R < 0){++z5;}
			else if(R < 0.2){++z6;}
			else if(R < 0.3){++z7;}
			else if(R < 0.4){++z8;}
			else if(R < 0.5){++z9;}
			else if(R < 4.5){++z10;}
			else if(R < 4.6){++z11;}
			else if(R < 4.7){++z12;}
			else if(R < 4.8){++z13;}
			else if(R < 4.9){++z15;}
			else if(R < 5){++z17;}
			else if(R < 19){++z19;}
			else if(R < 21){++z21;}
			else if(R < 23){++z23;}
			else{++zZ;
			*/
			/*
			x1=elements.getElement(i).getSx();
			y1=elements.getElement(i).getSy();
			z1=elements.getElement(i).getSz();

				if(0<x1 && x1<max && 0<y1 && y1<max && 0<z1 && z1<max){
				points.setElement(x1,y1,z1);
				}
			*/

		}
	}
	/*
	AnsiString text;

	text=" R<5 :  ";text+=z5;
	m1->Lines->Add(text);
	text=" R<6 :   ";text+=z6;
	m1->Lines->Add(text);
	text=" R<7 :  ";text+=z7;
	m1->Lines->Add(text);
	text=" R<8 :   ";text+=z8;
	m1->Lines->Add(text);
	text=" R<9 :   ";text+=z9;
	m1->Lines->Add(text);
	text=" R<10 :   ";text+=z10;
	m1->Lines->Add(text);
	text=" R<11 :   ";text+=z11;
	m1->Lines->Add(text);
	text=" R<12 :   ";text+=z12;
	m1->Lines->Add(text);
	text=" R<13 :   ";text+=z13;
	m1->Lines->Add(text);
	text=" R<15 :   ";text+=z15;
	m1->Lines->Add(text);
	text=" R<17 :   ";text+=z17;
	m1->Lines->Add(text);
	text=" R<19 :   ";text+=z19;
	m1->Lines->Add(text);
	text=" R<21 :   ";text+=z21;
	m1->Lines->Add(text);
	text=" R<23 :   ";text+=z23;
	m1->Lines->Add(text);
	text=" R>23 :   ";text+=zZ;
	m1->Lines->Add(text);
	*/
}


bool Delanouy3D::stosunekRdor(double stosunek) {
	bool flaga = false;
	double p, a, b, c, d, e, f, V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;
//double maxX = dx;
	int z5 = 0, z6 = 0, z7 = 0, z8 = 0, z9 = 0, z10 = 0, z11 = 0, z12 = 0, z13 = 0, z15 = 0, z17 = 0, z19 = 0, z21 = 0, z23 = 0, zZ = 0;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {


		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		//R=sqrt(elements.getElement(i).getR());

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);

		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );




		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;
		p = 0.5 * (a + d + e);

		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		r = (V * 3) / pol;

		R = R / r;

		if (R > stosunek) {
			flaga = true;
		}

		if (R < 5) {
			++z5;
		}
		else if (R < 6) {
			++z6;
		}
		else if (R < 7) {
			++z7;
		}
		else if (R < 8) {
			++z8;
		}
		else if (R < 9) {
			++z9;
		}
		else if (R < 10) {
			++z10;
		}
		else if (R < 11) {
			++z11;
		}
		else if (R < 12) {
			++z12;
		}
		else if (R < 13) {
			++z13;
		}
		else if (R < 15) {
			++z15;
		}
		else if (R < 17) {
			++z17;
		}
		else if (R < 19) {
			++z19;
		}
		else if (R < 21) {
			++z21;
		}
		else if (R < 23) {
			++z23;
		}
		else {
			++zZ;

			/*
			if(R < 0){++z5;}
			else if(R < 0.2){++z6;}
			else if(R < 0.3){++z7;}
			else if(R < 0.4){++z8;}
			else if(R < 0.5){++z9;}
			else if(R < 4.5){++z10;}
			else if(R < 4.6){++z11;}
			else if(R < 4.7){++z12;}
			else if(R < 4.8){++z13;}
			else if(R < 4.9){++z15;}
			else if(R < 5){++z17;}
			else if(R < 19){++z19;}
			else if(R < 21){++z21;}
			else if(R < 23){++z23;}
			else{++zZ;
			*/
			/*
			x1=elements.getElement(i).getSx();
			y1=elements.getElement(i).getSy();
			z1=elements.getElement(i).getSz();

				if(0<x1 && x1<max && 0<y1 && y1<max && 0<z1 && z1<max){
				points.setElement(x1,y1,z1);
				}
			*/

		}
	}
	/*
	AnsiString text;

	text=" R<5 :  ";text+=z5;
	m1->Lines->Add(text);
	text=" R<6 :   ";text+=z6;
	m1->Lines->Add(text);
	text=" R<7 :  ";text+=z7;
	m1->Lines->Add(text);
	text=" R<8 :   ";text+=z8;
	m1->Lines->Add(text);
	text=" R<9 :   ";text+=z9;
	m1->Lines->Add(text);
	text=" R<10 :   ";text+=z10;
	m1->Lines->Add(text);
	text=" R<11 :   ";text+=z11;
	m1->Lines->Add(text);
	text=" R<12 :   ";text+=z12;
	m1->Lines->Add(text);
	text=" R<13 :   ";text+=z13;
	m1->Lines->Add(text);
	text=" R<15 :   ";text+=z15;
	m1->Lines->Add(text);
	text=" R<17 :   ";text+=z17;
	m1->Lines->Add(text);
	text=" R<19 :   ";text+=z19;
	m1->Lines->Add(text);
	text=" R<21 :   ";text+=z21;
	m1->Lines->Add(text);
	text=" R<23 :   ";text+=z23;
	m1->Lines->Add(text);
	text=" R>23 :   ";text+=zZ;
	m1->Lines->Add(text);
	*/
	return flaga;
}

void Delanouy3D::stosunekRdorAktywneZiarna(IntList & nrAktywne) {

	double p, a, b, c, d, e, f, V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;
//double maxX = dx;
	int z5 = 0, z6 = 0, z7 = 0, z8 = 0, z9 = 0, z10 = 0, z11 = 0, z12 = 0, z13 = 0, z15 = 0, z17 = 0, z19 = 0, z21 = 0, z23 = 0, zZ = 0;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		if (!nrAktywne.sprawdzCzyJest(elements.getElement(i).getRodzajZiarna())) {


			p1 = elements.getElement(i).getP1();
			p2 = elements.getElement(i).getP2();
			p3 = elements.getElement(i).getP3();
			p4 = elements.getElement(i).getP4();
			//R=sqrt(elements.getElement(i).getR());

			x1 = points.getElement(p1).getX();
			x2 = points.getElement(p2).getX();
			x3 = points.getElement(p3).getX();
			x4 = points.getElement(p4).getX();
			y1 = points.getElement(p1).getY();
			y2 = points.getElement(p2).getY();
			y3 = points.getElement(p3).getY();
			y4 = points.getElement(p4).getY();
			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();
			z4 = points.getElement(p4).getZ();


			a2 = (-x3 + x1) * 2;
			b2 = (-y3 + y1) * 2;
			c2 = (-z3 + z1) * 2;
			l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

			a3 = (-x4 + x1) * 2;
			b3 = (-y4 + y1) * 2;
			c3 = (-z4 + z1) * 2;
			l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);

			a4 = (-x4 + x2) * 2;
			b4 = (-y4 + y2) * 2;
			c4 = (-z4 + z2) * 2;
			l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

			W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
			W = 1 / W;
			l2 *= W;
			l3 *= W;
			l4 *= W;

			X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
			Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
			Y = Y * -1;
			Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

			R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );




			a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
			b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
			c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
			d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
			e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
			f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

			pol = 0;
			p = 0.5 * (a + d + e);

			pol += sqrt(p * (p - a) * (p - d) * (p - e));
			p = 0.5 * (b + f + e);
			pol += sqrt(p * (p - b) * (p - f) * (p - e));
			p = 0.5 * (c + d + f);
			pol += sqrt(p * (p - c) * (p - d) * (p - f));
			p = 0.5 * (a + b + c);
			pol += sqrt(p * (p - a) * (p - b) * (p - c));

			V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
							(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

			r = (V * 3) / pol;

			R = R / r;

			if (R < 5) {
				++z5;
			}
			else if (R < 6) {
				++z6;
			}
			else if (R < 7) {
				++z7;
			}
			else if (R < 8) {
				++z8;
			}
			else if (R < 9) {
				++z9;
			}
			else if (R < 10) {
				++z10;
			}
			else if (R < 11) {
				++z11;
			}
			else if (R < 12) {
				++z12;
			}
			else if (R < 13) {
				++z13;
			}
			else if (R < 15) {
				++z15;
			}
			else if (R < 17) {
				++z17;
			}
			else if (R < 19) {
				++z19;
			}
			else if (R < 21) {
				++z21;
			}
			else if (R < 23) {
				++z23;
			}
			else {
				++zZ;

				/*
				if(R < 0){++z5;}
				else if(R < 0.2){++z6;}
				else if(R < 0.3){++z7;}
				else if(R < 0.4){++z8;}
				else if(R < 0.5){++z9;}
				else if(R < 4.5){++z10;}
				else if(R < 4.6){++z11;}
				else if(R < 4.7){++z12;}
				else if(R < 4.8){++z13;}
				else if(R < 4.9){++z15;}
				else if(R < 5){++z17;}
				else if(R < 19){++z19;}
				else if(R < 21){++z21;}
				else if(R < 23){++z23;}
				else{++zZ;
				*/
				/*
				x1=elements.getElement(i).getSx();
				y1=elements.getElement(i).getSy();
				z1=elements.getElement(i).getSz();

					if(0<x1 && x1<max && 0<y1 && y1<max && 0<z1 && z1<max){
					points.setElement(x1,y1,z1);
					}
				*/

			}
		}
	}
	/*
	AnsiString text;

	text=" R<5 :  ";text+=z5;
	m1->Lines->Add(text);
	text=" R<6 :   ";text+=z6;
	m1->Lines->Add(text);
	text=" R<7 :  ";text+=z7;
	m1->Lines->Add(text);
	text=" R<8 :   ";text+=z8;
	m1->Lines->Add(text);
	text=" R<9 :   ";text+=z9;
	m1->Lines->Add(text);
	text=" R<10 :   ";text+=z10;
	m1->Lines->Add(text);
	text=" R<11 :   ";text+=z11;
	m1->Lines->Add(text);
	text=" R<12 :   ";text+=z12;
	m1->Lines->Add(text);
	text=" R<13 :   ";text+=z13;
	m1->Lines->Add(text);
	text=" R<15 :   ";text+=z15;
	m1->Lines->Add(text);
	text=" R<17 :   ";text+=z17;
	m1->Lines->Add(text);
	text=" R<19 :   ";text+=z19;
	m1->Lines->Add(text);
	text=" R<21 :   ";text+=z21;
	m1->Lines->Add(text);
	text=" R<23 :   ";text+=z23;
	m1->Lines->Add(text);
	text=" R>23 :   ";text+=zZ;
	m1->Lines->Add(text);
	*/
}

void Delanouy3D::stosunekGamaKGMSHAktywneZiarna(IntList & nrAktywne) {

	double p, a, b, c, d, e, f, V, r, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
//double maxX = dx;
	int z95 = 0, z90 = 0, z85 = 0, z80 = 0, z75 = 0, z70 = 0, z65 = 0, z60 = 0, z55 = 0, z50 = 0, z45 = 0, z40 = 0, z35 = 0, z30 = 0, z25 = 0, z20 = 0, z15 = 0, z10 = 0, z05 = 0, zZ;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {

		if (!nrAktywne.sprawdzCzyJest(elements.getElement(i).getRodzajZiarna())) {


			p1 = elements.getElement(i).getP1();
			p2 = elements.getElement(i).getP2();
			p3 = elements.getElement(i).getP3();
			p4 = elements.getElement(i).getP4();
			//R=sqrt(elements.getElement(i).getR());

			x1 = points.getElement(p1).getX();
			x2 = points.getElement(p2).getX();
			x3 = points.getElement(p3).getX();
			x4 = points.getElement(p4).getX();
			y1 = points.getElement(p1).getY();
			y2 = points.getElement(p2).getY();
			y3 = points.getElement(p3).getY();
			y4 = points.getElement(p4).getY();
			z1 = points.getElement(p1).getZ();
			z2 = points.getElement(p2).getZ();
			z3 = points.getElement(p3).getZ();
			z4 = points.getElement(p4).getZ();

			a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
			b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
			c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
			d = sqrt((x4 - x1) * (x4 - x1) + (y4 - y1) * (y4 - y1) + (z4 - z1) * (z4 - z1));
			e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
			f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

			pol = 0;
			p = 0.5 * (a + d + e);

			pol += sqrt(p * (p - a) * (p - d) * (p - e));
			p = 0.5 * (b + f + e);
			pol += sqrt(p * (p - b) * (p - f) * (p - e));
			p = 0.5 * (c + d + f);
			pol += sqrt(p * (p - c) * (p - d) * (p - f));
			p = 0.5 * (a + b + c);
			pol += sqrt(p * (p - a) * (p - b) * (p - c));

			V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
							(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;



			double maxEdge = a;
			if (maxEdge < b) {
				maxEdge = b;
			}
			if (maxEdge < c) {
				maxEdge = c;
			}
			if (maxEdge < d) {
				maxEdge = d;
			}
			if (maxEdge < e) {
				maxEdge = e;
			}
			if (maxEdge < f) {
				maxEdge = f;
			}

			r = (V * 14.696936) / (pol * maxEdge);

			if (r > 0.95) {
				++z95;
			}
			else if (r > 0.90) {
				++z90;
			}
			else if (r > 0.85) {
				++z85;
			}
			else if (r > 0.80) {
				++z80;
			}
			else if (r > 0.75) {
				++z75;
			}
			else if (r > 0.7) {
				++z70;
			}
			else if (r > 0.65) {
				++z65;
			}
			else if (r > 0.6) {
				++z60;
			}
			else if (r > 0.55) {
				++z55;
			}
			else if (r > 0.5) {
				++z50;
			}
			else if (r > 0.45) {
				++z45;
			}
			else if (r > 0.4) {
				++z40;
			}
			else if (r > 0.35) {
				++z35;
			}
			else if (r > 0.3) {
				++z30;
			}
			else if (r > 0.25) {
				++z25;
			}
			else if (r > 0.2) {
				++z20;
			}
			else if (r > 0.15) {
				++z15;
			}
			else if (r > 0.1) {
				++z10;
			}
			else if (r > 0.05) {
				++z05;
			}
			else {
				++zZ;


			}
		}
	}
	/*
	AnsiString text;

	text=" >0.95 :  ";text+=z95;
	m1->Lines->Add(text);
	text=" >0.90 :   ";text+=z90;
	m1->Lines->Add(text);
	text=" >0.85 :  ";text+=z85;
	m1->Lines->Add(text);
	text=" >0.80 :   ";text+=z80;
	m1->Lines->Add(text);
	text=" >0.75 :   ";text+=z75;
	m1->Lines->Add(text);
	text=" >0.70 :   ";text+=z70;
	m1->Lines->Add(text);
	text=" >0.65 :   ";text+=z65;
	m1->Lines->Add(text);
	text=" >0.60 :   ";text+=z60;
	m1->Lines->Add(text);
	text=" >0.55 :   ";text+=z55;
	m1->Lines->Add(text);
	text=" >0.50 :   ";text+=z50;
	m1->Lines->Add(text);
	text=" >0.45 :   ";text+=z45;
	m1->Lines->Add(text);
	text=" >0.40 :   ";text+=z40;
	m1->Lines->Add(text);
	text=" >0.35 :   ";text+=z35;
	m1->Lines->Add(text);
	text=" >0.30 :   ";text+=z30;
	m1->Lines->Add(text);
	text=" >0.25 :   ";text+=z25;
	m1->Lines->Add(text);
	text=" >0.20 :   ";text+=z20;
	m1->Lines->Add(text);
	text=" >0.15 :   ";text+=z15;
	m1->Lines->Add(text);
	text=" >0.10 :   ";text+=z10;
	m1->Lines->Add(text);
	text=" >0.05 :   ";text+=z05;
	m1->Lines->Add(text);
	*/
}

void Delanouy3D::stosunekGamaKGMSH() {

	double p, a, b, c, d, e, f, V, r, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
//double maxX = dx;
	int z95 = 0, z90 = 0, z85 = 0, z80 = 0, z75 = 0, z70 = 0, z65 = 0, z60 = 0, z55 = 0, z50 = 0, z45 = 0, z40 = 0, z35 = 0, z30 = 0, z25 = 0, z20 = 0, z15 = 0, z10 = 0, z05 = 0, zZ;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {




		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		//R=sqrt(elements.getElement(i).getR());

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();

		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x4 - x1) * (x4 - x1) + (y4 - y1) * (y4 - y1) + (z4 - z1) * (z4 - z1));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;
		p = 0.5 * (a + d + e);

		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;



		double maxEdge = a;
		if (maxEdge < b) {
			maxEdge = b;
		}
		if (maxEdge < c) {
			maxEdge = c;
		}
		if (maxEdge < d) {
			maxEdge = d;
		}
		if (maxEdge < e) {
			maxEdge = e;
		}
		if (maxEdge < f) {
			maxEdge = f;
		}

		r = (V * 14.696936) / (pol * maxEdge);

		if (r > 0.95) {
			++z95;
		}
		else if (r > 0.90) {
			++z90;
		}
		else if (r > 0.85) {
			++z85;
		}
		else if (r > 0.80) {
			++z80;
		}
		else if (r > 0.75) {
			++z75;
		}
		else if (r > 0.7) {
			++z70;
		}
		else if (r > 0.65) {
			++z65;
		}
		else if (r > 0.6) {
			++z60;
		}
		else if (r > 0.55) {
			++z55;
		}
		else if (r > 0.5) {
			++z50;
		}
		else if (r > 0.45) {
			++z45;
		}
		else if (r > 0.4) {
			++z40;
		}
		else if (r > 0.35) {
			++z35;
		}
		else if (r > 0.3) {
			++z30;
		}
		else if (r > 0.25) {
			++z25;
		}
		else if (r > 0.2) {
			++z20;
		}
		else if (r > 0.15) {
			++z15;
		}
		else if (r > 0.1) {
			++z10;
		}
		else if (r > 0.05) {
			++z05;
		}
		else {
			++zZ;


		}

	}

	/*
	AnsiString text;

	text=" >0.95 :  ";text+=z95;
	m1->Lines->Add(text);
	text=" >0.90 :   ";text+=z90;
	m1->Lines->Add(text);
	text=" >0.85 :  ";text+=z85;
	m1->Lines->Add(text);
	text=" >0.80 :   ";text+=z80;
	m1->Lines->Add(text);
	text=" >0.75 :   ";text+=z75;
	m1->Lines->Add(text);
	text=" >0.70 :   ";text+=z70;
	m1->Lines->Add(text);
	text=" >0.65 :   ";text+=z65;
	m1->Lines->Add(text);
	text=" >0.60 :   ";text+=z60;
	m1->Lines->Add(text);
	text=" >0.55 :   ";text+=z55;
	m1->Lines->Add(text);
	text=" >0.50 :   ";text+=z50;
	m1->Lines->Add(text);
	text=" >0.45 :   ";text+=z45;
	m1->Lines->Add(text);
	text=" >0.40 :   ";text+=z40;
	m1->Lines->Add(text);
	text=" >0.35 :   ";text+=z35;
	m1->Lines->Add(text);
	text=" >0.30 :   ";text+=z30;
	m1->Lines->Add(text);
	text=" >0.25 :   ";text+=z25;
	m1->Lines->Add(text);
	text=" >0.20 :   ";text+=z20;
	m1->Lines->Add(text);
	text=" >0.15 :   ";text+=z15;
	m1->Lines->Add(text);
	text=" >0.10 :   ";text+=z10;
	m1->Lines->Add(text);
	text=" >0.05 :   ";text+=z05;
	m1->Lines->Add(text);
	*/
}


double Delanouy3D::stosunekRdorNajwiekszyDlaPunktuDODelanoya(int ktoryP, int OdElementu) {

	double p, a, b, c, d, e, f, V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;
//double maxX = dx;
	double maxR = 0;

	for (int i = OdElementu, ktoryE, ileE = elements.getIter(), p1, p2, p3, p4; i < ileE; ++i) {

		ktoryE = i;

		p1 = elements.getElement(ktoryE).getP1();
		p2 = elements.getElement(ktoryE).getP2();
		p3 = elements.getElement(ktoryE).getP3();
		p4 = elements.getElement(ktoryE).getP4();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);
		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);
		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );

		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;
		p = 0.5 * (a + d + e);

		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		if (V < 0) {
			return 999999999;
		}

		r = (V * 3) / pol;

		R = R / r;

		if (maxR < R) {
			maxR = R;
		}
	}


	return maxR;
}


double Delanouy3D::stosunekRdorNajwiekszyDlaPunktu(int ktoryP) {

	double V, r, R, pol;
	double Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, nx, ny, nz;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, X, Y, Z, W;
//double maxX = dx;
	double maxR = 0;
	double maxV = 0;

	for (int i = 0, ktoryE, ileE = laplas[ktoryP].getIter(), p1, p2, p3, p4; i < ileE; ++i) {

		ktoryE = laplas[ktoryP].getElement(i);

		p1 = elements.getElement(ktoryE).getP1();
		p2 = elements.getElement(ktoryE).getP2();
		p3 = elements.getElement(ktoryE).getP3();
		p4 = elements.getElement(ktoryE).getP4();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		a2 = (-x3 + x1) * 2;
		b2 = (-y3 + y1) * 2;
		c2 = (-z3 + z1) * 2;
		l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);
		a3 = (-x4 + x1) * 2;
		b3 = (-y4 + y1) * 2;
		c3 = (-z4 + z1) * 2;
		l3 = -(x4 * x4 - x1 * x1 + y4 * y4 - y1 * y1 + z4 * z4 - z1 * z1);
		a4 = (-x4 + x2) * 2;
		b4 = (-y4 + y2) * 2;
		c4 = (-z4 + z2) * 2;
		l4 = -(x4 * x4 - x2 * x2 + y4 * y4 - y2 * y2 + z4 * z4 - z2 * z2);

		W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;
		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

		R =  sqrt( (X - x4) * (X - x4) + (Y - y4) * (Y - y4) + (Z - z4) * (Z - z4) );

		pol = 0;
		Ax = x2 - x1;
		Bx = x3 - x1;
		Cx = x4 - x1;
		Ay = y2 - y1;
		By = y3 - y1;
		Cy = y4 - y1;
		Az = z2 - z1;
		Bz = z3 - z1;
		Cz = z4 - z1;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = Ay * Cz - Cy * Az;
		ny = Cx * Az - Ax * Cz;
		nz = Ax * Cy - Cx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		nx = By * Cz - Cy * Bz;
		ny = Cx * Bz - Bx * Cz;
		nz = Bx * Cy - Cx * By;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;

		V = (nx * Ax + ny * Ay + nz * Az) * 0.5;

		Ax = x2 - x4;
		Bx = x3 - x4;
		Ay = y2 - y4;
		By = y3 - y4;
		Az = z2 - z4;
		Bz = z3 - z4;

		nx = Ay * Bz - By * Az;
		ny = Bx * Az - Ax * Bz;
		nz = Ax * By - Bx * Ay;
		pol += sqrt(nx * nx + ny * ny + nz * nz) * 0.5;


		r = V / pol;

		if (V <= 0) {
			V = (V * -1) + 100000;
			if (maxV < V) {
				maxV = V;
			}
		}

		R = R / r;

		if (maxR < R) {
			maxR = R;
		}
	}

	if (maxV) {
		return maxV;
	}
	return maxR;

}

/*
double Delanouy3D::stosunekRdorNajwiekszyDlaPunktu(int ktoryP){

double p,a,b,c,d,e,f,V,r,R,pol;
double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
double a2,a3,a4,b2,b3,b4,c2,c3,c4,l2,l3,l4,X,Y,Z,W;
//double maxX = dx;
double maxR=0;
double maxV=0;

	for(int i=0,ktoryE,ileE=laplas[ktoryP].getIter(),p1,p2,p3,p4;i<ileE;++i){

		ktoryE=laplas[ktoryP].getElement(i);

		p1=elements.getElement(ktoryE).getP1();
		p2=elements.getElement(ktoryE).getP2();
		p3=elements.getElement(ktoryE).getP3();
		p4=elements.getElement(ktoryE).getP4();

		x1=points.getElement(p1).getX();
		x2=points.getElement(p2).getX();
		x3=points.getElement(p3).getX();
		x4=points.getElement(p4).getX();
		y1=points.getElement(p1).getY();
		y2=points.getElement(p2).getY();
		y3=points.getElement(p3).getY();
		y4=points.getElement(p4).getY();
		z1=points.getElement(p1).getZ();
		z2=points.getElement(p2).getZ();
		z3=points.getElement(p3).getZ();
		z4=points.getElement(p4).getZ();


		a2 = (-x3+x1)*2;b2 = (-y3+y1)*2;c2 = (-z3+z1)*2;l2 = -(x3*x3-x1*x1 + y3*y3-y1*y1 + z3*z3-z1*z1);
		a3 = (-x4+x1)*2;b3 = (-y4+y1)*2;c3 = (-z4+z1)*2;l3 = -(x4*x4-x1*x1 + y4*y4-y1*y1 + z4*z4-z1*z1);
		a4 = (-x4+x2)*2;b4 = (-y4+y2)*2;c4 = (-z4+z2)*2;l4 = -(x4*x4-x2*x2 + y4*y4-y2*y2 + z4*z4-z2*z2);

		W =  a2*b3*c4+a3*b4*c2+a4*b2*c3-c2*b3*a4-c3*b4*a2-c4*b2*a3;
		W=1/W;
		l2*=W;l3*=W;l4*=W;

		X = (b2*c3*l4+b3*c4*l2+b4*c2*l3-l2*c3*b4-l3*c4*b2-l4*c2*b3);
		Y = (a2*c3*l4+a3*c4*l2+a4*c2*l3-l2*c3*a4-l3*c4*a2-l4*c2*a3);Y=Y*-1;
		Z = (a2*b3*l4+a3*b4*l2+a4*b2*l3-l2*b3*a4-l3*b4*a2-l4*b2*a3);

		R =  sqrt( (X-x4)*(X-x4) + (Y-y4)*(Y-y4) + (Z-z4)*(Z-z4) );

		a=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
		b=sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2));
		c=sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3));
		d=sqrt((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4)+(z1-z4)*(z1-z4));
		e=sqrt((x4-x2)*(x4-x2)+(y4-y2)*(y4-y2)+(z4-z2)*(z4-z2));
		f=sqrt((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3)+(z4-z3)*(z4-z3));

		pol=0;
		p=0.5*(a+d+e);

		pol+=sqrt(p*(p-a)*(p-d)*(p-e));
		p=0.5*(b+f+e);
		pol+=sqrt(p*(p-b)*(p-f)*(p-e));
		p=0.5*(c+d+f);
		pol+=sqrt(p*(p-c)*(p-d)*(p-f));
		p=0.5*(a+b+c);
		pol+=sqrt(p*(p-a)*(p-b)*(p-c));

		V = ((x2-x1)*(y3-y1)*(z4-z1)+(x3-x1)*(y4-y1)*(z2-z1)+(x4-x1)*(y2-y1)*(z3-z1)-
		(z2-z1)*(y3-y1)*(x4-x1)-(z3-z1)*(y4-y1)*(x2-x1)-(z4-z1)*(y2-y1)*(x3-x1))/6;


		r = (V*3)/pol;

		if(V<=0){V = (V*-1)+100000; if(maxV<V){maxV=V;}}

		R= R/r;

		if(maxR<R){maxR=R;}
	}

if(maxV){return maxV;}
return maxR;

}
*/
/*
double Delanouy3D::naprawaElementow_Objetosci(int ktoryP){

double p,a,b,c,d,e,f,V,r,R,pol;
double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
double a2,a3,a4,b2,b3,b4,c2,c3,c4,l2,l3,l4,X,Y,Z,W;
double maxX = dx;
double maxR=0;

	for(int i=0,ktoryE,ileE=laplas[ktoryP].getIter(),p1,p2,p3,p4;i<ileE;++i){

		ktoryE=laplas[ktoryP].getElement(i);

		p1=elements.getElement(ktoryE).getP1();
		p2=elements.getElement(ktoryE).getP2();
		p3=elements.getElement(ktoryE).getP3();
		p4=elements.getElement(ktoryE).getP4();

		x1=points.getElement(p1).getX();
		x2=points.getElement(p2).getX();
		x3=points.getElement(p3).getX();
		x4=points.getElement(p4).getX();
		y1=points.getElement(p1).getY();
		y2=points.getElement(p2).getY();
		y3=points.getElement(p3).getY();
		y4=points.getElement(p4).getY();
		z1=points.getElement(p1).getZ();
		z2=points.getElement(p2).getZ();
		z3=points.getElement(p3).getZ();
		z4=points.getElement(p4).getZ();




		a2 = (-x3+x1)*2;b2 = (-y3+y1)*2;c2 = (-z3+z1)*2;l2 = -(x3*x3-x1*x1 + y3*y3-y1*y1 + z3*z3-z1*z1);
		a3 = (-x4+x1)*2;b3 = (-y4+y1)*2;c3 = (-z4+z1)*2;l3 = -(x4*x4-x1*x1 + y4*y4-y1*y1 + z4*z4-z1*z1);
		a4 = (-x4+x2)*2;b4 = (-y4+y2)*2;c4 = (-z4+z2)*2;l4 = -(x4*x4-x2*x2 + y4*y4-y2*y2 + z4*z4-z2*z2);

		W =  a2*b3*c4+a3*b4*c2+a4*b2*c3-c2*b3*a4-c3*b4*a2-c4*b2*a3;
		W=1/W;
		l2*=W;l3*=W;l4*=W;

		X = (b2*c3*l4+b3*c4*l2+b4*c2*l3-l2*c3*b4-l3*c4*b2-l4*c2*b3);
		Y = (a2*c3*l4+a3*c4*l2+a4*c2*l3-l2*c3*a4-l3*c4*a2-l4*c2*a3);Y=Y*-1;
		Z = (a2*b3*l4+a3*b4*l2+a4*b2*l3-l2*b3*a4-l3*b4*a2-l4*b2*a3);

		R =  sqrt( (X-x4)*(X-x4) + (Y-y4)*(Y-y4) + (Z-z4)*(Z-z4) );

		a=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
		b=sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2));
		c=sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3));
		d=sqrt((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4)+(z1-z4)*(z1-z4));
		e=sqrt((x4-x2)*(x4-x2)+(y4-y2)*(y4-y2)+(z4-z2)*(z4-z2));
		f=sqrt((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3)+(z4-z3)*(z4-z3));

		pol=0;
		p=0.5*(a+d+e);

		pol+=sqrt(p*(p-a)*(p-d)*(p-e));
		p=0.5*(b+f+e);
		pol+=sqrt(p*(p-b)*(p-f)*(p-e));
		p=0.5*(c+d+f);
		pol+=sqrt(p*(p-c)*(p-d)*(p-f));
		p=0.5*(a+b+c);
		pol+=sqrt(p*(p-a)*(p-b)*(p-c));

		V = ((x2-x1)*(y3-y1)*(z4-z1)+(x3-x1)*(y4-y1)*(z2-z1)+(x4-x1)*(y2-y1)*(z3-z1)-
		(z2-z1)*(y3-y1)*(x4-x1)-(z3-z1)*(y4-y1)*(x2-x1)-(z4-z1)*(y2-y1)*(x3-x1))/6;

		if(V<0){return 999999999;}

		r = (V*3)/pol;

		R= R/r;

		if(maxR<R){maxR=R;}
	}


return maxR;
}
*/

double Delanouy3D::stosunekGMSHNajwiekszyDlaPunktu(int ktoryP) {

	double p, a, b, c, d, e, f, V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double maxEdge;
//double maxX = dx;
	double maxR = 0;

	for (int i = 0, ktoryE, ileE = laplas[ktoryP].getIter(), p1, p2, p3, p4; i < ileE; ++i) {

		ktoryE = laplas[ktoryP].getElement(i);

		p1 = elements.getElement(ktoryE).getP1();
		p2 = elements.getElement(ktoryE).getP2();
		p3 = elements.getElement(ktoryE).getP3();
		p4 = elements.getElement(ktoryE).getP4();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();

		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;
		p = 0.5 * (a + d + e);

		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		if (V < 0) {
			return 999999999;
		}

		maxEdge = a;
		if (maxEdge < b) {
			maxEdge = b;
		}
		if (maxEdge < c) {
			maxEdge = c;
		}
		if (maxEdge < d) {
			maxEdge = d;
		}
		if (maxEdge < e) {
			maxEdge = e;
		}
		if (maxEdge < f) {
			maxEdge = f;
		}

		r = (V * 14.696936) / (pol * maxEdge);

		R = 1 - r;

		if (maxR < R) {
			maxR = R;
		}
	}


	return maxR;
}


double Delanouy3D::stosunekVDlaPunktu(int ktoryP) {

	double V;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
//double maxX = dx;
	double maxV = 0;
	double minV = 99999999;

	for (int i = 0, ktoryE, ileE = laplas[ktoryP].getIter(), p1, p2, p3, p4; i < ileE; ++i) {

		ktoryE = laplas[ktoryP].getElement(i);

		p1 = elements.getElement(ktoryE).getP1();
		p2 = elements.getElement(ktoryE).getP2();
		p3 = elements.getElement(ktoryE).getP3();
		p4 = elements.getElement(ktoryE).getP4();

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();


		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		if (V < 0) {
			return 999999999;
		}

		if (maxV < V) {
			maxV = V;
		}
		if (minV > V) {
			minV = V;
		}
	}


	return maxV - minV;
}


int Delanouy3D::poprawaElement(int model, double stosunekR_r, bool PoprawZPozaObszaru) {

	PunktList noweP(1000, 1000);
	IntList ktoreP(1000, 1000);
	IntList ktoreE(10000, 10000);

	int a;
//    double sx,sy,sz;



	stosunekRdorWyszukajElemnty(stosunekR_r);
	a = ustalNowePunkty(wybraneZleEl, PoprawZPozaObszaru, noweP);



	int ileP = points.getIter();

	for (int i = 0, ileE = noweP.getIter(); i < ileE; ++i) {
		ktoreP.setElement(ileP + i);
	}


	delanoyTablicaOPDodajPunkty(dx, dy, dz, false, noweP, 0.3, 3);

	switch (model) {

		case 1:
			optymalizacjaStosunekRdorWybranePunkty(0.1, ktoreP, 4);
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);

			break;
		case 2:
			sasiedniePunkty();
			wygladzanieLaplaceWybranePunkty(ktoreP);
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);
			break;
		case 3:
			sasiedniePunkty();
			for (int i = 0; i < 4; ++i) {
				wygladzanieLaplaceWagaWybranePunkty(2, 0.3, ktoreP, 1);
			}
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);
			break;

	}

	return a;
}



int Delanouy3D::poprawaElement_popZiarna(int model, double stosunekR_r, bool PoprawZPozaObszaru) {

	PunktList noweP(1000, 1000);
	IntList ktoreP(1000, 1000);
	IntList ktoreE(10000, 10000);

	int a;
	//double sx,sy,sz;
	stosunekRdorWyszukajElemnty(stosunekR_r);
	a = ustalNowePunkty(wybraneZleEl, PoprawZPozaObszaru, noweP);

	int ileP = points.getIter();

	for (int i = 0, ileE = noweP.getIter(); i < ileE; ++i) {
		ktoreP.setElement(ileP + i);
	}

	delanoyTablicaOPDodajPunkty(dx, dy, dz, false, noweP, 0.3, 3);

	switch (model) {
		case 1:
			optymalizacjaStosunekRdorWybranePunkty(0.1, ktoreP, 4);
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);
			break;
	}



	return a;
}

void Delanouy3D::znajdzElemetyPoDodaniuPunktow(IntList & ktoreP, IntList & ktoreE) {

	sasiednieElementy();
	ktoreE.ustawIter(0);

	for (int i = 0, ilePunktow = ktoreP.getIter(), p; i < ilePunktow; ++i) {

		p = ktoreP.getElement(i);

		for (int j = 0, ileElementow = laplas[p].getIter(); j < ileElementow; ++j) {

			ktoreE.setElement(laplas[p].getElement(j));

		}

	}
}



void Delanouy3D::poprawaElementWszystkieRdorPoprawaMiejscowa(int model, double stosunekR_r, bool PoprawZPozaObszaru) {

	int a = 1;
	int max = 0;

	while (max < 50) {

		if (a == 0) {
			break;
		}
		if (model == 1) {
			reMes();
		}
		a = poprawaElement(model, stosunekR_r, PoprawZPozaObszaru);
		++max;
	}


}

int Delanouy3D::poprawaElementMalejaco(int ** tabelementow, int model, double & stosunekR_r, double stosunekDocelowy, double & zmniejszO, bool PoprawZPozaObszaru) {

	PunktList noweP(1000, 1000);
	IntList ktoreP(1000, 1000);
	IntList ktoreE(4000, 4000);
	/////////////////////////////////////////////////

	int a, b;
	//double sx,sy,sz;
	//bool zwieksz;

	stosunekRdorWyszukajElemntyMalejacoSprawdz(stosunekR_r);
	b = wybraneZleEl.getIter();

	//okreslam kierunek szuaknia +/- stosunekR_r
	/*
	if(b<50){zwieksz=false;}
	else if(b>200){zwieksz=true;}

	if(zwieksz){
		while(1){

			stosunekRdorWyszukajElemntyMalejacoSprawdz(stosunekR_r);
			b=wybraneZleEl.getIter();

			if(b<50){       zmniejszO*=0.5;stosunekR_r -= zmniejszO;zmniejszO*=0.5;}
			else if(b>200){ stosunekR_r += zmniejszO;zmniejszO*=2;}
			else if(stosunekDocelowy>stosunekR_r){stosunekR_r=stosunekDocelowy;break;}
			else{break;}

		}
	}
	else{
		while(1){

			stosunekRdorWyszukajElemntyMalejacoSprawdz(stosunekR_r);
			b=wybraneZleEl.getIter();

			if(b<50){       stosunekR_r -= zmniejszO;zmniejszO*=2;}
			else if(b>200){ zmniejszO*=0.5;stosunekR_r += zmniejszO;zmniejszO*=0.5;}
			else if(stosunekDocelowy>stosunekR_r){stosunekR_r=stosunekDocelowy;break;}
			else{break;}

		}
	}

	*/

	a = ustalNowePunkty(wybraneZleEl, PoprawZPozaObszaru, noweP);

	int ileP = points.getIter();

	for (int i = 0, ileE = noweP.getIter(); i < ileE; ++i) {
		ktoreP.setElement(ileP + i);
	}


	delanoyTablicaOPDodajPunktyZTablicaPolozeniaElUstalJakoscEl(tabelementow, dx, dy, dz, false, noweP, 0.3, 3);

	switch (model) {

		case 1:
			optymalizacjaStosunekRdorWybranePunkty(0.1, ktoreP, 4);
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);

			break;
		case 2:
			sasiedniePunkty();
			wygladzanieLaplaceWybranePunkty(ktoreP);
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);
			break;
		case 3:
			sasiedniePunkty();
			for (int i = 0; i < 4; ++i) {
				wygladzanieLaplaceWagaWybranePunkty(2, 0.3, ktoreP, 1);
			}
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);
			break;

	}

	return a;
}

int Delanouy3D::poprawaElementPrzysp(int ** tabelementow, int model, double & stosunekR_r, bool PoprawZPozaObszaru) {

	PunktList noweP(1000, 1000);
	IntList ktoreP(1000, 1000);
	IntList ktoreE(4000, 4000);
	/////////////////////////////////////////////////

	int a, b;
	//double sx,sy,sz;
	//bool zwieksz;

	stosunekRdorWyszukajElemntyMalejacoSprawdz(stosunekR_r);
	b = wybraneZleEl.getIter();

	a = ustalNowePunkty(wybraneZleEl, PoprawZPozaObszaru, noweP);

	int ileP = points.getIter();

	for (int i = 0, ileE = noweP.getIter(); i < ileE; ++i) {
		ktoreP.setElement(ileP + i);
	}


	delanoyTablicaOPDodajPunktyZTablicaPolozeniaElUstalJakoscEl(tabelementow, dx, dy, dz, false, noweP, 0.3, 3);

	switch (model) {

		case 1:
			optymalizacjaStosunekRdorWybranePunkty(0.1, ktoreP, 4);
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);

			break;
		case 2:
			sasiedniePunkty();
			wygladzanieLaplaceWybranePunkty(ktoreP);
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);
			break;
		case 3:
			sasiedniePunkty();
			for (int i = 0; i < 4; ++i) {
				wygladzanieLaplaceWagaWybranePunkty(2, 0.3, ktoreP, 1);
			}
			znajdzElemetyPoDodaniuPunktow(ktoreP, ktoreE);
			uzupelnienieElementowWybranych(ktoreE);
			break;

	}

	return a;
}


void Delanouy3D::poprawaElementWszystkieRdorPoprawaMiejscowaMalejaco(int model, double stosunekR_r, bool PoprawZPozaObszaru) {


	int max = 0, a = 1, b = -1;
	double zmniejszanieO = 1;
	double obecnyStosunekR_r = 30;
	//double obecnyStosunekR_r=stosunekR_r;
	int ** tabElementow2;
	int wymiarZusunieciaTab = przygotojTabliceDoDelanoya(tabElementow2, 3, 0.3);
	double startR_r = 30;

	stosunekRdorWyszukajElemntyMalejacoWypelnij();



	/*
	while(max<500){

		//if(stosunekR_r>obecnyStosunekR_r){break;}
		if(a==0){break;}
		//if(model==1){reMes();}
		a=poprawaElementMalejaco(tabElementow2,model,obecnyStosunekR_r,stosunekR_r,zmniejszanieO,PoprawZPozaObszaru);
		++max;

		m1->Lines->Add(obecnyStosunekR_r);
		m1->Lines->Add(zmniejszanieO);
		m1->Lines->Add(a);
	}
	*/
	elements.getIter();


	while (stosunekR_r < obecnyStosunekR_r) {

		a = 1;
		max = 0;
		obecnyStosunekR_r = startR_r;

		while (max < 150) {

			if (a == 0) {
				break;
			}
			if (model == 1) {
				reMes();
			}
			a = poprawaElementMalejaco(tabElementow2, model, obecnyStosunekR_r, stosunekR_r, zmniejszanieO, PoprawZPozaObszaru);
			if (max == 0) {
				b = a;
			}
			++max;

		}

		//if(max < 149){break;}
		if (b == 0) {
			zmniejszanieO = 0.25;
		}
		if (b < 50) {
			zmniejszanieO = zmniejszanieO * 2;
		}
		else if (b > 200) {
			zmniejszanieO = zmniejszanieO / 2;
		}
		else if (b > 400) {
			zmniejszanieO = zmniejszanieO / 4;
		}

		startR_r -= zmniejszanieO;
		//m1->Lines->Add(startR_r);
		//m1->Lines->Add(zmniejszanieO);
	}





	//m1->Lines->Add("max");
	//m1->Lines->Add(max);

	for (int i = 0; i < wymiarZusunieciaTab; i++) {
		delete []tabElementow2[i];
	}
	delete []tabElementow2;
}


void Delanouy3D::poprawaElementWszystkieRdorPoprawaMiejscowaPrzysp(int model, double stosunekR_r, bool PoprawZPozaObszaru) {


	int max = 0, a = 5;
	int ** tabElementow2;
	int wymiarZusunieciaTab = przygotojTabliceDoDelanoya(tabElementow2, 3, 0.3);

	stosunekRdorWyszukajElemntyMalejacoWypelnij();

	while (max < 500) {

		if (a == 0) {
			break;
		}
		a = poprawaElementPrzysp(tabElementow2, model, stosunekR_r, PoprawZPozaObszaru);
		++max;

	}

	for (int i = 0; i < wymiarZusunieciaTab; i++) {
		delete []tabElementow2[i];
	}
	delete []tabElementow2;
}

void Delanouy3D::poprawaElementWszystkieRdorPoprawaMiejscowa_popZiarna(int model, double stosunekR_r, bool PoprawZPozaObszaru) {

	int startowyS = 20;
	IntList wzorG(10000, 10000);
//    int a=1;
	int max = 0;
	int iterP = points.getIter();
	PunktList nP(100, 100);

	dopasujElementDoZiarna(true);
	zapiszWzorGranicy(wzorG);

	for (; startowyS >= stosunekR_r; --startowyS) {
		int a = 1;
		while (max < 50) {

			if (a == 0) {
				break;
			}

			a = poprawaElement_popZiarna(model, startowyS, PoprawZPozaObszaru);
			++max;
			if (model == 1) {
				reMes();
			}

			poprawGranice(wzorG, nP);
			if (nP.getIter()) {

				delanoyTablicaOPDodajPunkty(dx, dy, dz, false, nP, 0.3, 3);

				dopasujElementDoZiarna(true);
				zapiszWzorGranicy(wzorG);
				iterP = points.getIter();
				a = 1;

			}
		}

		/*
		m1->Lines->Add("max");
		m1->Lines->Add(max);
		m1->Lines->Add(startowyS);
		*/
	}

}

/*
void Delanouy3D::poprawaElementWszystkieRdorPoprawaMiejscowa_popZiarna(int wymiar,int model,double stosunekR_r,bool PoprawZPozaObszaru){

	dopasujElementDoZiarna(true);
	zapiszWzorGranicy();
	int a=1;
	int max=0;
	PunktList nP(100,100);

	while(max<50){

		if(a==0){break;}

		a=poprawaElement_popZiarna(model,wymiar,stosunekR_r,PoprawZPozaObszaru);
		++max;
		if(model==1){reMes();}

		poprawGranice(nP);
		if(nP.getIter()){
			delanoyTablicaOPDodajPunkty(wymiar,wymiar,wymiar,false,nP,0.3,3);
			dopasujElementDoZiarna(true);
			zapiszWzorGranicy();
			a=1;
		}


	}

	m1->Lines->Add("max");
	m1->Lines->Add(max);

}
*/

void Delanouy3D::poprawaElementWszystkieRdor(int sposobZbieraniaWag, double stosunekR_r, bool PoprawZPozaObszaru, double waga, bool modelPoprawiania, int ileWag) {

	int a = 0;
	int b = 0;
	int krotnosc = 10;
	int max = 0;

	reMes();


	while (krotnosc && max < 50) {

		stosunekRdorWyszukajElemnty(stosunekR_r);
		//stosunekGamaKGMSHWyszukajElemnty(stosunekR_r);

		a = ustalNowePunkty(wybraneZleEl, PoprawZPozaObszaru);

		if (!a) {
			break;
		}

		else if (a == b) {
			--krotnosc;
		}
		else if (a != b) {
			krotnosc = 10;
		}
		b = a;
		++max;


		if (modelPoprawiania) {

			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWagaWybraneP(sposobZbieraniaWag, waga, points.getIter() - a, ileWag);
			}
			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWagaWybraneP(sposobZbieraniaWag, waga, points.getIter() - a, ileWag);
			}
			reMes();
		}

		else {

			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWaga(sposobZbieraniaWag, waga, ileWag);
			}
			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWaga(sposobZbieraniaWag, waga, ileWag);
			}
			reMes();

		}


	}

//m1->Lines->Add("max");
//m1->Lines->Add(max);

}

void Delanouy3D::poprawaElementWszystkieGMSH(int sposobZbieraniaWag, double stosunekR_r, bool PoprawZPozaObszaru, double waga, bool modelPoprawiania, int ileWag) {

	int a = 0;
	int b = 0;
	int krotnosc = 10;
	int max = 0;

	reMes();


	while (krotnosc && max < 50) {

		//stosunekRdorWyszukajElemnty(stosunekR_r);
		stosunekGamaKGMSHWyszukajElemnty(stosunekR_r);

		a = ustalNowePunkty(wybraneZleEl, PoprawZPozaObszaru);

		if (!a) {
			break;
		}

		else if (a == b) {
			--krotnosc;
		}
		else if (a != b) {
			krotnosc = 10;
		}
		b = a;
		++max;

		if (modelPoprawiania) {

			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWagaWybraneP(sposobZbieraniaWag, waga, points.getIter() - a, ileWag);
			}
			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWagaWybraneP(sposobZbieraniaWag, waga, points.getIter() - a, ileWag);
			}
			reMes();
		}

		else {

			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWaga(sposobZbieraniaWag, waga, ileWag);
			}
			reMes();
			sasiedniePunkty();
			for (int ff = 0; ff < 20; ++ff) {
				wygladzanieLaplaceWaga(sposobZbieraniaWag, waga, ileWag);
			}
			reMes();

		}


	}

//m1->Lines->Add("max");
//m1->Lines->Add(max);

}



void Delanouy3D::stosunekGamaKGMSHWyszukajElemnty(double stosunekR_r) {

	wybraneZleEl.ustawIter(0);
//zleElementy.ustawIter(0);
	double maxEdge;
	double p, a, b, c, d, e, f, V, r, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
//double maxX = dx;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {


		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();


		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();

		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;

		p = 0.5 * (a + d + e);
		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		maxEdge = a;
		if (maxEdge < b) {
			maxEdge = b;
		}
		if (maxEdge < c) {
			maxEdge = c;
		}
		if (maxEdge < d) {
			maxEdge = d;
		}
		if (maxEdge < e) {
			maxEdge = e;
		}
		if (maxEdge < f) {
			maxEdge = f;
		}

		r = (V * 14.696936) / (pol * maxEdge);


		if (r < stosunekR_r) {

			wybraneZleEl.setElement(i);

		}

	}

}

void Delanouy3D::stosunekRdorWyszukajElemnty(double stosunekR_r) {

	wybraneZleEl.ustawIter(0);
//zleElementy.ustawIter(0);

	double p, a, b, c, d, e, f, V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
//double maxX = dx;
	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {


		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		R = sqrt(elements.getElement(i).getR());

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();

		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;

		p = 0.5 * (a + d + e);
		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		r = (V * 3) / pol;
		R = R / r;

		if (R > stosunekR_r) {

			wybraneZleEl.setElement(i);

		}

	}

}


void Delanouy3D::stosunekRdorWyszukajElemntyMalejacoWypelnij() {

//zleElementy.ustawIter(0);

	double p, a, b, c, d, e, f, V, r, R, pol;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;


	for (int i = 0, p1, p2, p3, p4, ileE = elements.getIter(); i < ileE; ++i) {


		p1 = elements.getElement(i).getP1();
		p2 = elements.getElement(i).getP2();
		p3 = elements.getElement(i).getP3();
		p4 = elements.getElement(i).getP4();
		R = sqrt(elements.getElement(i).getR());

		x1 = points.getElement(p1).getX();
		x2 = points.getElement(p2).getX();
		x3 = points.getElement(p3).getX();
		x4 = points.getElement(p4).getX();
		y1 = points.getElement(p1).getY();
		y2 = points.getElement(p2).getY();
		y3 = points.getElement(p3).getY();
		y4 = points.getElement(p4).getY();
		z1 = points.getElement(p1).getZ();
		z2 = points.getElement(p2).getZ();
		z3 = points.getElement(p3).getZ();
		z4 = points.getElement(p4).getZ();

		a = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		b = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
		c = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3));
		d = sqrt((x1 - x4) * (x1 - x4) + (y1 - y4) * (y1 - y4) + (z1 - z4) * (z1 - z4));
		e = sqrt((x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2) + (z4 - z2) * (z4 - z2));
		f = sqrt((x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3));

		pol = 0;

		p = 0.5 * (a + d + e);
		pol += sqrt(p * (p - a) * (p - d) * (p - e));
		p = 0.5 * (b + f + e);
		pol += sqrt(p * (p - b) * (p - f) * (p - e));
		p = 0.5 * (c + d + f);
		pol += sqrt(p * (p - c) * (p - d) * (p - f));
		p = 0.5 * (a + b + c);
		pol += sqrt(p * (p - a) * (p - b) * (p - c));

		V = ((x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
						(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1)) / 6;

		r = (V * 3) / pol;
		R = R / r;

		elements.getElement(i).setStosunekR_r(R);

	}

}



void Delanouy3D::stosunekRdorWyszukajElemntyMalejacoSprawdz(double stosunekR_r) {

	wybraneZleEl.ustawIter(0);
	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		if (elements.getElement(i).getStosunekR_r() > stosunekR_r) {

			wybraneZleEl.setElement(i);

		}

	}

}

int Delanouy3D::ustalNowePunkty(IntList & zleElementy, bool PoprawZPozaObszaru) {

	int ileDodano = 0;
	int ileEe = elements.getIter();
	bool * zlyEFlaga = new bool[ileEe];
	for (int i = 0; i < ileEe; ++i) {
		zlyEFlaga[i] = false;
	}
	IntList delElementy(300, 200);

	Punkt temp;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	for (int i = 0, ileE = zleElementy.getIter(), ktory; i < ileE; ++i) {
		ktory = zleElementy.getElement(i);
		bool flaga = true;
		flaga = elementyDoDelnoyaPoprawaE(ktory, delElementy);


		for (int j = 0, ileEE = delElementy.getIter(); j < ileEE; ++j) {
			if (zlyEFlaga[delElementy.getElement(j)]) {
				flaga = false;
				break;
			}
		}


		if (flaga) {

			for (int j = 0, ileEE = delElementy.getIter(); j < ileEE; ++j) {
				zlyEFlaga[delElementy.getElement(j)] = true;
			}

			//najwazniejsze
			temp.setPunkt(ustalDodatkowyPunkt(delElementy));

			if (0 <= temp.getX() && temp.getX() <= maxX && 0 <= temp.getY() && temp.getY() <= maxY && 0 <= temp.getZ() && temp.getZ() <= maxZ) {
				points.setElement(temp);
				ileDodano++;
			}
			else {

				if (PoprawZPozaObszaru) {
					if (0 > temp.getX()) {
						temp.setX(0);
					}
					else if (maxX < temp.getX()) {
						temp.setX(maxX);
					}

					if (0 > temp.getY()) {
						temp.setY(0);
					}
					else if (maxY < temp.getY()) {
						temp.setY(maxY);
					}

					if (0 > temp.getZ()) {
						temp.setZ(0);
					}
					else if (maxZ < temp.getZ()) {
						temp.setZ(maxZ);
					}
					points.setElement(temp);
					ileDodano++;
				}

			}

		}

	}

	delete []zlyEFlaga;

	return ileDodano;
}

int Delanouy3D::ustalNowePunkty(IntList & zleElementy, bool PoprawZPozaObszaru, PunktList & dodaneP) {

	dodaneP.czysc(100, zleElementy.getIter() + 10);
	int ileDodano = 0;
	int ileEe = elements.getIter();
	bool * zlyEFlaga = new bool[ileEe];
	for (int i = 0; i < ileEe; ++i) {
		zlyEFlaga[i] = false;
	}
	IntList delElementy(300, 200);

	Punkt temp;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	for (int i = 0, ileE = zleElementy.getIter(), ktory; i < ileE; ++i) {
		ktory = zleElementy.getElement(i);
		bool flaga = true;
		flaga = elementyDoDelnoyaPoprawaE(ktory, delElementy);


		for (int j = 0, ileEE = delElementy.getIter(); j < ileEE; ++j) {
			if (zlyEFlaga[delElementy.getElement(j)]) {
				flaga = false;
				break;
			}
		}


		if (flaga) {

			for (int j = 0, ileEE = delElementy.getIter(); j < ileEE; ++j) {
				zlyEFlaga[delElementy.getElement(j)] = true;
			}

			//najwazniejsze
			temp.setPunkt(ustalDodatkowyPunkt(delElementy));

			if (0 <= temp.getX() && temp.getX() <= maxX && 0 <= temp.getY() && temp.getY() <= maxY && 0 <= temp.getZ() && temp.getZ() <= maxZ) {
				dodaneP.setElement(temp);
				ileDodano++;
			}
			else {

				if (PoprawZPozaObszaru) {
					if (0 > temp.getX()) {
						temp.setX(0);
					}
					else if (maxX < temp.getX()) {
						temp.setX(maxX);
					}

					if (0 > temp.getY()) {
						temp.setY(0);
					}
					else if (maxY < temp.getY()) {
						temp.setY(maxY);
					}

					if (0 > temp.getZ()) {
						temp.setZ(0);
					}
					else if (maxZ < temp.getZ()) {
						temp.setZ(maxZ);
					}
					dodaneP.setElement(temp);
					ileDodano++;
				}

			}

		}

	}

	delete []zlyEFlaga;

	return ileDodano;
}

Punkt Delanouy3D::ustalDodatkowyPunkt(IntList & delElementy) {

	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, tx, ty, tz, R;
	double maxX = dx;
	double maxY = dy;
	double maxZ = dz;

	int p1, p2, p3, p4;
	Punkt temp, dodaj, dodajS(-1, -1, -1);
	int wybrany = delElementy.getElement(0);
	dodaj.setPunkt(elements.getElement(wybrany).getSx(), elements.getElement(wybrany).getSy(), elements.getElement(wybrany).getSz(), 'a');

	int mm1, mm2, mm3, mm4, ee1, ee2, ee3, ee4;

	for (int j = 1, ktory, ileEE = delElementy.getIter(); j < ileEE; ++j) {
		ktory = delElementy.getElement(j);

		mm1 = 0;
		mm2 = 0;
		mm3 = 0;
		mm4 = 0;
		ee1 = 0, ee2 = 0, ee3 = 0, ee4 = 0;

		p1 = elements.getElement(ktory).getP1();
		p2 = elements.getElement(ktory).getP2();
		p3 = elements.getElement(ktory).getP3();
		p4 = elements.getElement(ktory).getP4();

		x1 = points.getElement(p1).getX();
		y1 = points.getElement(p1).getY();
		z1 = points.getElement(p1).getZ();

		x2 = points.getElement(p2).getX();
		y2 = points.getElement(p2).getY();
		z2 = points.getElement(p2).getZ();

		x3 = points.getElement(p3).getX();
		y3 = points.getElement(p3).getY();
		z3 = points.getElement(p3).getZ();

		x4 = points.getElement(p4).getX();
		y4 = points.getElement(p4).getY();
		z4 = points.getElement(p4).getZ();

		if (points.getElement(p1).getGranica() == 'g') {
			ee1 = 1;
		}
		if (points.getElement(p2).getGranica() == 'g') {
			ee2 = 1;
		}
		if (points.getElement(p3).getGranica() == 'g') {
			ee3 = 1;
		}
		if (points.getElement(p4).getGranica() == 'g') {
			ee4 = 1;
		}


		if (0 == x1 || x1 == maxX || 0 == y1 || y1 == maxY || 0 == z1 || z1 == maxZ) {
			mm1 = 1;
		}
		if (0 == x2 || x2 == maxX || 0 == y2 || y2 == maxY || 0 == z2 || z2 == maxZ) {
			mm2 = 1;
		}
		if (0 == x3 || x3 == maxX || 0 == y3 || y3 == maxY || 0 == z3 || z3 == maxZ) {
			mm3 = 1;
		}
		if (0 == x4 || x4 == maxX || 0 == y4 || y4 == maxY || 0 == z4 || z4 == maxZ) {
			mm4 = 1;
		}


		if (ee1 + ee2 + ee3 + ee4 > 2) {
			//m1->Lines->Add("Dodano na gr");

			if (ee4 + ee2 + ee1 == 3) {

				R = punktSciana3D(x4, y4, z4, x2, y2, z2, x1, y1, z1, temp);
				if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
					return temp;
				}

			}
			else if (ee4 + ee3 + ee2 == 3) {

				R = punktSciana3D(x4, y4, z4, x3, y3, z3, x2, y2, z2, temp);
				if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
					return temp;
				}

			}
			else if (ee4 + ee1 + ee3 == 3) {

				R = punktSciana3D(x4, y4, z4, x1, y1, z1, x3, y3, z3, temp);
				if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
					return temp;
				}

			}
			else if (ee1 + ee2 + ee3 == 3) {

				R = punktSciana3D(x1, y1, z1, x2, y2, z2, x3, y3, z3, temp);
				if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
					return temp;
				}

			}

		}
		else if (ee1 + ee2 + ee3 + ee4 > 1) {

			if (ee1 + ee2 == 2) {

				tx = (x1 + x2) * 0.5;
				ty = (y1 + y2) * 0.5;
				tz = (z1 + z2) * 0.5;

				if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
					//m1->Lines->Add("kra gr");
					temp.setPunkt(tx, ty, tz, 'g');
					return temp;
				}

			}
			else if (ee1 + ee3 == 2) {

				tx = (x1 + x3) * 0.5;
				ty = (y1 + y3) * 0.5;
				tz = (z1 + z3) * 0.5;

				if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
					//m1->Lines->Add("kra gr");
					temp.setPunkt(tx, ty, tz, 'g');
					return temp;
				}

			}
			else if (ee1 + ee4 == 2) {

				tx = (x1 + x4) * 0.5;
				ty = (y1 + y4) * 0.5;
				tz = (z1 + z4) * 0.5;

				if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
					//m1->Lines->Add("kra gr");
					temp.setPunkt(tx, ty, tz, 'g');
					return temp;
				}

			}
			else if (ee2 + ee3 == 2) {

				tx = (x2 + x3) * 0.5;
				ty = (y2 + y3) * 0.5;
				tz = (z2 + z3) * 0.5;

				if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x2) * (tx - x2) + (ty - y2) * (ty - y2) + (tz - z2) * (tz - z2))) {
					//m1->Lines->Add("kra gr");
					temp.setPunkt(tx, ty, tz, 'g');
					return temp;
				}

			}
			else if (ee2 + ee4 == 2) {

				tx = (x2 + x4) * 0.5;
				ty = (y2 + y4) * 0.5;
				tz = (z2 + z4) * 0.5;

				if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x2) * (tx - x2) + (ty - y2) * (ty - y2) + (tz - z2) * (tz - z2))) {
					//m1->Lines->Add("kra gr");
					temp.setPunkt(tx, ty, tz, 'g');
					return temp;

				}

			}
			else if (ee3 + ee4 == 2) {

				tx = (x3 + x4) * 0.5;
				ty = (y3 + y4) * 0.5;
				tz = (z3 + z4) * 0.5;

				if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x3) * (tx - x3) + (ty - y3) * (ty - y3) + (tz - z3) * (tz - z3))) {
					//m1->Lines->Add("kra gr");
					temp.setPunkt(tx, ty, tz, 'g');
					return temp;

				}

			}

		}



		if (mm1 + mm2 + mm3 + mm4 > 2) {

			if (mm4 + mm2 + mm1 == 3) {

				if (x4 == x2 && x2 == x1 && (x1 == 0 || x1 == maxX)) {
					R = punktSciana2D(0, 0, x1, y4, y2, y1, z4, z2, z1, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getX());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}
				else if (y4 == y2 && y2 == y1 && (y1 == 0 || y1 == maxY)) {
					R = punktSciana2D(0, 1, y1, x4, x2, x1, z4, z2, z1, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getY());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}
				else if (z4 == z2 && z2 == z1 && (z1 == 0 || z1 == maxZ)) {
					R = punktSciana2D(0, 2, z1, x4, x2, x1, y4, y2, y1, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getZ());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}

			}
			else if (mm4 + mm3 + mm2 == 3) {

				if (x4 == x2 && x2 == x3 && (x2 == 0 || x2 == maxX)) {
					R = punktSciana2D(1, 0, x2, y4, y2, y3, z4, z2, z3, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getX());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}
				else if (y4 == y2 && y2 == y3 && (y2 == 0 || y2 == maxY)) {
					R = punktSciana2D(1, 1, y2, x4, x2, x3, z4, z2, z3, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getY());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}
				else if (z4 == z2 && z2 == z3 && (z2 == 0 || z2 == maxZ)) {
					R = punktSciana2D(1, 2, z2, x4, x2, x3, y4, y2, y3, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getZ());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}

			}
			else if (mm4 + mm1 + mm3 == 3) {

				if (x4 == x3 && x3 == x1 && (x1 == 0 || x1 == maxX)) {
					R = punktSciana2D(2, 0, x1, y4, y3, y1, z4, z3, z1, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getX());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}
				else if (y4 == y3 && y3 == y1 && (y1 == 0 || y1 == maxY)) {
					R = punktSciana2D(2, 1, y1, x4, x3, x1, z4, z3, z1, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getY());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}
				else if (z4 == z3 && z3 == z1 && (z1 == 0 || z1 == maxZ)) {
					R = punktSciana2D(2, 2, z1, x4, x3, x1, y4, y3, y1, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getZ());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}

			}
			else if (mm1 + mm2 + mm3 == 3) {

				if (x3 == x2 && x2 == x1 && (x1 == 0 || x1 == maxX)) {
					R = punktSciana2D(3, 0, x1, y1, y2, y3, z1, z2, z3, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getX());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}
				}
				else if (y3 == y2 && y2 == y1 && (y1 == 0 || y1 == maxY)) {
					R = punktSciana2D(3, 1, y1, x1, x2, x3, z1, z2, z3, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getY());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}
				}
				else if (z3 == z2 && z2 == z1 && (z1 == 0 || z1 == maxZ)) {
					R = punktSciana2D(3, 2, z1, x1, x2, x3, y1, y2, y3, temp);
					//if(R!=-1){dodajS.setPunkt(temp.getX(),temp.getY(),temp.getZ());m1->Lines->Add(temp.getZ());}
					if (R != -1 && dodaj.odlegloscPunktu(temp) <= R) {
						dodajS.setPunkt(temp.getX(), temp.getY(), temp.getZ());
					}

				}

			}

		}
		if (mm1 + mm2 + mm3 + mm4 > 1) {

			if (mm1 + mm2 == 2) {
				if (x1 == x2 && (x1 == 0 || x1 == maxX)) {
					if ( (y1 == y2 && (y1 == 0 || y1 == maxY)) || (z1 == z2 && (z1 == 0 || z1 == maxZ)) ) {
						tx = (x1 + x2) * 0.5;
						ty = (y1 + y2) * 0.5;
						tz = (z1 + z2) * 0.5;

						if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
							//m1->Lines->Add("kra");
							temp.setPunkt(tx, ty, tz, 'a');
							return temp;

						}
					}
				}
				else if (y1 == y2 && z1 == z2 && (y1 == 0 || y1 == maxY) && (z1 == 0 || z1 == maxZ)) {
					tx = (x1 + x2) * 0.5;
					ty = (y1 + y2) * 0.5;
					tz = (z1 + z2) * 0.5;

					if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
						//m1->Lines->Add("kra");
						temp.setPunkt(tx, ty, tz, 'a');
						return temp;

					}
				}
			}
			else if (mm1 + mm3 == 2) {
				if (x1 == x3 && (x1 == 0 || x1 == maxX)) {
					if ( (y1 == y3 && (y1 == 0 || y1 == maxY)) || (z1 == z3 && (z1 == 0 || z1 == maxZ)) ) {
						tx = (x1 + x3) * 0.5;
						ty = (y1 + y3) * 0.5;
						tz = (z1 + z3) * 0.5;

						if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
							//m1->Lines->Add("kra");
							temp.setPunkt(tx, ty, tz, 'a');
							return temp;

						}
					}
				}
				else if (y1 == y3 && (y1 == 0 || y1 == maxY) && z1 == z3 && (z1 == 0 || z1 == maxZ)) {
					tx = (x1 + x3) * 0.5;
					ty = (y1 + y3) * 0.5;
					tz = (z1 + z3) * 0.5;

					if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
						//m1->Lines->Add("kra");
						temp.setPunkt(tx, ty, tz, 'a');
						return temp;

					}
				}
			}
			else if (mm1 + mm4 == 2) {
				if (x1 == x4 && (x1 == 0 || x1 == maxX)) {
					if ( (y1 == y4 && (y1 == 0 || y1 == maxY)) || (z1 == z4 && (z1 == 0 || z1 == maxZ)) ) {
						tx = (x1 + x4) * 0.5;
						ty = (y1 + y4) * 0.5;
						tz = (z1 + z4) * 0.5;

						if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
							//m1->Lines->Add("kra");
							temp.setPunkt(tx, ty, tz, 'a');
							return temp;

						}
					}
				}
				else if (y1 == y4 && (y1 == 0 || y1 == maxY) && z1 == z4 && (z1 == 0 || z1 == maxZ)) {
					tx = (x1 + x4) * 0.5;
					ty = (y1 + y4) * 0.5;
					tz = (z1 + z4) * 0.5;

					if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x1) * (tx - x1) + (ty - y1) * (ty - y1) + (tz - z1) * (tz - z1))) {
						//m1->Lines->Add("kra");
						temp.setPunkt(tx, ty, tz, 'a');
						return temp;

					}
				}
			}
			else if (mm2 + mm3 == 2) {
				if (x2 == x3 && (x2 == 0 || x2 == maxX)) {
					if ( (y2 == y3 && (y2 == 0 || y2 == maxY)) || (z2 == z3 && (z2 == 0 || z2 == maxZ)) ) {
						tx = (x2 + x3) * 0.5;
						ty = (y2 + y3) * 0.5;
						tz = (z2 + z3) * 0.5;

						if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x2) * (tx - x2) + (ty - y2) * (ty - y2) + (tz - z2) * (tz - z2))) {
							//m1->Lines->Add("kra");
							temp.setPunkt(tx, ty, tz, 'a');
							return temp;

						}
					}
				}
				else if (y2 == y3 && (y2 == 0 || y2 == maxY) && z2 == z3 && (z2 == 0 || z2 == maxZ)) {
					tx = (x2 + x3) * 0.5;
					ty = (y2 + y3) * 0.5;
					tz = (z2 + z3) * 0.5;

					if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x2) * (tx - x2) + (ty - y2) * (ty - y2) + (tz - z2) * (tz - z2))) {
						//m1->Lines->Add("kra");
						temp.setPunkt(tx, ty, tz, 'a');
						return temp;

					}
				}
			}
			else if (mm2 + mm4 == 2) {
				if (x2 == x4 && (x2 == 0 || x2 == maxX)) {
					if ( (y2 == y4 && (y2 == 0 || y2 == maxY)) || (z2 == z4 && (z2 == 0 || z2 == maxZ)) ) {
						tx = (x2 + x4) * 0.5;
						ty = (y2 + y4) * 0.5;
						tz = (z2 + z4) * 0.5;

						if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x2) * (tx - x2) + (ty - y2) * (ty - y2) + (tz - z2) * (tz - z2))) {
							//m1->Lines->Add("kra");
							temp.setPunkt(tx, ty, tz, 'a');
							return temp;

						}
					}
				}
				else if (y2 == y4 && (y2 == 0 || y2 == maxY) && z2 == z4 && (z2 == 0 || z2 == maxZ)) {
					tx = (x2 + x4) * 0.5;
					ty = (y2 + y4) * 0.5;
					tz = (z2 + z4) * 0.5;

					if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x2) * (tx - x2) + (ty - y2) * (ty - y2) + (tz - z2) * (tz - z2))) {
						//m1->Lines->Add("kra");
						temp.setPunkt(tx, ty, tz, 'a');
						return temp;

					}
				}
			}
			else if (mm3 + mm4 == 2) {
				if (x3 == x4 && (x3 == 0 || x3 == maxX)) {
					if ( (y3 == y4 && (y3 == 0 || y3 == maxY)) || (z3 == z4 && (z3 == 0 || z3 == maxZ)) ) {
						tx = (x3 + x4) * 0.5;
						ty = (y3 + y4) * 0.5;
						tz = (z3 + z4) * 0.5;

						if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x3) * (tx - x3) + (ty - y3) * (ty - y3) + (tz - z3) * (tz - z3))) {
							//m1->Lines->Add("kra");
							temp.setPunkt(tx, ty, tz, 'a');
							return temp;

						}
					}
				}
				else if (y3 == y4 && (y3 == 0 || y3 == maxY) && z3 == z4 && (z3 == 0 || z3 == maxX)) {
					tx = (x3 + x4) * 0.5;
					ty = (y3 + y4) * 0.5;
					tz = (z3 + z4) * 0.5;

					if (dodaj.odlegloscPunktu(tx, ty, tz) <= float((tx - x3) * (tx - x3) + (ty - y3) * (ty - y3) + (tz - z3) * (tz - z3))) {
						//m1->Lines->Add("kra");
						temp.setPunkt(tx, ty, tz, 'a');
						return temp;

					}
				}
			}


		}



	}

	if (dodajS.getX() != -1) {
		return dodajS;
	}
	return dodaj;
}


float Delanouy3D::punktSciana3D(const double & x1, const double & y1, const double & z1, const double & x2, const double & y2, const double & z2, const double & x3, const double & y3, const double & z3, Punkt & temp) {

	double a2, a3, a4, b2, b3, b4, c2, c3, c4, l2, l3, l4, W, X, Y, Z;

	a2 = (-x3 + x1) * 2;
	b2 = (-y3 + y1) * 2;
	c2 = (-z3 + z1) * 2;
	l2 = -(x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + z3 * z3 - z1 * z1);

	a3 = (-x2 + x1) * 2;
	b3 = (-y2 + y1) * 2;
	c3 = (-z2 + z1) * 2;
	l3 = -(x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + z2 * z2 - z1 * z1);

	a4 = y1 * z2 + y2 * z3 + y3 * z1 - z2 * y3 - z3 * y1 - z1 * y2;
	b4 = -(x1 * z2 + x2 * z3 + x3 * z1 - z2 * x3 - z3 * x1 - z1 * x2);
	c4 = x1 * y2 + x2 * y3 + x3 * y1 - y2 * x3 - y3 * x1 - y1 * x2;
	l4 = (x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - z1 * y2 * x3 - z2 * y3 * x1 - z3 * y1 * x2);

	W =  a2 * b3 * c4 + a3 * b4 * c2 + a4 * b2 * c3 - c2 * b3 * a4 - c3 * b4 * a2 - c4 * b2 * a3;


	if (W) {

		W = 1 / W;
		l2 *= W;
		l3 *= W;
		l4 *= W;

		X = (b2 * c3 * l4 + b3 * c4 * l2 + b4 * c2 * l3 - l2 * c3 * b4 - l3 * c4 * b2 - l4 * c2 * b3);
		Y = (a2 * c3 * l4 + a3 * c4 * l2 + a4 * c2 * l3 - l2 * c3 * a4 - l3 * c4 * a2 - l4 * c2 * a3);
		Y = Y * -1;
		Z = (a2 * b3 * l4 + a3 * b4 * l2 + a4 * b2 * l3 - l2 * b3 * a4 - l3 * b4 * a2 - l4 * b2 * a3);

	}
	else {

		X = (x1 + x2 + x3) / 3;
		Y = (y1 + y2 + y3) / 3;
		Z = (z1 + z2 + z3) / 3;

	}

	temp.setXYZ(X, Y, Z);
	temp.setGranica('g');

	return (X - x1) * (X - x1) + (Y - y1) * (Y - y1) + (Z - z1) * (Z - z1);

}

float Delanouy3D::punktSciana2D(int r, int rodzaj, const double & War0, const double & x1, const double & x2, const double & x3, const double & y1, const double & y2, const double & y3, Punkt & temp) {


	double a1, b1, c1, a2, b2, c2, X, Y;

	a1 = x1 - x2;
	b1 = y1 - y2;
	c1 = -(a1 * ((x1 + x2) * 0.5) + b1 * ((y1 + y2) * 0.5));

	a2 = x2 - x3;
	b2 = y2 - y3;
	c2 = -(a2 * ((x2 + x3) * 0.5) + b2 * ((y2 + y3) * 0.5));

	if (a1 * b2 != a2 * b1) {

		Y = (-a1 * c2 + a2 * c1) / (a1 * b2 - a2 * b1);
		X = (-c1 * b2 + c2 * b1) / (a1 * b2 - a2 * b1);



		switch (rodzaj) {

			case 0:
				temp.setPunkt(War0, X, Y, 'a');
				break;
			case 1:
				temp.setPunkt(X, War0, Y, 'a');
				break;
			case 2:
				temp.setPunkt(X, Y, War0, 'a');
				break;

		}

		return (X - x1) * (X - x1) + (Y - y1) * (Y - y1);
	}

	/*
	AnsiString a;

	m1->Lines->Add("error");
	a=x1;a+="-";a+=x2;a+="-";a+=x3;
	m1->Lines->Add(a);
	a=y1;a+="-";a+=y2;a+="-";a+=y3;
	m1->Lines->Add(a);
	a="TO!!!  :  ";a+=r;
	Lines->Add(a);
	*/

	return -1;
}

void Delanouy3D::zapiszNajVTrojkat() {

	double max = 0;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, wynik;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		x1 = points.getElement(elements.getElement(i).getP1()).getX();
		x2 = points.getElement(elements.getElement(i).getP2()).getX();
		x3 = points.getElement(elements.getElement(i).getP3()).getX();
		x4 = points.getElement(elements.getElement(i).getP4()).getX();

		y1 = points.getElement(elements.getElement(i).getP1()).getY();
		y2 = points.getElement(elements.getElement(i).getP2()).getY();
		y3 = points.getElement(elements.getElement(i).getP3()).getY();
		y4 = points.getElement(elements.getElement(i).getP4()).getY();

		z1 = points.getElement(elements.getElement(i).getP1()).getZ();
		z2 = points.getElement(elements.getElement(i).getP2()).getZ();
		z3 = points.getElement(elements.getElement(i).getP3()).getZ();
		z4 = points.getElement(elements.getElement(i).getP4()).getZ();



		wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
				(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1) ;

		if (max < wynik / 6) {
			max = wynik / 6;
		}

	}

	NajT = max;
}


void Delanouy3D::wyszukajDuzeTrojkaty(const double & stosunekR_r, IntList & zleElementy) {

	double naj = NajT * stosunekR_r;
	int ileZT = 0;
//double maxX=dx,wynikAll=0;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, wynik;

	for (int i = 0, ileE = elements.getIter(); i < ileE; i++) {

		x1 = points.getElement(elements.getElement(i).getP1()).getX();
		x2 = points.getElement(elements.getElement(i).getP2()).getX();
		x3 = points.getElement(elements.getElement(i).getP3()).getX();
		x4 = points.getElement(elements.getElement(i).getP4()).getX();

		y1 = points.getElement(elements.getElement(i).getP1()).getY();
		y2 = points.getElement(elements.getElement(i).getP2()).getY();
		y3 = points.getElement(elements.getElement(i).getP3()).getY();
		y4 = points.getElement(elements.getElement(i).getP4()).getY();

		z1 = points.getElement(elements.getElement(i).getP1()).getZ();
		z2 = points.getElement(elements.getElement(i).getP2()).getZ();
		z3 = points.getElement(elements.getElement(i).getP3()).getZ();
		z4 = points.getElement(elements.getElement(i).getP4()).getZ();

		wynik = (x2 - x1) * (y3 - y1) * (z4 - z1) + (x3 - x1) * (y4 - y1) * (z2 - z1) + (x4 - x1) * (y2 - y1) * (z3 - z1) -
				(z2 - z1) * (y3 - y1) * (x4 - x1) - (z3 - z1) * (y4 - y1) * (x2 - x1) - (z4 - z1) * (y2 - y1) * (x3 - x1);
		if (wynik < 0) {
			wynik *= -1;
		}
		if (naj < wynik / 6) {
			zleElementy.setElement(i);
			++ileZT;
		}
	}

//m1->Lines->Add(ileZT);

}


void Delanouy3D::rozbijajDuzeTrojkaty(const double & stosunekR_r, bool PoprawZPozaObszaru) {

	IntList zleElementy(2000, 2000);

	wyszukajDuzeTrojkaty(stosunekR_r, zleElementy);
	ustalNowePunkty(zleElementy, PoprawZPozaObszaru);

}
/*
void Delanouy3D::rysujZleEl(){

double maxX = dx;
double maxY = dy;
double maxZ = dz;

glPushMatrix();


glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);

	glBegin(GL_POINTS);
	glColor3f(0,0.5,1);

	for(int i,j=0,ileE=wybraneZleEl.getIter();j<ileE;++j){

		i=wybraneZleEl.getElement(j);
		glVertex3f(elements.getElement(i).getSx(),elements.getElement(i).getSy(),elements.getElement(i).getSz());

	}

	glEnd();


	glBegin(GL_LINES);
	glColor3f(1,0.5,0);


	for(int i,j=0,ileE=wybraneZleEl.getIter();j<ileE;++j){

		i=wybraneZleEl.getElement(j);

		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());


		glVertex3f(points.getElement(elements.getElement(i).getP1()).getX(),points.getElement(elements.getElement(i).getP1()).getY(), points.getElement(elements.getElement(i).getP1()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP2()).getX(),points.getElement(elements.getElement(i).getP2()).getY(), points.getElement(elements.getElement(i).getP2()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

		glVertex3f(points.getElement(elements.getElement(i).getP3()).getX(),points.getElement(elements.getElement(i).getP3()).getY(), points.getElement(elements.getElement(i).getP3()).getZ());
		glVertex3f(points.getElement(elements.getElement(i).getP4()).getX(),points.getElement(elements.getElement(i).getP4()).getY(), points.getElement(elements.getElement(i).getP4()).getZ());

	}



	glColor3f(1,1,1);
		glVertex3f(points.getElement(0).getX(),points.getElement(0).getY(), points.getElement(0).getZ());
		glVertex3f(points.getElement(1).getX(),points.getElement(1).getY(), points.getElement(1).getZ());
		glVertex3f(points.getElement(0).getX(),points.getElement(0).getY(), points.getElement(0).getZ());
		glVertex3f(points.getElement(4).getX(),points.getElement(4).getY(), points.getElement(4).getZ());
		glVertex3f(points.getElement(0).getX(),points.getElement(0).getY(), points.getElement(0).getZ());
		glVertex3f(points.getElement(2).getX(),points.getElement(2).getY(), points.getElement(2).getZ());

		glVertex3f(points.getElement(1).getX(),points.getElement(1).getY(), points.getElement(1).getZ());
		glVertex3f(points.getElement(3).getX(),points.getElement(3).getY(), points.getElement(3).getZ());
		glVertex3f(points.getElement(1).getX(),points.getElement(1).getY(), points.getElement(1).getZ());
		glVertex3f(points.getElement(5).getX(),points.getElement(5).getY(), points.getElement(5).getZ());

		glVertex3f(points.getElement(2).getX(),points.getElement(2).getY(), points.getElement(2).getZ());
		glVertex3f(points.getElement(3).getX(),points.getElement(3).getY(), points.getElement(3).getZ());
		glVertex3f(points.getElement(2).getX(),points.getElement(2).getY(), points.getElement(2).getZ());
		glVertex3f(points.getElement(6).getX(),points.getElement(6).getY(), points.getElement(6).getZ());

		glVertex3f(points.getElement(3).getX(),points.getElement(3).getY(), points.getElement(3).getZ());
		glVertex3f(points.getElement(7).getX(),points.getElement(7).getY(), points.getElement(7).getZ());

		glVertex3f(points.getElement(4).getX(),points.getElement(4).getY(), points.getElement(4).getZ());
		glVertex3f(points.getElement(6).getX(),points.getElement(6).getY(), points.getElement(6).getZ());
		glVertex3f(points.getElement(4).getX(),points.getElement(4).getY(), points.getElement(4).getZ());
		glVertex3f(points.getElement(5).getX(),points.getElement(5).getY(), points.getElement(5).getZ());

		glVertex3f(points.getElement(5).getX(),points.getElement(5).getY(), points.getElement(5).getZ());
		glVertex3f(points.getElement(7).getX(),points.getElement(7).getY(), points.getElement(7).getZ());

		glVertex3f(points.getElement(6).getX(),points.getElement(6).getY(), points.getElement(6).getZ());
		glVertex3f(points.getElement(7).getX(),points.getElement(7).getY(), points.getElement(7).getZ());


	glEnd();

glPopMatrix();
glPopMatrix();

}
*/

/*
void Delanouy3D::rysujPow(){


IntList elementy(3000,3000);
bool b1,b2,b3,b4;
bool *zazEl = new bool[elements.getIter()];
for(int i=0,ileE=elements.getIter();i<ileE;++i){zazEl[i]=true;}


	for(int i=0,ileE=elements.getIter(),p1,p2,p3,p4;i<ileE;++i){

		p1=elements.getElement(i).getP1();
		p2=elements.getElement(i).getP2();
		p3=elements.getElement(i).getP3();
		p4=elements.getElement(i).getP4();

		if(points.getElement(p1).getGranica()=='g'){b1=true;}else{b1=false;}
		if(points.getElement(p2).getGranica()=='g'){b2=true;}else{b2=false;}
		if(points.getElement(p3).getGranica()=='g'){b3=true;}else{b3=false;}
		if(points.getElement(p4).getGranica()=='g'){b4=true;}else{b4=false;}

		if(b4 && b2 && b1 && zazEl[elements.getElement(i).getE1()]){elementy.setElement(p4);elementy.setElement(p2);elementy.setElement(p1);zazEl[i]=false;}
		if(b4 && b3 && b2 && zazEl[elements.getElement(i).getE2()]){elementy.setElement(p4);elementy.setElement(p3);elementy.setElement(p2);zazEl[i]=false;}
		if(b4 && b1 && b3 && zazEl[elements.getElement(i).getE3()]){elementy.setElement(p4);elementy.setElement(p1);elementy.setElement(p3);zazEl[i]=false;}
		if(b1 && b2 && b3 && zazEl[elements.getElement(i).getE4()]){elementy.setElement(p1);elementy.setElement(p2);elementy.setElement(p3);zazEl[i]=false;}

	}

delete []zazEl;

int *pozPunktow = new int[points.getIter()];

	int iP=0;
	for(int i=0,ileP=points.getIter();i<ileP;++i){
		if(points.getElement(i).getGranica()=='g'){
			pozPunktow[i]=iP;++iP;
		}
	}

PunktList punkty(4000,4000);


	for(int i=0,j=0,ileP=points.getIter();i<ileP;++i){
		if(points.getElement(i).getGranica()=='g'){punkty.setElement(points.getElement(i));}
	}


double maxX = dx;
double maxY = dy;
double maxZ = dz;
glPushMatrix();


glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);


	glBegin(GL_LINES);
	glColor3f(1,0.5,0);

	for(int i=0,j=0,p1,p2,p3,ileE=elementy.getIter();i<ileE;i=i+3){

		p1=elementy.getElement(i);
		p2=elementy.getElement(i+1);
		p3=elementy.getElement(i+2);

		glVertex3f(punkty.getElement(pozPunktow[p1]).getX(),punkty.getElement(pozPunktow[p1]).getY(), punkty.getElement(pozPunktow[p1]).getZ());
		glVertex3f(punkty.getElement(pozPunktow[p2]).getX(),punkty.getElement(pozPunktow[p2]).getY(), punkty.getElement(pozPunktow[p2]).getZ());

		glVertex3f(punkty.getElement(pozPunktow[p2]).getX(),punkty.getElement(pozPunktow[p2]).getY(), punkty.getElement(pozPunktow[p2]).getZ());
		glVertex3f(punkty.getElement(pozPunktow[p3]).getX(),punkty.getElement(pozPunktow[p3]).getY(), punkty.getElement(pozPunktow[p3]).getZ());

		glVertex3f(punkty.getElement(pozPunktow[p3]).getX(),punkty.getElement(pozPunktow[p3]).getY(), punkty.getElement(pozPunktow[p3]).getZ());
		glVertex3f(punkty.getElement(pozPunktow[p1]).getX(),punkty.getElement(pozPunktow[p1]).getY(), punkty.getElement(pozPunktow[p1]).getZ());

	}

	glEnd();

glPopMatrix();
glPopMatrix();


delete []pozPunktow;
}

*/
/*
void Delanouy3D::rysujPunkty(int size){

glPointSize(size);
double maxX = dx;
double maxY = dy;
double maxZ = dz;


glPushMatrix();

glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);

glBegin(GL_POINTS);


for(int i=0;i<points.getIter();i++){

	if(points.getElement(i).getGranica()=='g'){glColor3f(1,0,0);}
	else if(points.getElement(i).getGranica()=='r'){glColor3f(1,1,0);}
	else{glColor3f(1,1,1);}

	glVertex3f(points.getElement(i).getX(),points.getElement(i).getY(),points.getElement(i).getZ());

}


glEnd();
glPopMatrix();
glPopMatrix();

}

*/
/*
void Delanouy3D::zapiszBMP(char *nazwa,void *handle){

	HWND BCBHandle = FindWindow("TAppBuilder", NULL);

			HDC dc = GetWindowDC(handle);
			TCanvas *ScreenCanvas = new TCanvas;
			Graphics::TBitmap *BMP = new Graphics::TBitmap;
			ScreenCanvas->Handle = dc;

			TRect rect = ScreenCanvas->ClipRect;
			rect.Right = rect.Right - rect.Left;
			rect.Bottom = rect.Bottom - rect.Top;
			rect.Top = 0;
			rect.Left =0;

			BMP->Width = rect.Right;
			BMP->Height= rect.Bottom;
			BMP->Canvas->CopyRect(rect, ScreenCanvas,ScreenCanvas->ClipRect);

			BMP->SaveToFile(nazwa);

			delete ScreenCanvas;
			delete BMP;


	   //	}

}

*/
void Delanouy3D::wczytaj3D(const char * nazwa) {

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);

	int ilePunktow;
	int ileElementow;

	wczytaj >> ilePunktow;
	wczytaj >> ileElementow;

	points.czysc(ilePunktow * 0.5, ilePunktow + 10);
	elements.czysc(ileElementow * 0.5, ileElementow + 10);

	double x, y, z;
	int p1, p2, p3, p4;

	for (int i = 0; i < ilePunktow; ++i) {

		wczytaj >> z;
		wczytaj >> x;
		wczytaj >> y;
		wczytaj >> z;

		points.setElement(x, y, z);

	}

	for (int i = 0; i < ileElementow; ++i) {

		wczytaj >> p1;
		wczytaj >> p1;
		wczytaj >> p2;
		wczytaj >> p3;
		wczytaj >> p4;

		elements.setElement(p1 - 1, p2 - 1, p3 - 1, p4 - 1);

	}
	wczytaj.close();
}

void Delanouy3D::wczytajPunktyPSS(const char * nazwa, int ktoreSerce) {


	double procent = 90;
	double c;
	double x, y, z;
	points.czysc(7000, 28000);
	elements.czysc(20000, 100000);

	double minX = 99999999, maxX = -99999999, tranX = 0;
	double minY = 99999999, maxY = -99999999, tranY = 0;
	double minZ = 99999999, maxZ = -99999999, tranZ = 0;
	double sX, sY, sZ;

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);
//6795   pk graniczne
//26084
//19021  pun
//99184  ele

//137048 pun
//741871 ele
//27463  pk graniczne

	int ilePunktow;
	int eleMentow;

	switch (ktoreSerce) {

		case 1:
			ilePunktow = 19021 ;
			eleMentow  = 99184;
			break;

		case 2:
			ilePunktow = 137048 ;
			eleMentow  = 741871;
			break;

		default:
			ilePunktow = 19021;
			eleMentow  = 99184;
	}

//int ilePunktow = 137048 ;
//int eleMentow  = 741871;

//wczytanie z pliku
	for (int i = 0; i < ilePunktow; ++i) {

		wczytaj >> c;
		wczytaj >> x;
		wczytaj >> y;
		wczytaj >> z;

		points.setElement(x, y, z, 'g');

		//wczytaj>>c;

		if (x < minX) {
			minX = x;
		}
		if (x > maxX) {
			maxX = x;
		}
		if (y < minY) {
			minY = y;
		}
		if (y > maxY) {
			maxY = y;
		}
		if (z < minZ) {
			minZ = z;
		}
		if (z > maxZ) {
			maxZ = z;
		}

	}

	for (int i = 0, p1, p2, p3, p4; i < eleMentow; ++i) {

		wczytaj >> p1;
		wczytaj >> p1;
		wczytaj >> p2;
		wczytaj >> p3;
		wczytaj >> p4;

		elements.setElement(p1 - 1, p2 - 1, p3 - 1, p4 - 1);

	}



	double skal;

	if ((minX - maxX) * (minX - maxX) > (minY - maxY) * (minY - maxY)) {
		if ((minX - maxX) * (minX - maxX) > (minZ - maxZ) * (minZ - maxZ)) { //x
			skal = procent / sqrt((minX - maxX) * (minX - maxX));
		}
		else { //z
			skal = procent / sqrt((minZ - maxZ) * (minZ - maxZ));
		}
	}
	else {
		if ((minY - maxY) * (minY - maxY) > (minZ - maxZ) * (minZ - maxZ)) { //y
			skal = procent / sqrt((minY - maxY) * (minY - maxY));
		}
		else { //z
			skal = procent / sqrt((minZ - maxZ) * (minZ - maxZ));
		}

	}

	for (int i = 0; i < ilePunktow; ++i) {
		// punkty.getElement(i).setPunkt(punkty.getElement(i).getX()+tranX,punkty.getElement(i).getY()+tranY,punkty.getElement(i).getZ()+tranZ);
		points.getElement(i) = points.getElement(i) * skal;
		//testG.setElement(punkty.getElement(i));
	}

	minX = 99999999;
	maxX = -99999999;
	tranX = 0;
	minY = 99999999;
	maxY = -99999999;
	tranY = 0;
	minZ = 99999999;
	maxZ = -99999999;
	tranZ = 0;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		x = points.getElement(i).getX();
		y = points.getElement(i).getY();
		z = points.getElement(i).getZ();

		if (x < minX) {
			minX = x;
		}
		if (x > maxX) {
			maxX = x;
		}
		if (y < minY) {
			minY = y;
		}
		if (y > maxY) {
			maxY = y;
		}
		if (z < minZ) {
			minZ = z;
		}
		if (z > maxZ) {
			maxZ = z;
		}

	}

	sX = (minX + maxX) * 0.5;
	sY = (minY + maxY) * 0.5;
	sZ = (minZ + maxZ) * 0.5;

	tranX = 50 - sX;
	tranY = 50 - sY;
	tranZ = 50 - sZ;


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		points.getElement(i).setPunkt(points.getElement(i).getX() + tranX, points.getElement(i).getY() + tranY, points.getElement(i).getZ() + tranZ, 'g');
	}
	//punkty.ustawIter(0);


	globalTranX = tranX;
	globalTranY = tranY;
	globalTranZ = tranZ;
	globalSkal = skal;

	//m1->Lines->Add(1/globalSkal);

	wczytaj.close();
}


void Delanouy3D::zwrocLiczby(double & a, double & b, double & c, string & liczba) {

	liczba.erase(0, 1);
	string l1, l2, l3;

	l1 = liczba.substr(0, 8);
	liczba.erase(0, 8);
	l2 = liczba.substr(0, 8);
	liczba.erase(0, 8);
	l3 = liczba.substr(0, 8);
	liczba.erase(0, 8);

	if (l1[5] == '-') {
		l1.insert(5, "e");
	}
	//if(l1[0]==' '){l1='0';}
	if (l2[5] == '-') {
		l2.insert(5, "e");
	}
	//if(l2[0]==' '){l2='0';}
	if (l3[5] == '-') {
		l3.insert(5, "e");
	}
	//if(l3[0]==' '){l3='0';}

	a = atof(l1.c_str());
	b = atof(l2.c_str());
	c = atof(l3.c_str());


}


/*
void Delanouy3D::wczytajPlikNas(const char *nazwa,int mnoznik,bool wczytajWarunki){

points.czysc(100000,100000);
elements.czysc(600000,600000);
IntList warunki(10,10);

if(wczytajWarunki){warunki.czysc(300000,300000);}

ifstream wczytaj(nazwa);
wczytaj.precision(20);


//int i=0;
string textInFile;
string testGRID="GRID";
string testCTETRA="CTETRA";
string testCTRIAX="CTRIAX";

int numer,temp,p1,p2,p3,p4;
double x,y,z;

while(!wczytaj.eof()){

	wczytaj>>textInFile;


	if(textInFile==testGRID){
		wczytaj>>numer;
		wczytaj>>x;wczytaj>>y;wczytaj>>z;
		wczytaj>>temp;

		points.setElement(x*mnoznik,y*mnoznik,z*mnoznik);

	}
	else if(textInFile==testCTETRA){
		wczytaj>>numer;
		wczytaj>>temp;
		wczytaj>>p1;wczytaj>>p2;wczytaj>>p3;wczytaj>>p4;

		elements.setElement(p1-1,p2-1,p3-1,p4-1);

	}
	else if(wczytajWarunki && textInFile==testCTRIAX){
		wczytaj>>numer;
		wczytaj>>temp;
		wczytaj>>p1;wczytaj>>p2;wczytaj>>p3;
		if(temp){
			warunki.setElement(temp);
			warunki.setElement(p1-1);
			warunki.setElement(p2-1);
			warunki.setElement(p3-1);
		}

	}
	else{break;}


}

wczytaj.close();

	if(wczytajWarunki){

		ustawWarunkiNaPodstawiePunktow(warunki);
		//m1->Lines->Add(warunki.getIter()*0.25);
	}

}

void Delanouy3D::ustawWarunkiNaPodstawiePunktow(IntList &warunkiFace){
//warunkiFace 0-war 1-p1 2-p2 3-p3

	uzupelnienieElementow();
	wyszukanieSasidnichElementowE();
	oznaczWezlyGraniczneNaPodstawieScianG();

	creatFaceAndEdge();
	sasiednieElementy();

	for(int i=0,war,p1,p2,p3,ileI=warunkiFace.getIter();i<ileI;i+=4){

		war=warunkiFace.getElement(i);

		p1=warunkiFace.getElement(i+1);
		p2=warunkiFace.getElement(i+2);
		p3=warunkiFace.getElement(i+3);

		for(int ii=0,e,f,ileII=laplas[p1].getIter();ii<ileII;++ii){

			e = laplas[p1].getElement(ii);
			f=elements.getElement(e).getNumerSciany(p1,p2,p3);

			if(f!=-1 && elements.getElement(e).sprawdzSciane(f,p1,p2,p3)){
				faces.getElement(mapFace[e].getElement(f)).setBC(war);
				break;

			}


		}

	}

}
*/


void Delanouy3D::wczytajPlikNas(const char * nazwa, int mnoznik, bool wczytajWarunki, double px, double py, double pz) {

	points.czysc(100000, 100000);
	elements.czysc(600000, 600000);
	IntList bc_connect(100000, 100000);
	IntList warunki(10, 10);

	if (wczytajWarunki) {
		warunki.czysc(300000, 300000);
	}

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);


//int i=0;
	string textInFile;
	const string testGRID = "GRID";
	const string testCTETRA = "CTETRA";
	const string testCTRIAX = "CTRIAX";
	const string testCTRIAX_BC = "CTRIAX_BC";

	bool bc_conn = false;
	int numer, temp, p1, p2, p3, p4, p5, p6, e1, e2;
	double x, y, z;

	while (!wczytaj.eof()) {

		textInFile = "ble";  //musi byc
		wczytaj >> textInFile;


		if (textInFile == testGRID) {
			wczytaj >> numer;
			wczytaj >> x;
			wczytaj >> y;
			wczytaj >> z;
			wczytaj >> temp;

			points.setElement(x * mnoznik + px, y * mnoznik + py, z * mnoznik + pz);

		}
		else if (textInFile == testCTETRA) {
			wczytaj >> numer;
			wczytaj >> temp;
			wczytaj >> p1;
			wczytaj >> p2;
			wczytaj >> p3;
			wczytaj >> p4;

			elements.setElement(p1 - 1, p2 - 1, p3 - 1, p4 - 1, temp);

		}
		else if (wczytajWarunki && textInFile == testCTRIAX) {
			wczytaj >> numer;
			wczytaj >> temp;
			wczytaj >> p1;
			wczytaj >> p2;
			wczytaj >> p3;
			if (temp) {
				warunki.setElement(temp);
				warunki.setElement(p1 - 1);
				warunki.setElement(p2 - 1);
				warunki.setElement(p3 - 1);
			}

		}
		else if (wczytajWarunki && textInFile == testCTRIAX_BC) {

			bc_conn = true;
			wczytaj >> e1;
			wczytaj >> p1;
			wczytaj >> p2;
			wczytaj >> p3;
			wczytaj >> e2;
			wczytaj >> p4;
			wczytaj >> p5;
			wczytaj >> p6;
			bc_connect.setElement(e1 - 1);
			bc_connect.setElement(p1 - 1);
			bc_connect.setElement(p2 - 1);
			bc_connect.setElement(p3 - 1);
			bc_connect.setElement(-(e2 - 1) - 2);
			bc_connect.setElement(p4 - 1);
			bc_connect.setElement(p5 - 1);
			bc_connect.setElement(p6 - 1);



			//std::cout<<"wspolzedne f0: "<<points.getElement(p1-1).getX()<<"  "<<points.getElement(p1-1).getY()<<"  "<<points.getElement(p1-1).getZ()<<" - "<<points.getElement(p2-1).getX()<<"  "<<points.getElement(p2-1).getY()<<"  "<<points.getElement(p2-1).getZ()<<" - "<<points.getElement(p3-1).getX()<<"  "<<points.getElement(p3-1).getY()<<"  "<<points.getElement(p3-1).getZ()<<endl;
			//std::cout<<"wspolzedne f1: "<<points.getElement(p4-1).getX()<<"  "<<points.getElement(p4-1).getY()<<"  "<<points.getElement(p4-1).getZ()<<" - "<<points.getElement(p5-1).getX()<<"  "<<points.getElement(p5-1).getY()<<"  "<<points.getElement(p5-1).getZ()<<" - "<<points.getElement(p-1).getX()<<"  "<<points.getElement(p6-1).getY()<<"  "<<points.getElement(p6-1).getZ()<<endl;





		}
		else {
			break;
		}


	}

	wczytaj.close();

	if (wczytajWarunki) {

		ustawWarunkiNaPodstawiePunktow(warunki);
		//m1->Lines->Add(warunki.getIter()*0.25);

		if (bc_conn) {
			for (int i = 0, ile = bc_connect.getIter(), nrF; i < ile; i += 8) {
				nrF = elements.getElement(bc_connect.getElement(i)).getNumerSciany(bc_connect.getElement(i + 1), bc_connect.getElement(i + 2), bc_connect.getElement(i + 3));
				elements.getElement(bc_connect.getElement(i)).setE(bc_connect.getElement(i + 4), nrF);

				//nrF = elements.getElement(bc_connect.getElement(i)).getNumerSciany(bc_connect.getElement(i+5),bc_connect.getElement(i+6),bc_connect.getElement(i+7));
				//elements.getElement(-bc_connect.getElement(i+4)-2).setE(-bc_connect.getElement(i)-2,nrF);

			}
			tablicaFC_connect();
		}
	}

}

void Delanouy3D::ustawWarunkiNaPodstawiePunktow(IntList & warunkiFace) {
//warunkiFace 0-war 1-p1 2-p2 3-p3

	uzupelnienieElementow();
	wyszukanieSasidnichElementowE();
	oznaczWezlyGraniczneNaPodstawieScianG();

	creatFaceAndEdge();
	sasiednieElementy();

	for (int i = 0, war, p1, p2, p3, ileI = warunkiFace.getIter(); i < ileI; i += 4) {

		war = warunkiFace.getElement(i);

		p1 = warunkiFace.getElement(i + 1);
		p2 = warunkiFace.getElement(i + 2);
		p3 = warunkiFace.getElement(i + 3);

		for (int ii = 0, e, f, ileII = laplas[p1].getIter(); ii < ileII; ++ii) {

			e = laplas[p1].getElement(ii);
			f = elements.getElement(e).getNumerSciany(p1, p2, p3);

			if (f != -1 && elements.getElement(e).sprawdzSciane(f, p1, p2, p3)) {
				faces.getElement(mapFace[e].getElement(f)).setBC(war);
				break;

			}


		}

	}

}

void Delanouy3D::wczytajPunktyPSSnast(const char * nazwa, int ktoreSerce) {


	double procent = 90;
	double c;
	double x, y, z;
	points.czysc(100000, 100000);
	elements.czysc(500000, 500000);

	double minX = 99999999, maxX = -99999999, tranX = 0;
	double minY = 99999999, maxY = -99999999, tranY = 0;
	double minZ = 99999999, maxZ = -99999999, tranZ = 0;
	double sX, sY, sZ;

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);


//wczytanie z pliku
	string textInFile;
	string liczby;
	string testGRID = "GRID";
	string testCTETRA = "CTETRA";

	string tempLiczby = "";
	int nastepny = 1, nn;
	int ileWezlowWczytuje = 0, podliczWezly = 0;

	for (int i = 0, p1, p2, p3, p4;; ++i) {

		wczytaj >> textInFile;

		if (textInFile == testGRID) {

			wczytaj >> nn;
			wczytaj >> liczby;

			if (liczby.size() != 25) {

				wczytaj.get((char *)tempLiczby.c_str(), 26 - liczby.size());
				liczby += tempLiczby;
				// m1->Lines->Add(liczby.c_str());
				tempLiczby = "";

			}

			wczytaj >> p1;
			wczytaj >> p1;

			if (nn == nastepny) {


				nastepny = nn + 1;
				zwrocLiczby(x, y, z, liczby);
				points.setElement(x, y, z);
				++ileWezlowWczytuje;
				++podliczWezly;

				if (x < minX) {
					minX = x;
				}
				if (x > maxX) {
					maxX = x;
				}
				if (y < minY) {
					minY = y;
				}
				if (y > maxY) {
					maxY = y;
				}
				if (z < minZ) {
					minZ = z;
				}
				if (z > maxZ) {
					maxZ = z;
				}


				//if(10<podliczWezly){m1->Lines->Add(ileWezlowWczytuje);podliczWezly=0;}
				//m1->Lines->Add(ileWezlowWczytuje);
			}
			else {

				for (int t = 0, ilet = nn - nastepny; t < ilet; ++t) {
					points.setElement(0, 0, 0);
				}


				nastepny = nn + 1;

				zwrocLiczby(x, y, z, liczby);
				points.setElement(x, y, z);
				++ileWezlowWczytuje;
				++podliczWezly;

				if (x < minX) {
					minX = x;
				}
				if (x > maxX) {
					maxX = x;
				}
				if (y < minY) {
					minY = y;
				}
				if (y > maxY) {
					maxY = y;
				}
				if (z < minZ) {
					minZ = z;
				}
				if (z > maxZ) {
					maxZ = z;
				}
				//m1->Lines->Add(ileWezlowWczytuje);


			}
		}
		else if (textInFile == testCTETRA) {


			wczytaj >> c;
			wczytaj >> c;
			wczytaj >> p1;
			wczytaj >> p2;
			wczytaj >> p3;
			wczytaj >> p4;

			elements.setElement(p1 - 1, p2 - 1, p3 - 1, p4 - 1);

		}
		else if (textInFile == "ENDDATA") { //m1->Lines->Add("ENDDATA");
			break;
		}
		else {
			break;
		}

	}

	/*
	m1->Lines->Add(minX);
	m1->Lines->Add(maxX);
	m1->Lines->Add(minY);
	m1->Lines->Add(maxY);
	m1->Lines->Add(minZ);
	m1->Lines->Add(maxZ);
	*/

	//m1->Lines->Add(ileWezlowWczytuje);

	double skal;

	if ((minX - maxX) * (minX - maxX) > (minY - maxY) * (minY - maxY)) {
		if ((minX - maxX) * (minX - maxX) > (minZ - maxZ) * (minZ - maxZ)) { //x
			skal = procent / sqrt((minX - maxX) * (minX - maxX));
		}
		else { //z
			skal = procent / sqrt((minZ - maxZ) * (minZ - maxZ));
		}
	}
	else {
		if ((minY - maxY) * (minY - maxY) > (minZ - maxZ) * (minZ - maxZ)) { //y
			skal = procent / sqrt((minY - maxY) * (minY - maxY));
		}
		else { //z
			skal = procent / sqrt((minZ - maxZ) * (minZ - maxZ));
		}

	}

	for (int i = 0, ilePunktow = points.getIter(); i < ilePunktow; ++i) {
		// punkty.getElement(i).setPunkt(punkty.getElement(i).getX()+tranX,punkty.getElement(i).getY()+tranY,punkty.getElement(i).getZ()+tranZ);
		points.getElement(i) = points.getElement(i) * skal;
		//testG.setElement(punkty.getElement(i));
	}

	minX = 99999999;
	maxX = -99999999;
	tranX = 0;
	minY = 99999999;
	maxY = -99999999;
	tranY = 0;
	minZ = 99999999;
	maxZ = -99999999;
	tranZ = 0;

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {

		x = points.getElement(i).getX();
		y = points.getElement(i).getY();
		z = points.getElement(i).getZ();

		if (x < minX) {
			minX = x;
		}
		if (x > maxX) {
			maxX = x;
		}
		if (y < minY) {
			minY = y;
		}
		if (y > maxY) {
			maxY = y;
		}
		if (z < minZ) {
			minZ = z;
		}
		if (z > maxZ) {
			maxZ = z;
		}

	}

	sX = (minX + maxX) * 0.5;
	sY = (minY + maxY) * 0.5;
	sZ = (minZ + maxZ) * 0.5;

	tranX = 50 - sX;
	tranY = 50 - sY;
	tranZ = 50 - sZ;


	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		points.getElement(i).setPunkt(points.getElement(i).getX() + tranX, points.getElement(i).getY() + tranY, points.getElement(i).getZ() + tranZ, 'g');
	}
	//punkty.ustawIter(0);


	globalTranX = tranX;
	globalTranY = tranY;
	globalTranZ = tranZ;
	globalSkal = skal;

	//m1->Lines->Add(1/globalSkal);

	wczytaj.close();
}

void Delanouy3D::zapisParaView(const char * nazwa) {

	ofstream zapis(nazwa);
	zapis.precision(20);

	int ileP = points.getIter();
	int ileE = elements.getIter();
	int ilePryzm = pryzmy.getIter();

	zapis << "# vtk DataFile Version 2.0" << endl;
	zapis << "Navier-Stokes: pdd_navstokes" << endl;
	zapis << "ASCII" << endl;
	zapis << "DATASET UNSTRUCTURED_GRID" << endl;
	zapis << "POINTS" << " " << ileP << " " << "double" << endl;

	for (int i = 0; i < ileP; i++) {
		zapis << points.getElement(i).getX() << "  " << points.getElement(i).getY() << "  " << points.getElement(i).getZ() << endl;
	}

	zapis << "CELLS" << " " << ilePryzm + ileE << " " << ilePryzm * 7 + ileE * 5 << endl;

	for (int i = 0; i < ileE; i++) {
		zapis << 4 << " " << elements.getElement(i).getP1() << " " << elements.getElement(i).getP2() << " " << elements.getElement(i).getP3() << " " << elements.getElement(i).getP4() << endl;
	}
	for (int i = 0; i < ilePryzm; i++) {
		zapis << 6 << " " << pryzmy.getElement(i).getP1() << " " << pryzmy.getElement(i).getP2() << " " << pryzmy.getElement(i).getP3() << " " << pryzmy.getElement(i).getP4() << " " << pryzmy.getElement(i).getP5() << " " << pryzmy.getElement(i).getP6() << endl;
	}

	zapis << "CELL_TYPES" << " " << ilePryzm + ileE << endl;

	for (int i = 0; i < ileE; i++) {
		zapis << 10 << endl;
	}
	for (int i = 0; i < ilePryzm; i++) {
		zapis << 13 << endl;
	}

	zapis << "POINT_DATA" << " " << ileP << endl;
	zapis << "SCALARS pressure double 1" << endl;
	zapis << "LOOKUP_TABLE default" << endl;

	for (int i = 0; i < ileP; i++) {
		zapis << 0 << endl;
	}

	zapis << "VECTORS velocity double" << endl;

	for (int i = 0; i < ileP; i++) {
		zapis << 0 << " " << 0 << " " << 0 << endl;
	}



	zapis.close();
}
