#include <cstdlib>
#include <modfem/mm_remesh/rozrost.h>



void Rozrost::BFSOnlyStep(double dlX, double dlY, double dlZ, double stosunekH, double dlWylou, double gestosc) {

	double Wstopien = dlWylou;
	double Hstopien = dlY / stosunekH;
	double dLkroku;

	pGranicy.ustawIter(0);

	//text->Lines->Add("Hstopien , Wstopien");
	//text->Lines->Add(Hstopien);
	//text->Lines->Add(Wstopien);

	//schodek
	for (int w = 0; w < dlZ + 1; ++w) {

		pGranicy.setElement(0, Hstopien, w, 'g');
		pGranicy.setElement(Wstopien, Hstopien, w, 'g');
		pGranicy.setElement(Wstopien, 0, w, 'g');

		dLkroku = Wstopien / (int(Wstopien * gestosc));

		for (int i = 1, ileP = Wstopien * gestosc; i < ileP; ++i) {
			pGranicy.setElement(i * dLkroku, Hstopien, w, 'g');
		}

		dLkroku = (Hstopien) / (int((Hstopien) * gestosc));

		for (int i = 1, ileP = (Hstopien) * gestosc; i < ileP; ++i) {
			pGranicy.setElement(Wstopien, 0 + i * dLkroku, w, 'g');
		}

	}

}

void Rozrost::granicaObecnePunkty(PunktList & p, bool tylkoG) {

	if (tylkoG) {
		int zlicz = 0;

		for (int i = 0, ileP = p.getIter(); i < ileP; ++i) {
			if (p.getElement(i).getGranica() == 'g') {
				zlicz++;
			}
		}

		ilePG = 0;
		usunGR.czysc(1000, 1000);
		numerZiarna = 0;
		pGranicy.czysc(1000, zlicz + 10);

		for (int i = 0, ileP = p.getIter(); i < ileP; ++i) {
			if (p.getElement(i).getGranica() == 'g') {
				pGranicy.setElement(p.getElement(i));
			}
		}

	}
	else {
		ilePG = 0;
		usunGR.czysc(1000, 1000);
		numerZiarna = 0;
		pGranicy.czysc(1000, p.getIter() + 10);

		for (int i = 0, ileP = p.getIter(); i < ileP; ++i) {
			pGranicy.setElement(p.getElement(i));
		}
	}

}

void Rozrost::warstwyZ2D(int ileWarstw) {

	double wielkoscX = pR->getWielkoscX();
	double wielkoscY = pR->getWielkoscY();
	double wielkoscZ = pR->getWielkoscZ();


	double wsp = double(1) / (ileWarstw - 1);

	pGranicy.czysc(pGranicy2D.getIter() + 10, ileWarstw * pGranicy2D.getIter() + 10);


	double x, y, z;
	bool flagaPoz00;
	char c;

	//text->Lines->Add(wsp);
	for (int i = 0, ileP = pGranicy2D.getIter(); i < ileP; i += 3) {

		x = pGranicy2D.getElement(i);
		y = pGranicy2D.getElement(i + 1);
		c = 'a';
		if (pGranicy2D.getElement(i + 2) == 1) {
			c = 'g';
		}


		flagaPoz00 = true;

		if (!((x == 0 || x == wielkoscX) && (y == 0 || y == wielkoscY))) {
			flagaPoz00 = false;
		}


		for (int i = 0, ileW = ileWarstw; i < ileW; ++i) {

			z = i * wsp;

			if ( !(flagaPoz00 && (z == 0 || z == wielkoscZ) ) ) {
				pGranicy.setElement(x, y, wielkoscZ * z, c);
			}

		}



		//pGranicy.setElement(x,y,wielkoscZ*0.4,'g');
		//pGranicy.setElement(x,y,wielkoscZ*0.425,'g');
		//pGranicy.setElement(x,y,wielkoscZ*0.45,'g');

		/*
		pGranicy.setElement(x,y,wielkosc*0.05,'g');
		if(!flagaPoz00){pGranicy.setElement(x,y,wielkoscZ*0,'g');}
		*/


	}

}
/*
void Rozrost::graniceAnizotropowe(TMemo *text){

int wielkosc=pR->getWielkosc();
int rozd=39;

rozd+=2;
int srodek=rozd*0.5;

double ii=0,jj=0,x,y;

pGranicy2D.czysc(rozd*rozd+10*2,rozd*rozd+10*2);

for(int i=0;i<rozd;i++){jj=0;

	for(int j=0;j<rozd;j++){


		if(i<srodek){x=double(i*i*0.04+i*0.2);}
		else{x=rozd-double(ii*ii*0.04+ii*0.2)-1;}

		if(j<srodek){y=double(j*j*0.04+j*0.2);}
		else{y=rozd - double(jj*jj*0.04+jj*0.2)-1;++jj;}

			pGranicy2D.setElement(x);
			pGranicy2D.setElement(y);

	}
	if(!(i<srodek)){ii++;}

}

pGranicy.zmianaKolejnosci();
}
*/


void Rozrost::graniceStrukturalna() {

	int wielkoscX = pR->getWielkoscX();
	int wielkoscY = pR->getWielkoscY();

	for (int i = 0; i <= wielkoscY; i += 5) {
		for (int j = 0; j <= wielkoscX; j += 5) {

			pGranicy2D.setElement(i);
			pGranicy2D.setElement(j);

		}
	}

}


void Rozrost::granicaTuba(double r, int ilPunkOkrag, int grubosc, double Sx, double Sy, double Sz) {

	double x, y;
//double x1,y1;
	double z;
//pGranicy.czysc(2000,1000);
	int rozZ = pR->getWielkoscZ() + 1;

	double iloscStopni;
	double pi = 3.1415926535;
	double stopien = 2 * pi / 360;
	double obrXY = stopien;
	//double obrXZ=stopien;
	int iLmalyOkrag = ilPunkOkrag;
	iloscStopni = 360 / ilPunkOkrag ;


	obrXY = 0;
	double noweR = r;
	//double noweROtoczenie=r-1;

	for (int i = 0; i < iLmalyOkrag; ++i) {

		z = 0;
		x = noweR * sin(obrXY); //*cos(obrXY);
		y = noweR * cos(obrXY); //*cos(obrXY);
		//x1 = noweROtoczenie*sin(obrXY);//*cos(obrXY);
		//y1 = noweROtoczenie*cos(obrXY);//*cos(obrXY);
		bool flaga = true;
		bool flaga2 = false;

		if (noweR + 0.01 > r) {
			flaga2 = true;
		}

		for (int k = 0, ileP = pGranicy.getIter(); k < ileP; k++) {
			if (0.7 >= pGranicy.getElement(k).odlegloscPunktuSqrt(x + Sx, y + Sy, z + Sz)) {
				flaga = false;
			}
		}

		if (flaga) {

			for (int w = 0; w < rozZ; w += grubosc) {

				pGranicy.setElement(x + Sx, y + Sy, w, 'g');
				//pGranicy.setElement(x1+Sx,y1+Sy,w,'a');

			}

		}

		obrXY += stopien * iloscStopni;

	}


}

void Rozrost::granicaZastawka(bool tylkoGranicaG, double r, int ilPunkOkrag, int grubosc, double Sx, double Sy, double Sz) {

	double x;
	double y;
	double z;
	pGranicy.czysc(2000, 1000);

	double iloscStopni;
	double pi = 3.1415926535;
	double stopien = 2 * pi / 360;
	double obrXY = stopien;
	double obrXZ = stopien;
	int iLmalyOkrag = ilPunkOkrag;
	iloscStopni = 360 / ilPunkOkrag ;

	obrXZ = 0;
	obrXY = 0;
	double noweR = r;
	for (int j = 0; j < r + 1; ++j) {


		obrXY = 0;

		for (int i = 0; i < iLmalyOkrag; ++i) {

			z = 0;
			x = noweR * sin(obrXY); //*cos(obrXY);
			y = noweR * cos(obrXY); //*cos(obrXY);
			bool flaga = true;
			bool flaga2 = false;

			if (noweR + 0.01 > r) {
				flaga2 = true;
			}

			for (int k = 0, ileP = pGranicy.getIter(); k < ileP; k++) {
				if (0.7 >= pGranicy.getElement(k).odlegloscPunktuSqrt(x + Sx, y + Sy, z + Sz)) {
					flaga = false;
				}
			}

			if (flaga) {
				double wsp = (grubosc - 1) * 0.5;
				if (tylkoGranicaG) {
					for (int ilW = 0; ilW < grubosc; ++ilW) {

						if (ilW == 0 || ilW == grubosc - 1 || flaga2) {
							pGranicy.setElement(x + Sx, y + Sy, z + Sz + ilW - wsp, 'g');
						}
						else {
							pGranicy.setElement(x + Sx, y + Sy, z + Sz + ilW - wsp, 'a');
						}
					}
				}
				else {
					for (int ilW = 0; ilW < grubosc; ++ilW) {

						pGranicy.setElement(x + Sx, y + Sy, Sz + ilW - wsp, 'g');

					}
				}
			}

			obrXY += stopien * iloscStopni;

		}
		--noweR;
		//obrXZ += stopien*iloscStopni;

	}


}

void Rozrost::granicaKuli(double r, int ilPunkOkrag, double Sx, double Sy, double Sz) {

	double x;
	double y;
	double z;
	pGranicy.czysc(2000, 1000);

	double iloscStopni;
	double pi = 3.1415926535;
	double stopien = 2 * pi / 360;
	double obrXY = stopien;
	double obrXZ = stopien;
	int iLmalyOkrag = ilPunkOkrag;
	iloscStopni = 360 / ilPunkOkrag ;

	obrXZ = 0;
	obrXY = 0;
	for (int j = 0; j < ilPunkOkrag; ++j) {

		obrXY = 0;

		for (int i = 0; i < iLmalyOkrag; ++i) {

			z = r * sin(obrXZ);
			x = r * cos(obrXZ) * cos(obrXY);
			y = r * cos(obrXZ) * sin(obrXY);
			bool flaga = true;

			for (int k = 0, ileP = pGranicy.getIter(); k < ileP; k++) {
				if (0.8 >= pGranicy.getElement(k).odlegloscPunktuSqrt(x + Sx, y + Sy, z + Sz)) {
					flaga = false;
				}
			}

			if (flaga) {
				pGranicy.setElement(x + Sx, y + Sy, z + Sz, 'g');
			}

			obrXY += stopien * iloscStopni;

		}

		obrXZ += stopien * iloscStopni;

	}


}

void Rozrost::sztucznyRozrost() {

	int rozX = pR->getWielkoscX();
	int rozY = pR->getWielkoscY();
	int rozZ = pR->getWielkoscZ();

	for (int z = 0, i; z < rozZ; z++) {
		for (int y = 0; y < rozY; y++) {
			for (int x = 0; x < rozX; x++) {

				i = x + y * rozX + z * rozX * rozY;

				if (z < 10) {
					pR->komorki[i].setRodzajKomorki(1);
				}
				else {
					pR->komorki[i].setRodzajKomorki(2);
				}

			}
		}
	}


//koniec
}

void Rozrost::pomniejszGranice(int promien, bool brzeg) {

	for (int i = 0, ilePG = pGranicy.getIter(); i < ilePG; i++) {
		pGranicy.getElement(i).setGranica('g');
	}

	for (int i = 0, ilePG = pGranicy.getIter(); i < ilePG; i++) {
		if (pGranicy.getElement(i).getGranica() == 'g') {
			for (int j = i + 1; j < ilePG; j++) {
				if (pGranicy.getElement(j).getGranica() == 'g') {
					if (promien > pGranicy.getElement(i).odlegloscPunktu(pGranicy.getElement(j))) {
						pGranicy.getElement(j).setGranica('a');
					}
				}
			}
		}
	}

	double wymiarPX =  pR->getWielkoscX();
	double wymiarPY =  pR->getWielkoscY();
	double wymiarPZ =  pR->getWielkoscZ();

	for (int i = 0; i < pGranicy.getIter(); i++) {

		if (!brzeg) {
			if (pGranicy.getElement(i).getX() == 0 || pGranicy.getElement(i).getX() == wymiarPX
					|| pGranicy.getElement(i).getY() == 0 || pGranicy.getElement(i).getY() == wymiarPY
					|| pGranicy.getElement(i).getZ() == 0 || pGranicy.getElement(i).getZ() == wymiarPZ) {
				pGranicy.getElement(i).setGranica('g');
			}
		}

		if (pGranicy.getElement(i).getGranica() == 'a') {
			pGranicy.usunElement(i--);
		}

	}
}



void Rozrost::graniceP2() {

	int wielkoscX = pR->getWielkoscX();
	int wielkoscY = pR->getWielkoscY();
//int wielkoscZ=pR->getWielkoscZ();

	pGranicy.czysc(1000, 1000);

	for (int i = 0; i <= wielkoscY; i += 2) {
		for (int j = 0; j <= wielkoscX; j += 2) {

			//pGranicy.setElement(j,i,wielkoscZ*0.5,'g');
			pGranicy.setElement(j, wielkoscY * 0.5, i, 'g');
			pGranicy.setElement(wielkoscX * 0.5, i, j, 'g');

		}
	}

}

void Rozrost::ustalGraniceRozrostu() {

	gR.czysc(100000, 200000);
	int wielkosc = pR->getWielkoscX() * pR->getWielkoscY() * pR->getWielkoscZ();
	for (int i = 0; i < wielkosc; i++) {

		if (pR->komorki[i].getRodzajKomorki()) {
			gR.setElement(i);
		}

	}
}


void Rozrost::rozrostLosujRange(int ileZiaren, int range) {
//pR->zerojKomorki();
	numerZiarna = 0;
	pR->resetZiaren(pR->ziarna.getIter());
	pR->resetKomorek();
//int ileZ =pR->ziarna.getIter();


//time_t start,koniec;
//start = clock();
	int dok  =  range;

	int WYMIARX = pR->getWielkoscX();
	int WYMIARY = pR->getWielkoscY();
	int WYMIARZ = pR->getWielkoscZ();

	int wymLosX = WYMIARX - 2;
	int wymLosY = WYMIARY - 2;
	int wymLosZ = WYMIARZ - 2;


	double dok1 = dok, dok2, dok3, dok4, dok5, dok6, dok7;
	bool flaga = false;

	int ileTabPX = WYMIARX + dok * 2;
	int ileTabPY = WYMIARY + dok * 2;
	int ileTabPZ = WYMIARZ + dok * 2;

	char ** TabP;
	TabP = new char * [ileTabPZ];
	for (int i(0); i < ileTabPZ; i++) {
		TabP[i] = new char[ileTabPX * ileTabPY];
	}

	for (int i = 0, ilePunkt = ileTabPX * ileTabPY; i < ileTabPZ; i++) {
		for (int j = 0; j < ilePunkt; j++) {
			TabP[i][j] = 'a';
		}
	}


	int x, y, z;
	int a = 0, b = 0;
	int los;

	for (int i = 0; i < ileZiaren; i++) {
		double superdok = dok * dok * dok;
		if (a > 500) {
			break;
		}

		x = 1 + rand() % wymLosX + dok;
		y = 1 + rand() % wymLosY + dok;
		z = 1 + rand() % wymLosZ + dok;

		dok6 = z - dok1;
		dok7 = z + dok1;
		dok2 = y - dok1;
		dok3 = y + dok1;
		dok4 = x - dok1;
		dok5 = x + dok1;

		flaga = true;

		for (int i0 = dok6; i0 < dok7; i0++) {
			for (int i1 = dok2; i1 < dok3; i1++) {
				for (int i2 = dok4; i2 < dok5; i2++) {
					if (TabP[i0][i2 + i1 * ileTabPX] == 'b') {

						if ((x - i2) * (x - i2) + (y - i1) * (y - i1) + (z - i0) * (z - i0) < superdok) {
							flaga = false;
							i--;
							a++;
							i1 = dok3;
							i0 = dok7;
							break;
						}

					}
				}

			}
		}

		if (flaga) {

			los = (x - dok) + (y - dok) * WYMIARX + (z - dok) * WYMIARX * WYMIARY;


			pR->komorki[los].setRodzajKomorki(numerZiarna);
			int num;
			num = pR->komorki[los].sasiad[4];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[10];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[12];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[13];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[15];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[21];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}

			numerZiarna++;
			pR->ziarna.setElement(numerZiarna, rand() % 101 * 0.01, rand() % 101 * 0.01, rand() % 101 * 0.01);

			TabP[z][y * ileTabPX + x] = 'b';
			b = a;
			a = 0;
		}

	}

//koniec = clock();


	for (int i = 0; i < ileTabPZ; i++) {
		delete []TabP[i];
	}
	delete []TabP;


//text->Lines->Add(double(koniec-start)/CLOCKS_PER_SEC);
//text->Lines->Add(pR->ziarna.getIter()-ileZ-1);
	ustalGraniceRozrostu();
}


void Rozrost::rozrostLosuj(int ileZiaren) {


	int los = 0;
	int wielkosc = pR->getWielkoscX() * pR->getWielkoscY() * pR->getWielkoscZ();
	for (int i = 0; i < ileZiaren; i++) {

		los = rand() * rand() % wielkosc ;

		if (pR->komorki[los].getRodzajKomorki() == 0) {

			numerZiarna++;

			pR->komorki[los].setRodzajKomorki(numerZiarna);

			int num;
			num = pR->komorki[los].sasiad[4];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[10];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[12];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[13];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[15];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}
			num = pR->komorki[los].sasiad[21];
			if (num != -1) {
				pR->komorki[num].setRodzajKomorki(numerZiarna);
			}

			pR->ziarna.setElement(numerZiarna, rand() % 101 * 0.01, rand() % 101 * 0.01, rand() % 101 * 0.01);

		}


	}

	ustalGraniceRozrostu();

}


void Rozrost::graniceP1Per() {

	pGranicy.czysc(100000, 100000);

	int rozX = pR->getWielkoscX();
	int rozY = pR->getWielkoscY();
	int rozZ = pR->getWielkoscZ();

	for (int z = 0, i = 0; z < rozZ; z++) {
		for (int y = 0; y < rozY; y++) {
			for (int x = 0, s1, s2, s3; x < rozX; x++) {

				i = x + y * rozX + z * rozX * rozY;

				s1 = pR->komorki[i].sasiad[13];
				s2 = pR->komorki[i].sasiad[15];
				s3 = pR->komorki[i].sasiad[21];

				int f0 = pR->komorki[i].getRodzajKomorki();
				int f1 = pR->komorki[s1].getRodzajKomorki();
				int f2 = pR->komorki[s2].getRodzajKomorki();
				int f3 = pR->komorki[s3].getRodzajKomorki();


				if (f0 != f1) {

					pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].getPolozenieSZ(), 'g');


					if (y == rozY - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');

					}

					if (y == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].getPolozenieSZ(), 'g');

					}

					if (z == rozZ - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

						if (y == rozY - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}
						else if (y == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						}
					}

					if (z == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[2], 'g');

						if (y == rozY - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						}
						else if (y == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[2], 'g');

						}
					}

				}

				if (f0 != f2) {

					pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');


					if (x == rozX - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');

					}

					if (x == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');

					}
					if (z == rozZ - 1) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}
						else if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}
					}

					if (z == 0) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						}
						else if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						}
					}

				}

				if (f0 != f3) {

					pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

					ilePG++;

					if (x == rozX - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

					}

					if (x == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

					}

					if (y == rozY - 1) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}
						else if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}

					}

					if (y == 0) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						}
						else if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						}

					}

				}


			}
		}
	}

//koniec
}


void Rozrost::graniceP1NiePer() {

	pGranicy.czysc(100000, 100000);

	int rozX = pR->getWielkoscX();
	int rozY = pR->getWielkoscY();
	int rozZ = pR->getWielkoscZ();

	for (int z = 0, i = 0; z < rozZ; z++) {
		for (int y = 0; y < rozY; y++) {
			for (int x = 0, s1, s2, s3; x < rozX; x++) {

				i = x + y * rozX + z * rozX * rozY;

				s1 = pR->komorki[i].sasiad[13];
				s2 = pR->komorki[i].sasiad[15];
				s3 = pR->komorki[i].sasiad[21];

				int f0 = pR->komorki[i].getRodzajKomorki();

				int f1, f2, f3;
				if (s1 == -1) {
					f1 = -1;
				}
				else {
					f1 = pR->komorki[s1].getRodzajKomorki();
				}
				if (s2 == -1) {
					f2 = -1;
				}
				else {
					f2 = pR->komorki[s2].getRodzajKomorki();
				}
				if (s3 == -1) {
					f3 = -1;
				}
				else {
					f3 = pR->komorki[s3].getRodzajKomorki();
				}


				if (f0 != f1 && x != rozX - 1) {

					pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].getPolozenieSZ(), 'g');

					if (y == rozY - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');

					}

					if (y == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].getPolozenieSZ(), 'g');

					}

					if (z == rozZ - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

						if (y == 0) {

							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						}
						else if (y == rozY - 1) {

							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}

					}

					if (z == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[2], 'g');

						if (y == 0) {

							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[2], 'g');

						}
						else if (y == rozY - 1) {

							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						}
					}


				}

				if (f0 != f2 && y != rozY - 1) {

					pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');


					if (x == rozX - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');

					}

					if (x == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].getPolozenieSZ(), 'g');

					}
					if (z == rozZ - 1) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}
						else if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}
					}

					if (z == 0) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						}
						else if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[2], 'g');

						}
					}


				}

				if (f0 != f3 && z != rozZ - 1) {

					pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

					ilePG++;

					if (x == rozX - 1) {
						pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

					}

					if (x == 0) {
						pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].getPolozenieSY(), pR->komorki[i].polozenie[5], 'g');

					}

					if (y == rozY - 1) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}

						if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[4], pR->komorki[i].polozenie[5], 'g');

						}
					}

					if (y == 0) {
						pGranicy.setElement(pR->komorki[i].getPolozenieSX(), pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						if (x == rozX - 1) {
							pGranicy.setElement(pR->komorki[i].polozenie[3], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						}

						if (x == 0) {
							pGranicy.setElement(pR->komorki[i].polozenie[0], pR->komorki[i].polozenie[1], pR->komorki[i].polozenie[5], 'g');

						}

					}

				}


			}
		}
	}

//koniec
}


/*
void Rozrost::graniceP(TMemo *text){

ilePG=0;

for(int i=0,pomoc=0;i<pR->pokazIloscKomorek();i++){

	if(pomoc !=pR->pokazRozdzielczoscPlanszy()-1){
		if(pR->komorki[i].rodzajKomorki()!=pR->komorki[i+1].rodzajKomorki()){

		pGranicy[ilePG].setPG(pR->komorki[i].polozenieKomorki().x2,pR->komorki[i].polozenieKomorki().py);
		pGranicy[ilePG].setSasiad(0,pR->komorki[i].rodzajKomorki());
		pGranicy[ilePG].setSasiad(1,pR->komorki[i+1].rodzajKomorki());
		ilePG++;

	   ///
			if(pomoc ==pR->pokazRozdzielczoscPlanszy()-2){
			pGranicy[ilePG].setPG(pR->komorki[i+1].polozenieKomorki().x2,pR->komorki[i+1].polozenieKomorki().py);

			ilePG++;
			}
			if(pomoc ==0){
			pGranicy[ilePG].setPG(pR->komorki[i].polozenieKomorki().x1,pR->komorki[i].polozenieKomorki().py);

			ilePG++;
			}
	   ///

		}
	pomoc++;
	}
	else{pomoc=0;}
}


for(int i=0;i<pR->pokazIloscKomorek()-pR->pokazRozdzielczoscPlanszy();i++){


		if(pR->komorki[i].rodzajKomorki()!=pR->komorki[i+pR->pokazRozdzielczoscPlanszy()].rodzajKomorki()){

		pGranicy[ilePG].setPG(pR->komorki[i].polozenieKomorki().px,pR->komorki[i].polozenieKomorki().y2);
		pGranicy[ilePG].setSasiad(2,pR->komorki[i].rodzajKomorki());
		pGranicy[ilePG].setSasiad(3,pR->komorki[i+pR->pokazRozdzielczoscPlanszy()].rodzajKomorki());
		ilePG++;

			if(i<pR->pokazRozdzielczoscPlanszy()){
			pGranicy[ilePG].setPG(pR->komorki[i].polozenieKomorki().px,pR->komorki[i].polozenieKomorki().y1);

			ilePG++;
			}
			if(i>pR->pokazIloscKomorek()-pR->pokazRozdzielczoscPlanszy()*2){
			pGranicy[ilePG].setPG(pR->komorki[i+1].polozenieKomorki().px,pR->komorki[i+1].polozenieKomorki().y2);

			ilePG++;
			}


		}

}


}
*/

/*
void Rozrost::graniceZiaren(TMemo *text){



for(int i=0;i<pR->pokazIloscKomorek()-1;i++){


		if(pR->komorki[i].rodzajKomorki()!=pR->komorki[i+1].rodzajKomorki()){

		pR->komorki[i].ustalGranice(true);

		}

}

for(int i=0;i<pR->pokazIloscKomorek()-pR->pokazRozdzielczoscPlanszy();i++){


		if(pR->komorki[i].rodzajKomorki()!=pR->komorki[i+pR->pokazRozdzielczoscPlanszy()].rodzajKomorki()){

		pR->komorki[i].ustalGranice(true);

		}

}



for(int i=0;i<pR->pokazIloscKomorek()-pR->pokazRozdzielczoscPlanszy();i++){

		if(pR->komorki[i].pobierzGranice()){

				if(pR->komorki[i+1+pR->pokazRozdzielczoscPlanszy()].pobierzGranice()){
				pR->komorki[i+1].ustalGranice(false);
				//pR->komorki[i+pR->pokazRozdzielczoscPlanszy()].ustalGranice(false);

				}

				if(pR->komorki[i-1+pR->pokazRozdzielczoscPlanszy()].pobierzGranice()){
			   // pR->komorki[i+pR->pokazRozdzielczoscPlanszy()].ustalGranice(false);
				pR->komorki[i-1].ustalGranice(false);
				}

		}

}


}
*/

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          1
//////////////////////////////////////////////////////////////////////////////////////

void Rozrost::rozrostZiaren() {

	int rozrostowa[27];
	int rodzaj;
	int iloscKomurek = pR->getWielkoscX() * pR->getWielkoscY() * pR->getWielkoscZ();

	for (int i = 0; i < iloscKomurek; i++) {

		rodzaj = pR->komorki[i].getRodzajKomorki();

		if (!pR->komorki[i].runda  && rodzaj != 0) {

			rozrostowa[26] = 0;

			for (int j = 0; j < 26; j++) {

				if (pR->komorki[pR->komorki[i].sasiad[j]].getRodzajKomorki() == rodzaj) {
					rozrostowa[26]++;
					rozrostowa[j] = 1;
				}

				else if (pR->komorki[pR->komorki[i].sasiad[j]].getRodzajKomorki()) {
					rozrostowa[j] = 2;
				}
				else {
					rozrostowa[j] = 0;
				}

			}

			for (int j = 0; j < 26; j++) {
				if (!rozrostowa[j]) {

					if (rozrostowa[26]*rozrostowa[26]*rozrostowa[26] > rand() % 1000000) {

						pR->komorki[pR->komorki[i].sasiad[j]].setRodzajKomorki(rodzaj);
						pR->komorki[pR->komorki[i].sasiad[j]].runda = 1;
					}

				}
			}

		}
	}

	for (int i = 0; i < iloscKomurek; i++) {
		pR->komorki[i].runda = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          1.1
//////////////////////////////////////////////////////////////////////////////////////

void Rozrost::rozrostZiarenGranicy(int wartoscPrawdo) {

	usunGR.ustawIter(0);
	int rozrostowa[27];
	int rodzaj;
	int ileElementow = gR.getIter();


	for (int i = 0, wybranyE; i < ileElementow; i++) {

		wybranyE = gR.getElement(i);
		rodzaj = pR->komorki[wybranyE].getRodzajKomorki();

		rozrostowa[26] = 0;

		for (int j = 0, sasiad; j < 26; j++) {
			sasiad = pR->komorki[wybranyE].sasiad[j];

			if (sasiad == -1) {
				rozrostowa[j] = 1;
			}
			else if (pR->komorki[sasiad].getRodzajKomorki() == rodzaj) {
				rozrostowa[26]++;
				rozrostowa[j] = 1;

			}
			else if (pR->komorki[sasiad].getRodzajKomorki()) {
				rozrostowa[j] = 2;
			}
			else {
				rozrostowa[j] = 0;
			}
		}

		//bool flaga=true;

		for (int j = 0; j < 26; j++) {
			if (!rozrostowa[j]) {

				if (rozrostowa[26]*rozrostowa[26]*rozrostowa[26] > rand() % wartoscPrawdo) {
					pR->komorki[pR->komorki[wybranyE].sasiad[j]].setRodzajKomorki(rodzaj);
					gR.setElement(pR->komorki[wybranyE].sasiad[j]);

				}
				//flaga = false;
			}
		}

		if (rozrostowa[26] == 26) {
			usunGR.setElement(i);
		}
	}

	for (int i = 0; i < usunGR.getIter(); i++) {
		gR.usunElement(usunGR.getElement(i));
	}


}
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          2
//////////////////////////////////////////////////////////////////////////////////////
/*
void Rozrost::rozrostOtoczenia(TMemo *text){

int rozrostowa[9];
int nr,rodzaj;

for(int i=0;i<pR->pokazIloscKomorek();i++){
rodzaj=pR->komorki[i].rodzajKomorki();
if(!pR->komorki[i].runda  && rodzaj==0){

rozrostowa[8]=0;

		for(int j=0;j<8;j++){

		if(pR->komorki[pR->komorki[i].sasiedziKomorki().s[j]].rodzajKomorki()==rodzaj){
		rozrostowa[8]++;
		rozrostowa[j]=1;

		}
		else{rozrostowa[j]=0;}

		}

		for(int j=0;j<8;j++){
		if(!rozrostowa[j]){

				if(rozrostowa[8]*rozrostowa[8]*rozrostowa[8]>rand()%5000){pR->komorki[pR->komorki[i].sasiedziKomorki().s[j]].ustawRodzajKomorki(rodzaj);
				pR->komorki[pR->komorki[i].sasiedziKomorki().s[j]].runda=1;
				}

		}
		}

}
}

for(int i=0;i<pR->pokazIloscKomorek();i++){pR->komorki[i].runda=0;}
}
*/
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          3
//////////////////////////////////////////////////////////////////////////////////////

void Rozrost::pomniejszenieZiarna(int wspolczynnik) {

	int rozrostowa[27];
	int rodzaj;
	int iloscKomurek = pR->getWielkoscX() * pR->getWielkoscY() * pR->getWielkoscZ();

	for (int i = 0; i < iloscKomurek; i++) {
		rodzaj = pR->komorki[i].getRodzajKomorki();
		if (!pR->komorki[i].runda  && rodzaj != 0) {

			rozrostowa[26] = 0;

			for (int j = 0; j < 26; j++) {

				if (pR->komorki[pR->komorki[i].sasiad[j]].getRodzajKomorki() == rodzaj) {
					rozrostowa[26]++;

				}
			}


			if (rozrostowa[26] < wspolczynnik) {

				pR->komorki[i].setRodzajKomorki(0);
				pR->komorki[i].runda = 1;


			}


		}
	}

	for (int i = 0; i < iloscKomurek; i++) {
		pR->komorki[i].runda = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          3.1
//////////////////////////////////////////////////////////////////////////////////////

void Rozrost::pomniejszenieZiarnaGranicy(int wspolczynnik) {

	usunGR.ustawIter(0);
	int rozrostowa[27];
	int rodzaj;

	int ileElementow = gR.getIter();

	for (int i = 0; i < ileElementow; i++) {
		pR->komorki[gR.getElement(i)].runda = 10;
	}

	for (int i = 0, wybranyE; i < ileElementow; i++) {

		wybranyE = gR.getElement(i);
		rodzaj = pR->komorki[wybranyE].getRodzajKomorki();

		rozrostowa[26] = 0;

		for (int j = 0, sasiad; j < 26; j++) {
			sasiad = pR->komorki[wybranyE].sasiad[j];
			if (sasiad != -1) {
				if (pR->komorki[sasiad].getRodzajKomorki() == rodzaj) {
					rozrostowa[26]++;

				}
			}
		}


		if (rozrostowa[26] < wspolczynnik) {

			pR->komorki[wybranyE].setRodzajKomorki(0);

			for (int j = 0, sasiad; j < 26; j++) {
				sasiad = pR->komorki[wybranyE].sasiad[j];
				if (sasiad != -1) {
					if (pR->komorki[sasiad].runda != 10) {

						gR.setElement(sasiad);
						pR->komorki[sasiad].runda = 10;

					}
				}
			}

			usunGR.setElement(i);

		}

	}

	for (int i = 0; i < gR.getIter(); i++) {
		pR->komorki[gR.getElement(i)].runda = 0;
	}
	for (int i = 0; i < usunGR.getIter(); i++) {
		gR.usunElement(usunGR.getElement(i));
	}

}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          4
//////////////////////////////////////////////////////////////////////////////////////
/*
void Rozrost::granice(TMemo *text,int wspolczynnik){

int rozrostowa[9];
int nr,rodzaj;

for(int i=0;i<pR->pokazIloscKomorek();i++){
pR->komorki[i].ustalGranice(true);
rodzaj=pR->komorki[i].rodzajKomorki();
if(!pR->komorki[i].runda  && rodzaj!=0){

rozrostowa[8]=0;

		for(int j=0;j<8;j++){

				if(pR->komorki[pR->komorki[i].sasiedziKomorki().s[j]].rodzajKomorki()==rodzaj){
				rozrostowa[8]++;

				}
		}


		if(rozrostowa[8]>wspolczynnik){

				pR->komorki[i].ustalGranice(false);
				pR->komorki[i].runda=1;


		}


}
}

for(int i=0;i<pR->pokazIloscKomorek();i++){pR->komorki[i].runda=0;}

}
*/

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          5
//////////////////////////////////////////////////////////////////////////////////////

void Rozrost::komplekosoweTworzenieZiaren(int wartoscPrawdo, int powtorzen, int ileRozrost, int ileOtocz, int RozrostOtoczenie) {


	for (int j = 0; j < 15; j++) {
		rozrostZiarenGranicy(1000000);
	}

	RozrostOtoczenie = 14;

	for (int i = 1; i < RozrostOtoczenie; --RozrostOtoczenie) {

		for (int i = 0; i < powtorzen; i++) {

			for (int j = 0; j < ileRozrost; j++) {
				rozrostZiarenGranicy(wartoscPrawdo);
			}

			for (int j = 0; j < ileOtocz; j++) {
				pomniejszenieZiarnaGranicy(RozrostOtoczenie);
			}

		}

		if (RozrostOtoczenie == 7) {
			powtorzen = powtorzen / 2 + 1;
		}

	}


}


/*
void Rozrost::rysujGranice(int polZ){

glPointSize(1);
glColor3f(1,1,1);

	int wielkoscX =  pR->getWielkoscX();
	int wielkoscY =  pR->getWielkoscY();
	int wielkoscZ =  pR->getWielkoscZ();
	int translateX = wielkoscX*0.5;
	int translateY = wielkoscY*0.5;
	int translateZ = wielkoscZ*0.5;

glPushMatrix();
glTranslatef(-wielkoscX,-wielkoscY,-wielkoscZ);
glPushMatrix();
glTranslatef(translateX,translateY,translateZ);

glBegin(GL_POINTS);

		for(int i=0,ilePG=pGranicy.getIter();i<ilePG;i++){

			if(pGranicy.getElement(i).getZ()==polZ){
				glVertex3f(pGranicy.getElement(i).getX(),pGranicy.getElement(i).getY(),pGranicy.getElement(i).getZ());
			}
			else if(polZ==-1){
				glVertex3f(pGranicy.getElement(i).getX(),pGranicy.getElement(i).getY(),pGranicy.getElement(i).getZ());
			}
		}

glEnd();

glPopMatrix();
glPopMatrix();


}
*/

/*
void Rozrost::rysujZiarnaGraniczne(){

glPointSize(1);

	int wielkoscX =  pR->getWielkoscX();
	int wielkoscY =  pR->getWielkoscY();
	int wielkoscZ =  pR->getWielkoscZ();
	int translateX = wielkoscX*0.5;
	int translateY = wielkoscY*0.5;
	int translateZ = wielkoscZ*0.5;

	Kolor barwa(1,0,0);

glPushMatrix();

glTranslatef(-wielkoscX,-wielkoscY,-wielkoscZ);

glPushMatrix();


glTranslatef(translateX,translateY,translateZ);

	for(int i=0,el;i<gR.getIter();i++){
		el=gR.getElement(i);
		if(pR->komorki[el].getRodzajKomorki()){
		barwa = pR->ziarna.getElement(pR->komorki[el].getRodzajKomorki()).pobierzKolor();
		glColor3f(barwa.r,barwa.g,barwa.b);

		//komorki[i].rysujKomorkaPowierzchnie();
		pR->komorki[el].rysujKomorkaPunkt();
		//komorki[i].rysujKomorkaKrawedzie();
		}

	}

glPopMatrix();
glPopMatrix();

}
*/

/*
void Rozrost::wczytajPunktyGraniczne2D(const char *nazwa,bool tylkoGraniczne,int sposobSkalowania){

	int wielkoscX =  pR->getWielkoscX();
	int wielkoscY =  pR->getWielkoscY();
	int wielkoscZ =  pR->getWielkoscZ();

ifstream wczytaj(nazwa);
wczytaj.precision(10);



string text;
double xW,yW;
int ileP;
double temp;
char granica;

wczytaj>>ileP;
pGranicy2D.czysc(ileP+10,ileP+10);
pGranicy.czysc(ileP+10,ileP+10);

if(tylkoGraniczne){

	for(int i=0;i<ileP;++i){

		wczytaj>>temp;
		wczytaj>>xW;
		wczytaj>>yW;
		wczytaj>>granica;
		wczytaj>>temp;

		if(granica == 'g'){

			pGranicy2D.setElement(xW);
			pGranicy2D.setElement(yW);
			pGranicy2D.setElement(1);

		}

	}
}
else{

	for(int i=0;i<ileP;++i){

		wczytaj>>temp;
		wczytaj>>xW;
		wczytaj>>yW;
		wczytaj>>granica;
		wczytaj>>temp;

		pGranicy2D.setElement(xW);
		pGranicy2D.setElement(yW);

		if(granica == 'g'){
			pGranicy2D.setElement(1);
		}
		else{pGranicy2D.setElement(0);}

	}


}

wczytaj.close();

	double maxY=0,mY,maxX=0,mX,minX=99999999,minY=99999999;
	for(int i=0,ileDp=pGranicy2D.getIter();i<ileDp;i+=3){

		if(pGranicy2D.getElement(i)>maxX){maxX=pGranicy2D.getElement(i);}
		if(pGranicy2D.getElement(i+1)>maxY){maxY=pGranicy2D.getElement(i+1);}

		if(pGranicy2D.getElement(i)<minX){minX=pGranicy2D.getElement(i);}
		if(pGranicy2D.getElement(i+1)<minY){minY=pGranicy2D.getElement(i+1);}


	}

	int dzielnik;
	if(1==sposobSkalowania){
		dzielnik=maxX;
		//m1->Lines->Add(" przez X x,y");
		//m1->Lines->Add(((maxX-minX)/dzielnik)*wielkoscX);
		//m1->Lines->Add(((maxY-minY)/dzielnik)*wielkoscX);
	}
	else{
		dzielnik=maxY;
		//m1->Lines->Add(" przez Y x,y");
		//m1->Lines->Add(((maxX-minX)/dzielnik)*wielkoscY);
		//m1->Lines->Add(((maxY-minY)/dzielnik)*wielkoscY);
	}

	for(int i=0,ileDp=pGranicy2D.getIter();i<ileDp;i+=3){

		mX=pGranicy2D.getElement(i)/dzielnik;
		pGranicy2D.podmienElement(i,mX*wielkoscX);

		mY=pGranicy2D.getElement(i+1)/dzielnik;
		pGranicy2D.podmienElement(i+1,mY*wielkoscY);

	}



	for(int i=0,ileDp=pGranicy2D.getIter();i<ileDp;i+=3){

		if(pGranicy2D.getElement(i+2)==1){
			pGranicy.setElement(pGranicy2D.getElement(i),pGranicy2D.getElement(i+1),0,'g');
		}
		else{
			pGranicy.setElement(pGranicy2D.getElement(i),pGranicy2D.getElement(i+1),0,'a');
		}

	}

}
*/

void Rozrost::wczytajPunktyGraniczne(const char * nazwa) {


	ifstream wczytaj(nazwa);
	wczytaj.precision(10);
	string text;

	double xW, yW, zW;
	char granicaW;
	int ileP;

	wczytaj >> ileP;

	pGranicy.czysc(ileP * 0.3, ileP * 0.3);

	for (int i = 0; i < ileP; i++) {

		wczytaj >> xW;
		wczytaj >> xW;
		wczytaj >> yW;
		wczytaj >> zW;
		wczytaj >> granicaW;


		if (granicaW == 'g') {
			pGranicy.setElement(xW, yW, zW, 'g');
		}
	}

	wczytaj.close();

}

void Rozrost::wczytajPunktyGranicznePSS(const char * nazwa) {


	double procent = 90;
	double c;
	double x, y, z;
	pGranicy.czysc(28000, 7000);

	double minX = 99999999, maxX = -99999999, tranX = 0;
	double minY = 99999999, maxY = -99999999, tranY = 0;
	double minZ = 99999999, maxZ = -99999999, tranZ = 0;
	double sX, sY, sZ;

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);
//6795
//26084
//19021
//gr16850
//wczytanie z pliku
	for (int i = 0; i < 19021; ++i) {

		wczytaj >> c;
		wczytaj >> x;
		wczytaj >> y;
		wczytaj >> z;

		pGranicy.setElement(x, y, z, 'g');

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

	for (int i = 0; i < 19021; ++i) {
		// punkty.getElement(i).setPunkt(punkty.getElement(i).getX()+tranX,punkty.getElement(i).getY()+tranY,punkty.getElement(i).getZ()+tranZ);
		pGranicy.getElement(i) = pGranicy.getElement(i) * skal;
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

	for (int i = 0, ileP = pGranicy.getIter(); i < ileP; ++i) {

		x = pGranicy.getElement(i).getX();
		y = pGranicy.getElement(i).getY();
		z = pGranicy.getElement(i).getZ();

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


	for (int i = 0, ileP = pGranicy.getIter(); i < ileP; ++i) {
		pGranicy.getElement(i).setPunkt(pGranicy.getElement(i).getX() + tranX, pGranicy.getElement(i).getY() + tranY, pGranicy.getElement(i).getZ() + tranZ, 'g');
	}
	//punkty.ustawIter(0);


}
