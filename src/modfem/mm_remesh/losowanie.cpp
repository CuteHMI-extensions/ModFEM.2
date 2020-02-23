#include <cstdlib>
#include <modfem/mm_remesh/losowanie.h>



double Losowanie::getWymiarMax() {
	if (wymiarX > wymiarY) {
		if (wymiarX > wymiarZ) {
			return wymiarX;
		}
		else {
			return wymiarZ;
		}

	}

	if (wymiarY > wymiarZ) {
		return wymiarY;
	}
	return wymiarZ;
}

void Losowanie::setPG(PunktList & pg) {
	testG.czysc(pg.getIter() + 100, pg.getIter() * 0.1);
	for (int i = 0, ileP = pg.getIter(); i < ileP; ++i) {
		testG.setElement(pg.getElement(i));
	}
	//for(int i=0,ileP=testG.getIter();i<ileP;++i){testG.getElement(i).setGranica('g');}
}

PunktList & Losowanie::zwrocPunktyZGranica() {

	//for(int i=0,ileP=testG.getIter();i<ileP;i++){
	//punkty.setElement(testG.getElement(i));
	// }

//punkty.losowyZmianaMiejsca();
	return punkty;

}

void Losowanie::dodajPunktyWWolnepola(char ** TabP, int ileTabPX, int ileTabPY, int ileTabPZ) {

	for (int i = 0; i < ileTabPZ; ++i) {
		for (int j = 0; j < ileTabPY; ++j) {
			for (int k = 0; k < ileTabPX; ++k) {

				if (TabP[i][j * ileTabPX + k] == 'a') {
					punkty.setElement(k, j, i);
				}

			}
		}
	}
}

void Losowanie::dodajPunktyPodSpawanie(int ileTabPX, int ileTabPY, int ileTabPZ) {

	int od_j = 80, do_j = 121;
	const int wartSkoku = 2;

	testG.czysc((do_j - od_j)*ileTabPZ * ileTabPX, ((do_j - od_j)*ileTabPZ * ileTabPX) * 0.1);

	for (int i = 0; i < ileTabPZ + 1; i += wartSkoku) {

		//for(int j=0;j<ileTabPY;++j){
		for (int j = od_j; j < do_j; j += wartSkoku) {

			for (int k = 0; k < ileTabPX + 1; k += wartSkoku) {

				testG.setElement(k, j, i, 'g');

				/*
				if(TabP[i][j*ileTabPX+k]=='a'){
					punkty.setElement(k,j,i);
				}
				*/

			}
		}
	}

}

void Losowanie::punktyPoczatkoweOctree(int podziel, int ileRazy) {



	int ileWierszyX = wymiarX / podziel;
	int ileWierszyY = wymiarY / podziel;
	int ileWierszyZ = wymiarZ / podziel;

	int ileKom = ileWierszyX * ileWierszyY * ileWierszyZ;

	bool * uzytyPunkt = new bool[punkty.getIter()];
	IntList komorkaPol(punkty.getIter() * 0.1, punkty.getIter() + 10);
	IntList * komorka = new IntList[ileKom];


	for (int i = 0; i < punkty.getIter(); ++i) {
		uzytyPunkt[i] = true;
	}


	for (int j = 0, xW, yW, zW, numerK; j < punkty.getIter(); ++j) {

		double x = punkty.getElement(j).getX();
		double y = punkty.getElement(j).getY();
		double z = punkty.getElement(j).getZ();

		if (x == wymiarX) {
			x = wymiarX - 1;
		}
		if (y == wymiarY) {
			y = wymiarY - 1;
		}
		if (z == wymiarZ) {
			z = wymiarZ - 1;
		}


		xW = x / podziel;
		yW = y / podziel;
		zW = z / podziel;

		numerK = zW * ileWierszyX * ileWierszyY + yW * ileWierszyX + xW;

		komorka[numerK].setElement(j);


	}

	IntList losowaKolejnosc(10, ileKom);

	for (int i = 0; i < ileKom; ++i) {
		losowaKolejnosc.setElement(i);
	}
	losowaKolejnosc.losowyZmianaMiejsca();

	for (int j = 0; j < ileRazy; ++j) {

		for (int i = 0, nr; i < ileKom; ++i) {

			nr = losowaKolejnosc.getElement(i);

			if (komorka[nr].getIter() > j) {
				komorkaPol.setElement(komorka[nr].getElement(j));
			}

		}

	}

//m1->Lines->Add(komorkaPol.getIter());

	for (int i = 0, ileP = komorkaPol.getIter(); i < ileP; ++i) {
		punkty.zamianaMiejsc(i, komorkaPol.getElement(i));
	}


	delete []komorka;
	delete []uzytyPunkt;

}

void Losowanie::skalujWymiar(double liczba) {

	wymiarX = liczba * wymiarX;
	wymiarY = liczba * wymiarY;
	wymiarZ = liczba * wymiarZ;

	for (int i = 0, ileP = punkty.getIter(); i < ileP; ++i) {
		punkty.getElement(i) = punkty.getElement(i) * liczba;
	}

}

void Losowanie::losujPromien(int wym, int ileP, int promien) {

	punkty.czysc(getWymiarMax()*getWymiarMax(), getWymiarMax()*getWymiarMax() + 100);

	int wymiar = wym + 1;
	double R = promien * promien;
	bool flaga;
	for (int i = 0, x, y, z, a = 0; i < ileP && a < 1000; i++) {

		x = rand() % wymiar;
		y = rand() % wymiar;
		z = rand() % wymiar;
		flaga = true;

		for (int j = 0, iP = punkty.getIter(); j < iP; j++) {

			if (R > punkty.getElement(j).odlegloscPunktu(x, y, z)) {
				i--;
				a++;
				flaga = false;
				break;
			}

		}

		if (flaga) {
			a = 0;
			punkty.setElement(x, y, z);
		}

	}

//m1->Lines->Add(punkty.getIter());

}

/*
void Losowanie::rysuj123(){

glBegin(GL_POINTS);

glColor3f(0,0,0);
int wymiar =100;
int ileTabP = wymiar+1;
char **TabP;
TabP = new char*[ileTabP];
for(int i(0);i<ileTabP;i++){TabP[i] = new char[ileTabP*ileTabP];}


for(int i=0;i<ileTabP;i++){
		for(int j=0;j<ileTabP*ileTabP;j++){TabP[i][j]='a';}
}

for(int i=0,ileP=testG.getIter();i<ileP;i++){

TabP[int(testG.getElement(i).getZ()+0.5)][int(testG.getElement(i).getX()+0.5)+int(testG.getElement(i).getY()+0.5)*ileTabP] = 'g';

}

		for(int i0=0;i0<ileTabP;i0++){
			for(int i1=0;i1<ileTabP;i1++){
				for(int i2=0;i2<ileTabP;i2++){

					if(TabP[i0][i2+i1*ileTabP] != 'a'){
					glVertex3f(i2,i1,i0);
					}

				}
			}
		}


glEnd();
}
*/

void Losowanie::losuj(bool losuj, int wymX, int wymY, int wymZ, int ileP, int promien, int odGranicy, int odPunktu, bool flagaPer, double ax, double b0, int rszuk) {


	pomoc.czysc(100, 310);

	wymiarX = wymX;
	wymiarY = wymY;
	wymiarZ = wymZ;

	startP[0].setPunkt(0, 0, 0);
	startP[1].setPunkt(wymiarX, 0, 0);
	startP[2].setPunkt(0, wymiarY, 0);
	startP[3].setPunkt(wymiarX, wymiarY, 0);
	startP[4].setPunkt(0, 0, wymiarY);
	startP[5].setPunkt(wymiarX, 0, wymiarZ);
	startP[6].setPunkt(0, wymiarY, wymiarZ);
	startP[7].setPunkt(wymiarX, wymiarY, wymiarZ);

	/*
	startP[0].setPunkt(0,0,0);
	startP[1].setPunkt(wymiar-1,0,0);
	startP[2].setPunkt(0,wymiar-1,0);
	startP[3].setPunkt(wymiar-1,wymiar-1,0);
	startP[4].setPunkt(0,0,wymiar-1);
	startP[5].setPunkt(wymiar-1,0,wymiar-1);
	startP[6].setPunkt(0,wymiar-1,wymiar-1);
	startP[7].setPunkt(wymiar-1,wymiar-1,wymiar-1);
	*/

	punkty.czysc(getWymiarMax()*getWymiarMax(), getWymiarMax()*getWymiarMax() + 100);


//time_t start,koniec;
//start = clock();
	int dok  =  promien;



//bool flaga=false;

	int ileTabPX = wymiarX + 1;
	int ileTabPY = wymiarY + 1;
	int ileTabPZ = wymiarZ + 1;

	char ** TabP;
	TabP = new char * [ileTabPZ];
	for (int i(0); i < ileTabPZ; i++) {
		TabP[i] = new char[ileTabPX * ileTabPY];
	}


//uzupelnianie macierzy
	/*
	TabP[0][0]='g';
	TabP[0][0]='g';
	TabP[0][0]='g';
	TabP[0][0]='g';

	TabP[0][0]='g';
	TabP[0][0]='g';
	TabP[0][0]='g';
	TabP[0][0]='g';
	*/

	for (int i = 0; i < ileTabPZ; i++) {
		for (int j = 0; j < ileTabPX * ileTabPY; j++) {
			TabP[i][j] = 'a';
		}
	}

	for (int i = 0, ilePP = testG.getIter(); i < ilePP; i++) {

		TabP[int(testG.getElement(i).getZ() + 0.5)][int(testG.getElement(i).getX() + 0.5) + int(testG.getElement(i).getY() + 0.5)*ileTabPX] = 'g';

	}

//uzupelnianie macierzy

	TabP[0][0] = 's';
	TabP[0][wymiarX] = 's';

	TabP[0][wymiarY * ileTabPX] = 's';
	TabP[0][wymiarY * ileTabPX + wymiarX] = 's';

	TabP[wymiarZ][0] = 's';
	TabP[wymiarZ][wymiarX] = 's';

	TabP[wymiarZ][wymiarY * ileTabPX] = 's';
	TabP[wymiarZ][wymiarY * ileTabPX + wymiarX] = 's';

//uzupelnianie macierzy

	if (losuj) {

		if (flagaPer) {
			losujKrawedziePer(dok, ileTabPX, ileTabPY, ileTabPZ, TabP, odGranicy, odPunktu, ax, b0, rszuk);
			losujPlaszczyznyPer(dok, ileTabPX, ileTabPY, ileTabPZ, TabP, odGranicy, odPunktu, ax, b0, rszuk);
		}
		else {
			losujKrawedzieNiePer(dok, ileTabPX, ileTabPY, ileTabPZ, TabP, odGranicy, odPunktu, ax, b0, rszuk);
			losujPlaszczyznyNiePer(dok, ileTabPX, ileTabPY, ileTabPZ, TabP, odGranicy, odPunktu, ax, b0, rszuk);
		}

		losujWObjetosci(dok, ileTabPX, ileTabPY, ileTabPZ, TabP, ileP, odGranicy, odPunktu, ax, b0, rszuk);
	}
	else {
		dodajPunktyWWolnepola(TabP, ileTabPX, ileTabPY, ileTabPZ);
	}


	punkty.losowyZmianaMiejsca();

	for (int i = 0, ilePP = testG.getIter(); i < ilePP; ++i) {
		punkty.setElement(testG.getElement(i));
	}



//koniec = clock();

	for (int i = 0; i < ileTabPZ; i++) {
		delete []TabP[i];
	}
	delete []TabP;


//m1->Lines->Add(double(koniec-start)/CLOCKS_PER_SEC);
//m1->Lines->Add(punkty.getIter());

//ustalGraniceRozrostu();

}

void Losowanie::losujPlaszczyznyNiePer(int dok, int ileTabPX, int ileTabPY, int ileTabPZ, char ** TabP, int odGranicy, int odPunktu, double ax, double b0, int rszuk) {

	double dok1 = dok, dok2, dok3, dok4, dok5, dok6, dok7;
	dok1 = dok + rszuk;
	bool flaga;
	int a = 0, b = 0;
	int ktory = 0;
//double superdok = dok*dok;
	double range0, range, max0, max;


	for (int i = 0, x = 0, y = 0, z = 0, rozX = ileTabPX - 1, rozY = ileTabPY - 1, rozZ = ileTabPZ - 1; a < 300; i++) {


		switch (ktory) {

			case 0:
				x = rand() % ileTabPX;
				y = rand() % ileTabPY;
				z = 0;
				break;
			case 1:
				y = rand() % ileTabPY;
				z = rand() % ileTabPZ;
				x = 0;
				break;
			case 2:
				z = rand() % ileTabPZ;
				x = rand() % ileTabPX;
				y = 0;
				break;
			case 3:
				x = rand() % ileTabPX;
				y = rand() % ileTabPY;
				z = rozZ;
				break;
			case 4:
				y = rand() % ileTabPY;
				z = rand() % ileTabPZ;
				x = rozX;
				break;
			case 5:
				z = rand() % ileTabPZ;
				x = rand() % ileTabPX;
				y = rozY;
				break;

		}


		dok2 = y - dok1;
		dok3 = y + dok1;
		dok4 = x - dok1;
		dok5 = x + dok1;
		dok6 = z - dok1;
		dok7 = z + dok1;

		if (dok2 < 0) {
			dok2 = 0;
		}
		if (dok4 < 0) {
			dok4 = 0;
		}
		if (dok6 < 0) {
			dok6 = 0;
		}
		if (dok3 > ileTabPY) {
			dok3 = ileTabPY;
		}
		if (dok5 > ileTabPX) {
			dok5 = ileTabPX;
		}
		if (dok7 > ileTabPZ) {
			dok7 = ileTabPZ;
		}

		range0 = dok * dok * dok;
		max0 = 99999999;

		flaga = true;

		for (int i0 = dok6; i0 < dok7; i0++) {
			for (int i1 = dok2; i1 < dok3; i1++) {
				for (int i2 = dok4; i2 < dok5; i2++) {

					if (TabP[i0][i2 + i1 * ileTabPX] != 'a') {
						max = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odPunktu
						if (max <= odPunktu) {
							flaga = false;
						}
						if (max < max0) {
							max0 = max;
						}
					}

					if (TabP[i0][i2 + i1 * ileTabPX] == 'g') {
						range = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odGranicy
						if (range <= odGranicy) {
							flaga = false;
						}
						if (range < range0) {
							range0 = range;
						}
					}

				}
			}
		}


		if (max0 < range0 * ax + b0 || flaga == false) {
			i--;
			a++;
		}
		else {

			switch (ktory) {

				case 0:

					punkty.setElement(x, y, 0, 'a');
					TabP[0][y * ileTabPX + x] = 'b';
					ktory = 1;
					break;

				case 1:

					punkty.setElement(0, y, z, 'a');
					TabP[z][y * ileTabPX] = 'b';
					ktory = 2;
					break;

				case 2:

					punkty.setElement(x, 0, z, 'a');
					TabP[z][x] = 'b';
					ktory = 3;
					break;

				case 3:

					punkty.setElement(x, y, rozZ, 'a');
					TabP[rozZ][y * ileTabPX + x] = 'b';
					ktory = 4;
					break;


				case 4:

					punkty.setElement(rozX, y, z, 'a');
					TabP[z][y * ileTabPX + rozX] = 'b';
					ktory = 5;
					break;

				case 5:

					punkty.setElement(x, rozY, z, 'a');
					TabP[z][rozY * ileTabPX + x] = 'b';
					ktory = 0;
					break;

			}

			b = a;
			a = 0;
		}

	}




}

void Losowanie::losujPlaszczyznyPer(int dok, int ileTabPX, int ileTabPY, int ileTabPZ, char ** TabP, int odGranicy, int odPunktu, double ax, double b0, int rszuk) {

	double dok1 = dok, dok2, dok3, dok4, dok5, dok6, dok7;
	dok1 = dok + rszuk;
	bool flaga;
	int a = 0, b = 0;
	int ktory = 0;
//double superdok = dok*dok;
	double range0, range, max0, max;


	for (int i = 0, x = 0, y = 0, z = 0, rozX = ileTabPX - 1, rozY = ileTabPY - 1, rozZ = ileTabPZ - 1; a < 300; i++) {


		switch (ktory) {

			case 0:
				x = rand() % ileTabPX;
				y = rand() % ileTabPY;
				z = 0;
				break;
			case 1:
				y = rand() % ileTabPY;
				z = rand() % ileTabPZ;
				x = 0;
				break;
			case 2:
				z = rand() % ileTabPZ;
				x = rand() % ileTabPX;
				y = 0;
				break;

		}

		dok2 = y - dok1;
		dok3 = y + dok1;
		dok4 = x - dok1;
		dok5 = x + dok1;
		dok6 = z - dok1;
		dok7 = z + dok1;

		if (dok2 < 0) {
			dok2 = 0;
		}
		if (dok4 < 0) {
			dok4 = 0;
		}
		if (dok6 < 0) {
			dok6 = 0;
		}
		if (dok3 > ileTabPY) {
			dok3 = ileTabPY;
		}
		if (dok5 > ileTabPX) {
			dok5 = ileTabPX;
		}
		if (dok7 > ileTabPZ) {
			dok7 = ileTabPZ;
		}

		range0 = dok * dok * dok;
		max0 = 99999999;

		flaga = true;

		for (int i0 = dok6; i0 < dok7; i0++) {
			for (int i1 = dok2; i1 < dok3; i1++) {
				for (int i2 = dok4; i2 < dok5; i2++) {

					if (TabP[i0][i2 + i1 * ileTabPX] != 'a') {
						max = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odPunktu
						if (max <= odPunktu) {
							flaga = false;
						}
						if (max < max0) {
							max0 = max;
						}
					}

					if (TabP[i0][i2 + i1 * ileTabPX] == 'g') {
						range = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odGranicy
						if (range <= odGranicy) {
							flaga = false;
						}
						if (range < range0) {
							range0 = range;
						}
					}

				}
			}
		}


		if (max0 < range0 * ax + b0 || flaga == false) {
			i--;
			a++;
		}
		else {

			switch (ktory) {

				case 0:

					punkty.setElement(x, y, 0, 'a');
					punkty.setElement(x, y, rozZ, 'a');

					TabP[0][y * ileTabPX + x] = 'b';
					TabP[rozZ][y * ileTabPX + x] = 'b';

					ktory = 1;
					break;

				case 1:

					punkty.setElement(0, y, z, 'a');
					punkty.setElement(rozX, y, z, 'a');

					TabP[z][y * ileTabPX] = 'b';
					TabP[z][y * ileTabPX + rozX] = 'b';

					ktory = 2;
					break;

				case 2:
					punkty.setElement(x, 0, z, 'a');
					punkty.setElement(x, rozY, z, 'a');

					TabP[z][x] = 'b';
					TabP[z][rozY * ileTabPX + x] = 'b';

					ktory = 0;
					break;
			}

			b = a;
			a = 0;
		}

	}

}



void Losowanie::losujKrawedziePer(int dok, int ileTabPX, int ileTabPY, int ileTabPZ, char ** TabP, int odGranicy, int odPunktu, double ax, double b0, int rszuk) {


	double dok1 = dok, dok2, dok3, dok4, dok5, dok6, dok7;
	dok1 = dok + rszuk;
	bool flaga;
	int a = 0, b = 0;
	int ktory = 0;
//double superdok = dok*dok;
	double range0, range, max0, max;

	for (int i = 0, x = 0, y = 0, z = 0, rozX = ileTabPX - 1, rozY = ileTabPY - 1, rozZ = ileTabPZ - 1; a < 100; i++) {


		switch (ktory) {

			case 0:
				x = rand() % ileTabPX;
				y = 0;
				z = 0;
				break;
			case 1:
				y = rand() % ileTabPY;
				z = 0;
				x = 0;
				break;
			case 2:
				z = rand() % ileTabPZ;
				x = 0;
				y = 0;
				break;

		}


		dok2 = y - dok1;
		dok3 = y + dok1;
		dok4 = x - dok1;
		dok5 = x + dok1;
		dok6 = z - dok1;
		dok7 = z + dok1;

		if (dok2 < 0) {
			dok2 = 0;
		}
		if (dok4 < 0) {
			dok4 = 0;
		}
		if (dok6 < 0) {
			dok6 = 0;
		}
		if (dok3 > ileTabPY) {
			dok3 = ileTabPY;
		}
		if (dok5 > ileTabPX) {
			dok5 = ileTabPX;
		}
		if (dok7 > ileTabPZ) {
			dok7 = ileTabPZ;
		}

		range0 = dok * dok * dok;
		max0 = 99999999;

		flaga = true;

		for (int i0 = dok6; i0 < dok7; i0++) {
			for (int i1 = dok2; i1 < dok3; i1++) {
				for (int i2 = dok4; i2 < dok5; i2++) {

					if (TabP[i0][i2 + i1 * ileTabPX] != 'a') {
						max = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odPunktu
						if (max <= odPunktu) {
							flaga = false;
						}
						if (max < max0) {
							max0 = max;
						}
					}

					if (TabP[i0][i2 + i1 * ileTabPX] == 'g') {
						range = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odGranicy
						if (range <= odGranicy) {
							flaga = false;
						}
						if (range < range0) {
							range0 = range;
						}
					}

				}
			}
		}


		if (max0 < range0 * ax + b0 || flaga == false) {
			i--;
			a++;
		}
		else {
			switch (ktory) {

				case 0:
					punkty.setElement(x, 0, 0, 'a');
					punkty.setElement(x, rozY, 0, 'a');
					punkty.setElement(x, 0, rozZ, 'a');
					punkty.setElement(x, rozY, rozZ, 'a');

					TabP[0][x] = 'b';
					TabP[0][rozY * ileTabPX + x] = 'b';
					TabP[rozZ][x] = 'b';
					TabP[rozZ][rozY * ileTabPX + x] = 'b';

					ktory = 1;
					break;

				case 1:

					punkty.setElement(0, y, 0, 'a');
					punkty.setElement(rozX, y, 0, 'a');
					punkty.setElement(0, y, rozZ, 'a');
					punkty.setElement(rozX, y, rozZ, 'a');

					TabP[0][y * ileTabPX] = 'b';
					TabP[0][y * ileTabPX + rozX] = 'b';
					TabP[rozZ][y * ileTabPX] = 'b';
					TabP[rozZ][y * ileTabPX + rozX] = 'b';

					ktory = 2;
					break;

				case 2:

					punkty.setElement(0, 0, z, 'a');
					punkty.setElement(0, rozY, z, 'a');
					punkty.setElement(rozX, 0, z, 'a');
					punkty.setElement(rozX, rozY, z, 'a');

					TabP[z][0] = 'b';
					TabP[z][y * ileTabPX] = 'b';
					TabP[z][rozX] = 'b';
					TabP[z][y * ileTabPX + rozX] = 'b';

					ktory = 0;
					break;
			}

			b = a;
			a = 0;
		}

	}

}

void Losowanie::losujKrawedzieNiePer(int dok, int ileTabPX, int ileTabPY, int ileTabPZ, char ** TabP, int odGranicy, int odPunktu, double ax, double b0, int rszuk) {


	double dok1 = dok, dok2, dok3, dok4, dok5, dok6, dok7;
	dok1 = dok + rszuk;
	bool flaga;
	int a = 0, b = 0;
	int ktory = 0;
//double superdok = dok*dok;
	double range0, range, max0, max;

	for (int i = 0, x = 0, y = 0, z = 0, rozX = ileTabPX - 1, rozY = ileTabPY - 1, rozZ = ileTabPZ - 1; a < 100; i++) {


		switch (ktory) {

			case 0:
				x = rand() % ileTabPX;
				y = 0;
				z = 0;
				break;
			case 1:
				y = rand() % ileTabPY;
				z = 0;
				x = 0;
				break;
			case 2:
				z = rand() % ileTabPZ;
				x = 0;
				y = 0;
				break;
			case 3:
				x = rand() % ileTabPX;
				y = rozY;
				z = 0;
				break;
			case 4:
				y = rand() % ileTabPY;
				z = rozZ;
				x = 0;
				break;
			case 5:
				z = rand() % ileTabPZ;
				x = rozX;
				y = 0;
				break;
			case 6:
				x = rand() % ileTabPX;
				y = 0;
				z = rozZ;
				break;
			case 7:
				y = rand() % ileTabPY;
				z = 0;
				x = rozX;
				break;
			case 8:
				z = rand() % ileTabPZ;
				x = 0;
				y = rozY;
				break;
			case 9:
				x = rand() % ileTabPX;
				y = rozY;
				z = rozZ;
				break;
			case 10:
				y = rand() % ileTabPY;
				z = rozZ;
				x = rozX;
				break;
			case 11:
				z = rand() % ileTabPZ;
				x = rozX;
				y = rozY;
				break;

		}

		dok2 = y - dok1;
		dok3 = y + dok1;
		dok4 = x - dok1;
		dok5 = x + dok1;
		dok6 = z - dok1;
		dok7 = z + dok1;

		if (dok2 < 0) {
			dok2 = 0;
		}
		if (dok4 < 0) {
			dok4 = 0;
		}
		if (dok6 < 0) {
			dok6 = 0;
		}
		if (dok3 > ileTabPY) {
			dok3 = ileTabPY;
		}
		if (dok5 > ileTabPX) {
			dok5 = ileTabPX;
		}
		if (dok7 > ileTabPZ) {
			dok7 = ileTabPZ;
		}

		range0 = dok * dok * dok;
		max0 = 99999999;

		flaga = true;

		for (int i0 = dok6; i0 < dok7; i0++) {
			for (int i1 = dok2; i1 < dok3; i1++) {
				for (int i2 = dok4; i2 < dok5; i2++) {

					if (TabP[i0][i2 + i1 * ileTabPX] != 'a') {
						max = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odPunktu
						if (max <= odPunktu) {
							flaga = false;
						}
						if (max < max0) {
							max0 = max;
						}
					}

					if (TabP[i0][i2 + i1 * ileTabPX] == 'g') {
						range = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odGranicy
						if (range <= odGranicy) {
							flaga = false;
						}
						if (range < range0) {
							range0 = range;
						}
					}

				}
			}
		}


		if (max0 < range0 * ax + b0 || flaga == false) {
			i--;
			a++;
		}
		else {
			switch (ktory) {

				case 0:
					punkty.setElement(x, 0, 0, 'a');
					TabP[0][x] = 'b';
					ktory = 1;
					break;

				case 1:
					punkty.setElement(0, y, 0, 'a');
					TabP[0][y * ileTabPX] = 'b';
					ktory = 2;
					break;

				case 2:
					punkty.setElement(0, 0, z, 'a');
					TabP[z][0] = 'b';
					ktory = 3;
					break;

				case 3:
					punkty.setElement(x, rozY, 0, 'a');
					TabP[0][rozY * ileTabPX + x] = 'b';
					ktory = 4;
					break;

				case 4:
					punkty.setElement(0, y, rozZ, 'a');
					TabP[rozZ][y * ileTabPX] = 'b';
					ktory = 5;
					break;

				case 5:
					punkty.setElement(rozX, 0, z, 'b');
					TabP[z][rozX] = 'a';
					ktory = 6;
					break;

				case 6:
					punkty.setElement(x, 0, rozZ, 'a');
					TabP[rozZ][x] = 'b';
					ktory = 7;
					break;

				case 7:
					punkty.setElement(rozX, y, 0, 'a');
					TabP[0][y * ileTabPX + rozX] = 'b';
					ktory = 8;
					break;

				case 8:
					punkty.setElement(0, rozY, z, 'a');
					TabP[z][y * ileTabPX] = 'b';
					ktory = 9;
					break;

				case 9:
					punkty.setElement(x, rozY, rozZ, 'a');
					TabP[rozZ][rozY * ileTabPX + x] = 'b';
					ktory = 10;
					break;

				case 10:
					punkty.setElement(rozX, y, rozZ, 'a');
					TabP[rozZ][y * ileTabPX + rozX] = 'b';
					ktory = 11;
					break;

				case 11:
					punkty.setElement(rozX, rozY, z, 'a');
					TabP[z][y * ileTabPX + rozX] = 'b';
					ktory = 0;
					break;
			}


			b = a;
			a = 0;
		}

	}

}


void Losowanie::losujWObjetosci(int dok, int ileTabPX, int ileTabPY, int ileTabPZ, char ** TabP, int ileP, int odGranicy, int odPunktu, double ax, double b0, int rszuk) {

	double dok1 = dok, dok2, dok3, dok4, dok5, dok6, dok7;
	dok1 = dok + rszuk;
	bool flaga;
	int a = 0, b = 0;

	double range0, range, max0, max;

	for (int i = 0, x, y, z; i < ileP && a < 200; i++) {


		x = rand() % ileTabPX;
		y = rand() % ileTabPY;
		z = rand() % ileTabPZ;

		dok2 = y - dok1;
		dok3 = y + dok1;
		dok4 = x - dok1;
		dok5 = x + dok1;
		dok6 = z - dok1;
		dok7 = z + dok1;

		if (dok2 < 0) {
			dok2 = 0;
		}
		if (dok4 < 0) {
			dok4 = 0;
		}
		if (dok6 < 0) {
			dok6 = 0;
		}
		if (dok3 > ileTabPY) {
			dok3 = ileTabPY;
		}
		if (dok5 > ileTabPX) {
			dok5 = ileTabPX;
		}
		if (dok7 > ileTabPZ) {
			dok7 = ileTabPZ;
		}

		range0 = dok * dok * dok;
		max0 = 99999999;

		flaga = true;

		for (int i0 = dok6; i0 < dok7; i0++) {
			for (int i1 = dok2; i1 < dok3; i1++) {
				for (int i2 = dok4; i2 < dok5; i2++) {

					if (TabP[i0][i2 + i1 * ileTabPX] != 'a') {
						max = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odPunktu
						if (max <= odPunktu) {
							flaga = false;
						}
						if (max < max0) {
							max0 = max;
						}
					}

					if (TabP[i0][i2 + i1 * ileTabPX] == 'g') {
						range = (i2 - x) * (i2 - x) + (i1 - y) * (i1 - y) + (i0 - z) * (i0 - z);
						//odGranicy
						if (range <= odGranicy) {
							flaga = false;
						}
						if (range < range0) {
							range0 = range;
						}
					}

				}
			}
		}


		if (max0 < range0 * ax + b0 || flaga == false) {
			i--;
			a++;
		}
		else {

			punkty.setElement(x, y, z, 'a');
			TabP[z][y * ileTabPX + x] = 'b';
			b = a;
			a = 0;

		}

	}

}


/*
void Losowanie::rysujPunktySzkieletu(){

glBegin(GL_POINTS);

glColor3f(0,0,0);

	for(int i=0;i<ilePSZ;i++){

	glVertex3f(punkty.getElement(i).getX(),punkty.getElement(i).getY(),punkty.getElement(i).getZ());

	}

glEnd();
}
*/

/*
void Losowanie::rysujPunkty(double maxX,double maxY,double maxZ){

glPushMatrix();

glTranslatef(-maxX,-maxY,-maxZ);

glPushMatrix();
glTranslatef(maxX*0.5,maxY*0.5,maxZ*0.5);

glBegin(GL_POINTS);
glColor3f(1,1,1);

for(int i=0;i<punkty.getIter();i++){

if(punkty.getElement(i).getGranica()=='g'){glColor3f(1,0,0);}
else{glColor3f(1,1,1);}
glVertex3f(punkty.getElement(i).getX(),punkty.getElement(i).getY(),punkty.getElement(i).getZ());

}

glColor3f(1,0,0);
for(int i=0;i<testG.getIter();i++){

glVertex3f(testG.getElement(i).getX(),testG.getElement(i).getY(),testG.getElement(i).getZ());

}

glColor3f(1,1,1);
glEnd();
glPopMatrix();
glPopMatrix();
}
*/


void Losowanie::testGranica(int rozdz) {

	testG.czysc(1000, 1000);

	for (int i = 0; i < wymiarY; i += rozdz) {
		for (int j = 0; j < wymiarX; j += rozdz) {

			testG.setElement(j, i, wymiarZ * 0.5);
			testG.setElement(j, wymiarY * 0.5, i);
			//testG.setElement(wymiar*0.5,i,j);

		}
	}

}

void Losowanie::punktyPoczatkoweSzkielet(int range, int prob) {

	double R = range * range;
	PunktList pS(30000, 30000);
	double x0, y0, z0, x1, y1, z1;
	bool flaga;

	for (int i = 0, ileP = punkty.getIter(); i < ileP; ++i) {
		if (!prob) {
			break;
		}
		--prob;
		flaga = true;
		x0 = punkty.getElement(i).getX();
		y0 = punkty.getElement(i).getY();
		z0 = punkty.getElement(i).getZ();

		for (int j = 0, ilePS = pS.getIter(); j < ilePS; j++) {

			x1 = pS.getElement(j).getX();
			y1 = pS.getElement(j).getY();
			z1 = pS.getElement(j).getZ();

			if ( (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0) < R) {
				flaga = false;
				break;
			}

		}

		if (flaga) {
			pS.setElement(punkty.getElement(i));
			punkty.usunElement(i);
			--i;
			--ileP;
		}

	}


	for (int i = 0, ileP = pS.getIter(); i < ileP; i++) {
		punkty.setElement(pS.getElement(i));
	}

	for (int i = 0, ileP = pS.getIter(), j = punkty.getIter() - 1; i < ileP; ++i, j--) {
		punkty.zamianaMiejsc(i, j);
	}

	ilePSZ = pS.getIter();
}

void Losowanie::wczytajPunktySztuczneSerce(const char * nazwa) {

	double procent = 90;
	double c;
	double x, y, z;
	testG.czysc(28000, 7000);
	punkty.ustawIter(0);
	double minX = 99999999, maxX = -99999999, tranX = 0;
	double minY = 99999999, maxY = -99999999, tranY = 0;
	double minZ = 99999999, maxZ = -99999999, tranZ = 0;
	double sX, sY, sZ;

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);
//6795
//26084

//wczytanie z pliku
	for (int i = 0; i < 6795; ++i) {

		wczytaj >> c;
		wczytaj >> x;
		wczytaj >> y;
		wczytaj >> z;

		punkty.setElement(x, y, z);

		wczytaj >> c;

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

	for (int i = 0; i < 6795; ++i) {
		// punkty.getElement(i).setPunkt(punkty.getElement(i).getX()+tranX,punkty.getElement(i).getY()+tranY,punkty.getElement(i).getZ()+tranZ);
		punkty.getElement(i) = punkty.getElement(i) * skal;
		testG.setElement(punkty.getElement(i));
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

	for (int i = 0, ileP = testG.getIter(); i < ileP; ++i) {

		x = testG.getElement(i).getX();
		y = testG.getElement(i).getY();
		z = testG.getElement(i).getZ();
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


	for (int i = 0, ileP = testG.getIter(); i < ileP; ++i) {
		testG.getElement(i).setPunkt(testG.getElement(i).getX() + tranX, testG.getElement(i).getY() + tranY, testG.getElement(i).getZ() + tranZ, true);
	}
	punkty.ustawIter(0);

}

void Losowanie::ileGranicznych() {

	int gra = 0;
	for (int i = 0; i < punkty.getIter(); ++i) {
		if (punkty.getElement(i).getGranica() == 'g') {
			gra++;
		}
	}
//m1->Lines->Add(gra);

}
