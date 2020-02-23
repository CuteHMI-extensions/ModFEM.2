#include <cstdlib>

#include <modfem/mm_remesh/wektor.h>



void StringList::setElement(string & a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		string * Lista1 = new string[max];
		for (int i = 0; iter < i; ++i) {
			Lista1[i] = Lista[i];
		}
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/* IntList , DoubleList , Punkt , Element , PunktList , ElementList                               */
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              IntList                                           */
////////////////////////////////////////////////////////////////////////////////////////////////////

void IntList::sortowanie() {

	int temp;
	for (int i = 1; i < iter; ++i) {
		for (int j = 0; j < iter - i; ++j) {
			if (Lista[j] < Lista[j + 1]) {
				temp = Lista[j + 1];
				Lista[j + 1] = Lista[j];
				Lista[j] = temp;
			}
		}
	}
}

int IntList::getElement2Kay(const int & key1, const int & key2) {


	for (int i = 0, ile = iter; i < ile; i += 3) {

		if (Lista[i] == key1) {
			if (Lista[i + 1] == key2) {
				return Lista[i + 2];
			}
		}

		if (Lista[i] == key2) {
			if (Lista[i + 1] == key1) {
				return Lista[i + 2];
			}
		}

	}
	return -10;

}

int IntList::funPodWarunki(const int & a1, const int & a2, int & war) {

	for (int i = 0, ile = iter; i < ile; i += 3) {

		if (Lista[i] == a1 && Lista[i + 1] == a2) {
			return Lista[i + 2];
		}
		if (Lista[i] == a2 && Lista[i + 1] == a1) {
			return Lista[i + 2];
		}

	}

	setElement(a1);
	setElement(a2);
	setElement(++war);
	return war;
}

void IntList::dodajIntList(IntList & a) {

	if (iter + a.iter + 1 > max) {
		for (int i = 0, ile = a.iter; i < ile; ++i) {

			Lista[iter++] = a.Lista[i];

		}
	}

	else {

		max = iter + a.iter + grow;
		int * Lista1 = new int[max];
		memcpy ( Lista1, Lista, sizeof(int)*iter );
		delete []Lista;
		Lista = Lista1;

		for (int i = 0, ile = a.iter; i < ile; ++i) {

			Lista[iter++] = a.Lista[i];

		}
	}

}

void IntList::losowyZmianaMiejsca() {

	int temp;

	for (int i = 0, p1, p2, ile = iter; i < ile; i++) {

		p1 = rand() % iter;
		p2 = rand() % iter;

		temp = Lista[p1];
		Lista[p1] = Lista[p2];
		Lista[p2] = temp;

	}

}
void IntList::operator=(IntList & a) {

	iter = a.iter;

	max = a.max;
	grow = a.grow;
	delete []Lista;
	int * Lista = new int[max];
	memcpy ( Lista, a.Lista, sizeof(int)*iter );

}

void IntList::zamienWyzerujIntListy(IntList & a, int grow, int max) {
	int * w, it, ma, gr;
	w = Lista;
	Lista = a.Lista;
	a.Lista = w;
	it = iter;
	ma = max;
	gr = grow;
	iter = a.iter;
	max = a.max;
	grow = a.grow;
	a.iter = iter;
	a.grow = gr;
	a.max = ma;
	a.czysc(grow, max);
}

void IntList::zamienIntListy(IntList & a) {
	int * w, it, ma, gr;
	w = Lista;
	Lista = a.Lista;
	a.Lista = w;

	it = iter;
	ma = max;
	gr = grow;
	iter = a.iter;
	max = a.max;
	grow = a.grow;
	a.iter = iter;
	a.grow = gr;
	a.max = ma;

}




void IntList::wczytajZPliku(const char * nazwa) {

	ifstream wczytaj(nazwa);

	int ileP;
	wczytaj >> ileP;

	czysc(0.2 * (ileP + 2), ileP + 2);

	for (int i = 0; i < ileP; ++i) {

		wczytaj >> Lista[i];
		iter++;
	}

	wczytaj.close();
}

void IntList::zapisDoPliku(const char * nazwa) {

	ofstream zapis(nazwa);

	zapis << iter << endl << endl;

	for (int i = 0; i < iter; i++) {

		zapis << Lista[i] << endl;

	}

	zapis.close();

}


bool IntList::sprCzyJestEleWList(const int & jaki) {
	for (int i = 0; i < iter; i++) {
		if (jaki == Lista[i]) {
			return false;
		}
	}
	return true;
}

bool IntList::sprawdzElement(const int & jaki) {

	if (jaki != -1) {

		for (int i = 0; i < iter; i++) {
			if (jaki == Lista[i]) {
				return false;
			}
		}

		return true;
	}

	return false;
}

void IntList::setElement(const int & a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		int * Lista1 = new int[max];
		memcpy ( Lista1, Lista, sizeof(int)*iter );
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              DoubleList                                        */
////////////////////////////////////////////////////////////////////////////////////////////////////


double DoubleList::getSrednia() {
	double srednia = 0;
	for (int i = 0; i < iter; ++i) {
		srednia += Lista[i];
	}
	return srednia / iter;

}

double DoubleList::getSuma() {

	double suma = 0;
	for (int i = 0; i < iter; ++i) {
		suma += Lista[i];
	}
	return suma;

}

void DoubleList::setElement(const double & a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		double * Lista1 = new double[max];
		memcpy ( Lista1, Lista, sizeof(double)*iter );
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Punkt                                             */
////////////////////////////////////////////////////////////////////////////////////////////////////




float Punkt::odlegloscPunktu(Punkt & a1) {
	return (x - a1.x) * (x - a1.x) + (y - a1.y) * (y - a1.y) + (z - a1.z) * (z - a1.z);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Element                                           */
////////////////////////////////////////////////////////////////////////////////////////////////////


void Element::getNumeryElSciany(int ktoraSciana, int punk[3]) {

	switch (ktoraSciana) {

		case 0:
			punk[0] = punkty[3];
			punk[1] = punkty[1];
			punk[2] = punkty[0];
			break;
		case 1:
			punk[0] = punkty[3];
			punk[1] = punkty[2];
			punk[2] = punkty[1];
			break;
		case 2:
			punk[0] = punkty[3];
			punk[1] = punkty[0];
			punk[2] = punkty[2];
			break;
		case 3:
			punk[0] = punkty[0];
			punk[1] = punkty[1];
			punk[2] = punkty[2];
			break;
	}

}

void Element::getNumeryElSciany(char ktoraSciana, int punk[3]) {

	switch (ktoraSciana) {

		case 'a':
			punk[0] = punkty[3];
			punk[1] = punkty[1];
			punk[2] = punkty[0];
			break;
		case 'b':
			punk[0] = punkty[3];
			punk[1] = punkty[2];
			punk[2] = punkty[1];
			break;
		case 'c':
			punk[0] = punkty[3];
			punk[1] = punkty[0];
			punk[2] = punkty[2];
			break;
		case 'd':
			punk[0] = punkty[0];
			punk[1] = punkty[1];
			punk[2] = punkty[2];
			break;
	}

}
/*
 void Element::getNumeryElSciany(char ktoraSciana,int punk[3]){

	switch(ktoraSciana){

		case 'a':punk[0] = punkty[3];punk[1] = punkty[0];punk[2] = punkty[1];break;
		case 'b':punk[0] = punkty[3];punk[1] = punkty[1];punk[2] = punkty[2];break;
		case 'c':punk[0] = punkty[3];punk[1] = punkty[2];punk[2] = punkty[0];break;
		case 'd':punk[0] = punkty[0];punk[1] = punkty[2];punk[2] = punkty[1];break;
	}

}
*/
bool Element::sprawdzCzyJestTakiPunkt(int nrPunktu) {
	if (nrPunktu == punkty[0] || nrPunktu == punkty[1] || nrPunktu == punkty[2] || nrPunktu == punkty[3]) {
		return true;
	}
	return false;
}

bool Element::sprawdzSciane(int numerS, int p1, int p2, int p3) {

	switch (numerS) {

		case 0:
			if ((p1 == punkty[3] || p1 == punkty[1] || p1 == punkty[0]) &&
					(p2 == punkty[3] || p2 == punkty[1] || p2 == punkty[0]) &&
					(p3 == punkty[3] || p3 == punkty[1] || p3 == punkty[0])) {
				return true;
			}
			break;

		case 1:
			if ((p1 == punkty[3] || p1 == punkty[2] || p1 == punkty[1]) &&
					(p2 == punkty[3] || p2 == punkty[2] || p2 == punkty[1]) &&
					(p3 == punkty[3] || p3 == punkty[2] || p3 == punkty[1])) {
				return true;
			}
			break;

		case 2:
			if ((p1 == punkty[3] || p1 == punkty[0] || p1 == punkty[2]) &&
					(p2 == punkty[3] || p2 == punkty[0] || p2 == punkty[2]) &&
					(p3 == punkty[3] || p3 == punkty[0] || p3 == punkty[2])) {
				return true;
			}
			break;

		case 3:
			if ((p1 == punkty[0] || p1 == punkty[1] || p1 == punkty[2]) &&
					(p2 == punkty[0] || p2 == punkty[1] || p2 == punkty[2]) &&
					(p3 == punkty[0] || p3 == punkty[1] || p3 == punkty[2])) {
				return true;
			}
			break;

	}

	return false;

}


int Element::getNumerSciany(int p1, int p2, int p3) {

	int sum = p1 + p2 + p3;
	if (sum == punkty[3] + punkty[1] + punkty[0]) {
		return 0;
	}
	else if (sum == punkty[3] + punkty[2] + punkty[1]) {
		return 1;
	}
	else if (sum == punkty[3] + punkty[0] + punkty[2]) {
		return 2;
	}
	else if (sum == punkty[0] + punkty[1] + punkty[2]) {
		return 3;
	}
	else {
		return -1;
	}

}

bool Element::podmienPunkt(const int & staryNumer, const int & nowyNumer) {

	if (punkty[0] == staryNumer) {
		punkty[0] = nowyNumer;
		return true;
	}
	else if (punkty[1] == staryNumer) {
		punkty[1] = nowyNumer;
		return true;
	}
	else if (punkty[2] == staryNumer) {
		punkty[2] = nowyNumer;
		return true;
	}
	else if (punkty[3] == staryNumer) {
		punkty[3] = nowyNumer;
		return true;
	}

	return false;

}


void Element::setElement(int p1, int p2, int p3, int p4, double Sx1, double Sy1, double Sz1, float R1, double Ax1, double Ay1, double Az1) {


	punkty[0] = p1;
	punkty[1] = p2;
	punkty[2] = p3;
	punkty[3] = p4;


	Sx = Sx1;
	Sy = Sy1;
	Sz = Sz1;
	R = R1;
	elementy[0] = -1;
	elementy[1] = -1;
	elementy[2] = -1;
	elementy[3] = -1;
	Ax = Ax1;
	Ay = Ay1;
	Az = Az1;

}

void Element::setElement(int p1, int p2, int p3, int p4, double Sx1, double Sy1, double Sz1, float R1, double Ax1, double Ay1, double Az1, float StosunekR_r) {


	punkty[0] = p1;
	punkty[1] = p2;
	punkty[2] = p3;
	punkty[3] = p4;


	Sx = Sx1;
	Sy = Sy1;
	Sz = Sz1;
	R = R1;
	elementy[0] = -1;
	elementy[1] = -1;
	elementy[2] = -1;
	elementy[3] = -1;
	Ax = Ax1;
	Ay = Ay1;
	Az = Az1;
	StosR_r = StosunekR_r;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Pryzmy                                            */
////////////////////////////////////////////////////////////////////////////////////////////////////


bool Pryzm::podmienPunkt(const int & staryNumer, const int & nowyNumer) {

	if (punkty[0] == staryNumer) {
		punkty[0] = nowyNumer;
		return true;
	}
	else if (punkty[1] == staryNumer) {
		punkty[1] = nowyNumer;
		return true;
	}
	else if (punkty[2] == staryNumer) {
		punkty[2] = nowyNumer;
		return true;
	}
	else if (punkty[3] == staryNumer) {
		punkty[3] = nowyNumer;
		return true;
	}
	else if (punkty[4] == staryNumer) {
		punkty[4] = nowyNumer;
		return true;
	}
	else if (punkty[5] == staryNumer) {
		punkty[5] = nowyNumer;
		return true;
	}

	return false;

}



////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              PunktyList                                        */
////////////////////////////////////////////////////////////////////////////////////////////////////


void PunktList::copyPunktList(PunktList & a) {

	iter = a.iter;
	grow = a.grow;
	max = a.max;

	delete []Lista;
	Lista = new Punkt[max];

	memcpy ( Lista, a.Lista, sizeof(Punkt) * (iter + 1) );

}

void PunktList::losowyZmianaMiejsca() {

	Punkt temp;

	for (int i = 0, p1, p2, ile = iter; i < ile; i++) {

		p1 = rand() % iter;
		p2 = rand() % iter;

		temp = Lista[p1];
		Lista[p1] = Lista[p2];
		Lista[p2] = temp;

	}

}

void PunktList::zmianaKolejnosci() {
	Punkt temp;

	for (int i = 0, ile = iter * 0.5; i < ile; i++) {

		temp = Lista[i + iter - 1 - i];
		Lista[i + iter - 1 - i] = Lista[i];
		Lista[i] = temp;

	}

}

void PunktList::setElement(Punkt & a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		Punkt * Lista1 = new Punkt[max + 1];
		memcpy ( Lista1, Lista, sizeof(Punkt)*iter );
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void PunktList::setElement(const double & x1, const double & y1, const double & z1) {
	if (max > iter) {
		Lista[iter].setPunkt(x1, y1, z1);
		iter++;
	}
	else {
		max = max + grow;
		Punkt * Lista1 = new Punkt[max + 1];
		memcpy ( Lista1, Lista, sizeof(Punkt)*iter );
		Lista1[iter].setPunkt(x1, y1, z1);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}



void PunktList::setElement(const double & x1, const double & y1, const double & z1, const char & g) {
	if (max > iter) {
		Lista[iter].setPunkt(x1, y1, z1, g);
		iter++;
	}
	else {
		max = max + grow;
		Punkt * Lista1 = new Punkt[max + 1];
		memcpy ( Lista1, Lista, sizeof(Punkt)*iter );
		Lista1[iter].setPunkt(x1, y1, z1, g);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

bool PunktList::setElementUniqat(const double & x1, const double & y1, const double & z1, const char & g) {

	Punkt temp;
	temp.setPunkt(x1, y1, z1, 'g');
	bool flaga = true;
	for (int i = 0; i < iter; ++i) {
		if (Lista[i] == temp) {
			flaga = false;
		}
	}

	if (flaga) {

		if (max > iter) {
			Lista[iter].setPunkt(x1, y1, z1, g);
			iter++;
		}
		else {
			max = max + grow;
			Punkt * Lista1 = new Punkt[max + 1];
			memcpy ( Lista1, Lista, sizeof(Punkt)*iter );
			delete []Lista;

			Lista = Lista1;
			Lista[iter].setPunkt(x1, y1, z1, g);
			iter++;
		}

		return true;
	}

	return false;

}


void PunktList::zapisDoPliku(const char * nazwa) {

	ofstream zapis(nazwa);
	zapis.precision(20);

	zapis << iter << endl << endl;

	for (int i = 0; i < iter; i++) {

		zapis << i + 1 << "  " << Lista[i].getX() << "  " << Lista[i].getY() << "  " << Lista[i].getZ() << "  " << Lista[i].getGranica() << endl;

	}

	zapis.close();
}

void PunktList::wczytajZPliku(const char * nazwa, bool tylkoGranica) {

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);
	string text;

	double xW, yW, zW;
	char granicaW;
	int ileP;

	wczytaj >> ileP;



	if (tylkoGranica) {

		czysc(int(ileP * 0.6), int(ileP * 0.6));

		for (int i = 0; i < ileP; i++) {
			wczytaj >> text;

			wczytaj >> xW;
			wczytaj >> yW;
			wczytaj >> zW;
			wczytaj >> granicaW;


			if (granicaW == 'g') {
				setElement(xW, yW, xW, 'g');
			}
		}

	}
	else {

		czysc(int(ileP * 1.1), int(ileP * 1.1));


		for (int i = 0; i < ileP; i++) {

			wczytaj >> xW;

			wczytaj >> xW;
			wczytaj >> yW;
			wczytaj >> zW;
			wczytaj >> granicaW;

			if (granicaW == 'g') {
				setElement(xW, yW, zW, 'g');
			}
			else {
				setElement(xW, yW, zW, 'a');
			}

		}

	}

	wczytaj.close();


}



void PunktList::wczytajZPlikuPrzerob(const char * nazwa) {

	ifstream wczytaj(nazwa);
	wczytaj.precision(20);
	string text;

	double xW, yW, zW, temp;

	while (1) {

		wczytaj >> text;
		if (text != "GRID") {
			break;
		}
		wczytaj >> xW;
		wczytaj >> xW;
		wczytaj >> yW;
		wczytaj >> zW;
		wczytaj >> temp;

		setElement(xW * 1.5434706, yW * 1.5434706, zW * 1.5434706);

	}

	wczytaj.close();

	ofstream zapis("nowePolozenieP.txt");
	zapis.precision(12);


	for (int i = 0; i < iter; i++) {

		zapis << "GRID" << "      " << i + 1 << "      " << Lista[i].getX() << "      " << Lista[i].getY() << "      " << Lista[i].getZ() << "  " << 0 << endl;

	}

	zapis.close();

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              ElementyList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////



void ElementList::ustawTempZtabB(int * tabBlock, int dl_tab) {


	for (int i = 0; i < iter; ++i) {
		Lista[i].setTemp(0);
	}

	for (int i = 0, nrZ; i < iter; ++i) {

		nrZ = Lista[i].getRodzajZiarna();

		for (int ii = 0; ii < dl_tab; ii += 2) {

			if (tabBlock[ii] == nrZ) {
				Lista[i].setTemp(tabBlock[ii + 1]);
				break;
			}

		}


	}

}

void ElementList::ustawTempnaZiarna() {

	for (int i = 0; i < iter; ++i) {
		Lista[i].ustawTempnaZiarna();
	}

}

int ElementList::getNumerSasSacianE1(int nrE) {

	int e = Lista[nrE].getE1();
	if (e < 0) {
		return -1;
	}
	//std::cout <<e<< " Method not implemented!: mmr_fa_neig " << nrE <<endl;
	return Lista[e].getNumerSciany(Lista[nrE].getP4(), Lista[nrE].getP2(), Lista[nrE].getP1());
}

int ElementList::getNumerSasSacianE2(int nrE) {
	int e = Lista[nrE].getE2();
	if (e < 0) {
		return -1;
	}
	return Lista[e].getNumerSciany(Lista[nrE].getP4(), Lista[nrE].getP3(), Lista[nrE].getP2());
}

int ElementList::getNumerSasSacianE3(int nrE) {
	int e = Lista[nrE].getE3();
	if (e < 0) {
		return -1;
	}
	return Lista[e].getNumerSciany(Lista[nrE].getP4(), Lista[nrE].getP1(), Lista[nrE].getP3());
}

int ElementList::getNumerSasSacianE4(int nrE) {
	int e = Lista[nrE].getE4();
	if (e < 0) {
		return -1;
	}
	return Lista[e].getNumerSciany(Lista[nrE].getP1(), Lista[nrE].getP2(), Lista[nrE].getP3());
}

int ElementList::getNumerSasSacianWybor(int nrE, int nrS) {

	int p1, p2, p3, p4;

	p1 = Lista[nrE].getP1();
	p2 = Lista[nrE].getP2();
	p3 = Lista[nrE].getP3();
	p4 = Lista[nrE].getP4();


	switch (nrS) {
		case 0: {
			return Lista[Lista[nrE].getE1()].getNumerSciany(p4, p2, p1);
		}
		break;
		case 1: {
			return Lista[Lista[nrE].getE2()].getNumerSciany(p4, p3, p2);
		}
		break;
		case 2: {
			return Lista[Lista[nrE].getE3()].getNumerSciany(p4, p1, p3);
		}
		break;
		case 3: {
			return Lista[Lista[nrE].getE4()].getNumerSciany(p1, p2, p3);
		}
		break;
	}

	return -1;
}


void ElementList::ustawSasiadow(const int & numerStary, const int & numerNowy) {

	for (int i = 0; i < iter; ++i) {

		if (Lista[i].getE1() == numerStary) {
			Lista[i].setE1(numerNowy);
		}
		if (Lista[i].getE2() == numerStary) {
			Lista[i].setE2(numerNowy);
		}
		if (Lista[i].getE3() == numerStary) {
			Lista[i].setE3(numerNowy);
		}
		if (Lista[i].getE4() == numerStary) {
			Lista[i].setE4(numerNowy);
		}

	}

}

void ElementList::podmienPunktWElementach(const int & staryNumer, const int & nowyNumer, IntList & elementy) {

	for (int i = 0, e, ileE = elementy.getIter(); i < ileE; ++i) {

		e = elementy.getElement(i);
		Lista[e].podmienPunkt(staryNumer, nowyNumer);

	}

}

void ElementList::ustawGraniceElementowZaz( IntList & nrZazZiaren) {



	for (int i = 0, e1, e2, e3, e4; i < iter; ++i) {

		if (nrZazZiaren.sprawdzCzyJest(Lista[i].getRodzajZiarna())) {

			e1 = Lista[i].getE1();
			e2 = Lista[i].getE2();
			e3 = Lista[i].getE3();
			e4 = Lista[i].getE4();

			if (e1 > -1 && !nrZazZiaren.sprawdzCzyJest(Lista[e1].getRodzajZiarna())) {
				Lista[i].setE1(-1);
			}
			if (e2 > -1 && !nrZazZiaren.sprawdzCzyJest(Lista[e2].getRodzajZiarna())) {
				Lista[i].setE2(-1);
			}
			if (e3 > -1 && !nrZazZiaren.sprawdzCzyJest(Lista[e3].getRodzajZiarna())) {
				Lista[i].setE3(-1);
			}
			if (e4 > -1 && !nrZazZiaren.sprawdzCzyJest(Lista[e4].getRodzajZiarna())) {
				Lista[i].setE4(-1);
			}


		}

	}

}

void ElementList::ustawGraniceElementowNieZaz( IntList & nrZazZiaren) {

	for (int i = 0, e1, e2, e3, e4; i < iter; ++i) {

		if (!nrZazZiaren.sprawdzCzyJest(Lista[i].getRodzajZiarna())) {

			e1 = Lista[i].getE1();
			e2 = Lista[i].getE2();
			e3 = Lista[i].getE3();
			e4 = Lista[i].getE4();

			if (e1 > -1 && nrZazZiaren.sprawdzCzyJest(Lista[e1].getRodzajZiarna())) {
				Lista[i].setE1(-1);
			}
			if (e2 > -1 && nrZazZiaren.sprawdzCzyJest(Lista[e2].getRodzajZiarna())) {
				Lista[i].setE2(-1);
			}
			if (e3 > -1 && nrZazZiaren.sprawdzCzyJest(Lista[e3].getRodzajZiarna())) {
				Lista[i].setE3(-1);
			}
			if (e4 > -1 && nrZazZiaren.sprawdzCzyJest(Lista[e4].getRodzajZiarna())) {
				Lista[i].setE4(-1);
			}


		}

	}

}



void ElementList::usunElementyNieZaz( IntList & nrZazZiaren) {

	for (int i = 0, e1, e2, e3, e4; i < iter; ++i) {
		if (!nrZazZiaren.sprawdzCzyJest(Lista[i].getRodzajZiarna())) {
			--iter;

			if (nrZazZiaren.sprawdzCzyJest(Lista[iter].getRodzajZiarna())) {

				Lista[i] = Lista[iter];

				e1 = Lista[i].getE1();
				e2 = Lista[i].getE2();
				e3 = Lista[i].getE3();
				e4 = Lista[i].getE4();

				if (e1 > -1) {
					if (iter == Lista[e1].getE1()) {
						Lista[e1].setE1(i);
					}
					else if (iter == Lista[e1].getE2()) {
						Lista[e1].setE2(i);
					}
					else if (iter == Lista[e1].getE3()) {
						Lista[e1].setE3(i);
					}
					else if (iter == Lista[e1].getE4()) {
						Lista[e1].setE4(i);
					}
				}
				if (e2 > -1) {
					if (iter == Lista[e2].getE1()) {
						Lista[e2].setE1(i);
					}
					else if (iter == Lista[e2].getE2()) {
						Lista[e2].setE2(i);
					}
					else if (iter == Lista[e2].getE3()) {
						Lista[e2].setE3(i);
					}
					else if (iter == Lista[e2].getE4()) {
						Lista[e2].setE4(i);
					}
				}
				if (e3 > -1) {
					if (iter == Lista[e3].getE1()) {
						Lista[e3].setE1(i);
					}
					else if (iter == Lista[e3].getE2()) {
						Lista[e3].setE2(i);
					}
					else if (iter == Lista[e3].getE3()) {
						Lista[e3].setE3(i);
					}
					else if (iter == Lista[e3].getE4()) {
						Lista[e3].setE4(i);
					}
				}
				if (e4 > -1) {
					if (iter == Lista[e4].getE1()) {
						Lista[e4].setE1(i);
					}
					else if (iter == Lista[e4].getE2()) {
						Lista[e4].setE2(i);
					}
					else if (iter == Lista[e4].getE3()) {
						Lista[e4].setE3(i);
					}
					else if (iter == Lista[e4].getE4()) {
						Lista[e4].setE4(i);
					}
				}

			}

			--i;

		}
	}

}


void ElementList::usunElementyZaz( IntList & nrZazZiaren) {

	for (int i = 0, e1, e2, e3, e4; i < iter; ++i) {
		if (nrZazZiaren.sprawdzCzyJest(Lista[i].getRodzajZiarna())) {

			Lista[i] = Lista[--iter];

			if (!nrZazZiaren.sprawdzCzyJest(Lista[iter].getRodzajZiarna())) {

				e1 = Lista[i].getE1();
				e2 = Lista[i].getE2();
				e3 = Lista[i].getE3();
				e4 = Lista[i].getE4();

				if (e1 > -1) {
					;
					if (iter == Lista[e1].getE1()) {
						Lista[e1].setE1(i);
					}
					else if (iter == Lista[e1].getE2()) {
						Lista[e1].setE2(i);
					}
					else if (iter == Lista[e1].getE3()) {
						Lista[e1].setE3(i);
					}
					else if (iter == Lista[e1].getE4()) {
						Lista[e1].setE4(i);
					}
				}
				if (e2 > -1) {
					if (iter == Lista[e2].getE1()) {
						Lista[e2].setE1(i);
					}
					else if (iter == Lista[e2].getE2()) {
						Lista[e2].setE2(i);
					}
					else if (iter == Lista[e2].getE3()) {
						Lista[e2].setE3(i);
					}
					else if (iter == Lista[e2].getE4()) {
						Lista[e2].setE4(i);
					}
				}
				if (e3 > -1) {
					if (iter == Lista[e3].getE1()) {
						Lista[e3].setE1(i);
					}
					else if (iter == Lista[e3].getE2()) {
						Lista[e3].setE2(i);
					}
					else if (iter == Lista[e3].getE3()) {
						Lista[e3].setE3(i);
					}
					else if (iter == Lista[e3].getE4()) {
						Lista[e3].setE4(i);
					}
				}
				if (e4 > -1) {
					if (iter == Lista[e4].getE1()) {
						Lista[e4].setE1(i);
					}
					else if (iter == Lista[e4].getE2()) {
						Lista[e4].setE2(i);
					}
					else if (iter == Lista[e4].getE3()) {
						Lista[e4].setE3(i);
					}
					else if (iter == Lista[e4].getE4()) {
						Lista[e4].setE4(i);
					}
				}

			}
			--i;

		}
	}

}





void ElementList::deleteElementPodmienS(int ktory) {

	int przes = --iter;



	if (przes > ktory) {
		int e1, e2, e3, e4;

		e1 = Lista[przes].getE1();
		e2 = Lista[przes].getE2();
		e3 = Lista[przes].getE3();
		e4 = Lista[przes].getE4();

		if (e1 > -1) {
			if (Lista[e1].getE1() == przes) {
				Lista[e1].setE1(ktory);
			}
			else if (Lista[e1].getE2() == przes) {
				Lista[e1].setE2(ktory);
			}
			else if (Lista[e1].getE3() == przes) {
				Lista[e1].setE3(ktory);
			}
			else if (Lista[e1].getE4() == przes) {
				Lista[e1].setE4(ktory);
			}
		}

		if (e2 > -1) {
			if (Lista[e2].getE1() == przes) {
				Lista[e2].setE1(ktory);
			}
			else if (Lista[e2].getE2() == przes) {
				Lista[e2].setE2(ktory);
			}
			else if (Lista[e2].getE3() == przes) {
				Lista[e2].setE3(ktory);
			}
			else if (Lista[e2].getE4() == przes) {
				Lista[e2].setE4(ktory);
			}
		}

		if (e3 > -1) {
			if (Lista[e3].getE1() == przes) {
				Lista[e3].setE1(ktory);
			}
			else if (Lista[e3].getE2() == przes) {
				Lista[e3].setE2(ktory);
			}
			else if (Lista[e3].getE3() == przes) {
				Lista[e3].setE3(ktory);
			}
			else if (Lista[e3].getE4() == przes) {
				Lista[e3].setE4(ktory);
			}
		}

		if (e4 > -1) {
			if (Lista[e4].getE1() == przes) {
				Lista[e4].setE1(ktory);
			}
			else if (Lista[e4].getE2() == przes) {
				Lista[e4].setE2(ktory);
			}
			else if (Lista[e4].getE3() == przes) {
				Lista[e4].setE3(ktory);
			}
			else if (Lista[e4].getE4() == przes) {
				Lista[e4].setE4(ktory);
			}
		}



		Lista[ktory] = Lista[przes];

	}

}

void ElementList::setElement(Element & a) {
	if (max > iter) {

		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		Element * Lista1 = new Element[max];
		memcpy ( Lista1, Lista, sizeof(Element)*iter );
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void ElementList::setElement(const int & p1, const int & p2, const int & p3, const int & p4) {
	if (max > iter) {
		Lista[iter].setElement(p1, p2, p3, p4);
		iter++;
	}
	else {
		max = max + grow;
		Element * Lista1 = new Element[max];
		memcpy ( Lista1, Lista, sizeof(Element)*iter );
		Lista1[iter].setElement(p1, p2, p3, p4);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void ElementList::setElement(const int & p1, const int & p2, const int & p3, const int & p4, const int & rodzZ) {
	if (max > iter) {
		Lista[iter].setElement(p1, p2, p3, p4, rodzZ);
		iter++;
	}
	else {
		max = max + grow;
		Element * Lista1 = new Element[max];
		memcpy ( Lista1, Lista, sizeof(Element)*iter );
		Lista1[iter].setElement(p1, p2, p3, p4, rodzZ);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void ElementList::setElement(const int & p1, const int & p2, const int & p3, const int & p4, const double & Sx1, const double & Sy1, const double & Sz1, const float & R1, const double & Ax1, const double & Ay1, const double & Az1) {
	if (max > iter) {
		Lista[iter].setElement(p1, p2, p3, p4, Sx1, Sy1, Sz1, R1, Ax1, Ay1, Az1);
		iter++;
	}
	else {
		max = max + grow;
		Element * Lista1 = new Element[max];
		memcpy ( Lista1, Lista, sizeof(Element)*iter );
		Lista1[iter].setElement(p1, p2, p3, p4, Sx1, Sy1, Sz1, R1, Ax1, Ay1, Az1);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void ElementList::setElement(const int & p1, const int & p2, const int & p3, const int & p4, const double & Sx1, const double & Sy1, const double & Sz1, const float & R1, const double & Ax1, const double & Ay1, const double & Az1, const float & R_r) {
	if (max > iter) {
		Lista[iter].setElement(p1, p2, p3, p4, Sx1, Sy1, Sz1, R1, Ax1, Ay1, Az1, R_r);
		iter++;
	}
	else {
		max = max + grow;
		Element * Lista1 = new Element[max];
		memcpy ( Lista1, Lista, sizeof(Element)*iter );
		Lista1[iter].setElement(p1, p2, p3, p4, Sx1, Sy1, Sz1, R1, Ax1, Ay1, Az1, R_r);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              PryzmList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////


void PryzmList::setElement(Pryzm & a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		Pryzm * Lista1 = new Pryzm[max];
		memcpy ( Lista1, Lista, sizeof(Pryzm)*iter );
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void PryzmList::setElement(const int & p1, const int & p2, const int & p3, const int & p4, const int & p5, const int & p6) {
	if (max > iter) {
		Lista[iter].setElement(p1, p2, p3, p4, p5, p6);
		iter++;
	}
	else {
		max = max + grow;
		Pryzm * Lista1 = new Pryzm[max];
		memcpy ( Lista1, Lista, sizeof(Pryzm)*iter );
		Lista1[iter].setElement(p1, p2, p3, p4, p5, p6);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                         FaceList                                               */
////////////////////////////////////////////////////////////////////////////////////////////////////


void FaceList::setElement(Face & a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		Face * Lista1 = new Face[max];
		memcpy ( Lista1, Lista, sizeof(Face)*iter );
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void FaceList::setElement(const int & ed1, const int & ed2, const int & ed3, const char & faceInEl1, const int & el) {
	if (max > iter) {
		Lista[iter].setFace(ed1, ed2, ed3, faceInEl1, el);
		iter++;
	}
	else {
		max = max + grow;
		Face * Lista1 = new Face[max];
		memcpy ( Lista1, Lista, sizeof(Face)*iter );
		Lista1[iter].setFace(ed1, ed2, ed3, faceInEl1, el);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void FaceList::setElement(const int & ed1, const int & ed2, const int & ed3, const char & faceInEl1, const int & el, const int & bc) {
	if (max > iter) {
		Lista[iter].setFace(ed1, ed2, ed3, faceInEl1, el, bc);
		iter++;
	}
	else {
		max = max + grow;
		Face * Lista1 = new Face[max];
		memcpy ( Lista1, Lista, sizeof(Face)*iter );
		Lista1[iter].setFace(ed1, ed2, ed3, faceInEl1, el, bc);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              EdgeList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////




void EdgeList::setElement(Edge & a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		Edge * Lista1 = new Edge[max];
		memcpy ( Lista1, Lista, sizeof(Edge)*iter );
		Lista1[iter] = a;
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

void EdgeList::setElement(const int & p1, const int & p2, const int & el1) {
	if (max > iter) {
		Lista[iter].setEdge(p1, p2, el1);
		iter++;
	}
	else {
		max = max + grow;
		Edge * Lista1 = new Edge[max];
		memcpy ( Lista1, Lista, sizeof(Edge)*iter );
		Lista1[iter].setEdge(p1, p2, el1);
		delete []Lista;

		Lista = Lista1;

		iter++;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              TemplateList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////


//template<typename TYP> void TLista<TYP>::setElement(TYP a){
//        if(max>iter){
//        Lista[iter]=a;
//        iter++;
//        }
//        else{
//        max = max+grow;
//        TYP *Lista1 = new TYP[max];
//        memcpy ( Lista1 , Lista , sizeof(TYP)*iter );
//        delete []Lista;

//        Lista = Lista1;
//        Lista[iter] = a;
//        iter++;
//        }
//}

/*
void test::rysuj(){

glBegin(GL_POINTS); //rysowanie prymitywu
glColor3f(0,0,0);

		for(int i=0;i<punkty.getIter();i+=2){
		glVertex3f(punkty.getElement(i),punkty.getElement(i+1),0);
		}

glEnd(); //koniec rysowania

}
*/


