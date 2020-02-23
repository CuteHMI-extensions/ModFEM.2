#include <modfem/mm_remesh/ziarno.h>

#define _USE_MATH_DEFINES
#include <math.h>

void ZiarnoList::obrocZiarno(int numerZiarna, PunktList & points, ElementList & elements, int kierunekOsi, double wartosc) {

	for (int i = 0; i < iter; ++i) {
		if (numerZiarna == Lista[i].pobierzRodzaj()) {
			numerZiarna = i;
			break;
		}
	}

	int ileP = points.getIter();
	bool * czyJestZaz = new bool[ileP];

	for (int i = 0; i < ileP; ++i) {
		czyJestZaz[i] = false;
	}

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		if (numerZiarna == elements.getElement(i).getRodzajZiarna()) {

			czyJestZaz[elements.getElement(i).getP1()] = true;
			czyJestZaz[elements.getElement(i).getP2()] = true;
			czyJestZaz[elements.getElement(i).getP3()] = true;
			czyJestZaz[elements.getElement(i).getP4()] = true;

		}

	}


	int zlicz = 0;

	for (int i = 0; i < ileP; ++i) {
		if (czyJestZaz[i]) {
			++zlicz;
		}
	}

	IntList punktyZiarna(zlicz * 0.1, zlicz + 10);

	double maxX = 0, minX = 999999999, maxY = 0, minY = 999999999, maxZ = 0, minZ = 999999999, x, y, z, x0, y0, z0, pX, pY, pZ;
	double alfa = wartosc * (M_PI / 180);


	for (int i = 0; i < ileP; ++i) {
		if (czyJestZaz[i]) {

			punktyZiarna.setElement(i);

			x = points.getElement(i).getX();
			y = points.getElement(i).getY();
			z = points.getElement(i).getZ();

			if (x > maxX) {
				maxX = x;
			}
			if (x < minX) {
				minX = x;
			}
			if (y > maxY) {
				maxY = y;
			}
			if (y < minY) {
				minY = y;
			}
			if (z > maxZ) {
				maxZ = z;
			}
			if (z < minZ) {
				minZ = z;
			}

		}
	}

	delete []czyJestZaz;
	x0 = (maxX + minX) * 0.5;
	y0 = (maxY + minY) * 0.5;
	z0 = (maxZ + minZ) * 0.5;

	switch (kierunekOsi) {

		case 1:

			for (int i = 0, p, ileP = punktyZiarna.getIter(); i < ileP; ++i) {

				p = punktyZiarna.getElement(i);
				x = points.getElement(p).getX();
				y = points.getElement(p).getY();
				z = points.getElement(p).getZ();

				pX = x0 + (x - x0) * cos(alfa) + (y - y0) * sin(alfa);
				pY = y0 + (y - y0) * cos(alfa) - (x - x0) * sin(alfa);

				points.getElement(p).setXYZ(pX, pY, z);

			}

			break;

		case 2:

			for (int i = 0, p, ileP = punktyZiarna.getIter(); i < ileP; ++i) {

				p = punktyZiarna.getElement(i);
				x = points.getElement(p).getX();
				y = points.getElement(p).getY();
				z = points.getElement(p).getZ();

				pX = x0 + (x - x0) * cos(alfa) + (z - z0) * sin(alfa);
				pZ = z0 + (z - z0) * cos(alfa) - (x - x0) * sin(alfa);

				points.getElement(p).setXYZ(pX, y, pZ);

			}

			break;

		case 3:

			for (int i = 0, p, ileP = punktyZiarna.getIter(); i < ileP; ++i) {

				p = punktyZiarna.getElement(i);
				x = points.getElement(p).getX();
				y = points.getElement(p).getY();
				z = points.getElement(p).getZ();

				pY = y0 + (y - y0) * cos(alfa) + (z - z0) * sin(alfa);
				pZ = z0 + (z - z0) * cos(alfa) - (y - y0) * sin(alfa);

				points.getElement(p).setXYZ(x, pY, pZ);

			}

			break;

	}
}

void ZiarnoList::przesunZiarno(int numerZiarna, PunktList & points, ElementList & elements, int kierunek, double wartosc) {

	for (int i = 0; i < iter; ++i) {
		if (numerZiarna == Lista[i].pobierzRodzaj()) {
			numerZiarna = i;
			break;
		}
	}

	int ileP = points.getIter();
	bool * czyJestZaz = new bool[ileP];

	for (int i = 0; i < ileP; ++i) {
		czyJestZaz[i] = false;
	}

	for (int i = 0, ileE = elements.getIter(); i < ileE; ++i) {

		if (numerZiarna == elements.getElement(i).getRodzajZiarna()) {

			czyJestZaz[elements.getElement(i).getP1()] = true;
			czyJestZaz[elements.getElement(i).getP2()] = true;
			czyJestZaz[elements.getElement(i).getP3()] = true;
			czyJestZaz[elements.getElement(i).getP4()] = true;

		}

	}

	double Px = 0, Py = 0, Pz = 0;

	switch (kierunek) {

		case 0:
			Px += wartosc;
			break;
		case 1:
			Px -= wartosc;
			break;
		case 2:
			Py += wartosc;
			break;
		case 3:
			Py -= wartosc;
			break;
		case 4:
			Pz += wartosc;
			break;
		case 5:
			Pz -= wartosc;
			break;

	}

	for (int i = 0; i < ileP; ++i) {

		if (czyJestZaz[i]) {
			points.getElement(i).przesunPunkt(Px, Py, Pz);
		}

	}


}


void ZiarnoList::setElement(int rodzaj, float r, float g, float b) {

	if (max > iter) {
		Lista[iter].ustalZiarno(rodzaj, r, g, b);
		iter++;
	}
	else {
		max = max + grow;
		Ziarno * Lista1 = new Ziarno[max];
		memcpy ( Lista1, Lista, sizeof(Ziarno)*iter );
		delete []Lista;

		Lista = Lista1;
		Lista[iter].ustalZiarno(rodzaj, r, g, b);
		iter++;
	}

}

void ZiarnoList::setElement(Ziarno a) {
	if (max > iter) {
		Lista[iter] = a;
		iter++;
	}
	else {
		max = max + grow;
		Ziarno * Lista1 = new Ziarno[max];
		memcpy ( Lista1, Lista, sizeof(Ziarno)*iter );
		delete []Lista;

		Lista = Lista1;
		Lista[iter] = a;
		iter++;
	}

}

/*
inline void ZiarnoList::identyfikacjaZiaren(){

int ileTemp=0;
Ziarno *temp;
temp = new Ziarno[ileZiaren];

for(int i=0;i<ileZiaren;i++){
temp[i].ustawRodzaj(-1);
}

for(int i=0;i<ileZiaren;i++){

ziarno[i].pobierzRodzaj();


		bool flaga=true;

		for(int j=0;j<ileTemp;j++){

		if(temp[j].pobierzRodzaj() == ziarno[i].pobierzRodzaj()){
		flaga=false;
		}

		}
		if(flaga){temp[ileTemp] = ziarno[i];ileTemp++;}


}

delete []ziarno;
ziarno = new Ziarno[ileTemp];

for(int j=0;j<ileTemp;j++){
ziarno[j] = temp[j];
ziarno[j].ustawRodzaj(j+1);
}

delete []temp;


}


ZbiorZiaren::ZbiorZiaren(int ile){

ileZiaren = ile;

ziarno = new Ziarno[ile];

for(int i=0;i<ile;i++){
ziarno[i].ustawRodzaj(0);
}

}
*/
