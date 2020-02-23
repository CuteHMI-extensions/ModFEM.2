#include <modfem/mm_remesh/figury.h>



int Kula::przesunKule(PunktList & points, int kierunek, double wartosc) {
	int ile = 0;
	double Px = 0, Py = 0, Pz = 0, odl;

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

	for (int i = 0, ileP = points.getIter(); i < ileP; ++i) {
		odl = points.getElement(i).odlegloscPunktuSqrt(x, y, z);
		if (odl < r || points.getElement(i).getGranica() == 'g') {
			++ile;
			points.getElement(i).przesunPunkt(Px, Py, Pz);
		}

	}
	x += Px;
	y += Py;
	z += Pz;
	return ile;
}
