#include <cstdlib>
#include <modfem/mm_remesh/plansza.h>

#include <modfem/mm_remesh/ziarno.h>

void Plansza::wczytajAC(const char * nazwa) {

	resetZiaren(wielkoscX * wielkoscY * wielkoscZ + 1000);
	int ile = wielkoscX * wielkoscY * wielkoscZ, a;

	ifstream wczytaj(nazwa);

	for (int i = 0; i < ile; ++i) {

		wczytaj >> a;
		wczytaj >> a;
		wczytaj >> a;
		wczytaj >> a;

		komorki[i].runda = 0;
		komorki[i].setRodzajKomorki(a);
	}


	wczytaj.close();

}

void Plansza::resetKomorek() {

	int ile = wielkoscX * wielkoscY * wielkoscZ;

	for (int i = 0; i < ile; ++i) {
		komorki[i].runda = 0;
		komorki[i].setRodzajKomorki(0);
	}

}

void Plansza::setPlansza(int wX, int wY, int wZ, bool periodyczne, bool wlaczAC) {

	rysujX = wX;
	rysujY = wY;
	rysujZ = wZ;

	warPer = periodyczne;
	wielkoscX = wX;
	wielkoscY = wY;
	wielkoscZ = wZ;

	if (!wlaczAC) {
		wielkoscX = 10;
		wielkoscY = 10;
		wielkoscZ = 10;
	}
	int wSred = (wielkoscX + wielkoscY + wielkoscZ) / 3;

	delete []komorki;

	komorki = new Komorka[wielkoscX * wielkoscY * wielkoscZ];

//int sasiad[26];
//int polozenie[6];

	for (int z = 0, iter = 0; z < wielkoscZ; z++) {

		for (int y = 0; y < wielkoscY; y++) {

			for (int x = 0; x < wielkoscX; x++) {

				if (periodyczne) {
					ustawSasiadowPolozenieKomorkiPer(komorki[iter], x, y, z);
				}
				else {
					ustawSasiadowPolozenieKomorkiNiePer(komorki[iter], x, y, z);
				}

				iter++;
			}
		}

	}


	resetKomorek();

	resetZiaren(wSred + 1000);
	niePokZiarn.czysc(wSred, wSred);


}


void Plansza::resetZiaren(int ileZ) {
	ziarna.czysc(ileZ, ileZ);


	for (int i = 0; i < ileZ; ++i) {

		ziarna.setElement(ileZ, rand() % 101 * 0.01, rand() % 101 * 0.01, rand() % 101 * 0.01);

	}

}

int Plansza::polNiePer(int ax, int ay, int az) {

	if (ax == -1 || ax == wielkoscX || ay == -1 || ay == wielkoscY || az == -1 || az == wielkoscZ ) {
		return -1;
	}

	return ax + ay * wielkoscX + az * wielkoscX * wielkoscY;

}



void Plansza::elementyDoZiaren(ElementList & elementy) {


	int ileZ = ziarna.getIter() - 1;

	for (int i = 0, x, y, z, ileE = elementy.getIter(), k1; i < ileE; i++) {

		x = elementy.getElement(i).getAx();
		y = elementy.getElement(i).getAy();
		z = elementy.getElement(i).getAz();

		/*
		if(x>=30){x=29;}
		if(y>=30){y=29;}
		if(z>=30){z=29;}
		*/

		k1 = (x + y * wielkoscX + z * wielkoscX * wielkoscY);;
		k1 = komorki[k1].getRodzajKomorki();
		if (ileZ < k1) {
			k1 = 0;
		}
		elementy.getElement(i).setRodzajZiarna(k1);

	}

}



void Plansza::ustawSasiadowPolozenieKomorkiPer(Komorka & k, int x, int y, int z) {


	k.sasiad[0] = polPerX(x - 1) + polPerY(y - 1) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;
	k.sasiad[1] = polPerX(x) + polPerY(y - 1) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;
	k.sasiad[2] = polPerX(x + 1) + polPerY(y - 1) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;

	k.sasiad[3] = polPerX(x - 1) + polPerY(y) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;
	k.sasiad[4] = polPerX(x) + polPerY(y) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;
	k.sasiad[5] = polPerX(x + 1) + polPerY(y) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;

	k.sasiad[6] = polPerX(x - 1) + polPerY(y + 1) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;
	k.sasiad[7] = polPerX(x) + polPerY(y + 1) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;
	k.sasiad[8] = polPerX(x + 1) + polPerY(y + 1) * wielkoscX + polPerZ(z - 1) * wielkoscX * wielkoscY;

	k.sasiad[9] = polPerX(x - 1) + polPerY(y - 1) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;
	k.sasiad[10] = polPerX(x) + polPerY(y - 1) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;
	k.sasiad[11] = polPerX(x + 1) + polPerY(y - 1) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;

	k.sasiad[12] = polPerX(x - 1) + polPerY(y) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;
	k.sasiad[13] = polPerX(x + 1) + polPerY(y) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;

	k.sasiad[14] = polPerX(x - 1) + polPerY(y + 1) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;
	k.sasiad[15] = polPerX(x) + polPerY(y + 1) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;
	k.sasiad[16] = polPerX(x + 1) + polPerY(y + 1) * wielkoscX + polPerZ(z) * wielkoscX * wielkoscY;

	k.sasiad[17] = polPerX(x - 1) + polPerY(y - 1) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;
	k.sasiad[18] = polPerX(x) + polPerY(y - 1) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;
	k.sasiad[19] = polPerX(x + 1) + polPerY(y - 1) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;

	k.sasiad[20] = polPerX(x - 1) + polPerY(y) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;
	k.sasiad[21] = polPerX(x) + polPerY(y) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;
	k.sasiad[22] = polPerX(x + 1) + polPerY(y) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;

	k.sasiad[23] = polPerX(x - 1) + polPerY(y + 1) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;
	k.sasiad[24] = polPerX(x) + polPerY(y + 1) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;
	k.sasiad[25] = polPerX(x + 1) + polPerY(y + 1) * wielkoscX + polPerZ(z + 1) * wielkoscX * wielkoscY;


	k.polozenie[0] = x;
	k.polozenie[1] = y;
	k.polozenie[2] = z;
	k.polozenie[3] = x + 1;
	k.polozenie[4] = y + 1;
	k.polozenie[5] = z + 1;

	k.setRodzajKomorki(0);

}

void Plansza::ustawSasiadowPolozenieKomorkiNiePer(Komorka & k, int x, int y, int z) {


	k.sasiad[0] = polNiePer((x - 1), (y - 1), (z - 1));
	k.sasiad[1] = polNiePer((x), (y - 1), (z - 1));
	k.sasiad[2] = polNiePer((x + 1), (y - 1), (z - 1));

	k.sasiad[3] = polNiePer((x - 1), (y), (z - 1));
	k.sasiad[4] = polNiePer((x), (y), (z - 1));
	k.sasiad[5] = polNiePer((x + 1), (y), (z - 1));

	k.sasiad[6] = polNiePer((x - 1), (y + 1), (z - 1));
	k.sasiad[7] = polNiePer((x), (y + 1), (z - 1));
	k.sasiad[8] = polNiePer((x + 1), (y + 1), (z - 1));

	k.sasiad[9] = polNiePer((x - 1), (y - 1), (z));
	k.sasiad[10] = polNiePer((x), (y - 1), (z));
	k.sasiad[11] = polNiePer((x + 1), (y - 1), (z));

	k.sasiad[12] = polNiePer((x - 1), (y), (z));
	k.sasiad[13] = polNiePer((x + 1), (y), (z));

	k.sasiad[14] = polNiePer((x - 1), (y + 1), (z));
	k.sasiad[15] = polNiePer((x), (y + 1), (z));
	k.sasiad[16] = polNiePer((x + 1), (y + 1), (z));

	k.sasiad[17] = polNiePer((x - 1), (y - 1), (z + 1));
	k.sasiad[18] = polNiePer((x), (y - 1), (z + 1));
	k.sasiad[19] = polNiePer((x + 1), (y - 1), (z + 1));

	k.sasiad[20] = polNiePer((x - 1), (y), (z + 1));
	k.sasiad[21] = polNiePer((x), (y), (z + 1));
	k.sasiad[22] = polNiePer((x + 1), (y), (z + 1));

	k.sasiad[23] = polNiePer((x - 1), (y + 1), (z + 1));
	k.sasiad[24] = polNiePer((x), (y + 1), (z + 1));
	k.sasiad[25] = polNiePer((x + 1), (y + 1), (z + 1));


	k.polozenie[0] = x;
	k.polozenie[1] = y;
	k.polozenie[2] = z;
	k.polozenie[3] = x + 1;
	k.polozenie[4] = y + 1;
	k.polozenie[5] = z + 1;

	k.setRodzajKomorki(0);

}



/*
void Plansza::rysujSasiadow(int ktory){

	komorki[ktory].rysujKomorkaKrawedzie();

	for(int i=0;i<26;i++){

	if(i<9){glColor3f(1,1,1);}
	else if(i<17){glColor3f(0,1,1);}
	else{glColor3f(0,1,0);}
	komorki[komorki[ktory].sasiad[i]].rysujKomorkaKrawedzie();

	}

}
*/

/*
void Plansza::rysujWszystkieKomorki(){

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	int ileZ = wielkoscX*wielkoscY*wielkoscZ;
	glColor4f(1,0,0,0.05);
	for(int i=0;i<ileZ;i++){

	komorki[i].rysujKomorkaKrawedzie();

	}

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

}
*/

/*
void Plansza::rysujZiarna(bool punktami,int polZ){

	int ileZ = wielkoscX*wielkoscY*wielkoscZ;

	int translateX =rysujX*0.5;
	int translateY =rysujY*0.5;
	int translateZ =rysujZ*0.5;

	Kolor barwa(1,0,0);

glPushMatrix();;

glTranslatef(-rysujX,-rysujY,-rysujZ);

glPushMatrix();


glTranslatef(translateX,translateY,translateZ);

	for(int i=0;i<ileZ;i++){

		if(komorki[i].getRodzajKomorki()){

			if(int(komorki[1].getPolozenieSZ())==polZ){

			barwa = ziarna.getElement(komorki[i].getRodzajKomorki()).pobierzKolor();
			glColor3f(barwa.r,barwa.g,barwa.b);
			if(punktami){komorki[i].rysujKomorkaPunkt();}
			else{komorki[i].rysujKomorkaPowierzchnie();}

			}
			else if(polZ==-1){

			barwa = ziarna.getElement(komorki[i].getRodzajKomorki()).pobierzKolor();
			glColor3f(barwa.r,barwa.g,barwa.b);
			if(punktami){komorki[i].rysujKomorkaPunkt();}
			else{komorki[i].rysujKomorkaPowierzchnie();}

			}

		//komorki[i].rysujKomorkaPunkt();
		//komorki[i].rysujKomorkaKrawedzie();
		}

	}

glPopMatrix();
glPopMatrix();

}
*/

/*
void Plansza::rysujElementyZiarn(PunktList &punkty,ElementList &elementy){



double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;



//double max= wielkosc;

	int translateX =rysujX*0.5;
	int translateY =rysujY*0.5;
	int translateZ =rysujZ*0.5;

Kolor barwa(1,0,0);

glPushMatrix();
glTranslatef(-rysujX,-rysujY,-rysujZ);
glPushMatrix();
glTranslatef(translateX,translateY,translateZ);


glBegin(GL_TRIANGLES);

	for(int i=0,ileE=elementy.getIter();i<ileE;i++){

	if(niePokZiarn.sprawdzElement(elementy.getElement(i).getRodzajZiarna())){

		barwa = ziarna.getElement(elementy.getElement(i).getRodzajZiarna()).pobierzKolor();
		glColor3f(barwa.r,barwa.g,barwa.b);

		x1 = punkty.getElement(elementy.getElement(i).getP1()).getX();
		x2 = punkty.getElement(elementy.getElement(i).getP2()).getX();
		x3 = punkty.getElement(elementy.getElement(i).getP3()).getX();
		x4 = punkty.getElement(elementy.getElement(i).getP4()).getX();

		y1 = punkty.getElement(elementy.getElement(i).getP1()).getY();
		y2 = punkty.getElement(elementy.getElement(i).getP2()).getY();
		y3 = punkty.getElement(elementy.getElement(i).getP3()).getY();
		y4 = punkty.getElement(elementy.getElement(i).getP4()).getY();

		z1 = punkty.getElement(elementy.getElement(i).getP1()).getZ();
		z2 = punkty.getElement(elementy.getElement(i).getP2()).getZ();
		z3 = punkty.getElement(elementy.getElement(i).getP3()).getZ();
		z4 = punkty.getElement(elementy.getElement(i).getP4()).getZ();

		glVertex3f(x1,y1,z1);glVertex3f(x2,y2,z2);glVertex3f(x3,y3,z3);
		glVertex3f(x2,y2,z2);glVertex3f(x4,y4,z4);glVertex3f(x3,y3,z3);
		glVertex3f(x1,y1,z1);glVertex3f(x3,y3,z3);glVertex3f(x4,y4,z4);
		glVertex3f(x1,y1,z1);glVertex3f(x4,y4,z4);glVertex3f(x2,y2,z2);

		///
			if((x1 == x2 && x2 == x3 && (x1==0 || x1 == max)) ||
			   (y1 == y2 && y2 == y3 && (y1==0 || y1 == max)) ||
			   (z1 == z2 && z2 == z3 && (z1==0 || z1 == max))
			){
			glVertex3f(x1,y1,z1);glVertex3f(x2,y2,z2);glVertex3f(x3,y3,z3);
			}

			if((x4 == x2 && x2 == x3 && (x4==0 || x4 == max)) ||
			   (y4 == y2 && y2 == y3 && (y4==0 || y4 == max)) ||
			   (z4 == z2 && z2 == z3 && (z4==0 || z4 == max))
			){
			glVertex3f(x2,y2,z2);glVertex3f(x4,y4,z4);glVertex3f(x3,y3,z3);
			}

			if((x1 == x4 && x4 == x3 && (x1==0 || x1 == max)) ||
			   (y1 == y4 && y4 == y3 && (y1==0 || y1 == max)) ||
			   (z1 == z4 && z4 == z3 && (z1==0 || z1 == max))
			){
			glVertex3f(x1,y1,z1);glVertex3f(x3,y3,z3);glVertex3f(x4,y4,z4);
			}

			if((x1 == x2 && x2 == x4 && (x1==0 || x1 == max)) ||
			   (y1 == y2 && y2 == y4 && (y1==0 || y1 == max)) ||
			   (z1 == z2 && z2 == z4 && (z1==0 || z1 == max))
			){
			glVertex3f(x1,y1,z1);glVertex3f(x4,y4,z4);glVertex3f(x2,y2,z2);
			}
		///

	}

	}

glEnd();


glPopMatrix();
glPopMatrix();



}
*/

/*
void Plansza::rysujElementyZiarnPrzekroj(PunktList &punkty,ElementList &elementy,int Px,int Py,int Pz,bool rysuj){

double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;


	int translateX =rysujX*0.5;
	int translateY =rysujY*0.5;
	int translateZ =rysujZ*0.5;

Kolor barwa(1,0,0);

glPushMatrix();
glTranslatef(-rysujX,-rysujY,-rysujZ);
glPushMatrix();
glTranslatef(translateX,translateY,translateZ);

glBegin(GL_TRIANGLES);

	for(int i=0,ileE=elementy.getIter();i<ileE;i++){

	if(niePokZiarn.sprCzyJestEleWList(elementy.getElement(i).getRodzajZiarna())){

		barwa = ziarna.getElement(elementy.getElement(i).getRodzajZiarna()).pobierzKolor();
		glColor3f(barwa.r,barwa.g,barwa.b);;

		x1 = punkty.getElement(elementy.getElement(i).getP1()).getX();
		x2 = punkty.getElement(elementy.getElement(i).getP2()).getX();
		x3 = punkty.getElement(elementy.getElement(i).getP3()).getX();
		x4 = punkty.getElement(elementy.getElement(i).getP4()).getX();

		y1 = punkty.getElement(elementy.getElement(i).getP1()).getY();
		y2 = punkty.getElement(elementy.getElement(i).getP2()).getY();
		y3 = punkty.getElement(elementy.getElement(i).getP3()).getY();
		y4 = punkty.getElement(elementy.getElement(i).getP4()).getY();

		z1 = punkty.getElement(elementy.getElement(i).getP1()).getZ();
		z2 = punkty.getElement(elementy.getElement(i).getP2()).getZ();
		z3 = punkty.getElement(elementy.getElement(i).getP3()).getZ();
		z4 = punkty.getElement(elementy.getElement(i).getP4()).getZ();

		if(elementy.getElement(i).getRodzajZiarna()==2 && rysuj){

			glVertex3f(x1,y1,z1);glVertex3f(x2,y2,z2);glVertex3f(x3,y3,z3);

			glVertex3f(x2,y2,z2);glVertex3f(x4,y4,z4);glVertex3f(x3,y3,z3);

			glVertex3f(x1,y1,z1);glVertex3f(x3,y3,z3);glVertex3f(x4,y4,z4);

			glVertex3f(x1,y1,z1);glVertex3f(x4,y4,z4);glVertex3f(x2,y2,z2);

		}
		else if(x1<Px && x2<Px && x3<Px && x4<Px &&
		   y1<Py && y2<Py && y3<Py && y4<Py &&
		   z1<Pz && z2<Pz && z3<Pz && z4<Pz){

			glVertex3f(x1,y1,z1);glVertex3f(x2,y2,z2);glVertex3f(x3,y3,z3);

			glVertex3f(x2,y2,z2);glVertex3f(x4,y4,z4);glVertex3f(x3,y3,z3);

			glVertex3f(x1,y1,z1);glVertex3f(x3,y3,z3);glVertex3f(x4,y4,z4);

			glVertex3f(x1,y1,z1);glVertex3f(x4,y4,z4);glVertex3f(x2,y2,z2);
		}

	}

	}

glEnd();

glPopMatrix();
glPopMatrix();
}
*/

