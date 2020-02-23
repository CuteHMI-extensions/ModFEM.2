#ifndef PLANSZA_H_
#define PLANSZA_H_

#include "ziarno.h"

class Komorka{

int rodzajKomorki;

public:

int runda;
int sasiad[26];
int polozenie[6];

//void setKomorka(int* s,int* p){for(int i=0;i<26;i++){rodzajKomorki=-1;sasiad[i]=s[i];}for(int i=0;i<6;i++){polozenie[i]=p[i];}}

void rysujKomorkaKrawedzie();
void rysujKomorkaPowierzchnie();
void rysujKomorkaPunkt();

double getPolozenieSX(){return (polozenie[3]+polozenie[0])*0.5;}
double getPolozenieSY(){return (polozenie[4]+polozenie[1])*0.5;}
double getPolozenieSZ(){return (polozenie[5]+polozenie[2])*0.5;}

void setRodzajKomorki(int rodzaj){ rodzajKomorki = rodzaj ;}
int getRodzajKomorki(){return rodzajKomorki;}


};



/*
void Komorka::rysujKomorkaPunkt(){

glPointSize(2.5);

glBegin(GL_POINTS);

glVertex3f(polozenie[0],polozenie[1],polozenie[2]);

glEnd();

}
*/

/*
void Komorka::rysujKomorkaPowierzchnie(){

glBegin(GL_QUADS);

glVertex3f(polozenie[0],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[2]);
glVertex3f(polozenie[0],polozenie[4],polozenie[2]);

glVertex3f(polozenie[0],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[1],polozenie[5]);
glVertex3f(polozenie[0],polozenie[1],polozenie[5]);

glVertex3f(polozenie[0],polozenie[1],polozenie[2]);
glVertex3f(polozenie[0],polozenie[4],polozenie[2]);
glVertex3f(polozenie[0],polozenie[4],polozenie[5]);
glVertex3f(polozenie[0],polozenie[1],polozenie[5]);




glVertex3f(polozenie[0],polozenie[1],polozenie[5]);
glVertex3f(polozenie[3],polozenie[1],polozenie[5]);
glVertex3f(polozenie[3],polozenie[4],polozenie[5]);
glVertex3f(polozenie[0],polozenie[4],polozenie[5]);

glVertex3f(polozenie[0],polozenie[4],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[5]);
glVertex3f(polozenie[0],polozenie[4],polozenie[5]);

glVertex3f(polozenie[3],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[5]);
glVertex3f(polozenie[3],polozenie[1],polozenie[5]);

glEnd();

}
*/

/*
void Komorka::rysujKomorkaKrawedzie(){

glBegin(GL_LINES);

glVertex3f(polozenie[0],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[1],polozenie[2]);

glVertex3f(polozenie[0],polozenie[4],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[2]);

glVertex3f(polozenie[0],polozenie[1],polozenie[5]);
glVertex3f(polozenie[3],polozenie[1],polozenie[5]);

glVertex3f(polozenie[0],polozenie[4],polozenie[5]);
glVertex3f(polozenie[3],polozenie[4],polozenie[5]);


glVertex3f(polozenie[0],polozenie[1],polozenie[2]);
glVertex3f(polozenie[0],polozenie[1],polozenie[5]);

glVertex3f(polozenie[0],polozenie[4],polozenie[2]);
glVertex3f(polozenie[0],polozenie[4],polozenie[5]);

glVertex3f(polozenie[3],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[1],polozenie[5]);

glVertex3f(polozenie[3],polozenie[4],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[5]);


glVertex3f(polozenie[0],polozenie[1],polozenie[2]);
glVertex3f(polozenie[0],polozenie[4],polozenie[2]);

glVertex3f(polozenie[0],polozenie[1],polozenie[5]);
glVertex3f(polozenie[0],polozenie[4],polozenie[5]);

glVertex3f(polozenie[3],polozenie[1],polozenie[2]);
glVertex3f(polozenie[3],polozenie[4],polozenie[2]);

glVertex3f(polozenie[3],polozenie[1],polozenie[5]);
glVertex3f(polozenie[3],polozenie[4],polozenie[5]);

glEnd();

}
*/

class Plansza{

int wielkoscX,wielkoscY,wielkoscZ;
int rysujX,rysujY,rysujZ;
bool warPer;

int *wyszukajSasiadow(int x,int y,int z);
void ustawSasiadowPolozenieKomorkiPer(Komorka &Komorka,int x,int y,int z);
void ustawSasiadowPolozenieKomorkiNiePer(Komorka &k,int x,int y,int z);
int polPerX(int a){if(a==-1){return wielkoscX-1;}else if(a==wielkoscX){return 0;}return a;}
int polPerY(int a){if(a==-1){return wielkoscY-1;}else if(a==wielkoscY){return 0;}return a;}
int polPerZ(int a){if(a==-1){return wielkoscZ-1;}else if(a==wielkoscZ){return 0;}return a;}

int polNiePer(int ax,int ay,int az);



public:

//do usuniecia

IntList niePokZiarn;

//

ZiarnoList ziarna;
Komorka* komorki;


//TMemo *m1;

//void setMemo(TMemo*a){m1=a;}

Plansza(int wX,int wY,int wZ,bool periodyczne){setPlansza(wX,wY,wZ,periodyczne,true);}

int getWielkoscX(){return wielkoscX;}
int getWielkoscY(){return wielkoscY;}
int getWielkoscZ(){return wielkoscZ;}

bool getWarPer(){return warPer;};
void setPlansza(int wX,int wY,int wZ,bool periodyczne,bool wlaczAC);
void wczytajAC(const char *nazwa);


void rysujSasiadow(int ktory);
void rysujWszystkieKomorki();
void rysujZiarna(bool punktami,int polZ);
void rysujZiarnaGraniczne();
void rysujElementyZiarn(PunktList &punkty,ElementList &elementy);
void rysujElementyZiarnPrzekroj(PunktList &punkty,ElementList &elementy,int x,int y,int z,bool rysowac);
void elementyDoZiaren(ElementList &elementy);
void resetZiaren(int ileZ);
void resetKomorek();

void przesunZiarno(int numerZiarna,PunktList &points,ElementList &elements,int kierunek,double wartosc){
    ziarna.przesunZiarno(numerZiarna,points,elements,kierunek,wartosc);
}

void obrocZiarno(int numerZiarna,PunktList &points,ElementList &elements,int kierunek,double wartosc){
    ziarna.obrocZiarno(numerZiarna,points,elements,kierunek,wartosc);
}

~Plansza(){delete []komorki;}

};


#endif // PLANSZA_H_



