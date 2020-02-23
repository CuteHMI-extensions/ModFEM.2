#ifndef ZIARNO_H_
#define ZIARNO_H_

#include <math.h>
#include <cstdlib>

#include "wektor.h"

class Kolor{

public:
float r,g,b;

Kolor(float R,float G,float B ){r=R;g=G;b=B;}
Kolor(){r=0;g=0;b=0;}

};

class Ziarno{

int rodzaj;
Kolor rgb;

public:

void ustalZiarno(int rodz, float R , float G , float B ){rodzaj = rodz;rgb.r = R; rgb.g=G; rgb.b=B;}
void ustawRodzaj(int rodz){rodzaj = rodz;}
void ustawKolor(Kolor k){ rgb.r = k.r; rgb.g=k.g; rgb.b=k.b;}
void ustawKolor( float R , float G , float B ){ rgb.r = R; rgb.g=G; rgb.b=B;}

int pobierzRodzaj(){return rodzaj;}
Kolor &pobierzKolor(){return rgb;}

float pobierzKolorR(){return rgb.r;}
float pobierzKolorG(){return rgb.g;}
float pobierzKolorB(){return rgb.b;}



};

class ZiarnoList{

int iter,grow,max;

Ziarno *Lista;

public:

ZiarnoList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new Ziarno[max];}
ZiarnoList(){iter=0;max=100;grow=100;Lista = new Ziarno[max];}

int getIter(){return iter;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}

void setElement(Ziarno a);
void setElement(int rodzaj,float r,float g,float b);

Ziarno &getElement(int ktory){return Lista[ktory];}
Ziarno &getLastElement(){return Lista[iter-1];}
Ziarno &getLastElement(int odKonca){return Lista[iter-(1+odKonca)];}

void czysc(int grow1,int max1){grow=grow1;max=max1; delete []Lista; Lista = new Ziarno[max];iter=0;}
void ustawIter(int iter1){iter = iter1;}
void zmienIterO(int iter1){iter += iter1;}
void identyfikacjaZiaren();

void przesunZiarno(int numerZiarna,PunktList &points,ElementList &elements,int kierunek,double wartosc);
void obrocZiarno(int numerZiarna,PunktList &points,ElementList &elements,int kierunekOsi,double wartosc);

~ZiarnoList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}
};



#endif // ZIARNO_H_
