#ifndef FIGURY_H_
#define FIGURY_H_

#include "wektor.h"

class Kula{

double x,y,z,r;

public:

Kula(double xX,double yY,double zZ,double rR){x=xX;r=rR;z=zZ;y=yY;}
Kula(){x=0,y=0,z=0,r=0;}

void ustawKule(double xX,double yY,double zZ,double rR){x=xX;r=rR;z=zZ;y=yY;}
int przesunKule(PunktList &points,int kierunek,double wartosc);

void zwrocPolozenie(Punkt &temp){temp.setPunkt(x,y,z);}

};

#endif // FIGURY_H_
