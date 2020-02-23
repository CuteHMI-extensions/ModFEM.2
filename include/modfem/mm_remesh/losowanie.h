#ifndef LOSOWANIE_H_
#define LOSOWANIE_H_

#include "wektor.h"
#include <cstdlib>

class Losowanie{
int ilePSZ;

int wymiarX,wymiarY,wymiarZ;
Punkt startP[8];
IntList pomoc;


PunktList testG;
PunktList punkty;

void wybierzKrawedz(int ktora);

void losujPlaszczyznyPer(int dok,int ileTabPX,int ileTabPY,int ileTabPZ,char **TabP,int odGranicy,int odPunktu,double ax,double b0, int rszuk);
void losujPlaszczyznyNiePer(int dok,int ileTabPX,int ileTabPY,int ileTabPZ,char **TabP,int odGranicy,int odPunktu,double ax,double b0, int rszuk);
void losujKrawedziePer(int dok,int ileTabPX,int ileTabPY,int ileTabPZ,char **TabP,int odGranicy,int odPunktu,double ax,double b0, int rszuk);
void losujKrawedzieNiePer(int dok,int ileTabPX,int ileTabPY,int ileTabPZ,char **TabP,int odGranicy,int odPunktu,double ax,double b0, int rszuk);

void losujWObjetosci(int dok,int ileTabPX,int ileTabPY,int ileTabPZ,char **TabP,int ilePunktow,int odGranicy,int odPunktu,double ax,double b0, int rszuk);
void dodajPunktyWWolnepola(char **TabP,int ileTabPX,int ileTabPY,int ileTabPZ);

void dodajPunktyPodSpawanie(char **TabP,int ileTabPX,int ileTabPY,int ileTabPZ);

public:

void setPG(PunktList &pg);
void skalujWymiar(double liczba);

int getIlePSZ(){return ilePSZ;}
//TMemo *m1;

void punktyPoczatkoweSzkielet(int range,int prob);
void punktyPoczatkoweOctree(int podziel,int ileRazy);

void losowanie(){;}
void testGranica(int rozdz);
//void setMemo(TMemo *m){m1=m;}

void losuj(bool losuj,int wymX,int wymY,int wymZ,int ileP,int rozdz,int odGranicy,int odPunktu,bool flagaPer,double ax,double b0, int rszuk);
void losujPromien(int wym,int ileP,int promien);
void losujMieszanie(int wym,int ileP,int promien);

int getIlePunktow(){return punkty.getIter();}
double getWymiarX(){return wymiarX;}
double getWymiarY(){return wymiarY;}
double getWymiarZ(){return wymiarZ;}
double getWymiarMax();

PunktList & zwrocPunkty(){return punkty;}
PunktList & zwrocPunktyZGranica();
void rysujPunkty(double maxX,double maxY,double maxZ);
void rysujPunktySzkieletu();
void rysuj123();

void wczytajPunktySztuczneSerce(const char *nazwa);
void ileGranicznych();
void dodajPunktyPodSpawanie(int ileTabPX,int ileTabPY,int ileTabPZ);

};



#endif //LOSOWANIE_H_
