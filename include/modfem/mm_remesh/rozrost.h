#ifndef ROZROST_H_
#define ROZROST_H_

#include "plansza.h"
#include "wektor.h"

class Rozrost{

int numerZiarna;
int ilePG;
public:

//IntList aktualnyBrzeg;
Plansza * pR;
PunktList pGranicy;
DoubleList pGranicy2D;

//granicaRozrostu
IntList gR;
IntList usunGR;

Rozrost(Plansza * p){ilePG=0;usunGR.czysc(1000,1000);pR = p;numerZiarna=0;pGranicy.czysc(1000,1000);}

void reset(Plansza * p){ilePG=0;usunGR.czysc(1000,1000);pR = p;numerZiarna=0;pGranicy.czysc(1000,1000);}
void rozrostLosuj(int ileZiaren);
void rozrostLosujRange(int ileZiaren,int range);

void rozrostZiaren();
void rozrostZiarenGranicy(int wartoscPrawdo);
void rozrostOtoczenia();

void pomniejszenieZiarnaGranicy(int wspolczynnik);
void pomniejszenieZiarna(int wspolczynnik);
void granice(int wspolczynnik);
void graniceZiaren();
void graniceP1Per();
void graniceP1NiePer();



void komplekosoweTworzenieZiaren(int wartoscPrawdo,int powtorzen,int ileRozrost,int ileOtocz,int RozrostOtoczenie);
void pomniejszGranice(int promien,bool brzeg);

void ustalGraniceRozrostu();
void rysujGranice(int polZ);
void rysujZiarnaGraniczne();
int getilePG(){return ilePG;}
Plansza* getPlansza(){return pR;}
void wczytajPunktyGraniczne(const char *nazwa);
void wczytajPunktyGraniczne2D(const char *nazwa,bool tylkoGraniczneA,int sposobSkalowania);
void wczytajPunktyGranicznePSS(const char *nazwa);

//wlasne granice
void sztucznyRozrost();
void graniceP2();
void granicaKuli(double r,int ilPunkOkrag,double Sx,double Sy,double Sz);
void granicaZastawka(bool tylkoGranicaG,double r,int ilPunkOkrag,int grubosc,double Sx,double Sy,double Sz);
void granicaTuba(double r,int ilPunkOkrag,int grubosc,double Sx,double Sy,double Sz);
void BFSOnlyStep(double dlX,double dlY,double dlZ,double stosunekH,double dlWylou,double gestosc);

void graniceStrukturalna();
void graniceAnizotropowe();
void warstwyZ2D(int ileWarstw);

void granicaObecnePunkty(PunktList &p, bool tylkoG);


};


#endif //ROZROST_H_
