#ifndef WEKTOR_H_
#define WEKTOR_H_

#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <iostream>


using std::string;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::cout;



class StringList{

int iter,grow,max;

string *Lista;

public:


StringList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new string[max];}
StringList(){iter=0;max=10;grow=5;Lista = new string[max];}

void setElement(string &a);
int setElementUniqatReturnPosition_plus1(string &a){for(int i=0;i<iter;i++){if(a == Lista[i]){return i+1;}}setElement(a);return iter;}
int setElementUniqatReturnPosition(string &a){for(int i=0;i<iter;i++){if(a == Lista[i]){return i;}}setElement(a);return iter-1;}

string& getElement(int ktory){return Lista[ktory];}
string& getLastElement(){return Lista[iter-1];}
string& getLastElement(int ktory_od_konca){return Lista[iter-(ktory_od_konca+1)];}
void czysc(int grow1,int max1){grow=grow1;max=max1; if(Lista!=NULL){delete []Lista; Lista = new string[max];iter=0;}}
void ustawIter(int iter1){iter = iter1;}
void usunElement(const int &ktory){Lista[ktory]=Lista[--iter];}
~StringList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}

};


////////////////////////////////////////////////////////////////////////////////////////////////////
/* IntList , DoubleList , Punkt , Element , PunktList , ElementList                               */
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              IntList                                           */
////////////////////////////////////////////////////////////////////////////////////////////////////

        
class IntList{

int iter,grow,max;

int *Lista;
public:

IntList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new int[max];}
IntList(){iter=0;max=10;grow=5;Lista = new int[max];}

int getIter(){return iter;}
int getM(){return max;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}

int funPodWarunki(const int &a1,const int &a2,int &war);

int getElementKey(const int key){for(int i=0;i<iter;i+=2){if(Lista[i]==key){return Lista[i+1];}}return -10;}
int getElement2Kay(const int &key1,const int &key2);
int uniqatPozycja(const int &jaki){for(int i=0;i<iter;++i){if(Lista[i]==jaki){return i;}}return (-1);} //jesli nie znajdzie zwraca -1

void setElement(const int &a);
bool setElementUniqat(const int &a){for(int i=0;i<iter;i++){if(a == Lista[i]){return false;}}setElement(a);return true;}
bool setElementUniqatOddNumbers(const int &a){for(int i=1;i<iter;i+=2){if(a == Lista[i]){return false;}}setElement(a);return true;}
bool setElementUniqatEvenNumbers(const int &a){for(int i=0;i<iter;i+=2){if(a == Lista[i]){return false;}}setElement(a);return true;}
void podmienElement(const int &ktory,const int &a){if(ktory<iter)Lista[ktory]=a;}
int getElement(int ktory){return Lista[ktory];}
int getLastElement(){return Lista[iter-1];}
void czysc(int grow1,int max1){grow=grow1;max=max1; if(Lista!=NULL){delete []Lista; Lista = new int[max];iter=0;}}
void ustawIter(int iter1){iter = iter1;}
void usunElement(const int &ktory){Lista[ktory]=Lista[--iter];}
void zamienMiejscami(const int poz_a,const int poz_b){int temp =Lista[poz_b]; Lista[poz_b] = Lista[poz_a];Lista[poz_b] = temp;}
~IntList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}

bool sprawdzElement(const int &jaki); //jesli jest zwraca true dodatkowo sprawdza czy jest -1(false)

bool sprCzyJestEleWList(const int &jaki);//jesli jest zwraca false
bool sprawdzCzyJest(const int &jaki){for(int i=0;i<iter;++i){if(Lista[i]==jaki){return true;}}return false;} //jesli jest zwraca true

void zapisDoPliku(const char *nazwa);
void wczytajZPliku(const char *nazwa);
void zamienIntListy(IntList &a);
void zamienWyzerujIntListy(IntList &a,int grow,int max);
void operator=(IntList &a);
void losowyZmianaMiejsca();
void dodajIntList(IntList &a);
void sortowanie();

};

////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              DoubleList                                        */
////////////////////////////////////////////////////////////////////////////////////////////////////

class DoubleList{

int iter,grow,max;

double *Lista;
public:

DoubleList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new double[max];}
DoubleList(){iter=0;max=10;grow=2;Lista = new double[max];}

int getIter(){return iter;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}

void setElement(const double &a);
void podmienElement(const int &ktory,const double &a){if(ktory<iter)Lista[ktory]=a;}
double getElement(int ktory){return Lista[ktory];}
double getLastElement(){return Lista[iter-1];}
void czysc(const int &grow1,const int &max1){grow=grow1;max=max1; if(Lista!=NULL){delete []Lista; Lista = new double[max];iter=0;}}
void ustawIter(int iter1){iter = iter1;}
void podzielPrzezL(const double &a){for(int i=0;i<iter;i++){Lista[i]/=a;}}
double getSrednia();
double getSuma();


~DoubleList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}
};


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Punkt                                             */
////////////////////////////////////////////////////////////////////////////////////////////////////

class Punkt{

char granica;
double x,y,z;
double vx,vy,vz;

public:

//IntList sasiad;

//void wypiszWsp(TMemo *m1){AnsiString t;t="[ ";t+=x;t+=" -- ";t+=y;t+=" -- ";t+=z;t+=" ]";m1->Lines->Add(t);}

Punkt(){;}
Punkt(double x1,double y1,double z1){x=x1;y=y1;z=z1;granica='a';}

void setPunktG(const double &x1,const double &y1,const double &z1){x=x1;y=y1;z=z1;}
void setPunkt(Punkt a){x=a.getX();y=a.getY();z=a.getZ();granica = a.getGranica();}
void setPunkt(const double &x1,const double &y1,const double &z1){x=x1;y=y1;z=z1;granica='a';}
void setPunkt(const double &x1,const double &y1,const double &z1,const char &g){x=x1;y=y1;z=z1;granica=g;}
void setX(const double &x1){x=x1;}
void setY(const double &y1){y=y1;}
void setZ(const double &z1){z=z1;}
void setXYZ(const double &x1,const double &y1,const double &z1){x=x1;y=y1;z=z1;}


void setGranica(const char &a){granica=a;}
void przesunPunkt(double xX,double yY,double zZ){x+=xX;y+=yY;z+=zZ;}

char getGranica(){return granica;}
double getX(){return x;}
double getY(){return y;}
double getZ(){return z;}
void setVecNorm(double v1,double v2,double v3){vx=v1;vy=v2;vz=v3;}
void getVecNorm(double *tab){tab[0]=vx;tab[1]=vy;tab[2]=vz;}

float odlegloscPunktu(Punkt &a1);
float odlegloscPunktuSqrt(Punkt &a1){
return sqrt((x-a1.x)*(x-a1.x))+((y-a1.y)*(y-a1.y))+((z-a1.z)*(z-a1.z));
}


float odlegloscPunktu(const double &x1,const double &y1,const double &z1){return (x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1);}
float odlegloscPunktuSqrt(const double &x1,const double &y1,const double &z1){return sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1));}

Punkt & operator*(const double &liczba){x*=liczba;y*=liczba;z*=liczba; return *this;}
void operator=(Punkt &a){x=a.getX();y=a.getY();z=a.getZ();granica = a.getGranica();}

//bool operator==(Punkt &a)const{return (x==a.getX()&&y==a.getY()&&z==a.getZ()&&granica&&a.getGranica()) ;}
bool operator==(Punkt &a)const{return (x==a.getX()&&y==a.getY()&&z==a.getZ()) ;}
bool operator!=(Punkt &a)const{return !(*this==a);}



};


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Element                                           */
////////////////////////////////////////////////////////////////////////////////////////////////////

class Element{


int rodzajZiarna,temp;
char f_Fluid;
int punkty[4];
int elementy[4];
double Sx,Sy,Sz;
float R,StosR_r;
double Ax,Ay,Az;
double uk_val[1];//double uk_val[1]; ALE



public:

Element(){rodzajZiarna=-1;elementy[0]=-1;elementy[1]=-1;elementy[2]=-1;elementy[3]=-1;Ax=-1;Ay=-1;Az=-1;Sx=-1;Sy=-1;Sz=-1;R=-1;StosR_r=999999;punkty[0]=-1;punkty[1]=-1;punkty[2]=-1;punkty[3]=-1;}
Element(int p1,int p2,int p3,int p4){punkty[0]=p1;punkty[1]=p2;punkty[2]=p3;punkty[3]=p4;elementy[0]=-1;elementy[1]=-1;elementy[2]=-1;elementy[3]=-1;rodzajZiarna=-1;Ax=-1;Ay=-1;Az=-1;Sx=-1;Sy=-1;Sz=-1;R=-1;StosR_r=999999;}
Element(int p1,int p2,int p3,int p4,double Sx1,double Sy1,double Sz1,float R1,double Ax1,double Ay1,double Az1){punkty[0]=p1;punkty[1]=p2;punkty[2]=p3;punkty[3]=p4;Sx=Sx1;Sy=Sy1;Sz=Sz1;R=R1;elementy[0]=-1;elementy[1]=-1;elementy[2]=-1;elementy[3]=-1;Ax=Ax1;Ay=Ay1;Az=Az1;rodzajZiarna=-1;StosR_r=999999;}

void ustawTempnaZiarna(){temp=rodzajZiarna;}

void setTemp(int temp1){temp=temp1;}
int getTemp(){return temp;}

void get_uk_val(int id_gauss,double *uk){id_gauss*=3;uk[0]=uk_val[id_gauss+0];uk[1]=uk_val[id_gauss+1];uk[2]=uk_val[id_gauss+2];}
void set_uk_val(int id_gauss,double *uk){id_gauss*=3;uk_val[id_gauss+0]=uk[0];uk_val[id_gauss+1]=uk[1];uk_val[id_gauss+2]=uk[2];}
void set_uk_val_zero(){for(int i=0;i<12;++i){uk_val[i]=0;}}
void setElement(int p1,int p2,int p3,int p4){punkty[0]=p1;punkty[1]=p2;punkty[2]=p3;punkty[3]=p4;elementy[0]=-1;elementy[1]=-1;elementy[2]=-1;elementy[3]=-1;}
void setElement(int p1,int p2,int p3,int p4,int rodzZ){punkty[0]=p1;punkty[1]=p2;punkty[2]=p3;punkty[3]=p4;elementy[0]=-1;elementy[1]=-1;elementy[2]=-1;elementy[3]=-1;rodzajZiarna=rodzZ;}
void setElement(int p1,int p2,int p3,int p4,double Sx1,double Sy1,double Sz1,float R1,double Ax1,double Ay1,double Az1);
void setElement(int p1,int p2,int p3,int p4,double Sx1,double Sy1,double Sz1,float R1,double Ax1,double Ay1,double Az1,float StosunekR_r);

void setPSO(double Sx1,double Sy1,double Sz1,float R1){Sx=Sx1;Sy=Sy1;Sz=Sz1;R=R1;}
void setSSE(double Ax1,double Ay1,double Az1){Ax=Ax1;Ay=Ay1;Az=Az1;}
void setElementySasiad(int e1,int e2,int e3,int e4){elementy[0]=e1;elementy[1]=e2;elementy[2]=e3;elementy[3]=e4;}

void setRodzajZiarna(int a){rodzajZiarna= a;}
int getRodzajZiarna(){return rodzajZiarna;}

float getStosunekR_r(){return StosR_r;}
void setStosunekR_r(float s){StosR_r=s;}
void setFluid(char a){f_Fluid = a;}
int isFluid(){if(f_Fluid=='f'){return 1;}return 0;}

int getP1(){return punkty[0];}
int getP2(){return punkty[1];}
int getP3(){return punkty[2];}
int getP4(){return punkty[3];}

/*
int getP1(){return punkty[1];}
int getP2(){return punkty[2];}
int getP3(){return punkty[0];}
int getP4(){return punkty[3];}
*/
int getP(int ktory){return punkty[ktory];}
int* getPunkty(){return punkty;}

int getE1(){return elementy[0];}
int getE2(){return elementy[1];}
int getE3(){return elementy[2];}
int getE4(){return elementy[3];}
int getE(int ktory){return elementy[ktory];}

void setE1(int e){elementy[0]=e;}
void setE2(int e){elementy[1]=e;}
void setE3(int e){elementy[2]=e;}
void setE4(int e){elementy[3]=e;}
void setE(int e,int ktory){elementy[ktory]=e;}

/* */
void setP1(int p){punkty[0]=p;}
void setP2(int p){punkty[1]=p;}
void setP3(int p){punkty[2]=p;}
void setP4(int p){punkty[3]=p;}

void setAx(double Ax1){Ax=Ax1;}
void setAy(double Ay1){Ay=Ay1;}
void setAz(double Az1){Az=Az1;}

int* getElementy(){return elementy;}

double getSx(){return Sx;}
double getSy(){return Sy;}
double getSz(){return Sz;}

float getR(){return R;}
double getAx(){return Ax;}
double getAy(){return Ay;}
double getAz(){return Az;}

bool podmienPunkt(const int &ktoryNr,const int &);
int getNumerSciany(int p1,int p2,int p3);
bool sprawdzSciane(int numerS,int p1,int p2,int p3);
bool sprawdzCzyJestTakiPunkt(int nrPunktu);
void getNumeryElSciany(int ktoraSciana,int punk[3]);
void getNumeryElSciany(char ktoraSciana,int punk[3]);

};


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Pryzmy                                            */
////////////////////////////////////////////////////////////////////////////////////////////////////

class Pryzm{


int rodzajZiarna;
int punkty[6];

public:

Pryzm(){rodzajZiarna=-1;;}
Pryzm(int p1,int p2,int p3,int p4,int p5,int p6){punkty[0]=p1;punkty[1]=p2;punkty[2]=p3;punkty[3]=p4;punkty[4]=p5;punkty[5]=p6;}


void setElement(int p1,int p2,int p3,int p4,int p5,int p6){punkty[0]=p1;punkty[1]=p2;punkty[2]=p3;punkty[3]=p4;punkty[4]=p5;punkty[5]=p6;}


void setRodzajZiarna(int a){rodzajZiarna= a;}
int getRodzajZiarna(){return rodzajZiarna;}

int getP1(){return punkty[0];}
int getP2(){return punkty[1];}
int getP3(){return punkty[2];}
int getP4(){return punkty[3];}
int getP5(){return punkty[4];}
int getP6(){return punkty[5];}
int* getPunkty(){return punkty;}


/* */
void setP1(int p){punkty[0]=p;}
void setP2(int p){punkty[1]=p;}
void setP3(int p){punkty[2]=p;}
void setP4(int p){punkty[3]=p;}
void setP5(int p){punkty[4]=p;}
void setP6(int p){punkty[5]=p;}


bool podmienPunkt(const int &ktoryNr,const int &nowyNumer);

};


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              PunktyList                                        */
////////////////////////////////////////////////////////////////////////////////////////////////////

class PunktList{

int iter,grow,max;

Punkt *Lista;
public:

PunktList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new Punkt[max+1];}
PunktList(){iter=0;max=10;grow=5;Lista = new Punkt[max+1];}

int getIter(){return iter;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}

void setElement(const double &x1,const double &y1,const double &z1);
void setElement(const double &x1,const double &y1,const double &z1,const char &g);
bool setElementUniqat(const double &x1,const double &y1,const double &z1,const char &g);
void setElement(Punkt &a);
void zamianaMiejsc(const int &a,const int &b){Lista[iter] = Lista[b];Lista[b]=Lista[a];Lista[a]=Lista[iter];}

void podmienElement(int ktory,Punkt &a){if(ktory<iter){Lista[ktory].setPunkt(a.getX(),a.getY(),a.getZ());}}
void podmienElement(int ktory,const double &x1,const double &y1,const double &z1){if(ktory<iter){Lista[ktory].setPunkt(x1,y1,z1);}}

Punkt &getElement(int ktory){return Lista[ktory];}
Punkt &getLastElement(){return Lista[iter-1];}
Punkt &getLastElement(int odKonca){return Lista[iter-(1+odKonca)];}
void zmianaKolejnosci();
void losowyZmianaMiejsca();
                                    
void czysc(int grow1,int max1){grow=grow1;max=max1; if(Lista!=NULL){delete []Lista; Lista = new Punkt[max];iter=0;}}
void ustawIter(int iter1){iter = iter1;}
void zmienIterO(int iter1){iter += iter1;}
void usunElement(int ktory){Lista[ktory]=Lista[--iter];}
PunktList & operator*(const double &liczba){for(int i=0;i<iter;i++){Lista[i]=Lista[i]*liczba;} return *this;}

void zapisDoPliku(const char *nazwa);
void wczytajZPliku(const char *nazwa,bool tylkoGranica);
void wczytajZPlikuPrzerob(const char *nazwa);

void copyPunktList(PunktList &a);
bool porownajPunkty(int p1,int p2){if (Lista[p1]==Lista[p2]){return true;}else{return false;};}

~PunktList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}
};



////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              ElementyList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////

class ElementList{

int iter,grow,max;
Element *Lista;

public:



ElementList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new Element[max];}
ElementList(){iter=0;max=10;grow=2;Lista = new Element[max];}

int getIter(){return iter;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}

void setElement(const int &p1,const int &p2,const int &p3,const int &p4);
void setElement(const int &p1,const int &p2,const int &p3,const int &p4,const int &rodzZ);
void setElement(const int &p1,const int &p2,const int &p3,const int &p4,const double &Sx1,const double &Sy1,const double &Sz1,const float &R1,const double &Ax1,const double &Ay1,const double &Az1);
void setElement(const int &p1,const int &p2,const int &p3,const int &p4,const double &Sx1,const double &Sy1,const double &Sz1,const float &R1,const double &Ax1,const double &Ay1,const double &Az1,const float &R_r);
void setElement(Element &a);

void podmienElement(int ktory,Element a){if(ktory<iter)Lista[ktory] = a;}
void podmienElement(int ktory,int p1,int p2,int p3,int p4){if(ktory<iter){Lista[ktory].setElement(p1,p2,p3,p4);}}
void podmienElement(int ktory,int p1,int p2,int p3,int p4,double sx,double sy,double sz,float r,double Ax1,double Ay1,double Az1){if(ktory<iter)Lista[ktory].setElement(p1,p2,p3,p4,sx,sy,sz,r,Ax1,Ay1,Az1);}
void podmienElement(int ktory,int na){if(ktory<iter)Lista[ktory] = Lista[na];}

Element& getElement(int ktory){return Lista[ktory];}
Element& getLastElement(){return Lista[iter-1];}
Element& getLastElement(int odKonca){return Lista[iter-(1+odKonca)];}

void czysc(int grow1,int max1){grow=grow1;max=max1;; delete []Lista; Lista = new Element[max];iter=0;}
void ustaw_ukval_zero(){for(int i=0;i<iter;i++){Lista[i].set_uk_val_zero();}}
void ustawIter(int iter1){iter = iter1;}
void zmienIterO(int iter1){iter += iter1;}
void deleteElement(int ktory){Lista[ktory]=Lista[--iter];}
void deleteElementPodmienS(int ktory);

//wycinanie ziaren
void usunElementyZaz( IntList &nrZazZiaren);
void usunElementyNieZaz( IntList &nrZazZiaren);
void ustawGraniceElementowZaz( IntList &nrZazZiaren);
void ustawGraniceElementowNieZaz( IntList &nrZazZiaren);
void podmienPunktWElementach(const int &staryNumer,const int &nowyNumer,IntList &elementy);
void ustawSasiadow(const int &numerStary,const int &numerNowy);


int getNumerSasSacianE1(int nrE);
int getNumerSasSacianE2(int nrE);
int getNumerSasSacianE3(int nrE);
int getNumerSasSacianE4(int nrE);
int getNumerSasSacianWybor(int nrE,int nrS);

void ustawTempnaZiarna();
void ustawTempZtabB(int *tabBlock,int dl_tab);

~ElementList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}

};


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              PryzmList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////

class PryzmList{

int iter,grow,max;
Pryzm *Lista;

public:



PryzmList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new Pryzm[max];}
PryzmList(){iter=0;max=10;grow=2;Lista = new Pryzm[max];}

int getIter(){return iter;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}

void setElement(const int &p1,const int &p2,const int &p3,const int &p4,const int &p5,const int &p6);
void setElement(Pryzm &a);

void podmienElement(int ktory,Pryzm a){if(ktory<iter)Lista[ktory] = a;}
void podmienElement(int ktory,int p1,int p2,int p3,int p4,int p5,int p6){if(ktory<iter){Lista[ktory].setElement(p1,p2,p3,p4,p5,p6);}}
void podmienElement(int ktory,int na){if(ktory<iter)Lista[ktory] = Lista[na];}

Pryzm& getElement(int ktory){return Lista[ktory];}
Pryzm& getLastElement(){return Lista[iter-1];}
Pryzm& getLastElement(int odKonca){return Lista[iter-(1+odKonca)];}

void czysc(int grow1,int max1){grow=grow1;max=max1;; delete []Lista; Lista = new Pryzm[max];iter=0;}
void ustawIter(int iter1){iter = iter1;}
void zmienIterO(int iter1){iter += iter1;}
void deleteElement(int ktory){Lista[ktory]=Lista[--iter];}



~PryzmList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}


};





////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Face                                             */
////////////////////////////////////////////////////////////////////////////////////////////////////

class Face{

char faceInEl;
int inEl;
int e1,e2,e3;
int bc;

public:

Face(){;}
Face(const int &ed1,const int &ed2,const int &ed3,const char &faceInEl1,const int &inEl1){e1=ed1;e2=ed2;e3=ed3;faceInEl=faceInEl1;inEl=inEl1;}

void setFace(const Face &a){e1=a.e1;e2=a.e2;e3=a.e3;faceInEl = a.faceInEl;inEl=a.inEl;bc=0;}
void setFace(const int &ed1,const int &ed2,const int &ed3,const char &faceInEl1,const int &inEl1){e1=ed1;e2=ed2;e3=ed3;faceInEl=faceInEl1;inEl=inEl1;bc=0;}
void setFace(const int &ed1,const int &ed2,const int &ed3,const char &faceInEl1,const int &inEl1,const int &bc1){e1=ed1;e2=ed2;e3=ed3;faceInEl=faceInEl1;inEl=inEl1;bc=bc1;}

void setFaceAndEl(const char &a,const int &e){faceInEl=a;inEl=e;}
void setE1(const int &ed1){e1=ed1;}
void setE2(const int &ed2){e2=ed2;}
void setE3(const int &ed3){e3=ed3;}
void setInEl(const int &e){inEl=e;}
void setAllE(const int &ed1,const int &ed2,const int &ed3){e1=ed1;e2=ed2;e3=ed3;}
void setFaceInEl(const char &a){faceInEl=a;}
void setBC(int bc1){bc=bc1;}

char getFaceInEl(){return faceInEl;}
int getEd1(){return e1;}
int getEd2(){return e2;}
int getEd3(){return e3;}
int getInEl(){return inEl;}
int getBC(){return bc;}


void operator=(const Face &a){e1=a.e1;e2=a.e2;e3=a.e3;faceInEl = a.faceInEl;inEl=a.inEl;bc=a.bc;}
bool operator==(const Face &a)const{return (e1==a.e1&&e2==a.e2&&e3==a.e3&&faceInEl==a.faceInEl&&inEl==a.inEl) ;}
bool operator!=(const Face &a)const{return !(*this==a);}

};



////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                         FaceList                                               */
////////////////////////////////////////////////////////////////////////////////////////////////////

class FaceList{

int iter,grow,max;
Face *Lista;

public:

FaceList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new Face[max];}
FaceList(){iter=0;max=10;grow=2;Lista = new Face[max];}

int getIter(){return iter;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}


void setElement(const int &ed1,const int &ed2,const int &ed3,const char &faceInEl1,const int &el);
void setElement(const int &ed1,const int &ed2,const int &ed3,const char &faceInEl1,const int &el,const int &bc);
void setElement(Face &a);

void podmienElement(int ktory,Face a){if(ktory<iter)Lista[ktory] = a;}
void podmienElement(int ktory,const int &ed1,const int &ed2,const int &ed3,const char &faceInEl1,const int &el){if(ktory<iter){Lista[ktory].setFace(ed1,ed2,ed3,faceInEl1,el);}}
void podmienElement(int ktory,int na){if(ktory<iter)Lista[ktory] = Lista[na];}

Face& getElement(int ktory){return Lista[ktory];}
Face& getLastElement(){return Lista[iter-1];}
Face& getLastElement(int odKonca){return Lista[iter-(1+odKonca)];}

void czysc(int grow1,int max1){grow=grow1;max=max1;; delete []Lista; Lista = new Face[max];iter=0;}
void ustawIter(int iter1){iter = iter1;}
void zmienIterO(int iter1){iter += iter1;}
void deleteElement(int ktory){Lista[ktory]=Lista[--iter];}



~FaceList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}


};





////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              Edge                                             */
////////////////////////////////////////////////////////////////////////////////////////////////////

class Edge{

int p1,p2,el;

public:

Edge(){;}
Edge(const int &pun1,const int &pun2){p1=pun1;p2=pun2;}

void setEdge(const Edge &a){p1=a.p1;p2=a.p2;el=a.el;}
void setEdge(const int &pun1,const int &pun2,const int &el1){p1=pun1;p2=pun2;el=el1;}

void setP1(const int &pun1){p1=pun1;}
void setP2(const int &pun2){p2=pun2;}
void setEl(const int &el1){el=el1;}

int getP1(){return p1;}
int getP2(){return p2;}
int getEl(){return el;}

void operator=(const Edge &a){p1=a.p1;p2=a.p2;el=a.el;}
bool operator==(const Edge &a)const{return (p1==a.p1&&p2==a.p2&&el==a.el);}
bool operator!=(const Edge &a)const{return !(*this==a);}

};


////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              EdgeList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////

class EdgeList{

int iter,grow,max;
Edge *Lista;

public:

EdgeList(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new Edge[max];}
EdgeList(){iter=0;max=10;grow=2;Lista = new Edge[max];}

int getIter(){return iter;}
void setGrow(int grow1){if(grow>1){grow=grow1;}}

void setElement(const int &p1,const int &p2,const int &el1);
void setElement(Edge &a);

void podmienElement(int ktory,Edge a){if(ktory<iter)Lista[ktory] = a;}
void podmienElement(int ktory,const int &p1,const int &p2,const int &el1){if(ktory<iter){Lista[ktory].setEdge(p1,p2,el1);}}
void podmienElement(int ktory,int na){if(ktory<iter)Lista[ktory] = Lista[na];}

Edge& getElement(int ktory){return Lista[ktory];}
Edge& getLastElement(){return Lista[iter-1];}
Edge& getLastElement(int odKonca){return Lista[iter-(1+odKonca)];}

void czysc(int grow1,int max1){grow=grow1;max=max1;; delete []Lista; Lista = new Edge[max];iter=0;}
void ustawIter(int iter1){iter = iter1;}
void zmienIterO(int iter1){iter += iter1;}
void deleteElement(int ktory){Lista[ktory]=Lista[--iter];}



~EdgeList(){if(Lista!=NULL){delete []Lista;Lista = NULL;}}


};



////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                              TemplateList                                      */
////////////////////////////////////////////////////////////////////////////////////////////////////

//template <typename TYP> class TLista{

//int iter,grow,max;

//TYP *Lista;
//public:

//TLista(int grow1,int max1){iter=0;grow=grow1;max=max1;Lista = new TYP[max];}
//TLista(){iter=0;max=10;grow=2;Lista = new TYP[max];}

//int getIter(){return iter;}
//void setGrow(int grow1){if(grow>1){grow=grow1;}}

//void setElement(TYP a);
//void podmienElement(int ktory,TYP a){if(ktory<iter)Lista[ktory]=a;}
//TYP getElement(int ktory){return Lista[ktory];}
//TYP getLastElement(){return Lista[iter-1];}

//void czysc(int grow1,int max1){grow=grow1;max=max1; if(Lista!=NULL){delete []Lista; Lista = new TYP[max];iter=0;}}
//void ustawIter(int iter1){iter = iter1;}

//~TLista(){if(*Lista!=0){delete []Lista;*Lista = 0;}}

//};


class test{

public:
IntList punkty;

void rysuj();
void zaladujP(int ile);
test(){punkty.setGrow(15000);;}

};




#endif // WEKTOR_H_
