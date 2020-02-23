#ifndef DELANOUY_3D_H_
#define DELANOUY_3D_H_

#include <cmath>

#include "wektor.h"


class Delanouy3D{

//do  usuniecia
double globalTranX,globalTranY,globalTranZ,globalSkal;
//
int flagaLaplasa;

 
double dx,dy,dz;
double NajT;
int schowekLW;
int ilePunktowCzworosciany;

void ustalSSE();
int maximum4(const double &m1,const double &m2,const double &m3,const double &m4);
int segreguj4(double *m,int *kol);

IntList pZiarno;
IntList nrZazZiaren;
IntList wzorGranicy;



public:
IntList fluidBlock,fluidGroup,uppermostBlock;
//



//do usuniecia
int do_kontaktu[192];

IntList wybraneZleEl;
void rysujZleEl();


//======================
PryzmList pryzmy;
ElementList elements;
PunktList points;
//TMemo *m1;
IntList *laplas;
//chwilowe

FaceList faces;
EdgeList edges;
IntList *mapFace;
int *facebc_connect;
bool flag_bc_connect;


void count_bc_surface(FILE* Interactive_output);
void pobierzWarunkiZFacow(IntList &warunkiBrzegowe,int *mapaPunktow1);
void pobierzWarunkiZFacow(IntList &warunkiBrzegowe);
double dlugosc_do_sciany(int punktID,int elID);
double srednieRamiePunktu(int punktID,double mnoznik);
void optymalizacjaStosunekRdorCalosc_zmienna_dl_kroku(double dlKroku,bool bezGranicy);

void smoothing(double l_step,int repetition);
void optymalizacjaStosunekRdorCalosc(double dlKroku,  bool bezGranicy);
void optymalizacjaStosunekRdorNajlepszayWynikCalosc(double dlKroku,  bool bezGranicy);

void info_bc_connect();
void zmien_do_kontaktu(int pozycja);
void mmr_dopasuj_punkty(int face_id,int fa_connect,int *el_id,int zmien);
int get_facebc_connect(int ktory){if(faces.getIter()>ktory && flag_bc_connect){return facebc_connect[ktory];}else{return -1;}}
void creatFaceAndEdge();
void czyscFaceAndEdge();
void czyscFace();
void czyscEdge();
void stosunekRdor(FILE *out);

//void setMemo(TMemo *m){m1=m;}
void zapiszaAktualnaLibczeW(){schowekLW=points.getIter();}
void przestawIterPunktowNaSchowek(){points.ustawIter(schowekLW);}
void zapisParaView(const char *nazwa);

~Delanouy3D();

Delanouy3D(double xdl,double ydl,double zdl){

facebc_connect = new int[1];facebc_connect[0]=-1;flag_bc_connect=false;
flagaLaplasa=0;
laplas = new IntList[8];
mapFace = new IntList[1];
Start(xdl,ydl,zdl,8);
wybraneZleEl.czysc(100000,100000);

for(int i=0;i<192;++i){do_kontaktu[i]=0;}



}


ElementList &zwrocRefElementow(){return elements;}
PunktList &zwrocRefPunktow(){return points;}

int rysuj(bool ujemne);
void rysujPow();
void rysujPunkty(int size);
int rysujWybraneZiarna(IntList &niePokZiarn);
int rysujWybraneZiarnaPrzekroj(IntList &niePokZiarn,int x,int y,int z,bool calosc);
void sprawdzSasiadow();          // sprawdz zgodnosc sasiednich scian
void sprawdzSasiednieElementy(); //              -||-
double V_objetosc(FILE *out);
double V_objetoscP();
double V_objetoscT(FILE *out);
double V_objetoscTPopraw();
double V_objetoscT(int ktory);
bool V_objetoscP(int ktory);


//optymalny delanouy
int przygotojTabliceDoDelanoya(int ** &tabElementow2,int WETE2,double WEws2); //zwraca wartosc Z potrzeba do usuniecia tablicy
int poprawaElementMalejaco(int ** tabElementow2,int model,double &stostunekObecny,double stosunekDocelowy,double &zmniejszO,bool PoprawZPozaObszaru);
int poprawaElementPrzysp(int ** tabElementow2,int model,double &stostunekObecny,bool PoprawZPozaObszaru);


void Start(double xdl,double ydl,double zdl);
void Start(double xdl,double ydl,double zdl,int ileP);
void ustawWIelkoscObszaru(double xdl,double ydl,double zdl){dx=xdl;dy=ydl;dz=zdl;}
void delanoyTablicaOP(double xdl,double ydl,double zdl,bool reMes,PunktList &Tab,double ws,int te,bool dokladneWyszukanie);
void delanoyTablicaOPDodajPunkty(double xdl,double ydl,double zdl,bool reMes,PunktList &Tab,double ws,int te);
void delanoyTablicaOPDodajPunktyZTablicaPolozeniaEl(int **tabelementow,double xdl,double ydl,double zdl,bool reMes,PunktList &Tab,double ws,int te);
void delanoyTablicaOPDodajPunktyZTablicaPolozeniaElUstalJakoscEl(int **tabelementow,double xdl,double ydl,double zdl,bool reMes,PunktList &Tab,double ws,int te);

void delanoyTablicaDodajPunktyWOparciuOTabliceElementow(PunktList &Tab,IntList &zElementu);

bool sprawdzElement(int jaki,IntList &lista);
bool elementyDoDelnoya(int wybrany,double x,double y,double z,IntList &delElementy);
bool elementyDoDelnoyaZeSprObj(int wybrany,double x,double y,double z,IntList &delElementy);
bool sprawdzCzyNalezyDoE(int nrE,const double &x,const double &y,const double &z);

int wyszukajElementOP(const int &wybranyTE,const double &x,const double &y,const double &z, IntList &zlyElement, IntList &pomoc);
int wyszukajElement(const double &x,const double &y,const double &z);
int wyszukajElementPoprawaOP(const double &x,const double &y,const double &z, IntList &zlyElement);
int wyszukajElementGwiazda(const double &x,const double &y,const double &z,IntList &sprawdz,const int &wielkoscG);
int sprawdzenieElementow(const double &x,const double &y,const double &z,IntList &sprawEle,const int &od);


//dopasowanie elementow do ziarna

void dodajPunktPopElement();
void dodajPunktWSrodkuElementow();
void dodajPunktWSrodkuNajwiekszejPow(IntList &zleElementy);
int dopasujElementDoZiarna(bool pokazWynik);
void dopasujElementDoZiarnaNaPodstawieFunkcji(double wartoscKata,double rZasat,double rTube);
int dopasujElementDoZiarnaPunktyStartowe(PunktList &startP);

void poprawaDopasowaniaZiaren(const int &dokladnoscDop,const int &rodzaj);
void poprawaDopasowaniaZiarenNajObjAryt(const int &dokladnoscDop);
void poprawaDopasowaniaZiarenNajKrawAryt(const int &dokladnoscDop);

void powierzchnieScian(double wynik[4],const int &p1,const int &p2,const int &p3,const int &p4);
void objetoscSrodAryt(double wynik[4],const int &p1,const int &p2,const int &p3,const int &p4);
void krawedzieSrodAryt(double wynik[4],const int &p1,const int &p2,const int &p3,const int &p4);


void poprawaDopasowaniaZiaren2Sas(IntList &zleE);
void oznaczWezlyNaPow_b();
void odznaczWezlyNaPow_b();


//wygraldzanie ReMes

void reMes();
void reMes(double wX,double wY,double wZ);
void sasiedniePunkty();
void sasiedniePunktyPryzm();

void sasiednieElementy();
void wygladzanieLaplace();
void wygladzanieLaplaceWaga(int sposobZbieraniaWag,double waga,int ileWag);

void wygladzanieLaplaceWOparciuOWybranePunktu(IntList &punk,int ilePow);
void wygladzanieLaplaceWagaWOparciuOWybranePunktu(IntList &punk,int sposobZbieraniaWag,double waga,int ileWag,int ilePow);

// poprawa podczas przesuniecia
void wygladzanieLaplaceWagaWybraneP(int sposobZbieraniaWag,double waga,int punktyOD,int ileWag);
void wygladzanieLaplaceWagaWybranePunkty(int sposobZbieraniaWag,double waga,IntList &punk,int ileWag);
void wygladzanieLaplaceWybranePunkty(IntList &punk);
void wyszukajPunktuDoWygladzenia(IntList &punk,int numerZiarna,int ilePok);
void optymalizacjaStosunekRdor(double dlKroku);
void optymalizacjaStosunekGMSH(double dlKroku);
void optymalizacjaV(double dlKroku);
void optymalizacjaStosunekRdorWybranePunkty(double dlKroku,IntList &wybraneP,int ileRazy);
void optymalizacjaStosunekRdorDODelanoya(double dlKroku,int odElementu);
void optymalizacjaWCeluUzyskaniaPrawidlowejV(double dlKroku,int numerP);



double stosunekRdorNajlepszayWynikDlaPunktu(int ktoryPunkt);
double stosunekRdorNajlepszayWynikDlaPunktu1(int ktoryPunkt);


void optymalizacjaVWybranePunkty(double dlKroku,IntList &wybraneP,int ileRazy);


void zbierzWagi(int nrSposobuZbieraniaWag,int ktoryP,DoubleList &wagi,int ileWag){

switch(nrSposobuZbieraniaWag){
    case 1:zbierzWagiOdleglosciKrawedzi(ktoryP,wagi);break;
    case 2:zbierzWagiOdleglosciOdPunktowWagi(ktoryP,wagi,ileWag);break;
    case 3:zbierzWagiOdleglosciKrawedziPoprawne(ktoryP,wagi,ileWag);break;

}

}

//nr 1
void zbierzWagiOdleglosciKrawedzi(int ktoryP,DoubleList &wagi);
//nr 2
void zbierzWagiOdleglosciOdPunktowWagi(int ktoryP,DoubleList &wagi,int ileWag);
//nr 3
void zbierzWagiOdleglosciKrawedziPoprawne(int ktoryP,DoubleList &wagi,int ileWag);

//regeneraziaZIaren

void zapiszWzorGranicy();
void zapiszWzorGranicy(IntList &wzorGranicy);
void poprawGranice();

void usunieciePunktowKoloG();
void poprawGraniceUsuniecieP();

void poprawGranice(PunktList &noweP);
void poprawGranice(IntList &wzorGranicy,PunktList &noweP);



void ZapiszDoPlikuAbaqus(const char *nazwa);
void ZapiszDoPlikuNas(const char *nazwa,double dzielnikZ,bool wartosciDo1,bool scianyProstopa);
void ZapiszDoPlikuNasPryzm(const char *nazwa,double dzielnikZ,bool wartosciDo1,bool scianyProstopa);
void ZapiszDoPlikuNasZPrzesunieciem(const char *nazwa,double dzielnikZ,bool scianyProstopa,int transX,int transY,int treansZ,double globalnyDzielnik,int *tab);
void ZapiszDoPlikuNasWOparciuOPowTetra(const char *nazwa,double dzielnikZ,bool wartosciDo1,bool scianyProstopa);
void ZapiszDoPlikuNasWOparciuOPowHybrydPSS(const char *nazwa,double dzielnikZ,bool wartosciDo1,bool scianyProstopa);
void ZapiszDoPlikuNasWOparciuOPowTetraKanal(const char *nazwa,int wartoscMaxDo);
void ZapiszDoPlikuNasWOparciuOPowHybrydKanal(const char *nazwa,int wartoscMaxDo);
void ZapiszDoPlikuNasZPrzesunieciem_BC_Face(const char *nazwa,double dzielnikZ,double transX,double transY,double treansZ,double globalnyDzielnik);    

void ZapiszDoPlikuNasWOparciuOPowHybrydDowolny(const char *nazwa,double dzielnikZ,bool wartosciDo1,int warunki);

void ZapiszDoPlikuSiatkaPow(const char *nazwa);
void zapiszBMP(const char *nazwa,void *handle);

//poprawa siatki

void znajdzElemetyPoDodaniuPunktow(IntList &ktoreP,IntList &ktoreE);
int poprawaElement(int model,double stosunekR_r,bool PoprawZPozaObszaru);

int poprawaElement_popZiarna(int model,double stosunekR_r,bool PoprawZPozaObszaru);
void poprawaElementWszystkieRdor(int sposobZbieraniaWag,double stosunekR_r,bool PoprawZPozaObszaru,double waga,bool modelPop,int ileWag);
void poprawaElementWszystkieRdorPoprawaMiejscowa(int model,double stosunekR_r,bool PoprawZPozaObszaru);
void poprawaElementWszystkieRdorPoprawaMiejscowaMalejaco(int model,double stosunekR_r,bool PoprawZPozaObszaru);
void poprawaElementWszystkieRdorPoprawaMiejscowaPrzysp(int model,double stosunekR_r,bool PoprawZPozaObszaru);

void poprawaElementWszystkieRdorPoprawaMiejscowa_popZiarna(int model,double stosunekR_r,bool PoprawZPozaObszaru);


void poprawaElementWszystkieGMSH(int sposobZbieraniaWag,double stosunekR_r,bool PoprawZPozaObszaru,double waga,bool modelPop,int ileWag);

void stosunekRdorCout();
void stosunekRdorZElementow();

bool stosunekRdor(double flaga);
void stosunekRdorAktywneZiarna(IntList &nrAktywne);
double stosunekRdorWybranyElement(int ktoryElement);

double naprawaElementow_Objetosci(int ktoryPunkt);
double stosunekRdorNajwiekszyDlaPunktu(int ktoryPunkt);
double stosunekRdorNajwiekszyDlaPunktuDODelanoya(int ktoryPunkt,int odElementu);

double stosunekGMSHNajwiekszyDlaPunktu(int ktoryPunkt);
double stosunekVDlaPunktu(int ktoryPunkt);
void stosunekGamaKGMSH();
void stosunekGamaKGMSHAktywneZiarna(IntList &nrAktywne);

void stosunekRdorWyszukajElemnty(double stosunekR_r);
void stosunekRdorWyszukajElemntyMalejacoSprawdz(double stosunekR_r);
void stosunekRdorWyszukajElemntyMalejacoWypelnij();
void stosunekGamaKGMSHWyszukajElemnty(double stosunekR_r);

int ustalNowePunkty(IntList &zleElementy,bool PoprawZPozaObszaru);
int ustalNowePunkty(IntList &zleElementy,bool PoprawZPozaObszaru,PunktList &dodanePunkty);
bool elementyDoDelnoyaPoprawaE(int wybrany,IntList &delElementy);
Punkt ustalDodatkowyPunkt(IntList &delElementy);
//przypadek 2D
float punktSciana2D(int ro1,int rodzaj,const double &War0,const double &x1,const double &x2,const double &x3,const double &y1,const double &y2,const double &y3,Punkt &temp);
float punktSciana3D(const double &x1,const double &y1,const double &z1,const double &x2,const double &y2,const double &z2,const double &x3,const double &y3,const double &z3,Punkt &temp);


void pokazMinMaxV();

void rozbijajDuzeTrojkaty(const double &stosunekR_r,bool PoprawZPozaObszaru);
void wyszukajDuzeTrojkaty(const double &stosunekR_r,IntList &zleElementy);
void zapiszNajVTrojkat();




//usuwanie ziaren

void usunZaznaczoneEl( IntList &nrZazZiaren){elements.ustawGraniceElementowNieZaz(nrZazZiaren);elements.usunElementyZaz(nrZazZiaren);}
void usunNieZaznaczoneEl( IntList &nrZazZiaren){elements.ustawGraniceElementowZaz(nrZazZiaren);elements.usunElementyNieZaz(nrZazZiaren);}

void usunWybraneElementy(bool zaznaczone,bool usun_true_oznacz_false,IntList &nrZ);
void usunPunktyNiepotrzebne();
void zaznaczElementyZPPocz();

//z 2D do 3D

void wczytaj3D(const char *nazwa);

// zageszczenie Siatki Na granicy
void zageszczenieSiatkiNaGranicy();
void zageszczenieSiatkiPowierzchni();
void zapiszDoPlikuZageszczeniePunktowNaPow(bool zastosujGlobPrzesuniecia,int ileRazy,const char *nazwa);

//PSS
//void wczytajPlikNas(const char *nazwa,int mnoznik,bool wczytajWarunki);
void wczytajAnsys1(const char *nazwa,double powieksz,double transX,double transY,double transZ);
void wczytajPlikNas(const char *nazwa,int mnoznik,bool warunki,double px,double py,double pz);

void wczytajPlikNas(const char *nazwa,int mnoznik,bool wczytajWarunki);
void ustawWarunkiNaPodstawiePunktow(IntList &warunkiFace);
void wczytajPunktyPSS(const char *nazwa,int ktoreSerce);
void wczytajPunktyPSSnast(const char *nazwa,int ktoreSerce);
void zwrocLiczby(double &a,double &b,double &c,string &liczba);

void oznaczWezlyGraniczneNaPodstSasP();
void oznaczWezlyGraniczneNaPodstawieScianG();
void wyszukanieSasidnichElementowE();

void wyszukanieSasidniegoElementu(int nrEl,int numerSciany);
void wyszukanieSasidniegoElementu(int nrEl);

// zageszczenieSiatki PSS

void dodajNowePunktyNaKrawedziach();
void sprawdzWaznoscPunktow();
void uzupelnienieElementow();
void uzupelnienieElementuWybranego(int ktory);
void uzupelnienieElementowWybranych(IntList &elWybrane);


double ileWezlowG();

//Nowa warstwaPryzmatyczna
void tworzeniePryzm(double gruboscWarstwyPryzmatycznej, int ileoscWarstw,int procentowaWielkoscW,double wagaWygladzania, double doklOptyWygladzania);
void tworzeniePryzmWygladzanie(double grubosc,int procentowaWielkoscW);
void tworzeniePryzmNormalna(double grubosc,int ileWarstw,bool flagaJakosci);

void wyszukajScianyNaPryzmyNormalna(IntList &scianyPryzm);
int ileScianyNaPryzmy();

//void wyszukajNormalnePowierzchniDlaPunktowGranicznych(PunktList &wektoryNormalne, bool flagaJakosci);
void wyszukajNormalnePowierzchniDlaPunktowGranicznych();
void wyszukajWlementyScianGranicznychZPunktem(int numerWezla,IntList &elementySas,DoubleList &polaEleSas);
void wyliczWektorNormalnyDlaPunkty(Punkt &szukWektor,int numerWezla,IntList &elementySas,DoubleList &poleEleSas,PunktList &sprawdz,bool flagaJakosci);
void wyliczWektorNormalny(double &x1_zwWyn,double &y1_zwWyn,double &z1_zwWyn);


// warstwa przyscienna

	void uwtorzWartstwePrzysciennaDlaZiarna(int nrZ,double dlOdsuniecia);

	void znajdzWszystkieScianyGraniczne(IntList &sciany);
	void znajdzScianyGraniczneDlaZiarna(IntList &sciany,int nrZiarna,bool czyKrawPrzes);
	void obliczGradiety(IntList &sciany,double dlOdsuniecia);
	void zmianaWartGrad(IntList &sciany,double dlOdsuniecia);
	void wyznaczOdleglosc(IntList &sciany,double dlOdsuniecia);
	void tworzenieWarstwyPrzysciennej(IntList &sciany);
	void rysujCzworosciany();
	//void czyscWarstwePrzyscienna(){elementyCz.czysc(401,401);punktyCzworoscianow.czysc(100,100);}

    void oznaczWszytkieWezlyGraniczneNaBrzegach();

	//chwilowe
	IntList sciany;

// warstwa przysciena struct

//void tworzeniePryzm(double grubosc,int ileWarstw,bool zeSkalowaniem);
void wyszukajScianyNaPryzmy(IntList &scianyPryzm);
int rysujPryzmy(bool ujemne);
void rysujPryzme(int ktora);
void skalujIPrzesunCzworosciany(double gruboscWarstwy);
void skalujCzworosciany(double gX,double gY,double gZ);
void zapiszPunktyGraniczne();
void wczytajPunktyGraniczne();

double powSciany(int p1,int p2,int p3);
//rozpoznawanieP³aszczyzn

void rozpoznawaniePlaszczyzn();
void wczytajPrzerob(const char *nazwa);

//skasuj

//void wypiszPunkty(){for(int i=0,ileP=points.getIter();i<ileP;i++){points.getElement(i).wypiszWsp(m1);}}
void z();
//

//PSS wyszukajSciany
void wyszukajPowZew_start(int rodzajGranicy);
void wyszukajPowZew_wygladScianyGranicznej(int rodzajGranicy,char **mapE_F);
void wyszukajPowZew_wyszkiwaniePow(char **mapE_F,IntList &pkScian,bool *pointFlag,char obecnyNumer);
void wyszukajPowZew_zmienKoloryWezlow(char **mapE_F);
double* ustawaElementyNaSrodku(double stosunekWielkoscPlanszy);
int wyszukajEdgesInEl(int nrEl,int *nrEdges);
void wyszukajElementsWithEdge(int nrEd,int *nrElements);


double *zworcNormalnaDoFace(int p1,int p2,int p3); //dlugosc wektora przed dormalizacja , vn_x , vn_y , vn_z 

//ruch

//spawanie
double ruchSpawanie(double obecnyKrok, double od_start,double minPoprawy,double doX,double doY,double dl,double zmiejszaPrzes,double limit,double szerokosc);
bool powSinus(double doX,double doY,double dl,double zmiejszaPrzes,double limit,double szerokosc);
void maxZgranica();
void maxZgranica1();
double minZwezglyGraniczne();

//bfs
void ruchPrzestrzen(int ileWarstw,int obecny_krok, int od_krok,int ileKrok,double minPoprawy,double px0,double py0,double pz0,double px1,double py1,double pz1,double endX,double endY,double endZ);
void initPrzestrzenRuch(double px0,double py0,double pz0,double px1,double py1,double pz1);
void tworzenieGranicy2D_XY(int rodzaj);  
void wyszukajPunktuDoWygladzeniaNaPodstawiePunktow(IntList &punk,IntList &punktyZaz,int ilePok);
void optymalizacjaStosunekRdorNajlepszayWynikCaloscDlaWybranychPunktow(IntList &wybranePunkty,double dlKroku, bool bezGranicy);
IntList punktyPR;

void kopiujSiatke(Delanouy3D &a);
void polozeniePunktowWybranegoEl(int el,int localP,double *coor);
void siatkaStatyczna(double dl_x,double dl_y,double dl_z,double il_w_x,double il_w_y,double il_w_z);

//warunki abaqus
void warunkiNaPodstawieZiaren();
void warunkiNaPodstawieZiaren_pop_in_out(int numer);
void warunkiNaPodstawieZiaren_reset(int numer);
void tablicaFC_connect();
void getFacePoints(int face_id,int *points);
void podzielNaObjekty(int *tabBlock,int dl_tab);
void rozdzielSasiednieMaterialy(bool block);
void warunkiNaPodstawieDanych(int *matMap_in,int dl_mat,int *blockMap_in,int dl_block,int *war,int dl_war);

int materialNumber();
int materialIDs(int *tab,int l_tab);
void zapiszWarunkiPoczatkoweNaPodstawieDanych(const char *nazwa,double *tempB,int dl_tempB,int* Block,int dl_Block);
void wczytajWarunkiBrzegowe(int *warTab);
void zapiszWarunkiBrzegowe(int *warTab);
void creatFaceAndEdgeRememberBC();
int getlowestGroupID();
void oznaczPunktyNiepotrzebne();

void zmienWarunkeSciana(double r_max,double r_min,int zWar,int naWar);

};





#endif //DELANOUY_3D_H_





