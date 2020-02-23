#ifndef _TICCO_MESH_3D_H_
#define _TICCO_MESH_3D_H_


#include "texty.h"
#include "wektor.h"
#include "plansza.h"
#include "rozrost.h"
#include "losowanie.h"
#include "delanouy3d.h"
#include "figury.h"

class TicooMesh3D{


time_t losowa;
double wielkoscTabX;
double wielkoscTabY;
double wielkoscTabZ;
double wielkoscTabMax;


Plansza plansza;
Rozrost rozrost;
Kula kula;
IntList Trasa;

Losowanie losowane;
Delanouy3D siatka,siatka1;


//oznaczenie plyny


//remesh
int aktualneIDElementu;


public:

bool empty();
void zmienWarunkeSciana(double r_max,double r_min,int zWar,int naWar);
//void mmr_set_fluid(int el_id,char a);
//int mmr_is_fluid(int el_id);

void count_bc_surface(FILE *out);
void mesh_quality_stats(FILE *out);
//implementacja interfejsu mms_remesh
void utr_smoothing(double l_step,int repetition);
void mmr_get_vec_norm_at_node(int node,double *vec_norm);
int mmr_node_status(int Node);
int mmr_grups_number();
int mmr_grups_ids(int *tab,int l_tab);
void mmr_split_into_blocks_add_contact(const char *workdir,int *tabMat,int l_mat,int *tabBlock,int l_block,int *war,int l_war,double *tempBlock,int l_tempBlock,int *groups_fluid,int l_groups_fluid);
int mmr_get_fa_el_bc_connect(int face_id,int *el_id);
int mmr_el_groupID(int idEl);
void mmr_el_node_coor(int Nel, int* Nodes, double* Node_coor);
int mmr_el_eq_neig(int El, int* Neig, int* Neig_sides);
void mmr_node_coor(int Node,double *Coor);
int mmr_el_faces(int El,int *Face, int *Orient);
int mmr_el_edges(int El,int *Edges);
int mmr_edge_elems(int Ed,int *Edge_elems);
void mmr_edge_nodes(int Ed,int *Edge_nodes);
void mmr_fa_node_coor(int Fa,int *Nodes,double *Node_coor);
void mmr_fa_edges(int Fa,int *Fa_edges,int *Ed_orient);
double mmr_el_hsize(int El);
void mmr_fa_area(int Fa,double *Area,double *Vec_norm);
void getElwithFace(int Fa,int *El,int *Neig_sides);
//void mmr_set_new_and_get_old_uk_val(int idEl,int idGP,double *new_uk_val,double *old_uk_val);
void mmr_get_coor_from_motion_element(int idEl,int idLP,double *coor,int flagaSiatki);    
//flagaSiatki = 0 - s0
//              1 - s1
void copyMesh(int flaga);
//flaga = 0 - s1->s0
//        1 - s0->s1

void creatEdgeFace();
void czyscEdgeFace();
void czyscEdge();
void czyscFace();

void delElement(int el){siatka.elements.deleteElement(el);}
void delPoints(int nP){siatka.points.usunElement(nP);}
//void delFace(int nP){siatka.faces.usunElement(nP);}
//void delEdge(int nP){siatka.edges.usunElement(nP);}

//

TicooMesh3D(int wielkoscTab);
int getNumberElements() {return siatka.elements.getIter();}
int getNumberPoints() {return siatka.points.getIter();}
int getNumberFaces(){return siatka.faces.getIter();}
int getNumberEdges(){return siatka.edges.getIter();}
double getWielkoscTabX(){return wielkoscTabX;}
double getWielkoscTabY(){return wielkoscTabY;}
double getWielkoscTabZ(){return wielkoscTabZ;}
double getWielkoscTabMax(){return wielkoscTabMax;}


//plansz
void wczytajAC(const char *nazwa);
void setPlansza(double wX,double wY,double wZ,bool periodyczne,bool wlaczAC);
void elementyDoZiarnaAC();
void elementyDoZiarnaPunktyGraniczne();

//rozrost ziarna
void rozrostLosuj(int iloscZiaren);
void rozrostLosujRange(int iloscZiaren,int range);
void rozrostZiarenGranicy(int wspPrawdop);
void pomniejszenieZiarnaGranicy(int wspPomniejszaniaZiaren);
void rozrostZiarnaZOtoczeniem(int wspPrawdop,int wpsOtoczenia,int ileRozrost,int ileOtocz,int ilePetla);
void komplekosoweTworzenieZiaren(int wspPrawdop,int ilePetla,int ileRozrost,int ileOtocz,int wspPomniejszaniaZiaren);
void granica();

//losuj
void losowanieDodatkowychPunktow(bool pobierzPunkty,bool losuj,int ileP,int rozdz,int odG,int odP,bool flagaPer,double ax,double b0,int rSzuk);
void ustawPunktySpaw1();
void zmianaKolejnosciPunktow(int ileRazy,int wielkoscKomurki,int ileNaKomurke);

//siatka
void delanouy(bool dokladneWyszukanie);
void reMes();
void wygladzanieReMES(int ileRazy,double waga,bool wyglZWaga);
double objetoscT(FILE *out);
void sasiedniePunkty();
void sasiednieElementy();
void wygladzanieLaplaceWaga(int ileR,double waga);
void wygladzanieLaplace(int ileR);
void optymalizacjaStosunekRdor(int ile,double dokladnosc);
void poprawWszystkieElFAST(double docelowyStosunekRdor,bool narozneTez);
void conwersjaDoTicooMesh3DZapisu();
int getFace_bc(int nrFace){return siatka.faces.getElement(nrFace).getBC();}
void setFace_bc(int nrFace,int bc){siatka.faces.getElement(nrFace).setBC(bc);}
void init_all_change(int a);
//operacje na siatkach //skalowanie przesuwanie itp.

void siatkaNaSrodek(double powiekszPrzestrzen);


//zapis
//void wczytajNAS(const char *nazwa,int mnoznik);
void wczytajAbaqus(const char *nazwa,int mnoznik,double px,double py,double pz);
void wczytajNAS(const char *nazwa,int mnoznik,double px,double py,double pz);
void zapiszDoPlikuNas(const char *nazwa,double podzielZ, bool brick,double tX,double tY,double tZ,double podzielAll);

void zapiszAbqus(const char *nazwa);
void wczytajPunktyGraniczne(const char *nazwa);
void zapiszPunkty(const char *nazwa);
void zapisParaView(const char *nazwa);
void wczytajPunktyPSS(const char *nazwa);
void wczytajPunktyPSSnast(const char *nazwa);
void zapiszPSS1DoPlikuNasWOparciuOPowTetra(const char *nazwa);
void zapiszPSS1DoPlikuNasWOparciuOPowHybrydPSS(const char *nazwa);
void ZapiszDoPlikuNasZPrzesunieciem_BC_Face(const char *nazwa);
  
//ruch
double ruchSpawanie(double obecnyKrok, double krok_start,double minPoprawy,double doX,double doY,double dl,double zmiejszaPrzes,double limit,double szerokosc);
void ruchPrzestrzen(int ileWarstw,int obecny_krok, int od_krok,int ileKrok,double minPoprawy,double px0,double py0,double pz0,double px1,double py1,double pz1,double endX,double endY,double endZ);


//tworzenie siatki
void createCube(const char *nazwa,int node_x,int node_y,int node_z,double size_x,double size_y,double size_z,double divide,int *warunki);      
};

#endif //_TICCO_MESH_3D_H_
