
/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


// child t4 number = [refType(0-2)][whichFace(0-4)][whichT4?]
/*static const int   offspringAtFace[3][4][4]= {
  {{0,1,2,4},
   {0,1,3,5},
   {0,2,3,6},
   {1,2,3,7}},

  {{0,1,2,4},
   {0,1,3,7},
   {0,2,3,5},
   {1,2,3,6}},

  {{0,1,2,4},
   {0,1,3,5},
   {0,2,3,7},
   {1,2,3,6}},
   };
*/
// sons-faces vertices indexes in array V[10]
// [son number(1-4)][ref kind(0-2)(eRef_67,49,58)][face vertex(0-2)]
// sons 1-4 are immune to ref kind
static const unsigned int faceVrtsIdx[9][3][3] = {
  {{UNKNOWN,UNKNOWN,UNKNOWN},
   {UNKNOWN,UNKNOWN,UNKNOWN},
   {UNKNOWN,UNKNOWN,UNKNOWN}},// son 0 is edge!
  {{4,5,7},{4,5,7},{4,5,7}}, // face son 1
  {{4,6,8},{4,6,8},{4,6,8}}, // face son 2
  {{5,6,9},{5,6,9},{5,6,9}}, // face son 3
  {{7,8,9},{7,8,9},{7,8,9}}, // face son 4
  {{4,6,7},{4,5,9},{4,5,8}}, // face son 5
  {{5,6,7},{4,6,9},{5,6,8}}, // face son 6
  {{6,7,8},{4,7,9},{5,7,8}}, // face son 7
  {{6,7,9},{4,8,9},{5,8,9}}  // face son 8
};

// sons-elems vertices indexes in array v[10]
// [son number(0-7)][ref kind(0-2)][elem vertex(0-3)]
// elems 1-4 are immune to ref kind(ref type)
static const int elemVrtsIdx[8][3][4] = {
  {{0,4,5,7},{0,4,5,7},{0,4,5,7}}, // elem son 0, type as parent
  {{4,1,6,8},{4,1,6,8},{4,1,6,8}}, // elem son 1, type as parent
  {{5,6,2,9},{5,6,2,9},{5,6,2,9}}, // elem son 2, type as parent
  {{7,8,9,3},{7,8,9,3},{7,8,9,3}}, // elem son 3, type as parent
  {{4,6,5,7},{4,6,5,9},{4,6,5,8}}, // elem son 4, type II
  {{4,6,7,8},{4,5,7,9},{4,5,7,8}}, // elem son 5, type I
  {{5,7,6,9},{4,8,6,9},{5,8,6,9}}, // elem son 6, type II
  {{6,7,8,9},{4,7,8,9},{5,7,8,9}}  // elem son 7,
};

// numbers (0-7) of element-sons neighbouring with internal face
// [face son nr(1-8)][ref type(0-2)][neig (0-1)]
// as you can see faces 5-9 are immune to ref type
static const unsigned int innerFaceNeigs[9][3][2] = {
  {{UNKNOWN,UNKNOWN},{UNKNOWN,UNKNOWN},{UNKNOWN,UNKNOWN}}, // son 0 is edge!!
  {{0,4},{0,5},{0,5}}, // face 1
  {{1,5},{1,6},{1,4}}, // face 2
  {{2,6},{2,4},{2,6}}, // face 3
  {{3,7},{3,7},{3,7}}, // face 4
  {{4,5},{4,5},{4,5}}, // face 5
  {{4,6},{4,6},{4,6}}, // face 6
  {{5,7},{5,7},{5,7}}, // face 7
  {{6,7},{6,7},{6,7}}, // face 8
};


// child elements faces orientations (auto-generated: see below)
// [subRefKind][parent type][child no]
//static const BYTE t4ChildType[3][EL_T4_ALLIN][8];
// Auto-generating t4childType.
/*
  for(int kind(0);kind<3;++kind) {
  for(int parTyp(0); parTyp < EL_T4_ALLIN; ++parTyp) {
  for(int childNo(0); childNo<8;++childNo){
  t4childtype[kind][parTyp][childNo]=-1;
  }
  }
  }       
*/
/**  @} */
