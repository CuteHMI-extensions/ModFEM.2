#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

//#include <vld.h>	//Visual Leak Detector

//PARALLELISM
//#include<omp.h>

#include <limits>

//#ifndef _DEBUG
//#define NDEBUG 1
//#endif

#include <cassert>
#include <sstream>

#include <modfem/uth_log.h>

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */


//#define PRINT_LISTS 1

#ifndef IN
#define IN
#endif

#ifndef OUT
#define OUT
#endif

#ifndef NULL
#define NULL 0
#endif

//#if (__cplusplus <= 199711L)
//    #ifndef NULL
//        #define NULL NULL
//#warning "NULL defined as NULL"
//#endif
//#endif

typedef unsigned char BYTE;
#ifdef _WIN32
#ifndef __int8_t_defined
# define __int8_t_defined
typedef signed char		int8_t;
#endif


#ifndef __uint8_t_defined
# define __uint8_t_defined
typedef unsigned char		uint8_t;
#endif

#else

#include <stdint.h>
#endif

typedef BYTE Tsmall;
typedef int	Tind;
typedef unsigned int uTind;
typedef double Tval;
typedef uTind ID;
typedef long int Tlong;
typedef unsigned long int uTlong;

// FIXED error: narrowing conversion of ‘4294967295u’ from ‘unsigned int’ to ‘int’ inside { } [-Wnarrowing] //
const Tind UNKNOWN_INT = std::numeric_limits<Tind>::min();
// -------------------------------------------------------------------------------------------------------- //
const uTind UNKNOWN = std::numeric_limits<uTind>::max();
const uTind FIRST = 1;
const int X = 0;
const int Y = 1;
const int Z = 2;
const double SMALL = 1e-10;

enum eFlagIdx
{
	GROUP_ID = 0,
	B_COND = 0,
	F_TYPE = 1,
	REFINEMENT = 1,
	EL_TYPE = 2
};

enum eFaceFlag
{
	F_IN = -1,	// value is important
	F_OUT = 1,	//
	F_NORMAL = 0, // 000
	F_FLIPPED = 4 // 100
};
// eElemFlag is bit field. At each bit position where is 1
// normal vector for this face in element is pointing inside elem.
// If at position is 0 vector is pointing outside elem.
enum eElemFlag
{	// for tetrahedron (t4)
	EL_TYPEI = 5, //0101,
	EL_TYPEII = 12, //1100,
	EL_TYPE0 = 1, //0001,
	EL_TYPE1 = 2, //0010,
	EL_TYPE2 = 4, //0100,
	EL_TYPE3 = 8, //1000,
	EL_T4_ALLIN = 15, //1111
	EL_TYPE4 = 16, //10000
	EL_TYPE00 = 0, //0000
	// for prisms
	EL_TYPE_PI = 9, //01001
	EL_TYPE_PII = 25 //11001
};

enum eMMC
{
	ACTIVE = 1,
	INACTIVE = -1,
	REF_ISO = 8,
	UNREFINED = 0
};

enum eKind
{
	eVertex,
	eEdge,
	eFace,
	eElement
};

// Here we define way of handling complicated boundary geometry.
template <typename T> class FirstGeometryModule;
typedef FirstGeometryModule < double > GeometryModule;

template<class T>
inline void safeDelete (T *& p)
{
	if ( p != NULL )  {
		delete p;
		p = NULL;
	}
}

template<class T>
inline void safeDeleteArray (T *& p)
{
	if ( p != NULL )  {
		delete [] p;
		p = NULL;
	}
}

extern std::ostringstream & mmv_out_stream;


//template< class T1, class T2, bool T1bigger = (sizeof(T1) > sizeof(T2)) >
//struct static_max;
//
//template< class T1, class T2 >
//struct static_max<T1,T2,true>
//{
//   typedef T1 type;
//};

//template< class T1, class T2 >
//struct static_max<T1,T2,false>
//{
//    typedef T2 type;
//};

//template< class T1, class T2, bool T1bigger = (sizeof(T1) > sizeof(T2)) >
//struct static_min;

//template< class T1, class T2 >
//struct static_min<T1,T2,false>
//{
//   typedef T1 type;
//};

//template< class T1, class T2 >
//struct static_min<T1,T2,true>
//{
//    typedef T2 type;
//};


/** @}*/


#endif // COMMON_H_INCLUDED
