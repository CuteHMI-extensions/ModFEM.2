/*
 * -- SuperLU MT routine (version 2.1) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * These macros define which machine will be used.
 *
 * Modified:  March 20, 2013  version 2.1
 */

#ifndef __SUPERLU_MACHINES /* allow multiple inclusions */
#define __SUPERLU_MACHINES

#define PTHREAD         0
#define OPENMP		1
#define SGI	        2
#define ORIGIN	        3
#define DEC	        4
#define CRAY_PVP	5
#define SUN             6

#ifdef __PTHREAD
#define MACH PTHREAD
#endif

#ifdef __OPENMP
#define MACH OPENMP
#endif

#ifdef __SOLARIS
#define MACH SUN 
#endif

#ifdef __SGI
#define MACH SGI 
#endif

#ifdef __ORIGIN
#define MACH ORIGIN 
#endif

#ifdef __DEC
#define MACH DEC 
#endif

#ifdef __CRAY
#define MACH CRAY_PVP 
#endif

///
/// HACK/FIX:
/// Because of wrong results from multithreaded SuperLU solver,
/// OPENMP is set here, but disabaled when called.
/// So there is __NO__ multithreaded solver enabled,
/// regardless from flag set here.
#define MACH OPENMP


#endif /* __SUPERLU_MACHINES */
