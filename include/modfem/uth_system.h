#ifndef _uth_system_
#define _uth_system_
/******************************************************************************
File uth_system.h - interface with features directly dependend on OS.

Time measurments:
  time_C - (C standard procedure) to return time in seconds from some date

  time_init   - to initiate time measurments
  time_clock  - to return wall clock time from initialization
  time_CPU    - to return CPU  time from initialization
  time_print  - to print CPU and wall clock time from initialization
******************************************************************************/
#ifdef __cplusplus
extern "C"
{
#endif

/**--------------------------------------------------------
  time_init   - to initiate time measurments
---------------------------------------------------------*/
extern void time_init();

/**--------------------------------------------------------
  time_C - (C standard procedure) to return time in seconds from some date
---------------------------------------------------------*/
extern double time_C();

/**--------------------------------------------------------
  time_clock  - to return wall clock time from initialization
---------------------------------------------------------*/
extern double time_clock();

/**--------------------------------------------------------
  time_CPU    - to return CPU  time from initialization
 ---------------------------------------------------------*/
extern double time_CPU();

/**--------------------------------------------------------
  time_print  - to print CPU and wall clock time from initialization
 ---------------------------------------------------------*/
extern void time_print();


#ifdef __cplusplus
}
#endif

#endif
