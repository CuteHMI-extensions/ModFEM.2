/************************************************************************
File tmh_intf.h - generic(?) interface to thread management routines

Contains declarations of constants and interface routines:
  tmr_init_multithreading - to initialize (multi)thread management


------------------------------  			
History:        
	02.2013 - Krzysztof Banas, initial version

*************************************************************************/

#ifndef _tmh_intf_
#define _tmh_intf_

#include"pdh_intf.h"

#ifdef __cplusplus
extern "C"{
#endif

/** Constants */

#define TMC_PRINT_NOT 0
#define TMC_PRINT_ERRORS 1
#define TMC_PRINT_INFO 2
#define TMC_PRINT_ALLINFO 3


/*** FUNCTION DECLARATIONS - headers for external functions ***/

/**--------------------------------------------------------
  tmr_init_multithreading - to initialize (multi)thread management
---------------------------------------------------------*/
extern int tmr_init_multithreading(
  char* Work_dir,
  int Argc, 
  char **Argv,
  FILE *Interactive_input,
  FILE *Interactive_output,
  int Control,  // not used !
  int Monitor
 );


#ifdef __cplusplus
}
#endif 

#endif



