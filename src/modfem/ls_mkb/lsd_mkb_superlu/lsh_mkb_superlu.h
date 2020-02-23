
#ifndef _lsh_mkb_superlu_
#define _lsh_mkb_superlu_

// API for all mkb solvers with some constants...
#include <modfem/ls_mkb/lsh_mkb_intf.h>


/*** Data types ***/

/* definition of lst_mkb_solvers - data type for multi-level iterative solver */
typedef struct {

  int solver_id;           /* solver_id */
  int SM_and_LV_id;

  int monitor;


  int *crs_row;
  int *crs_col;
  double *crs_val;
  double  *rhs;
  int offset;
  int free_flag; //flag - 0 if not needed delete crs matrices

} lst_mkb_superlu_solvers;

/*** Data types ***/

/* GLOBAL VARIABLES */
//extern int lsv_mkb_superlu_cur_solver_id;   /* ID of the current solver */
extern lst_mkb_superlu_solvers lsv_mkb_superlu_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */



#endif
