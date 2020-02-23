
#ifndef _lsh_mkb_mumps_
#define _lsh_mkb_mumps_

// API for all mkb solvers with some constants...
#include <modfem/ls_mkb/lsh_mkb_intf.h>


/*** Data types ***/

/* definition of lst_mkb_solvers - data type for multi-level iterative solver */
typedef struct {

  int solver_id;           /* solver_id */
  int SM_and_LV_id;

  int monitor;


  int *coo_row;
  int *coo_col;
  double *coo_val;
  double  *rhs;
  int offset;
  int nnz;

} lst_mkb_mumps_solvers;

/*** Data types ***/

/* GLOBAL VARIABLES */
//extern int lsv_mkb_mumps_cur_solver_id;   /* ID of the current solver */
extern lst_mkb_mumps_solvers lsv_mkb_mumps_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */



#endif
