
#ifndef _lsh_mkb_pardiso_
#define _lsh_mkb_pardiso_

#include "mkl.h"

// API for all mkb solvers with some constants...
#include <modfem/ls_mkb/lsh_mkb_intf.h>


/*** DATA TYPES ***/
typedef struct{
  int mtype;  // Matrix type
  int maxfct; // Maximum number of numerical factorizations.
  int mnum;   // Which factorization to use.
  int msglvl; // Print statistical information in file

  int iparm1; //use solver defaults?
  int iparm2; //fill-in reducing ordering.
  int iparm3; //currently is not used.
  int iparm4; //preconditioned CGS.
  int iparm5; //user permutation.
  int iparm6; //write solution on x.

  int iparm8; //iterative refinement step.
  int iparm9; //must be set to 0.
  int iparm10; //pivoting perturbation.
  int iparm11; //scaling vectors.
  int iparm12; //must be set to 0.
  int iparm13; //improved accuracy using (non-)symmetric weighted matchings.

  int iparm18; //-1
  int iparm19; //MFlops of factorization.

  int iparm21; //pivoting for symmetric indefinite matrices.

  int iparm27; //matrix checker. default is 0
  int iparm28; //single or double precision of MKB_PARDISO.
  int iparm60; //version of MKB_PARDISO.


} lst_mkb_pardiso_config;

/* definition of lst_mkb_solvers - data type for multi-level iterative solver */
typedef struct {

  int solver_id;           /* solver_id */
  int SM_and_LV_id;

  int monitor;

  lst_mkb_pardiso_config pardiso_config; //PARDISO control data
  bool isPardisoConfigInitializedext;
  int pt[64];
  int *iparm;

  int *crs_row;
  int *crs_col;
  double *crs_val;
  double  *rhs;
  int Ndof;
  int offset;
  int analysis;
  int free_flag; //flag - 0 if not needed delete crs matrices

} lst_mkb_pardiso_solvers;

/*** Data types ***/

/* GLOBAL VARIABLES */
//extern int lsv_mkb_pardiso_cur_solver_id;   /* ID of the current solver */
extern lst_mkb_pardiso_solvers lsv_mkb_pardiso_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */

#endif
