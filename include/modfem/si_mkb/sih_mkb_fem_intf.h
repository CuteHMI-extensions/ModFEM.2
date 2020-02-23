/*************************************************************
sih_mkb_fem_intf.h - functions from the FEM code called by the block
   multigrid preconditioned Krylow methods iterative solver module:
 fem_proj_sol_lev - to L2 project solution dofs between mesh levels
 fem_vec_norm - to compute a norm of global vector in parallel
 fem_sc_prod - to compute a scalar product of two global vectors
 fem_exchange_dofs - to exchange dofs between processors
*************************************************************/


#ifndef _sih_mkb_fem_intf_
#define _sih_mkb_fem_intf_


#ifdef __cplusplus
extern "C" {
#endif

/**--------------------------------------------------------
  fem_proj_sol_lev - to L2 project solution dofs between mesh levels
---------------------------------------------------------*/
extern int fem_proj_sol_lev( /* returns: >=0 - success; <0 - error code*/
	int Solver_id,        /* in: solver data structure to be used */
	int Subsystem_id,        /* in: solver data structure to be used */
	int Ilev_from,    /* in: level number to project from */
        double* Vec_from, /* in: vector of values to project */
	int Ilev_to,      /* in: level number to project to */
        double* Vec_to    /* out: vector of projected values */
	);

/**--------------------------------------------------------
  fem_vec_norm - to compute a norm of global vector in parallel
---------------------------------------------------------*/
extern double fem_vec_norm( /* returns: L2 norm of global Vector */
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  int Nrdof,            /* in: number of vector components */
  double* Vector        /* in: local part of global Vector */
  );

/**--------------------------------------------------------
  fem_sc_prod - to compute a scalar product of two global vectors 
---------------------------------------------------------*/
extern double fem_sc_prod( 
                       /* retruns: scalar product of Vector1 and Vector2 */
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  int Nrdof,           /* in: number of vector components */
  double* Vector1,     /* in: local part of global Vector */
  double* Vector2      /* in: local part of global Vector */
  );

/**--------------------------------------------------------
  fem_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
extern int fem_exchange_dofs(
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  double* Vec_dofs  /* in: vector of dofs to be exchanged */
  );

#ifdef __cplusplus
}
#endif

#endif
