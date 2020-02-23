// lsh_mkb_core_fem_intf.h - interface required by mkb_core 
//   implementations are usually provided by some adapter like e.g.
//   the ones in lsd_mkb or lsd_ns_supg directories, that implement
//   an interface known to FEM code - in particular in ModFEM
//   the FEM part of code interacting with the solver is a solver
//   interface module with the current practice of one-to-one
//   correspondence between interface and adapter modules (e.g.
//   sid_mkb for lsd_mkd, sid_ns_supg for lsd_ns_supg etc.)

#ifndef _lsh_mkb_core_fem_intf_
#define _lsh_mkb_core_fem_intf_

#include <stdio.h>

/**--------------------------------------------------------
  lsr_mkb_core_fem_proj_sol_lev - to L2 project solution dofs between mesh levels
---------------------------------------------------------*/
int lsr_mkb_core_fem_proj_sol_lev( /* returns: 1 - success; <=0 - error code*/
	int Solver_id,        /* in: solver data structure to be used */
        int Subsystem_id,  /* in: subsystem data structure to be used */
	int Ilev_from,    /* in: level number to project from */
        double* Vec_from, /* in: vector of values to project */
	int Ilev_to,      /* in: level number to project to */
        double* Vec_to    /* out: vector of projected values */
  );

/**--------------------------------------------------------
  lsr_mkb_core_fem_vec_norm - to compute a norm of global vector in parallel
---------------------------------------------------------*/
extern double lsr_mkb_core_fem_vec_norm( /* returns: L2 norm of global Vector */
  int Solver_id,      /* in: pointer to solver data structure to be passed
	                       to data structure dependent routines */
  int Subsystem_id,  /* in: subsystem data structure to be used */
  int Level_id,	/* in: index of the mesh (level) */
  int Nrdof,            /* in: number of vector components */
  double* Vector        /* in: local part of global Vector */
  );

/**--------------------------------------------------------
  lsr_mkb_core_fem_sc_prod - to compute a scalar product of two global vectors 
---------------------------------------------------------*/
extern double lsr_mkb_core_fem_sc_prod( 
                       /* returns: scalar product of Vector1 and Vector2 */
  int Solver_id,      /* in: pointer to solver data structure to be passed
	                       to data structure dependent routines */
  int Subsystem_id,  /* in: subsystem data structure to be used */
  int Level_id,	/* in: index of the mesh (level) */
  int Nrdof,           /* in: number of vector components */
  double* Vector1,     /* in: local part of global Vector */
  double* Vector2      /* in: local part of global Vector */
  );

/**--------------------------------------------------------
  lsr_mkb_core_fem_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
extern int lsr_mkb_core_fem_exchange_dofs(
  int Solver_id,      /* in: pointer to solver data structure to be passed
	                       to data structure dependent routines */
  int Subsystem_id,  /* in: subsystem data structure to be used */
  int Level_id,	/* in: index of the mesh (level) */
  double* Vec_dofs  /* in: vector of dofs to be exchanged */
  );

#endif
