/*************************************************************
File contains procedures:
  lsr_mkb_core_fem_proj_sol_lev - to L2 project solution dofs between mesh levels
  lsr_mkb_core_fem_vec_norm - to compute a norm of global vector in parallel
  lsr_mkb_core_fem_sc_prod - to compute a scalar product of two global vectors
  lsr_mkb_core_fem_exchange_dofs - to exchange dofs between processors

In the current implementation the library requires to be linked with
the fem code providing the functions called in this file and declared
in file "sid_mkb/sih_mkb_fem_intf.h":
  fem_proj_sol_lev - to L2 project solution dofs between mesh levels
  fem_vec_norm - to compute a norm of global vector (possibly in parallel)
  fem_sc_prod - to compute a scalar product of two global vectors ( -||- )
  fem_exchange_dofs - to exchange dofs between processors
These functions are usually defined in solver interface modules, sid_....
(usually the file sid_xxx/sis_xxx_fem_intf.c)

It can be made more elegant by defining a data structure with function
pointers for the above functions. The pointers will be filled with
values during the execution.

History:
	05.2013 - Krzysztof Banas, initial version

*************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

/* functions headers of implemented interface */
#include "lsd_mkb_core/lsh_mkb_core_fem_intf.h"

// core module data structure
#include "lsd_mkb_core/lsh_mkb_core.h"


/* interface with the FEM code - routines called by the solver */
#include <modfem/si_mkb/sih_mkb_fem_intf.h>

/* LAPACK procedures */
#include <modfem/lin_alg_intf.h>


/*---------------------------------------------------------
  lsr_mkb_core_fem_proj_sol_lev - to L2 project solution dofs between mesh levels
---------------------------------------------------------*/
int lsr_mkb_core_fem_proj_sol_lev( /* returns: 1 - success; <=0 - error code*/
	int Solver_id,        /* in: solver data structure to be used */
	int Subsystem_id,        /* in: solver data structure to be used */
	int Ilev_from,    /* in: level number to project from */
		double* Vec_from, /* in: vector of values to project */
	int Ilev_to,      /* in: level number to project to */
		double* Vec_to    /* out: vector of projected values */
	)
{

  if(Ilev_from != Ilev_to){
	fem_proj_sol_lev(Solver_id,Subsystem_id,Ilev_from,Vec_from,Ilev_to,Vec_to);
  }

  return(1);
}

/*---------------------------------------------------------
  lsr_mkb_core_fem_vec_norm - to compute a norm of global vector in parallel
---------------------------------------------------------*/
double lsr_mkb_core_fem_vec_norm( /* returns: L2 norm of global Vector */
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,
  int Level_id,
  int Nrdof,            /* in: number of vector components */
  double* Vector        /* in: local part of global Vector */
  )
{

  int IONE=1; double vec_norm;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kbw
  printf("In lsr_mkb_core_fem_vec_norm: Solver_id %d, Subsystem_id %d, Level_id %d\n",
	 Solver_id, Subsystem_id, Level_id);
//kew*/

#ifdef PARALLEL
/*kbw
{
#pragma omp critical(printing)
  {
	int i;
	printf("Vector before parallel vector norm\n");
	for(i=0;i<Nrdof;i++) printf("%20.10lf",Vector[i]);
	printf("\n");
	//getchar();
  }
	  }
/*kew*/

	vec_norm = fem_vec_norm(Solver_id, Subsystem_id, Level_id, Nrdof, Vector);

#else
/*kbw
{
#pragma omp critical(printing)
  {
	int i;
	printf("Vector before sequential vector norm\n");
	for(i=0;i<Nrdof;i++) printf("%20.10lf",Vector[i]);
	printf("\n");
	//getchar();
  }
	  }
/*kew*/

	vec_norm = dnrm2_(&Nrdof, Vector, &IONE);

#endif

/*kbw
{
#pragma omp critical(printing)
  {
	printf("Computed norm %lf\n", vec_norm);
	getchar();
  }
	  }
/*kew*/


  return(vec_norm);

}


/*---------------------------------------------------------
  lsr_mkb_core_fem_sc_prod - to compute a scalar product of two global vectors
---------------------------------------------------------*/
double lsr_mkb_core_fem_sc_prod(
			   /* returns: scalar product of Vector1 and Vector2 */
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,
  int Level_id,
  int Nrdof,           /* in: number of vector components */
  double* Vector1,     /* in: local part of global Vector */
  double* Vector2      /* in: local part of global Vector */
  )
{

  int IONE=1; double sc_prod;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL

	sc_prod = fem_sc_prod(Solver_id,Subsystem_id,Level_id,Nrdof,Vector1,Vector2);

#else

	sc_prod = ddot_(&Nrdof, Vector1, &IONE, Vector2, &IONE);

#endif

  return(sc_prod);
}


/*---------------------------------------------------------
  lsr_mkb_core_fem_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
int lsr_mkb_core_fem_exchange_dofs(
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,
  int Level_id,
  double* Vec_dofs  /* in: vector of dofs to be exchanged */
)
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL

	fem_exchange_dofs(Solver_id, Subsystem_id, Level_id, Vec_dofs);

#endif

  return(1);
}

