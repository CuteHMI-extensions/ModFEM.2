/*************************************************************
lss_mkb_intf.c - file contains definitions of procedures:

- solver management
  lsr_mkb_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
  lsr_mkb_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
  lsr_mkb_destroy - to destroy a particular instance of the solver

- SM and LV management (using lah_intf.h API)
  lsr_mkb_create_matrix - to allocate space for a global system matrix
  lsr_mkb_clear_matrix - to initialize the structure of system matrix
  lsr_mkb_fill_assembly_table_int_ent - to create a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
  lsr_mkb_assemble_local_stiff_mat_with_table - to assemble entries to the global
						   stiffness matrix and the global load vector using the
						   provided local stiffness matrix, load vector and
						   the proper part of the global assembly table
  lsr_mkb_assemble_local_sm - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  lsr_mkb_free_matrix - to free space for solver data structure (system
					   matrix, preconditioner and possibly other)

//////////////////// utility to create data structures on GPU ///////////////
  lsr_create_solver_structures_accel
//////////////////// utility to get CRS parameters and pointers ///////////////
  lsr_get_crs_data

- preconditioner management (for iterative solvers)
  lsr_mkb_create_precon - to create preconditioner
  lsr_mkb_fill_precon - to prepare preconditioner
  (freeing preconditioner is in lsr_mkb_free_matrix - may be should be separated...)


- multigrid utility
  lsr_mkb_get_pdeg_coarse - to get enforced pdeg for the coarse mesh

History:
		08.2013 - Krzysztof Banas, initial version

*************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* provided interface of the solver - headers of routines defined in this file */
/* plus interface of three routines of direct solvers acting through MKB interface */
#include <modfem/ls_mkb/lsh_mkb_intf.h>

/* internal information for the mkb_core solver implementation module */
#include "./lsd_mkb_core/lsh_mkb_core.h"


// interface of linear algebra supporting package
#include <modfem/ls_mkb/lah_intf.h>


// special implementation for ns_supg - as switcheable extension
//#pragma warning "NS_SUPG_MKB_EXTENSIONS off in this version."

//#define NS_SUPG_MKB_EXTENSIONS
//#define AMG_MKB_EXTENSIONS

//TODO: is there a need for different ns supg than amg?
#ifdef NS_SUPG_MKB_EXTENSIONS
//#include "amg_ext/amg_ext.h"
#include "amg_ext/amg_ext.h"
//#include "./lsd_ns_supg_ext/lsh_ns_supg_ext_intf.h"
#endif

#ifdef AMG_MKB_EXTENSIONS
//#include "amg_ext/lad_amg/lah_amg_interface.h"
#include "amg_ext/amg_ext.h"
#endif


/* GLOBAL VARIABLES */
int   lsv_mkb_nr_solvers=0;        /* the number of solvers in the problem */
int   lsv_mkb_cur_solver_id;                /* ID of the current problem */
lst_mkb_solvers lsv_mkb_solver[LSC_MAX_NUM_SOLV];        /* array of solvers */

/*** Procedures ***/

#ifdef __cplusplus
extern "C" {
#endif


/*------------------------------------------------------------
  lsr_mkb_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
------------------------------------------------------------*/
int lsr_mkb_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_type, // type of solver (as defined in problem input file)
  int Parallel,   /* in: parameter specifying sequential (LSC_SEQUENTIAL)*/
				  /*     or parallel (LSC_PARALLEL) execution */
  int* Max_num_levels_p,  /* in: requested number of levels for multigrid: */
				  /*      1 - enforce single level solver */
				  /*     >1 - enforce multi-level solver */
				  /* out: actual number of levels !!! */
  int* Storage_type, /* in: requested storage type (NOT IMPLEMENTED YET) */
					  /* out: actual storage type */
  char* Filename,  /* in: name of the file with control parameters */
  int Max_iter, /* maximal number of iterations, -1 for values from Filename */
  int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
  double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */
  int Monitoring_level /* Level of output, -1 for Filename */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* increase the counter for solvers */
  lsv_mkb_nr_solvers++;

  /* set the current solver ID */
  lsv_mkb_cur_solver_id = lsv_mkb_nr_solvers-1;

  lsv_mkb_solver[lsv_mkb_cur_solver_id].solver_id = lsv_mkb_cur_solver_id;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].solver_type = Solver_type;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].parallel = Parallel;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_init: Solver_type %d, Solver id %d, Parallel %d, Nr_levels %d\n",
	 Solver_type, lsv_mkb_cur_solver_id, Parallel, *Max_num_levels_p);
  printf("Filename %s\n", Filename);
  printf("Monitor %d, Nr_iter %d, Conv_type %d, Conv_meas %15.12lf\n",
	 Monitoring_level, Max_iter, Error_type, Error_tolerance);
/*kew*/

  if(Solver_type == 0){

	// direct solver - single level
	lsv_mkb_solver[lsv_mkb_cur_solver_id].nr_levels = 1;
	// direct solver - storage decided by the solver
	lsv_mkb_solver[lsv_mkb_cur_solver_id].level[0].storage_type = LSC_STORAGE_UNDEFINED;
	// we keep direct solver interface as simple as possible (Filename can be null)
	lsr_mkb_direct_init(lsv_mkb_cur_solver_id, Filename, Monitoring_level);

  }
  // Solver_types 0..100 are reserved for mkb (including ns_supg_extensions)
  else if(Solver_type>0 && Solver_type<100){

	// mkb solver (including ns_supg_extensions)
	int storage_type = lsr_mkb_core_init(lsv_mkb_cur_solver_id, Solver_type,
					 Parallel,  Max_num_levels_p,
					 Filename, Max_iter, Error_type, Error_tolerance,
					 Monitoring_level);

	lsv_mkb_solver[lsv_mkb_cur_solver_id].nr_levels = *Max_num_levels_p;
	int ilev;
	for(ilev=0; ilev<lsv_mkb_solver[lsv_mkb_cur_solver_id].nr_levels; ilev++){
	  lsv_mkb_solver[lsv_mkb_cur_solver_id].level[ilev].storage_type = storage_type;
	}

  }
  else if(Solver_type == MULTI_GRID_AMG){

	int storage_type = lsr_mkb_core_init(lsv_mkb_cur_solver_id, Solver_type,
					 Parallel,  Max_num_levels_p,
					 Filename, Max_iter, Error_type, Error_tolerance,
					 Monitoring_level);
	*Max_num_levels_p = 1;
	lsv_mkb_solver[lsv_mkb_cur_solver_id].nr_levels = *Max_num_levels_p;
	lsv_mkb_solver[lsv_mkb_cur_solver_id].level[0].storage_type = storage_type;
  }
  else{

	// place for amg or anything else....

	printf("Unknown solver type %d in lsr_mkb_init. Exiting.\n", Solver_type);
	exit(-1);

  }

  *Max_num_levels_p = lsv_mkb_solver[lsv_mkb_cur_solver_id].nr_levels;
  *Storage_type = lsv_mkb_solver[lsv_mkb_cur_solver_id].level[0].storage_type;

  return(lsv_mkb_cur_solver_id);
}




/*---------------------------------------------------------
  lsr_mkb_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
---------------------------------------------------------*/
int lsr_mkb_solve( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
						/* <0 - error code */
	int Solver_id,      /* in: solver ID */
	int Comp_type,  /* in: indicator for the scope of computations: */
					/*   LSC_SOLVE - solve the system */
					/*   LSC_RESOLVE - resolve for the new right hand side */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* B,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int* Nr_iter, 	/* in:	the maximum iterations to be performed */
			/* out:	actual number of iterations performed */
	double* Toler, 	/* in:	tolerance level for chosen measure */
			/* out:	the final value of convergence measure */
	int Monitor,	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data, 3 - 2+iteration data */
	double* Conv_rate /* out: convergence rate */
	)
{

  int info;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_solve: Solver_type %d, Solver id %d, Comp_type %d\n",
	 solver_type, Solver_id, Comp_type);
/*kew*/

  if(solver_type == 0){

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[0].SM_and_LV_id;
	int nrdofgl = lsv_mkb_solver[Solver_id].level[0].nrdofgl;

	// direct solver
	// we keep direct solver interface as simple as possible
	info = lsr_mkb_direct_solve(Solver_id, Comp_type, SM_and_LV_id, nrdofgl, X, B, Monitor);

  }
  else if(solver_type == MKB_CORE_NS_SUPG_SOLVER){

#ifdef NS_SUPG_MKB_EXTENSIONS
//    int info = lsr_ns_supg_ext_solve(Solver_id, Comp_type, Ini_zero, X, B,
//				     Nr_iter, Toler, Monitor, Conv_rate);
//	  NS_SUPG_MKB_CORE_EXTENSIONS
	  lsr_mkb_core_solve(Solver_id, Comp_type,
			  &lsv_mkb_solver[Solver_id].level[0].SM_and_LV_id,
			  &lsv_mkb_solver[Solver_id].level[0].nrdofgl,
			  Ini_zero, X, B, Nr_iter, Toler, Monitor, Conv_rate);
	if(info<0){
	  // error handling
	}
#endif

  }
  // solver_types 1..100 are reserved for mkb_core
  else if(solver_type>0 && solver_type<100){

	int nr_levels = lsv_mkb_solver[Solver_id].nr_levels;

	// if single level
	if(nr_levels==1){

	  int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[0].SM_and_LV_id;
	  int nrdofgl = lsv_mkb_solver[Solver_id].level[0].nrdofgl;

/*kbw
  printf("\nbefore calling lsr_mkb_core_solve: Solver id %d, Comp_type %d, nr_levels %d\n",
	 Solver_id, Comp_type, nr_levels);
  printf("SM_and_LV_id\tNdof\n");
  printf("%d\t\t%d\n", SM_and_LV_id, nrdofgl);
/*kew*/

	  lsr_mkb_core_solve(Solver_id, Comp_type, &SM_and_LV_id, &nrdofgl, Ini_zero,  X, B,
			 Nr_iter, Toler, Monitor, Conv_rate);

	}
	else{ // if geometric multigrid - many levels managed by solver

	  int l_matrix_id[LSC_MAX_NUM_LEV];
	  int l_nrdofgl[LSC_MAX_NUM_LEV];
	  int ilev;
	  for(ilev=0;ilev<nr_levels;ilev++){
	l_matrix_id[ilev] = lsv_mkb_solver[Solver_id].level[ilev].SM_and_LV_id;
	l_nrdofgl[ilev] = lsv_mkb_solver[Solver_id].level[ilev].nrdofgl;
	  }

/*kbw
  printf("\nbefore calling lsr_mkb_core_solve: Solver id %d, Comp_type %d, nr_levels %d\n",
	 Solver_id, Comp_type, nr_levels);
  printf("SM_and_LV_id\tNdof\n");
  for(ilev=0;ilev<nr_levels;ilev++){

	printf("%d\t\t%d\n",
	   l_matrix_id[ilev], l_nrdofgl[ilev]);
  }
/*kew*/

	  // we create a list of matrix IDs
	  lsr_mkb_core_solve(Solver_id, Comp_type, l_matrix_id, l_nrdofgl, Ini_zero,  X, B,
			 Nr_iter, Toler, Monitor, Conv_rate);

	} // end if geometric multigrid

  } // end if classic mkb (nr_subsystems==1)
  else if(solver_type == MULTI_GRID_AMG){
	  int nr_levels = lsv_mkb_core_solver[Solver_id].subsystem[0].nr_level;
	  int l_matrix_id[LSC_MAX_NUM_LEV];
	  int l_nrdofgl[LSC_MAX_NUM_LEV];
	  int ilev;
	  for(ilev=0;ilev<nr_levels;ilev++){
		  l_matrix_id[ilev] = ilev;
		  l_nrdofgl[ilev] = 0;
	  }
	  l_nrdofgl[0] = lsv_mkb_solver[Solver_id].level[0].nrdofgl;
	  l_nrdofgl[ilev-1] = lsv_mkb_solver[Solver_id].level[0].nrdofgl;
	  lsv_mkb_solver[Solver_id].level[nr_levels-1].nrdofgl = lsv_mkb_solver[Solver_id].level[0].nrdofgl;
	lsr_mkb_core_solve(Solver_id, Comp_type, l_matrix_id, l_nrdofgl, Ini_zero,  X, B,
			 Nr_iter, Toler, Monitor, Conv_rate);

  }
  else{
	// place for or anything else....

	printf("Unknown solver type %d in lsr_mkb_solve. Exiting.\n", solver_type);
	exit(-1);
  }

  #ifdef NS_SUPG_MKB_EXTENSIONS
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){
	  double norm = lsr_ns_supg_ext_compres(Solver_id, X, lsv_mkb_solver[Solver_id].level[0].nrdofgl);
	  printf("\nAfter solving a system of linear equations - norm of residuum = %20.15lf\n",norm);
	  return 0;
  }
  #endif

#ifndef PARALLEL
/*ok_kbw*/
  // check solution
  int nr_levels = lsv_mkb_solver[Solver_id].nr_levels;
  int matrix_id = lsv_mkb_solver[Solver_id].level[nr_levels-1].SM_and_LV_id;
  int nrdofgl = lsv_mkb_solver[Solver_id].level[nr_levels-1].nrdofgl;
  // just compute the residual
  double *vtemp;
  vtemp = (double *) malloc(nrdofgl*sizeof(double));
  lar_compute_residual ( matrix_id, 1, 0, nrdofgl, X, NULL, vtemp);
  int idof; double res_norm=0.0; double res_norm_max = 0.0;
  for(idof = 0; idof < nrdofgl; idof++){
	double temp = vtemp[idof]*vtemp[idof];
	res_norm += temp;
	if(res_norm_max < temp) res_norm_max = temp;
  }
  free(vtemp);
  printf("\nAfter solving a system of %d linear equations - norm of residuum: \n\tL2 = %20.15lf, MAX = %20.15lf\n",
	 nrdofgl, sqrt(res_norm), sqrt(res_norm_max));
//*kew*/
#endif

  return(0);

}


/*---------------------------------------------------------
  lsr_mkb_destroy - to destroy a particular instance of the solver
  (should be done in reverse order wrt creation)
---------------------------------------------------------*/
int lsr_mkb_destroy(
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
  )
{

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_destroy: Solver_type %d, Solver id %d\n",
	 solver_type, Solver_id);
/*kew*/

  int info=1;

  /* set the current solver ID */
  lsv_mkb_cur_solver_id = Solver_id;

  if(solver_type == 0){

	// direct solver
	// we keep direct solver interface as simple as possible
	info = lsr_mkb_direct_destroy(lsv_mkb_cur_solver_id);

  }
  // solver_types 0..100 are reserved for mkb_core (including ns_supg_extensions)
  else if(solver_type>0 && solver_type<100){

	// mkb_core solver (including ns_supg_extensions)
	info = lsr_mkb_core_destroy(lsv_mkb_cur_solver_id);

  }
  else if(solver_type == MULTI_GRID_AMG){


  }
  else{

	// place for anything else....

	printf("Unknown solver type %d in lsr_mkb_init. Exiting.\n", solver_type);
	exit(-1);

  }

  /* decrease the counter for solvers */
  lsv_mkb_nr_solvers--;

  /* set the current solver ID */
  if(lsv_mkb_cur_solver_id == lsv_mkb_nr_solvers) lsv_mkb_cur_solver_id = lsv_mkb_nr_solvers-1;
  else{
	printf("Solver destroyed %d (from %d) is not the last solver in lsr_mkb_destroy!!!\n",
	   Solver_id, lsv_mkb_nr_solvers);
	exit(0);
  }

  return(info);
}


/*---------------------------------------------------------
  lsr_mkb_create_matrix - to allocate space for a global system matrix
---------------------------------------------------------*/
int lsr_mkb_create_matrix(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,    /* in: solver ID (used to identify the subproblem) */
  int Level_id,     /* in: level ID */
  int Storage_type, /* in: enforced storage type; if -1 (DEFAULT_STORAGE) */
					/*     storage type is decided based on Block_size */
  int Nrblocks,     /* in: number of DOF blocks */
  int Nrdof_glob,   /* in: total number of DOFs */
  int Block_size,   /* in: size of SM blocks (-1 - non-constant size) */
  int Max_sm_size,  /* in: maximal size of the stiffness matrix */
  int* Nrdofbl,	    /* in: list of numbers of dofs in a block */
  int* Posglob,	    /* in: list of global numbers of first dof */
  int* Nroffbl,	    /* in: list of numbers of off diagonal blocks */
  int** L_offbl	    /* in: list of lists of off diagonal blocks */
  )
{

  int i, j;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_create_matrix: solver %d, level %d, storage_type %d\n",
	 Solver_id, Level_id, Storage_type);
  printf("nrblocks %d, max_sm_size %d, nrdof_glob %d, block_size %d\n",
	 Nrblocks, Max_sm_size, Nrdof_glob, Block_size);
//*kew*/
/*kbw
// !!! offset 1 numbering of blocks moved to lad_block !!!
  printf("block_id,\tnrdofbl,\tposglob,\tnroffbl\tneighbors\n");
  for(i=0; i<Nrblocks; i++){
	printf("%d\t\t%d\t\t%d\t\t%d\t",i,Nrdofbl[i], Posglob[i], Nroffbl[i]);
	for(j=0;j<Nroffbl[i]; j++){
	  printf("%10d",L_offbl[i][j]);
	}
	printf("\n");
  }
//*kew*/

  int info;

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){


	#ifdef NS_SUPG_MKB_EXTENSIONS
		// entry to ns_supg extensions
	//	  mkb_level = &lsv_mkb_solver[Solver_id].subsystem[0].level[Level_id];
	//      mkb_level->Nrdofgl = Nrdof_glob;
	//      mkb_level->Nr_dof_blocks = Nrblocks;
		info = lsr_ns_supg_ext_create_matrix(Solver_id, Level_id, Nrblocks,
						 Nrdof_glob, Max_sm_size,
						 Nrdofbl, Posglob,
						 Nroffbl, L_offbl);
		lst_mkb_levels* mkb_level = &lsv_mkb_solver[Solver_id].level[Level_id];
		mkb_level->nrdofgl = Nrdof_glob;
	#endif

  } // end if ns_supg extensions
  // solver_type 0 - direct solver
  // solver_types 1..100 are reserved for mkb_core (including ns_supg_extensions)
  else if((solver_type>=0 && solver_type<100) || solver_type == MULTI_GRID_AMG){

	if(solver_type == 0 && Level_id != 0){
	  printf("Wrong level %d != 0 for direct solver through MKB at create_matrix! Exiting.\n", Level_id );
	  exit(-1);
	}
	if(solver_type == 0)lsr_mkb_direct_create(Solver_id, NULL, 0);

	lst_mkb_levels* mkb_level = &lsv_mkb_solver[Solver_id].level[Level_id];
	mkb_level->nrdofgl = Nrdof_glob;
	//mkb_level->nr_dof_blocks = Nrblocks;
	//mkb_level->block_size = Block_size;
	if(Storage_type!=LSC_STORAGE_UNDEFINED){
	  mkb_level->storage_type = Storage_type; // enforced storage type
	}

	int matrix_id = lar_allocate_SM_and_LV(mkb_level->storage_type, solver_type,
					   Nrdof_glob, Max_sm_size, Nrblocks,
					   Block_size, Nrdofbl, Posglob, Nroffbl, L_offbl);

	lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id = matrix_id;

	if(matrix_id < 0){
	  printf("error in allocating SM and LV (%d) in lsr_mkb_create_matrix. Exiting.\n",
		 matrix_id);
	  exit(-1);
	}

/*kbw
printf("Allocated SM and LV %d (%d), Nrdof_glob %d (%d)\n",
	   matrix_id, lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id,
	   Nrdof_glob, lsv_mkb_solver[Solver_id].level[Level_id].nrdofgl);
//*kew*/

  } // end if standard mkb (direct or iterative)
  else{

	// place for amg or anything else....

	printf("Unknown solver type %d in lsr_mkb_create_matrix. Exiting.\n", solver_type);
	exit(-1);

  }

  return(info);
}


/*---------------------------------------------------------
  lsr_mkb_clear_matrix - to initialize  system matrix
---------------------------------------------------------*/
int lsr_mkb_clear_matrix(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id,    /* in: level ID */
  int Comp_type          /* in: indicator for the scope of computations: */
						 /*   LSC_SOLVE - solve the system */
						 /*   LSC_RESOLVE - resolve for the new rhs vector */
  )
{

  int info = 0;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_clear_matrix: solver %d, level %d, comp_type %d\n",
	 Solver_id, Level_id, Comp_type);
/*kew*/

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;
  // branch if ns_supg extensions
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){


#ifdef NS_SUPG_MKB_EXTENSIONS
	// entry to ns_supg extensions
	info = lsr_ns_supg_ext_clear_matrix(Solver_id, Level_id, Comp_type);
#endif

  } // end if ns_supg extensions
  else if(solver_type==MULTI_GRID_AMG){
	  info = lar_initialize_SM_and_LV(lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id, LAC_SCOPE_SM_AND_LV);
  }
  else{

	if(solver_type == 0 && Level_id != 0){
	  printf("Wrong level %d != 0 for direct solver through MKB at clear_matrix! Exiting.\n",Level_id );
	  exit(-1);
	}

	int scope_type;
	if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
	else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
	else scope_type=LAC_SCOPE_SM_AND_LV; // default

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
	info = lar_initialize_SM_and_LV(SM_and_LV_id, scope_type);


  } // end if classic mkb


  return(info);
}

/**-----------------------------------------------------------
  lsr_mkb_fill_assembly_table_int_ent - to create a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
extern int lsr_mkb_fill_assembly_table_int_ent(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int* Assembly_table_int_ent /* part of the global assembly table */
  )
{

  int info=0;

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*kbw
  printf("\nEntering lsr_mkb_fill_assembly_table_int_ent: Solver_type %d, Solver id %d\n",
	 solver_type, Solver_id);
/*kew*/

  // branch if ns_supg extensions
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){

#ifdef NS_SUPG_MKB_EXTENSIONS
	// entry to ns_supg extensions
	// ???
#endif

  } // end if ns_supg extensions
  else{

	if(solver_type == 0 && Level_id != 0){
	  printf("Wrong level %d != 0 for direct solver through MKB! Exiting.\n", Level_id);
	  exit(-1);
	}

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
	info = lar_fill_assembly_table_int_ent( SM_and_LV_id,
						Nr_dof_bl,  L_bl_id, L_bl_nrdof,
						Assembly_table_int_ent);


  } // end if not ns_supg extensions

  return(info);
}



/*------------------------------------------------------------
  lsr_mkb_assemble_local_stiff_mat_with_table - to assemble entries to the
						   global stiffness matrix and the global load vector using
						   the provided local stiffness matrix and load vector
------------------------------------------------------------*/
int lsr_mkb_assemble_local_stiff_mat_with_table(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator for the scope of computations: */
						 /*   LSC_SOLVE - solve the system */
						 /*   LSC_RESOLVE - resolve for the new rhs vector */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int *Assembly_table_int_ent, /* part of the global assembly table */
  int* Local_to_global_int_ent, /* part of the global local_to_global table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
  )
{

  int info;

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*kbw
  printf("\nEntering lsr_mkb_assemble_local_stiff_mat_with_table: Solver_type %d, Solver id %d\n",
	 solver_type, Solver_id);
/*kew*/

  // branch if ns_supg extensions
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){

#ifdef NS_SUPG_MKB_EXTENSIONS
	// entry to ns_supg extensions
	// ???
#endif

  } // end if ns_supg extensions
  else{

	if(solver_type == 0 && Level_id != 0){
	  printf("Wrong level %d != 0 for direct solver through MKB! Exiting.\n", Level_id);
	  exit(-1);
	}

	int scope_type;
	if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
	else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
	else scope_type=LAC_SCOPE_SM_AND_LV; // default

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
	info = lar_assemble_SM_and_LV_with_table( SM_and_LV_id, scope_type,
						  Nr_dof_bl, Assembly_table_int_ent,
						  Local_to_global_int_ent,
						  Stiff_mat, Rhs_vect, Rewr_dofs);

  } // end if classic mkb not ns_supg extensions

  return(info);
}


/*------------------------------------------------------------
  lsr_mkb_assemble_local_sm - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
------------------------------------------------------------*/
int lsr_mkb_assemble_local_sm(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator for the scope of computations: */
						 /*   LSC_SOLVE - solve the system */
						 /*   LSC_RESOLVE - resolve for the new rhs vector */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
  )
{

  int info=0;

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*kbw
  printf("\nEntering lsr_mkb_assemble_local_stiff_mat: Solver_type %d, Solver id %d\n",
	 solver_type, Solver_id);
/*kew*/
/*kbw
  {
	int nrdofall=0; int i;
	printf("\nEntering lsr_mkb_assemble_local_sm: solver %d, level %d, comp_type %d, rewr_dofs %c\n",
	   Solver_id, Level_id, Comp_type, *Rewr_dofs);
	printf("nr_dof_bl %d\n", Nr_dof_bl);
	printf("block_index,\tblock_id,\tnrdofbl\n");
	for(i=0; i<Nr_dof_bl; i++){
	  printf("%d\t\t%d\t\t%d\t\t%d\n",i,L_bl_id[i], L_bl_nrdof[i]);
	  nrdofall += L_bl_nrdof[i];
	}
	printf("SM\n");
	for(i=0; i<nrdofall*nrdofall; i++) printf("%10.6lf",Stiff_mat[i]);
	printf("\n");
	printf("LV\n");
	for(i=0; i<nrdofall; i++) printf("%10.6lf",Rhs_vect[i]);
	printf("\n");
  }
/*kew*/

  // branch if ns_supg extensions
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){
	//TODO: check if the right place
#ifdef NS_SUPG_MKB_EXTENSIONS
	// entry to ns_supg extensions
	info = lsr_ns_supg_ext_assemble_local_sm(Solver_id, Level_id, Comp_type,
						 Nr_dof_bl, L_bl_id, L_bl_nrdof,
						 Stiff_mat, Rhs_vect, Rewr_dofs);
#endif

  } // end if ns_supg extensions
  else{

	if(solver_type == 0 && Level_id != 0){
	  printf("Wrong level %d != 0 for direct solver through MKB! Exiting.\n",Level_id );
	  exit(-1);
	}

	int scope_type;
	if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
	else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
	else scope_type=LAC_SCOPE_SM_AND_LV; // default

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
	info = lar_assemble_SM_and_LV( SM_and_LV_id, scope_type,
				   Nr_dof_bl, L_bl_id, L_bl_nrdof,
				   Stiff_mat, Rhs_vect, Rewr_dofs);



  } // end if not ns_supg extensions

  return(info);
}

//////////////////// utility to create data structures on GPU ///////////////
int lsr_create_solver_structures_accel( /* returns: >0 - success code, <0 - error code */
  int Solver_id,  /* in: solver ID */
  int Level_id    /* in: level ID */
					){


  // go down to lad_... modules and pass created arrays to tmr routine
  /* tmr_create_solver_structures_accel_crs(  */
  /* 					 int Nrdof_glob, */
  /* 					 int nnz, */
  /* 					 int* crs_col_ind, */
  /* 					 int* crs_row_ptr, */
  /* 					 double* crs_val, */
  /* 					 double* rhs */
  /* 					 itd, itp ); */

  int info=0;

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

  // branch if ns_supg extensions
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){
	// not yet on GPU....
  }
  else{

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
	info = lar_create_solver_structures_accel(SM_and_LV_id);

  }

  return(info);

}


//////////////////// utility to get CRS parameters and pointers ///////////////
int lsr_get_crs_data( /* returns: >0 - success code, <0 - error code */
  int Solver_id,  /* in: solver ID */
  int Level_id,    /* in: level ID */
  int* Nrdof_glob_p, /* output pointers: nrdof, nnz etc. - self explanatory */
  int* Nnz_p,
  int** Crs_col_ind_p,
  int** Crs_row_ptr_p,
  double** Crs_val_p,
  double** Rhs_p
			 )
{
  int info=0;

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

  // branch if ns_supg extensions
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){
	// not yet on GPU....
  }
  else{

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
	info = lar_get_crs_data(SM_and_LV_id, Nrdof_glob_p, Nnz_p,
				Crs_col_ind_p, Crs_row_ptr_p, Crs_val_p, Rhs_p );

  }

  return(info);
}

/*---------------------------------------------------------
  lsr_mkb_free_matrix - to free space for solver data structure (system
					   matrix, preconditioner and possibly other)
---------------------------------------------------------*/
int lsr_mkb_free_matrix(
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_free_matrix: Solver_type %d, Solver id %d\n",
	 solver_type, Solver_id);
/*kew*/

  if(solver_type == 0){

	// direct solver
	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[0].SM_and_LV_id;
	lar_free_SM_and_LV(SM_and_LV_id);
	lsr_mkb_direct_free(SM_and_LV_id);

  }
  else if(solver_type == MKB_CORE_NS_SUPG_SOLVER){

#ifdef NS_SUPG_MKB_EXTENSIONS
	// entry to ns_supg extensions
	lsr_ns_supg_ext_free_matrix(Solver_id);
#endif

  }
  else{

	/* in a big loop over mesh levels */
	int ilev;
	for(ilev=0;ilev<lsv_mkb_solver[Solver_id].nr_levels;ilev++){

	  int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[ilev].SM_and_LV_id;

	  lsr_mkb_core_destroy_precon(Solver_id, ilev, SM_and_LV_id);

	  lar_free_SM_and_LV(SM_and_LV_id);

	}

  } // end if classic mkb

  return(0);
}


/*---------------------------------------------------------
lsr_mkb_create_precon - to create preconditioner
---------------------------------------------------------*/
int lsr_mkb_create_precon( /* returns:   >0 number of diagonal blocks */
						  /*	       <=0 - error */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id           /* in: level ID */
  )
{

  int info=1;

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_create_precon: Solver_type %d, Solver id %d\n",
	 solver_type, Solver_id);
/*kew*/

  // solver_types 1..100 are reserved for mkb_core (including ns_supg_extensions)
  if(solver_type>0 && solver_type<100){

	int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;

	lsr_mkb_core_create_precon(Solver_id, Level_id, SM_and_LV_id);


  } // end if classic mkb_core (NOT one of direct solvers through MKB interface)
  else if(solver_type==0){
	// do nothing - direct solvers do not use preconditioning...
  }
  else{

	// place for amg or anything else....

	printf("Wrong solver type %d in lsr_mkb_create_precon. Exiting.\n", solver_type);
	//exit(-1);

  }

  return(info);
}

/*---------------------------------------------------------
  lsr_mkb_fill_precon - to prepare preconditioner by factorizing the stiffness
					matrix, either only diagonal blocks or block ILU(0)
---------------------------------------------------------*/
int lsr_mkb_fill_precon(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id           /* in: level ID */
  )
{

  int info = 0;

/*kbw
  printf("before calling lar_fill_preconditioner: solver %d, level %d\n",
	 Solver_id, Level_id);
/*kew*/

  int solver_type = lsv_mkb_solver[Solver_id].solver_type;

/*ok_kbw*/
  printf("\nEntering lsr_mkb_fill_precon: Solver_type %d, Solver id %d\n",
	 solver_type, Solver_id);
/*kew*/

  // solver_types 1..100 are reserved for mkb_core (including ns_supg_extensions)
  if(solver_type>0 && solver_type<100){

  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){

		#ifdef NS_SUPG_MKB_EXTENSIONS
		// entry to ns_supg extensions
		info = lsr_ns_supg_ext_fill_precon(Solver_id, Level_id);
		#endif

	} // end if ns_supg extensions
	else{
		int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
		info = lsr_mkb_core_fill_precon(Solver_id, Level_id, SM_and_LV_id );
	}

  } // end if classic mkb_core not ns_supg extensions
  else if(solver_type==0){
	// do nothing - direct solvers do not use preconditioning...
  }

	else if(solver_type==MULTI_GRID_AMG)
	{
		int SM_and_LV_id = lsv_mkb_solver[Solver_id].level[Level_id].SM_and_LV_id;
		info = lar_fill_preconditioner( SM_and_LV_id );
		//TODO: fill all the levels nr dof?
		lsv_mkb_solver[lsv_mkb_cur_solver_id].level[lsv_mkb_solver[lsv_mkb_cur_solver_id].nr_levels-1].nrdofgl =
				lsv_mkb_solver[lsv_mkb_cur_solver_id].level[0].nrdofgl;

	}
  return(info);
}



///////////////////////////// GEOMETRIC MULTIGRID UTILITY //////////////////////

/*---------------------------------------------------------
  lsr_mkb_get_pdeg_coarse - to get enforced pdeg for the coarse mesh
---------------------------------------------------------*/

int lsr_mkb_get_pdeg_coarse( // returns: enforced pdeg for the coarse mesh
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id /* in: level number */
  )
{

  if(lsv_mkb_solver[Solver_id].nr_levels >1){
	return lsr_mkb_core_get_pdeg_coarse(Solver_id, Level_id);
  }

  return(1);
}

#ifdef __cplusplus
}
#endif
