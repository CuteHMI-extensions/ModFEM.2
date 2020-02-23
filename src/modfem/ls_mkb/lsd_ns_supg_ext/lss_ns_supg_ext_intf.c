/*************************************************************
lsh_ns_supg_ext_intf.h - interface of the extensions of the mkb solver of
   linear equations for ns_supg problem module in ModFEM. Similar to mkb,
   the extension uses la... modules for storage of the stiffness matrix and basic
   linear algebra operations (see interface include/lah_intf.h ) and
   employs Block versions of standard iterative methods (Jacobi, Gauss-Seidel,
   additive Schwarz, multiplicative Schwarz) and ILU(0) for preconditioning.
   The extension is switched on by providing special input file
   with the first line: NS_SUPG_SOLVER_DATA ,
   otherwise standard mkb operations are performed.

   The file contains the provided interface with the headers of routines
   called from the FEM code

   The actual solution phase (lsr_ns_supg_ext_solve procedure) has a parameter
   that controls whether the solver is a distributed memory parallel solver
   or a standard sequential solver. In the first case the solver calls
   the finite element external part to perform three operations:
	   fem_vec_norm - to compute a global (inter-processor) vector norm
	   fem_sc_prod - to compute a global (inter-processor) scalar product
	   fem_exchange_dofs - to exchange the values of degrees of freedom
					  between subdomains (in overlapping Schwarz manner)
   For multigrid solution the solver calls the finite element part with
	   fem_proj_sol_lev - to project solution between levels (grids)
   The required interface of the solver containing headers of the
   procedures above is defined in the file "lsh_ns_supg_fem_intf.h"

Contents (declarations of the following routines):
  lsr_ns_supg_ext_init - to create a new solver instance, read its control
						 parameters and initialize its data structure
  lsr_ns_supg_ext_create_matrix - to allocate space for a global system matrix
  lsr_ns_supg_ext_create_precon - to create preconditioner blocks corresponding
								  to small subdomains of neighboring elements
  lsr_ns_supg_ext_clear_matrix - to initialize data structure of system matrix
  lsr_ns_supg_ext_assemble_local_sm - to assemble entries to the global stiffness
							matrix and the global load vector using the provided
							local stiffness matrix and load vector
  lsr_ns_supg_ext_show_matrix - to show matrix structure using graphic library
  lsr_ns_supg_ext_fill_precon - to prepare special preconditioner for ns_supg
  lsr_ns_supg_ext_solve - to solve a system of equations, given previously
						  constructed system matrix and preconditioner
  lsr_ns_supg_ext_free_matrix - to free space for solver data structure (system
								matrix, preconditioner and possibly other)

History:
		08.2015 - Krzysztof Banas, initial version

*************************************************************/

   // The vectors in ns_supg extension are divided into two parts:
   // velocities (3 DOFs per DOF entity) and pressure (1 DOF per DOF entity)
   // In each part the internal DOFs are grouped first and than overlap DOFs
   // follow to enable execution in PARALLEL (MPI) mode (this requires
   // the interception of callbacks to FEM code and replacing standard
   // operations with separate operations for velocity DOFs and pressure DOFs
   // For an example vector u its velocity and pressure parts are denoted
   // by u_v and u_p. The induced decomposition of system matrix A into parts
   // Avv, Avp, Apv and App follows naturally.

   // As a consequence, in ns_supg extensions there are 4, 5 or 6 subsystems:
   // 0. original (may be omitted)
   // 1. Avv - velocity-velocity submatrix
   // 2. Avp - velocity-pressure submatrix
   // 3. Apv - pressure-velocity submatrix
   // 4. App - pressure-pressure submatrix
   // 5. S - Schur complement matrix, S = App - Apv * Avv^{-1} * Avp
   // How to create Schur complement?
   // Given information on non-zero structure of components Apv i Avp
   // the structure of Schur component can be determined
   // then, using standard operations matrix-matrix product can be performed, i.e.
// for(i=0;i<N;i++){
//   for(j=0;j<N;j++){
//     if(S(i,j) non-zero){ // this check is matrix format dependent
//       S(i,j) = inner_product( Apv_row * {Avv^{-1}*Avp}_column )

// standard inner product of two sparse vectors in CRS is:
// (two vectors v and w, with v_nnz and w_nnz non-zeros, with nnz_pos(.,index)
// indicating position of non-zero entry in v or w)
// v_i=0; w_i=0; product=0.0;
// while(v_i < v_nnz && w_i < w_nnz){
//    if( nnz_pos(v, v_i) < nnz_pos(w, w_i) ){ v_i++; }
//    if( nnz_pos(v, v_i) > nnz_pos(w, w_i) ){ w_i++; }
//    if( nnz_pos(v, v_i) == nnz_pos(w, w_i) ) {
//      product += v[nnz_pos(v, v_i)] * w[nnz_pos(w, w_i)]; v_i++; w_i++;
//    }
// }

// If Apv and Avp are stored in CRS the algorithm may look as follows:
// 1. A_temp = Avv^{-1} * Avp
// 2. for(j=0;j<N;j++) {
//   3. rewrite j-th column of A_temp to CRS
//   4. calculate S(i,j) for all i for which non-zero entries exist in S
// }

   // Original is used for matrix vector product (first step of compreres,
   // i.e. compres) - or not used at all
   // Avv is used for preconditioning velocities (it can employ different
   // preconditioners from mkb)
   // Apv and Avp are used for creating Schur complement matrix:
   // S = App - Apv * Avv^{-1} * Avp and possess preconditioners obtained as
   // products with Avv^{-1}
   // App is used for creating S only
   // S posses diferent non-zero structure than its components! (but this structure
   // can be computed by e.g. algorithm for inner product of two CRS vectors, with
   // calculations of inner product exchanged with indication of non-zero entry in S


   // the lad_... modules may need some extensions to support ns_supg extenions:
   // lar_mm_product_with_left_block_diagonal
   // lar_mm_product_with_right_block_diagonal
   // lar_mm_Schur_product (uses special information on non-zero pattern
   //                       of Schur complement matrix)
   // (test for pattern - rewrite to full, perform multiplication, check pattern)
   //

   // for the rest of operations the existing interface of lad_... modules
   // can be used, e.g.:
   // lar_compute_residual for spmv products with original submatrices,
   // lar_compute_preconditioned_residual for products with matrices
   // stored in preconditioner arrays (such as e.g.  Mpv = Apv * Avv^{-1}
   //  or Mvp = Avv^{-1} * Avp)


   // The preconditioner algorithm in terms of lar_ procedures looks as follows:

   // Version 1 (with existing interface)
   //
   // input vector u^k = { u_v^k, u_p^k } (both include internal and overlap DOFs)
   //
   // 1. u_v_1 = lar_compute_preconditioned_residual(Avv, u_v^k as B(RHS), X=0)
   // 1a. exchange of u_v_1 DOFs
   // 2. u_p_1 = lar_compute_residual(Apv, u_v_1 as X, B=0)
   // 2a. exchange of u_p_1 DOFs
   // 3. u_p_2 = u_p^k - u_p_1
   // 4. solve: u_p^{k+1} = S^{-1} * u_p_2 (direct-global?, multigrid?)
   // 4a. u_p^{k+1} made unique for all subdomains (exchange of u_p DOFs?)
   // 5. u_v_2 = lar_compute_residual(Avp, u_p^{k+1} as X, B=0)
   // 5a. exchange of u_v_2 DOFs
   // 6. u_v_1 = u_v^k - u_v_2
   // 7. u_v^{k+1} = lar_compute_preconditioned_residual(Avv, u_v_1 as B(RHS), X=0)
   // 7a. exchange of u_v^{k+1} DOFs
   //
   // output vector: u^{k+1} = {  u_v^{k+1}, u_p^{k+1} }

   // Version 2 (with extended interface and preconditioner matrix:
   //            Mpv = Apv * Avv^{-1} )
   //
   // input vector u^k = { u_v^k, u_p^k } (both include internal and overlap DOFs)
   //
   // 1. u_p_1 = lar_compute_preconditioned_residual(Mpv, u_v^k as B(RHS), X=0)
   // 1a. exchange of u_p_1 DOFs
   // 2. u_p_2 = u_p^k - u_p_1
   // 3. solve: u_p^{k+1} = S^{-1} * u_p_2 (direct-global?, multigrid?)
   // 3a. u_p^{k+1} made unique for all subdomains (exchange of u_p DOFs?)
   // 4. u_v_2 = lar_compute_residual(Dvp, u_p^{k+1} as X, B=0)
   // 4a. exchange of u_v_2 DOFs
   // 5. u_v_1 = u_v^k - u_v_2
   // 6. u_v^{k+1} = lar_compute_preconditioned_residual(Avv, u_v_1 as B(RHS), X=0)
   // 6a. exchange of u_v^{k+1} DOFs
   //
   // output vector: u^{k+1} = {  u_v^{k+1}, u_p^{k+1} }

   // Version 3 (with extended interface and preconditioner matrices:
   //            Mpv = Apv * Avv^{-1} and Mvp = Avv^{-1} * Avp )
   //
   // input vector u^k = { u_v^k, u_p^k } (both include internal and overlap DOFs)
   //
   // 1. u_p_1 = lar_compute_preconditioned_residual(Mpv, u_v^k as B(RHS), X=0)
   // 1a. exchange of u_p_1 DOFs
   // 2. u_p_2 = u_p^k - u_p_1
   // 3. solve: u_p^{k+1} = S^{-1} * u_p_2 (direct-global?, multigrid?)
   // 3a. u_p^{k+1} made unique for all subdomains (exchange of u_p DOFs?)
   // 4. u_v_1 = lar_compute_preconditioned_residual(Mvp, u_p^{k+1} as B(RHS), X=0)
   // 4a. exchange of u_v_1 DOFs
   // 5. u_v_2 = lar_compute_preconditioned_residual(Avv, u_v^k as B(RHS), X=0)
   // 5a. exchange of u_v_2 DOFs
   // 6. u_v^{k+1} = u_v_2 - u_v_1
   //
   // output vector: u^{k+1} = {  u_v^{k+1}, u_p^{k+1} }



#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* provided interface of the solver - headers of routines defined in this file */
#include <modfem/ls_mkb/lsh_mkb_intf.h>

/* internal information for the solver module */
#include "../lsd_mkb_core/lsh_mkb_core.h"

// interface of linear algebra supporting package
#include <modfem/ls_mkb/lah_intf.h>

#ifdef __cplusplus
extern "C" {
#endif

/**-----------------------------------------------------------
   lsr_ns_supg_ext_init - to create a new solver instance, read its control
   parameters and initialize its data structure
------------------------------------------------------------*/
int lsr_ns_supg_ext_init( /* returns: >0 - solver ID, <0 - error code */
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	int Parallel,      /* parameter specifying sequential (LSC_SEQUENTIAL) */
	/* or parallel (LSC_PARALLEL) execution */
	int* Max_num_levels_p,  /* in: number of levels for multigrid: */
	/*     1 - enforce single level solver */
	/*     >1 - enforce the number of levels */
	/* out: actual number of levels !!! */
	char* Filename,  /* in: name of the file with control parameters */
	int Max_iter, /* maximal number of iterations, -1 for values from Filename */
	int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
	double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */
	int Monitoring_level /* Level of output, -1 for Filename */
			  )
{


  // Read or assign control data for all subsystems

  // initiate local variable
  int nr_level = *Max_num_levels_p;

  /*kbw*/
  printf("\nIn lsr_ns_supg_ext_init: Solver id %d, Parallel %d, Nr_levels %d, Filename %s\n",
	 lsv_mkb_cur_solver_id, Parallel, nr_level, Filename);
  printf("Monitor %d, Nr_iter %d, Conv_type %d, Conv_meas %15.12lf\n",
	 Monitoring_level, Max_iter, Error_type, Error_tolerance);
  /*kew*/

  // To begin with, we assign default data to all subsystems (they can be read form file)

  // Version 1a. 4 subsystems for SM + additional subsystem? for Schur complement
  lsv_mkb_solver[lsv_mkb_cur_solver_id].solver_id = lsv_mkb_cur_solver_id;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].parallel = Parallel;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[0].nr_level = 1;

  lsv_mkb_solver[lsv_mkb_cur_solver_id].nr_subsystems = 6;

  // we initiate four subsystems for four parts of SM
  int isubs;
  for (isubs = 1; isubs <= 4; isubs++){

	lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
	lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

	lst_mkb_levels *mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);

	mkb_level->Solver = GMRES;
	if (Max_iter < 0) mkb_level->Max_iter = 100;
	else mkb_level->Max_iter = Max_iter;
	if (Error_type < 0) mkb_level->Conv_type = REL_RES_INI;
	else  mkb_level->Conv_type = Error_type;
	if (Error_tolerance < 0.0) mkb_level->Conv_meas = 1e-6;
	else mkb_level->Conv_meas = Error_tolerance;
	if (Monitoring_level < 0) mkb_level->Monitor = LSC_SILENT;
	else mkb_level->Monitor = Monitoring_level;

	/* default solver is single level GMRES with 50 Krylov vectors */
	mkb_level->GMRES_type = STANDARD_GMRES;
	mkb_level->Krylbas = 50;
	mkb_level->Pdeg = -1;
	mkb_level->Precon = BLOCK_GS; /* preconditioning by a single iteration */
	mkb_level->Nr_prec = 1;       /* of block Gauss-Seidel */
	mkb_level->Block_type = 1; /* single DOF structure subdomains */
	mkb_level->Nr_pre_smooth = 1;
	mkb_level->Nr_post_smooth = 0;

  } // the end of loop over subsystems


  /*ok_kbw*/
  //if(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[0].level[nr_level-1].Monitor>LSC_ERRORS){

  printf("\nInitiated data for solver (id = %d): %d mesh levels with parameters:\n",
	 lsv_mkb_cur_solver_id,
	 lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[0].nr_level);

  for (isubs = 1; isubs <= 4; isubs++){

	lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
	lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;


	lst_mkb_levels *mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);


	printf("\nSubsystem: %d\n", isubs);
	printf("Solver type %d, GMRES type  %d, Krylbas %d, Pdeg %d\n",
	   mkb_level->Solver, mkb_level->GMRES_type, mkb_level->Krylbas,
	   mkb_level->Pdeg);
	printf("Preconditioner %d, Number of sweeps %d, Block type/size %d\n",
	   mkb_level->Precon, mkb_level->Nr_prec, mkb_level->Block_type);
	printf("Max_iter %d, Conv_type %d, Conv_meas %20.15lf\n",
	   mkb_level->Max_iter, mkb_level->Conv_type, mkb_level->Conv_meas);
	printf("Monitor %d, Nr_pre_smmoth %d, Nr_post_smooth %d \n",
	   mkb_level->Monitor,
	   mkb_level->Nr_pre_smooth, mkb_level->Nr_post_smooth);
  }
  //}
  /*kew*/


  *Max_num_levels_p = 1; // single level for the time being


  return 0;

}

/**--------------------------------------------------------
   lsr_ns_supg_ext_create_matrix - to allocate space for a global system matrix
---------------------------------------------------------*/
int lsr_ns_supg_ext_create_matrix(
		/* returns: >=0 - success code, <0 - error code */
		int Solver_id,   /* in: solver ID (used to identify the subproblem) */
		int Level_id,    /* in: level ID */
		int Nrblocks,    /* in: number of DOF blocks */
		int Nrdof_glob,  /* in: total number of DOFs */
		int Max_sm_size, /* in: maximal size of the stiffness matrix */
		int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
		int* Posglob,	   /* in: list of global numbers of first dof */
		int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
		int** L_offbl	   /* in: list of lists of off diagonal blocks */
		)
{

  int info, idof, ineig, isubs, ibl;
  lst_mkb_levels *mkb_level;

  int block_size_input = Nrdof_glob / Nrblocks;

  /*kbw*/
  printf("Entering lsr_ns_supg_ext_create_matrix\n");
  /*kew*/

  // checking
#ifdef DEBUG_SIM
  if(Level_id != 0){
	printf("lsr_ns_supg_ext_create_matrix called with Level_id %d. Exiting.\n",
	   Level_id);
	exit(-1);
  }
  if(Nrdof_glob%Nrblocks != 0 || block_size_input != 4){
	printf("Attempt to initialize ns_supg_ext solver data structure \n");
	printf("with non-constant or wrong block_size %d. Exiting!\n", block_size_input);
	printf("Nrdof_glob %d, Nrblocks %d\n", Nrdof_glob, Nrblocks);
	exit(-1);
  }
  // we check whether input matrices satisfy assumptions used in simplifications
  for (ibl = 0; ibl< Nrblocks; ibl++){
	if(Posglob[ibl] != 4*ibl){
	  printf("Blocks not in standard order in lsr_ns_supg_ext_create_matrix! Exiting.\n");
	  printf("Block %d, Nrdofbl %d, Posglob %d\n", ibl, Nrdofbl[ibl], Posglob[ibl]);
	  exit(-1);
	}
	if(Nrdofbl[ibl]!=4){
	  printf("Blocks not with 4 DOFs in lsr_ns_supg_ext_create_matrix! Exiting.\n");
	  printf("Block %d, Nrdofbl %d, Posglob %d\n", ibl, Nrdofbl[ibl], Posglob[ibl]);
	  exit(-1);
	}
  }
#endif

  // we have to alocate local versions of arrays since after splitting we have
  // twice as much blocks as before - velocity blocks and pressure blocks are separated
  int nrblocks_ns_supg_ext = 2 * Nrblocks;

  // we allocate tables and split velocities and pressure
  int* nrdofbl_loc = (int *)malloc(nrblocks_ns_supg_ext*sizeof(int));
  int* posglob_loc = (int *)malloc(nrblocks_ns_supg_ext*sizeof(int));
  int* nroffbl_loc = (int *)malloc(nrblocks_ns_supg_ext*sizeof(int));
  int** l_offbl_loc = (int **)malloc(nrblocks_ns_supg_ext*sizeof(int*));
  for (ibl = 0; ibl < Nrblocks; ibl++){

	// for each subsystem we need one set of neighboring blocks, either velocities
	// or pressure, never both, hence Nroffbl[ibl] is not changed
	// but we add one more block for splitted diagonal blocks (part of it becomes off-diagonal)
	l_offbl_loc[ibl] = (int *)malloc((Nroffbl[ibl] + 1)*sizeof(int));
	l_offbl_loc[Nrblocks + ibl] = (int *)malloc((Nroffbl[ibl] + 1)*sizeof(int));

  }

  /* set the current solver ID */
  lsv_mkb_cur_solver_id = Solver_id;


  // we use the same trick as for distributed execution - for each subsystem we assume
  // the existence of DOF blocks without associated rows of SM
  // hence for each subsystem we can use the whole RHS vector (if needed), but
  // we can manipulate the locations of ghost and not ghost blocks




  // subsystem 1 - Avv
  isubs = 1;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);


  // for subsystem 1 we can use standard approach - it operates on velocities only
  // numbering of blocks remains unchanged, only nrdofbl and posglob change


  // original block_size = Nrdof_glob / Nrblocks - for ns_supg_ext should be == 4
  // we change Nrdof_glob
  mkb_level->Nrdofgl = Nrblocks * 3; // 3 is block_size for velocities
  mkb_level->Nr_dof_blocks = Nrblocks;

  // we fill local arrays, but only for velocities-velocities coupling (hence
  // we use only  Nrblocks blocks
  for (ibl = 0; ibl < Nrblocks; ibl++){

	nrdofbl_loc[ibl] = 3;
	posglob_loc[ibl] = 3 * ibl;
	nroffbl_loc[ibl] = Nroffbl[ibl];
	for (ineig = 0; ineig < Nroffbl[ibl]; ineig++){
	  l_offbl_loc[ibl][ineig] = L_offbl[ibl][ineig];
	}


  }

/*kbw

  printf("Allocating SM for subsystem %d: \n", isubs);
  printf("Max_sm_size %d, Nrblocks %d, Nrdofgl %d\n",
  Max_sm_size, Nrblocks, mkb_level->Nrdofgl);
  printf("Block_type %d, Precon %d\n",
  mkb_level->Block_type, mkb_level->Precon);

  for(ibl = 0; ibl < Nrblocks; ibl++){
  printf("Block %d, Nrdofbl %d, Posglob %d, Nroffbl %d, L_offbl:\n",
  ibl, nrdofbl_loc[ibl], posglob_loc[ibl], nroffbl_loc[ibl]);
  for(ineig = 0; ineig < nroffbl_loc[ibl]; ineig++){
  printf("%5d", l_offbl_loc[ibl][ineig]);
  }
  printf("\n");
  }
  printf("Entering lar_allocate_SM_and_LV\n");

//*kew*/

  mkb_level->SM_and_LV_id = lar_allocate_SM_and_LV(Max_sm_size, Nrblocks,
				   mkb_level->Nrdofgl, mkb_level->Block_size,
				   nrdofbl_loc, posglob_loc, nroffbl_loc, l_offbl_loc,
				   mkb_level->Block_type, mkb_level->Precon);


  if (mkb_level->SM_and_LV_id >= 0) {
	info = 0;
	/*kbw
	  printf("OK!\n");
	  //*kew*/

  }
  else {

	printf("Error allocating SM for subsystem %d. Exiting. \n", isubs);
	exit(-1);

  }


  // subsystem 4 - App
  isubs = 4;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);


  // two approaches possible:
  // 1. use the whole vector (with pressures at the end)
  // not used yet

  // 2. use only pressure part (with suitable management of input/outpur vectors)

  // original block_size = Nrdof_glob / Nrblocks - for ns_supg_ext should be == 4
  // we change Nrdof_glob
  mkb_level->Nrdofgl = Nrblocks; // block_size for pressure equal 1
  mkb_level->Nr_dof_blocks = Nrblocks;

  // we fill local arrays, but only for pressures-pressures coupling (hence
  // we use only  Nrblocks blocks - we start with 0 not with Nrblocks!!!)
  for (ibl = 0; ibl < Nrblocks; ibl++){

	nrdofbl_loc[ibl] = 1;
	posglob_loc[ibl] = ibl;
	// the same as for velocities
	nroffbl_loc[ibl] = Nroffbl[ibl];
	for (ineig = 0; ineig < Nroffbl[ibl]; ineig++){
	  l_offbl_loc[ibl][ineig] = L_offbl[ibl][ineig];
	}


  }

/*kbw

	printf("Allocating SM for subsystem %d: \n", isubs);
	printf("Max_sm_size %d, Nrblocks %d, Nrdofgl %d\n",
	Max_sm_size, Nrblocks, mkb_level->Nrdofgl);
	printf("Block_type %d, Precon %d\n",
	mkb_level->Block_type, mkb_level->Precon);

	for(ibl = 0; ibl < Nrblocks; ibl++){
	printf("Block %d, Nrdofbl %d, Posglob %d, Nroffbl %d, L_offbl:\n",
	ibl, nrdofbl_loc[ibl], posglob_loc[ibl], nroffbl_loc[ibl]);
	for(ineig = 0; ineig < nroffbl_loc[ibl]; ineig++){
	printf("%5d", l_offbl_loc[ibl][ineig]);
	}
	printf("\n");
	}
	printf("Entering lar_allocate_SM_and_LV\n");

//*kew*/


  mkb_level->SM_and_LV_id = lar_allocate_SM_and_LV(Max_sm_size, Nrblocks,
				   mkb_level->Nrdofgl, mkb_level->Block_size,
				   nrdofbl_loc, posglob_loc, nroffbl_loc, l_offbl_loc,
				   mkb_level->Block_type, mkb_level->Precon);


  if (mkb_level->SM_and_LV_id >= 0) {
	info = 0;
	/*kbw
	  printf("OK!\n");
	  //*kew*/

  }
  else {

	printf("Error allocating SM for subsystem %d. Exiting. \n", isubs);
	exit(-1);

  }


  // subsystem 3 - Apv
  isubs = 3;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);


  // similar to subsystem 2 but the role of velocities and pressure is interchanged
  // original block_size = Nrdof_glob / Nrblocks - for ns_supg_ext should be == 4
  // we do not change Nrdof_glob - we use the whole RHS vector
  mkb_level->Nrdofgl = Nrblocks * 4;
  mkb_level->Nr_dof_blocks = 2 * Nrblocks; // twice as much blocks

  // we fill local arrays, now for the whole vector, with splitted velocities and pressures
  for (ibl = 0; ibl < Nrblocks; ibl++){

	// velocity blocks
	nrdofbl_loc[ibl] = 3;
	posglob_loc[ibl] = 3 * ibl;
	// for subsystem 3 we multiply pressure blocks by velocity blocks - the latter have no rows
	nroffbl_loc[ibl] = -1;

	// pressure blocks
	nrdofbl_loc[ibl + Nrblocks] = 1;
	posglob_loc[ibl + Nrblocks] = 3 * Nrblocks + ibl;
	// for subsystem 3 we multiply pressure blocks by velocity blocks
	nroffbl_loc[ibl + Nrblocks] = Nroffbl[ibl];
	for (ineig = 0; ineig < Nroffbl[ibl]; ineig++){
	  l_offbl_loc[ibl + Nrblocks][ineig] = L_offbl[ibl][ineig];
	}

	// each pressure block has one more neighbour - its corresponding velocity block
	//  (the velocity block from the Dia array - now not on the diagonal!!!)
	nroffbl_loc[ibl + Nrblocks]++;
	l_offbl_loc[ibl + Nrblocks][Nroffbl[ibl]] = ibl;


  }

/*kbw

	printf("Allocating SM for subsystem %d: \n", isubs);
	printf("Max_sm_size %d, Nrblocks %d, Nrdofgl %d\n",
	Max_sm_size, Nrblocks, mkb_level->Nrdofgl);
	printf("Block_type %d, Precon %d\n",
	mkb_level->Block_type, mkb_level->Precon);

	for(ibl = 0; ibl < 2*Nrblocks; ibl++){
	printf("Block %d, Nrdofbl %d, Posglob %d, Nroffbl %d, L_offbl:\n",
	ibl, nrdofbl_loc[ibl], posglob_loc[ibl], nroffbl_loc[ibl]);
	for(ineig = 0; ineig < nroffbl_loc[ibl]; ineig++){
	printf("%5d", l_offbl_loc[ibl][ineig]);
	}
	printf("\n");
	}
	printf("Entering lar_allocate_SM_and_LV\n");

//*kew*/


  mkb_level->SM_and_LV_id = lar_allocate_SM_and_LV(Max_sm_size, 2 * Nrblocks,
				   mkb_level->Nrdofgl, mkb_level->Block_size,
				   nrdofbl_loc, posglob_loc, nroffbl_loc, l_offbl_loc,
				   mkb_level->Block_type, mkb_level->Precon);


  if (mkb_level->SM_and_LV_id >= 0) {
	info = 0;
	/*kbw
	  printf("OK!\n");
	  //*kew*/

  }
  else {

	printf("Error allocating SM for subsystem %d. Exiting. \n", isubs);
	exit(-1);

  }


  // subsystem 2 - Avp
  isubs = 2;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);


  // for subsystem 2 we use special approach - blocks associated with stiffness matrix rows
  // correspond to velocities (3 dofs) but rows are multiplied by vectors associated
  // with pressure (1 dof). Hence we need to use the whole vector of unknowns
  // with velocity dofs at the beginning and pressure dofs at the end
  // (we do not store LV in lad_ module, we manage RHS vector in lsr_ns_supg_ext_intf.c)

  // original block_size = Nrdof_glob / Nrblocks - for ns_supg_ext should be == 4
  // we do not change Nrdof_glob - we use the whole RHS vector
  mkb_level->Nrdofgl = Nrblocks * 4;
  mkb_level->Nr_dof_blocks = 2 * Nrblocks; // twice as much blocks

  // we fill local arrays, now for the whole vector, with splitted velocities and pressures
  for (ibl = 0; ibl < Nrblocks; ibl++){

	// velocity blocks
	nrdofbl_loc[ibl] = 3;
	posglob_loc[ibl] = 3 * ibl;
	// for subsystem 2 we multiply velocity blocks by pressure blocks
	nroffbl_loc[ibl] = Nroffbl[ibl];
	for (ineig = 0; ineig < Nroffbl[ibl]; ineig++){
	  l_offbl_loc[ibl][ineig] = Nrblocks + L_offbl[ibl][ineig];
	}

	// each velocity block has one more neighbour - its corresponding pressure block
	//  (the pressure block from the diagonal - now not on the diagonal!!!)
	nroffbl_loc[ibl]++;
	l_offbl_loc[ibl][Nroffbl[ibl]] = Nrblocks + ibl;

	// pressure blocks
	nrdofbl_loc[ibl + Nrblocks] = 1;
	posglob_loc[ibl + Nrblocks] = 3 * Nrblocks + ibl;
	// for subsystem 2 we multiply velocity blocks by pressure blocks - the latter have no rows
	nroffbl_loc[ibl + Nrblocks] = -1;


  }

  /*kbw

	printf("Allocating SM for subsystem %d: \n", isubs);
	printf("Max_sm_size %d, Nrblocks %d, Nrdofgl %d\n",
	Max_sm_size, Nrblocks, mkb_level->Nrdofgl);
	printf("Block_type %d, Precon %d\n",
	mkb_level->Block_type, mkb_level->Precon);

	for(ibl = 0; ibl < 2*Nrblocks; ibl++){
	printf("Block %d, Nrdofbl %d, Posglob %d, Nroffbl %d, L_offbl:\n",
	ibl, nrdofbl_loc[ibl], posglob_loc[ibl], nroffbl_loc[ibl]);
	for(ineig = 0; ineig < nroffbl_loc[ibl]; ineig++){
	printf("%5d", l_offbl_loc[ibl][ineig]);
	}
	printf("\n");
	}
	printf("Entering lar_allocate_SM_and_LV\n");

	//*kew*/


  mkb_level->SM_and_LV_id = lar_allocate_SM_and_LV(Max_sm_size, 2 * Nrblocks,
				   mkb_level->Nrdofgl, mkb_level->Block_size,
				   nrdofbl_loc, posglob_loc, nroffbl_loc, l_offbl_loc,
				   mkb_level->Block_type, mkb_level->Precon);


  if (mkb_level->SM_and_LV_id >= 0) {

	/*kbw
	  printf("OK!\n");
	  //*kew*/

  }
  else {

	printf("Error allocating SM for subsystem %d. Exiting. \n", isubs);
	exit(-1);

  }

  free(nrdofbl_loc);
  free(posglob_loc);
  free(nroffbl_loc);
  for (ibl = 0; ibl < 2 * Nrblocks; ibl++) free(l_offbl_loc[ibl]);
  free(l_offbl_loc);


  return(info);

}


/**--------------------------------------------------------
   lsr_ns_supg_ext_create_precon - to create preconditioner blocks corresponding
   to small subdomains of neighboring elements
---------------------------------------------------------*/
int lsr_ns_supg_ext_create_precon( /* returns:   >0 number of diagonal blocks */
					/*	       <=0 - error */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id           /* in: level ID */
	)
{

  // Create block data structures for utilized preconditioners
  return 0;
}

/**--------------------------------------------------------
   lsr_ns_supg_ext_clear_matrix - to initialize block structure of system matrix
---------------------------------------------------------*/
int lsr_ns_supg_ext_clear_matrix(
	/* returns: >=0 - success code, <0 - error code */
	 int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	 int Level_id,    /* in: level ID */
	 int Comp_type    /* in: indicator for the scope of computations: */
	 /*   NS_SUPG_EXT_SOLVE - solve the system */
	 /*   NS_SUPG_EXT_RESOLVE - resolve for the new rhs vector */
				 )
{

  /*kbw*/
  printf("Entering lsr_ns_supg_ext_clear_matrix\n");
  /*kew*/

#ifdef DEBUG_SIM
  if(Level_id != 0){
	printf("lsr_ns_supg_ext_create_matrix called with Level_id %d. Exiting.\n",
	   Level_id);
	exit(-1);
  }
#endif

  int isubs;
  for (isubs = 1; isubs <= 4; isubs++){

	lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
	lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = Level_id;

	// Clear matrices of subsystems
	int scope_type;
	if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
	else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
	else scope_type=LAC_SCOPE_SM_AND_LV; // default
	int SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[Level_id].SM_and_LV_id;
	int info = lar_initialize_SM_and_LV(SM_and_LV_id, scope_type);

  }

  return(0);
}

/**-----------------------------------------------------------
   lsr_ns_supg_ext_assemble_local_sm - to assemble entries to the global stiffness matrix
   and the global load vector using the provided local
   stiffness matrix and load vector
------------------------------------------------------------*/
int lsr_ns_supg_ext_assemble_local_sm(
		/* returns: >=0 - success code, <0 - error code */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id,          /* in: level ID */
	int Comp_type,         /* in: indicator for the scope of computations: */
	/*   NS_SUPG_EXT_SOLVE - solve the system */
	/*   NS_SUPG_EXT_RESOLVE - resolve for the new rhs vector */
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

  // checking
#ifdef DEBUG_SIM
  if(Level_id != 0){
	printf("lsr_ns_supg_ext_create_matrix called with Level_id %d. Exiting.\n",
	   Level_id);
	exit(-1);
  }
#endif

  /*kbw
	{
	printf("In lsr_ns_supg_ext_assemble_local_sm: Solver_id %d, Level_id %d, Comp_type %d, Nr_dof_bl %d\n",
	Solver_id, Level_id, Comp_type, Nr_dof_bl);
	int i,j,ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
	pli = 0; nrdof=0;
	for(ibl=0;ibl<Nr_dof_bl; ibl++) nrdof+=L_bl_nrdof[ibl];
	for(ibl=0;ibl<Nr_dof_bl; ibl++){
	printf("bl_id %d, bl_nrdof %d\n",
	L_bl_id[ibl],L_bl_nrdof[ibl]);
	nri = L_bl_nrdof[ibl];
	plj=0;
	for(jbl=0;jbl<Nr_dof_bl;jbl++){
	//printf("Stiff_mat transposed! (blocks %d:%d -> group of rows %d, group of columns %d)\n",
	//	 jbl,ibl,jbl,ibl);
	//nrj = L_bl_nrdof[jbl];
	//for(i=0;i<nri;i++){
	//  jaux = plj+(pli+i)*nrdof;
	//  for(j=0;j<nrj;j++){
	//    printf("%20.15lf",stiff_mat[jaux+j]);
	//  }
	//  printf("\n");
	//}
	printf("Stiff_mat  (blocks %d:%d -> group of rows %d, group of columns %d)\n",
	jbl,ibl,jbl,ibl);
	nrj = L_bl_nrdof[jbl];
	for(j=0;j<nrj;j++){
	for(i=0;i<nri;i++){
	jaux = plj+(pli+i)*nrdof;
	printf("%20.15lf",Stiff_mat[jaux+j]);
	}
	printf("\n");
	}
	plj += nrj;
	}
	printf("Rhs_vect (block %d -> group of rows %d):\n", ibl, ibl);
	for(i=0;i<nri;i++){
	printf("%20.15lf\n",Rhs_vect[pli+i]);
	}
	printf("\n");
	pli += nri;
	}
	getchar();
	//}
	}
	/*kew*/

  /* set the current solver ID */
  lsv_mkb_cur_solver_id = Solver_id;

  // we allocate tables to split velocities and pressure
  int* l_bl_id_loc = (int *)malloc(2 * Nr_dof_bl*sizeof(int));
  int* l_bl_nrdof_loc = (int *)malloc(2 * Nr_dof_bl*sizeof(int));
  int nr_dofs_loc = 2 * Nr_dof_bl * 4; // we assume large blocksize = 4
  double *stiff_mat_loc = (double *)malloc(nr_dofs_loc*nr_dofs_loc*sizeof(double));
  double *rhs_vect_loc = (double *)malloc(nr_dofs_loc*sizeof(double));

  // Assemble one element stiffness matrix into four subsystems
  int isubs;
  int i, j, ibl, jbl, pli, plj, nri, nrj, nrdof, jaux;
  int pli_loc, plj_loc, nri_loc, nrj_loc, nrdof_loc, jaux_loc;
  lst_mkb_levels *mkb_level;

  // subsystem 1 - Avv
  isubs = 1;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);

  // initialization
  memset(stiff_mat_loc, 0, nr_dofs_loc*nr_dofs_loc*sizeof(double));
  memset(rhs_vect_loc, 0, nr_dofs_loc*sizeof(double));

  // block numbering is not changed, but 3x3 blocks are cut out of 4x4 blocks
  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	l_bl_id_loc[ibl] = L_bl_id[ibl];
	l_bl_nrdof_loc[ibl] = 3;
  }

  nrdof = 0;
  for (ibl = 0; ibl < Nr_dof_bl; ibl++) nrdof += L_bl_nrdof[ibl];
  nrdof_loc = Nr_dof_bl * 3;
  pli = 0; // velocity blocks start at the beginning
  pli_loc = 0;
  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	nri = L_bl_nrdof[ibl];
	nri_loc = 3;
	plj = 0;
	plj_loc = 0;
	for (jbl = 0; jbl < Nr_dof_bl; jbl++){
	  nrj = L_bl_nrdof[jbl];
	  nrj_loc = 3;
	  for (j = 0; j < nrj_loc; j++){
	for (i = 0; i < nri_loc; i++){
	  jaux = plj + (pli + i)*nrdof;
	  jaux_loc = plj_loc + (pli_loc + i)*nrdof_loc;
	  stiff_mat_loc[jaux_loc + j] = Stiff_mat[jaux + j];
	}
	  }
	  plj += nrj;
	  plj_loc += nrj_loc;
	}
	for (i = 0; i < nri_loc; i++){
	  rhs_vect_loc[pli_loc + i] = Rhs_vect[pli + i];
	}
	pli += nri;
	pli_loc += nri_loc;
  }

  //printf("Entering lar_assemble_SM_and_LV for Avv\n");

  int scope_type;
  if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
  else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
  else scope_type=LAC_SCOPE_SM_AND_LV; // default
  int SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[Level_id].SM_and_LV_id;
  int info = lar_assemble_SM_and_LV(SM_and_LV_id, scope_type,
					Nr_dof_bl, l_bl_id_loc, l_bl_nrdof_loc,
					stiff_mat_loc, rhs_vect_loc, Rewr_dofs);


  // subsystem 4 - App
  isubs = 4;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);

  // initialization
  memset(stiff_mat_loc, 0, nr_dofs_loc*nr_dofs_loc*sizeof(double));
  memset(rhs_vect_loc, 0, nr_dofs_loc*sizeof(double));

  // block numbering is not changed, but 1x1 blocks are cut out of 4x4 blocks
  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	l_bl_id_loc[ibl] = L_bl_id[ibl];
	l_bl_nrdof_loc[ibl] = 1;
  }

  nrdof = 0;
  for (ibl = 0; ibl < Nr_dof_bl; ibl++) nrdof += L_bl_nrdof[ibl];
  nrdof_loc = Nr_dof_bl;
  pli = 3; // pressure blocks start at 3
  pli_loc = 0;
  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	nri = L_bl_nrdof[ibl];
	nri_loc = 1;
	plj = 3;
	plj_loc = 0;
	for (jbl = 0; jbl < Nr_dof_bl; jbl++){
	  nrj = L_bl_nrdof[jbl];
	  nrj_loc = 1;
	  for (j = 0; j < nrj_loc; j++){
	for (i = 0; i < nri_loc; i++){
	  jaux = plj + (pli + i)*nrdof;
	  jaux_loc = plj_loc + (pli_loc + i)*nrdof_loc;
	  stiff_mat_loc[jaux_loc + j] = Stiff_mat[jaux + j];
	}
	  }
	  plj += nrj;
	  plj_loc += nrj_loc;
	}
	for (i = 0; i < nri_loc; i++){
	  rhs_vect_loc[pli_loc + i] = Rhs_vect[pli + i];
	}
	pli += nri;
	pli_loc += nri_loc;
  }

  //printf("Entering lar_assemble_SM_and_LV for App\n");
  int scope_type;
  if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
  else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
  else scope_type=LAC_SCOPE_SM_AND_LV; // default

  SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[Level_id].SM_and_LV_id;
  info = lar_assemble_SM_and_LV(SM_and_LV_id, scope_type,
				Nr_dof_bl, l_bl_id_loc, l_bl_nrdof_loc,
				stiff_mat_loc, rhs_vect_loc, Rewr_dofs);


  // subsystem 3 - Apv
  isubs = 3;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);

  // initialization
  memset(stiff_mat_loc, 0, nr_dofs_loc*nr_dofs_loc*sizeof(double));
  memset(rhs_vect_loc, 0, nr_dofs_loc*sizeof(double));

  // twice as much blocks as for Avv and App
  int nr_dof_bl_loc = 2 * Nr_dof_bl;
  int nr_dof_bl_glob = mkb_level->Nr_dof_blocks;

  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	// velocity blocks
	l_bl_id_loc[ibl] = L_bl_id[ibl];
	l_bl_nrdof_loc[ibl] = 3;
	// pressure blocks
	l_bl_id_loc[ibl + Nr_dof_bl] = nr_dof_bl_glob / 2 + L_bl_id[ibl];
	l_bl_nrdof_loc[ibl + Nr_dof_bl] = 1;
  }

  /*kbw
	for(i = 0; i < 2*Nr_dof_bl; i++){
	printf("Block %d (%d), Nrdofbl %d:\n",
	i, l_bl_id_loc[i], l_bl_nrdof_loc[i]);
	}
	//*kew*/

  // in lar_assemble_SM_and_LV it is assumed that we supply full square matrices!!!!!!!!!
  int nrdof_i = 4 * Nr_dof_bl; // there are nrdof_i columns (horizontal i direction)
  int nrdof_j = 4 * Nr_dof_bl; // there are nrdof_j rows (vertical j direction)
  // int nrdof_loc = nrdof_i*nrdof_j; // total length of stiff_mat_loc
  nrdof = 4 * Nr_dof_bl; // standard size for original stiff_mat (nrdofxnrdof)
  pli = 0; // pressure-velocity blocks start at 0 for ibl, i (horizontal) direction
  pli_loc = 0;
  // we assemble only non-zero values in pressure rows (pressure-velocity blocks),
  // but take into account (include in stiff_mat_loc) zeros in pressure-pressure blocks
  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	nri = L_bl_nrdof[ibl];
	nri_loc = 3;
	plj = 3; // pressure-velocity blocks start at 3 for jbl, j (vertical) direction
	// we must provide the full SM with zeros for velocity-velocity blocks,
	// hence pressure-velocity blocks start at 3*Nr_dof_bl
	plj_loc = 3 * Nr_dof_bl;
	for (jbl = 0; jbl < Nr_dof_bl; jbl++){
	  nrj = L_bl_nrdof[jbl];
	  nrj_loc = 1;
	  for (j = 0; j < nrj_loc; j++){
	for (i = 0; i < nri_loc; i++){
	  jaux = plj + (pli + i)*nrdof;
	  jaux_loc = plj_loc + (pli_loc + i)*nrdof_j;
	  stiff_mat_loc[jaux_loc + j] = Stiff_mat[jaux + j];
	  /*kbw
		printf("rewriting SM entry %d to sm_loc %d: %20.15lf -> %20.15lf\n",
		jaux+j, jaux_loc+j, Stiff_mat[jaux+j], stiff_mat_loc[jaux_loc+j]);

		/*kew*/

	}
	  }
	  plj += nrj;
	  plj_loc += nrj_loc;
	}
	for (i = 0; i < nri_loc; i++){
	  rhs_vect_loc[pli_loc + i] = Rhs_vect[pli + i];
	}
	pli += nri;
	pli_loc += nri_loc;
  }

  //printf("Entering lar_assemble_SM_and_LV for Apv\n");
  int scope_type;
  if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
  else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
  else scope_type=LAC_SCOPE_SM_AND_LV; // default

  SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[Level_id].SM_and_LV_id;
  info = lar_assemble_SM_and_LV(SM_and_LV_id, scope_type,
				nr_dof_bl_loc, l_bl_id_loc, l_bl_nrdof_loc,
				stiff_mat_loc, rhs_vect_loc, Rewr_dofs);


  // subsystem 2 - Avp
  isubs = 2;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].cur_subsystem = isubs;
  lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].cur_level = 0;

  mkb_level = &(lsv_mkb_solver[lsv_mkb_cur_solver_id].subsystem[isubs].level[0]);

  // initialization
  memset(stiff_mat_loc, 0, nr_dofs_loc*nr_dofs_loc*sizeof(double));
  memset(rhs_vect_loc, 0, nr_dofs_loc*sizeof(double));

  // twice as much blocks as for Avv and App
  nr_dof_bl_loc = 2 * Nr_dof_bl;
  nr_dof_bl_glob = mkb_level->Nr_dof_blocks;

  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	// velocity blocks
	l_bl_id_loc[ibl] = L_bl_id[ibl];
	l_bl_nrdof_loc[ibl] = 3;
	// pressure blocks
	l_bl_id_loc[ibl + Nr_dof_bl] = nr_dof_bl_glob / 2 + L_bl_id[ibl];
	l_bl_nrdof_loc[ibl + Nr_dof_bl] = 1;
  }

  /*kbw
	for(i = 0; i < 2*Nr_dof_bl; i++){
	printf("Block %d (%d), Nrdofbl %d:\n",
	i, l_bl_id_loc[i], l_bl_nrdof_loc[i]);
	}
	//*kew*/

  // in lar_assemble_SM_and_LV it is assumed that we supply full square matrices!!!!!!!!!
  nrdof_i = 4 * Nr_dof_bl; // there are nrdof_i columns (horizontal i direction)
  nrdof_j = 4 * Nr_dof_bl; // there are nrdof_j rows (vertical j direction)
  // int nrdof_loc = nrdof_i*nrdof_j; // total length of stiff_mat_loc
  nrdof = 4 * Nr_dof_bl; // standard size for original stiff_mat (nrdofxnrdof)
  pli = 3; // velocity-pressure blocks start at 3 for ibl, i (horizontal) direction
  // we assemble only non-zero values in velocity rows (velocity-pressure blocks)
  // but we must provide the full SM with zeros for velocity-velocity blocks,
  // hence velocity-pressure blocks start at 3*Nr_dof_bl
  pli_loc = 3 * Nr_dof_bl;
  for (ibl = 0; ibl < Nr_dof_bl; ibl++){
	nri = L_bl_nrdof[ibl];
	nri_loc = 1;
	plj = 0; // velocity-pressure blocks start at 0 for jbl, j (vertical) direction
	plj_loc = 0;
	for (jbl = 0; jbl < Nr_dof_bl; jbl++){
	  nrj = L_bl_nrdof[jbl];
	  nrj_loc = 3;
	  for (j = 0; j < nrj_loc; j++){
	for (i = 0; i < nri_loc; i++){
	  jaux = plj + (pli + i)*nrdof;
	  jaux_loc = plj_loc + (pli_loc + i)*nrdof_j;
	  stiff_mat_loc[jaux_loc + j] = Stiff_mat[jaux + j];
	  /*kbw
		printf("rewriting SM entry %d to sm_loc %d: %20.15lf -> %20.15lf\n",
		jaux+j, jaux_loc+j, Stiff_mat[jaux+j], stiff_mat_loc[jaux_loc+j]);

		/*kew*/

	}
	  }
	  plj += nrj;
	  plj_loc += nrj_loc;
	}
	for (i = 0; i < nri_loc; i++){
	  rhs_vect_loc[pli_loc + i] = Rhs_vect[pli + i];
	}
	pli += nri;
	pli_loc += nri_loc;
  }

  //printf("Entering lar_assemble_SM_and_LV for Avp\n");
  int scope_type;
  if(Comp_type==LSC_SOLVE) scope_type=LAC_SCOPE_SM_AND_LV;
  else if(Comp_type==LSC_RESOLVE) scope_type=LAC_SCOPE_LV;
  else scope_type=LAC_SCOPE_SM_AND_LV; // default

  SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[Level_id].SM_and_LV_id;
  info = lar_assemble_SM_and_LV(SM_and_LV_id, scope_type,
				nr_dof_bl_loc, l_bl_id_loc, l_bl_nrdof_loc,
				stiff_mat_loc, rhs_vect_loc, Rewr_dofs);





  free(l_bl_id_loc);
  free(l_bl_nrdof_loc);
  free(stiff_mat_loc);
  free(rhs_vect_loc);

  return 1;

}



int lsr_ns_supg_ext_show_matrix(
				/* returns: >=0 - success code, <0 - error code */
	  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	  int Level_id          /* in: level ID */

				)
{
  return 0;
}


/**--------------------------------------------------------
lsr_ns_supg_ext_fill_precon - to prepare preconditioner by factorizing the stiffness matrix,
either only diagonal blocks or block ILU(0)
---------------------------------------------------------*/
int lsr_ns_supg_ext_fill_precon(
		/* returns: >=0 - success code, <0 - error code */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id           /* in: level ID */
				)
{

  // Create preconditioner for Avv

  // Create Schur complement matrix

  // Possibly create other preconditioner matrices (with extended
  // lad_... interface procedures)
  return 0;
}

/**--------------------------------------------------------
   lsr_ns_supg_ext_solve - to solve a system of equations, given previously constructed
   system matrix, preconditioner
---------------------------------------------------------*/
  int lsr_ns_supg_ext_solve( /* returns: convergence indicator: */
				/* 1 - convergence */
				/* 0 - noconvergence */
				/* <0 - error code */
	int Solver_id,  /* in: solver ID */
	int Ndof, 	/* in: 	the number of degrees of freedom */
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

  // Call GMRES
  return 0;
}

/**--------------------------------------------------------
lsr_ns_supg_ext_free_matrix - to free space for solver data structure (system
							  matrix, preconditioner and possibly other)
---------------------------------------------------------*/
int lsr_ns_supg_ext_free_matrix(
		int Solver_id   /* in: solver ID (used to identify the subproblem) */
				)
{

  // Free data structures in reverse order (1, 4, 3, 2)
  int isubs;   lst_mkb_levels *mkb_level;   int SM_and_LV_id;

  isubs = 2;
  mkb_level = &lsv_mkb_solver[Solver_id].subsystem[isubs].level[0];
  SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[0].SM_and_LV_id;
  //if(mkb_level->Precon != NO_PRECON) lar_free_preconditioner(SM_and_LV_id);
  lar_free_SM_and_LV(SM_and_LV_id);

  isubs = 3;
  mkb_level = &lsv_mkb_solver[Solver_id].subsystem[isubs].level[0];
  SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[0].SM_and_LV_id;
  //if(mkb_level->Precon != NO_PRECON) lar_free_preconditioner(SM_and_LV_id);
  lar_free_SM_and_LV(SM_and_LV_id);

  isubs = 4;
  mkb_level = &lsv_mkb_solver[Solver_id].subsystem[isubs].level[0];
  SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[0].SM_and_LV_id;
  //if(mkb_level->Precon != NO_PRECON) lar_free_preconditioner(SM_and_LV_id);
  lar_free_SM_and_LV(SM_and_LV_id);

  isubs = 1;
  mkb_level = &lsv_mkb_solver[Solver_id].subsystem[isubs].level[0];
  SM_and_LV_id = lsv_mkb_solver[Solver_id].subsystem[isubs].level[0].SM_and_LV_id;
  //if(mkb_level->Precon != NO_PRECON) lar_free_preconditioner(SM_and_LV_id);
  lar_free_SM_and_LV(SM_and_LV_id);



  return 0;
}

/**--------------------------------------------------------
   lsr_ns_supg_ext_compreres - to compute the residual of the left
	  preconditioned system of equations, v = M^-1 * ( b - Ax )
	  ( used also to compute the product v = M^-1 * Ax)
---------------------------------------------------------*/
extern void lsr_ns_supg_ext_compreres(
	int Solver_id,      /* in: pointer to solver data structure to be passed
				   to data structure dependent routines */
	int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Control,	/* in: indicator whether to compute residual (1)
			   or matrix-vector product (0)  */
	int Ini_zero,	/* in: flag for zero input vector X */
	int Ndof, 	/* in: number of unknowns (components of X and V) */
	double* X, 	/* in: input vector */
	double* B,	/* in: the rhs vector, if NULL take rhs */
	/*     from block data structure */
	double* V 	/* out: output vector, V = M^-1 * ( B - A*X ) */
					  )
{

  // Perform preconditioner algorithm for ns_supg extensions

}

#ifdef __cplusplus
}
#endif
