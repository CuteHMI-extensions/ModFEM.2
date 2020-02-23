/*************************************************************
File contains procedures:

- solver management routines
lsr_mkb_core_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
lsr_mkb_core_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
lsr_mkb_core_destroy - to delete a solver instance

- preconditioner management routines
lsr_mkb_core_create_precon - to create preconditioner
lsr_mkb_core_fill_precon - to prepare preconditioner by factorizing the stiffness
					matrix, either only diagonal blocks or block ILU(0)
lsr_mkb_core_destroy_precon - to free preconditioner data structure

- solver algorithms
lsr_mkb_core_gmres - to solve a system of linear equations Ax = b
	using the left preconditioned GMRES method

lsr_mkb_core_comp_norm_rhs - to compute the norm of the preconditioned rhs
				   v = || M^-1 *  b ||

lsr_mkb_core_standard - to solve the problem by standard iterations
			  x_out ~= A^-1 * b, || x^k - x^k-1 ||_max < Toler

lsr_mkb_core_vcycle - to perform one V-cycle of multigrid (as multi-level
			smoother): x_out = x_in + M^-1 * ( b - A * x_in )

lsr_mkb_core_smooth - to perform one smoothing step using different algorithms
			x_out = x_in + M^-1 * ( b - A * x_in )

lsr_mkb_core_precon - to perform preconditioning using different algorithms
			x_out = M^-1 * b

lsr_mkb_core_compreres - to compute the residual of the left preconditioned
	system of equations, v = M^-1 * ( b - Ax )
		( used also to compute the product v = M^-1 * Ax)

lsr_mkb_core_solve_coarse - to launch solver for the coarse grid problem

----------
Core solver procedure:

lsr_mkb_core_gmres	- to solve a system of linear equations Ax = b
	using the left preconditioned GMRES method as defined in
	"Templates for the Solution of Linear Systems: Building
	Blocks for Iterative Methods", Barrett, Berry, Chan,
	Demmel, Donato, Dongarra, Eijkhout, Pozo, Romine,
	and van der Vorst, SIAM Publications, 1994

The version is a general purpose procedure, independent
of the problem solved and the structure of the system matrix.
Neither the system matrix A, nor the preconditioner M, nor the
RHS vector b do not appear in lsr_mkb_core_gmres.


Required routines:

(provided internally in this file:)
lsr_mkb_core_compreres - to compute the residual of the left preconditioned
	system of equations, V = M^-1 * ( B - A*X ) (Control=1),
	used also to perform preconditioned matrix-vector
	multiplication, V = M^-1 * A * X  (Control=0)

(external:)
lsr_mkb_core_fem_proj_sol_lev - to project solution between mesh levels
lsr_mkb_core_fem_exchange_dofs - to exchange dofs between subdomains
lsr_mkb_core_fem_vec_norm - to perform parallel vector norm
lsr_mkb_core_fem_sc_prod - to perform parallel scalar product

The algorithms use la.... modules for storing system matrices, preconditioners
and right hand side vectors. The structure for each such triple is identified
by single ID (not all components of a triple must be present). la.... modules
provide also basic operations on matrices such as (detailed specification in
lah_intf.h file in lsd_mkb_core/include directory and implementations in
different lad_.... subdirectories of lsd_mkb_core):

At the end of file some independent utilities:
lsr_mkb_core_util_dvector - to allocate a double vector: name[0..ncom-1]:
lsr_mkb_core_util_ivector - to allocate an integer vector: name[0..ncom-1]:
lsr_mkb_core_util_imatrix - to allocate an integer matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
lsr_mkb_core_util_dmatrix - to allocate a double matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
lsr_mkb_core_util_chk_list - to check whether a number is on the list
lsr_mkb_core_util_put_list - to put Num on the list List with length Ll
lsr_mkb_core_util_d_zero - to zero a double vector
lsr_mkb_core_util_i_zero - to zero an integer vector
lsr_mkb_core_util_sort - to heap-sort an array
lsr_mkb_core_util_dgetrf - quasi-LU decomposition of a matrix
lsr_mkb_core_util_dgetrs - to perform forward reduction and back substitution
	of the RHS vector for solving a system of linear equations

lsr_mkb_core_get_pdeg_coarse - to get enforced pdeg for the coarse mesh

For all procedures the following BLAS/LAPACK procedures are required:
	dnrm2, daxpy, dcopy, dgemv, dtrsv, dscal, ddot,
		dgetrf, dgetrs, dgetri, dgemm, drot, drotg

History:
		05.2001 - Krzysztof Banas, initial version

*************************************************************/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#include "./lsh_mkb_core.h"
#include "./lsh_mkb_core_fem_intf.h"
#include <modfem/ls_mkb/lah_intf.h>

/* LAPACK and BLAS procedures */
#include <modfem/lin_alg_intf.h>

//#define NS_SUPG_MKB_CORE_EXTENSIONS
//#define AMG_MKB_CORE_EXTENSIONS

#ifdef NS_SUPG_MKB_CORE_EXTENSIONS
#include "../amg_ext/amg_ext.h"
//#include "./lsd_ns_supg_ext/lsh_ns_supg_ext_intf.h"
#endif


#ifdef AMG_MKB_CORE_EXTENSIONS
//#include "amg_ext/lad_amg/lah_amg_interface.h"
#include "../amg_ext/amg_ext.h"
#endif


/* minimal number of dofs to make use of DGEMV reasonable */
#define MIN_DOF_DGEMV 10000
#define SMALL 1.0e-15 /* approximation of round-off error */
#define TOL 1.0e-9    /* default tolerance for different iterative processes */
//  const double SMALL = 1.0e-15;


/* GLOBAL VARIABLES */
int lsv_mkb_core_cur_solver_id;   /* ID of the current solver */
lst_mkb_core_solvers lsv_mkb_core_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */


/**-----------------------------------------------------------
  lsr_mkb_core_init - to create a new solver instance, read its control parameters
			 and initialize its data structure
------------------------------------------------------------*/
extern int lsr_mkb_core_init( /* returns: storage type for matrices */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Solver_type, // type of solver (as defined in problem input file)
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

/* local variables */
  lst_mkb_core_levels *mkb_core_level, temp_level; /* mesh levels */
  int solver, gmres_type, krylbas;
  int i, ilev, nr_level;

/* auxiliary variables */
  FILE *fp=NULL;
  char keyword[20];


/*++++++++++++++++ executable statements ++++++++++++++++*/

// initiate global variable
  lsv_mkb_core_cur_solver_id = Solver_id;

  // initiate local variable
  if(*Max_num_levels_p<=0) *Max_num_levels_p = 1;
  nr_level = *Max_num_levels_p;


/*ok_kbw*/
  printf("\nEntering lsr_mkb_core_init: Solver_type %d, Solver id %d, Parallel %d, Nr_levels %d\n",
	 Solver_type, lsv_mkb_core_cur_solver_id, Parallel, nr_level);
  printf("Filename %s\n", Filename);
  printf("Monitor %d, Nr_iter %d, Conv_type %d, Conv_meas %15.12lf\n",
	 Monitoring_level, Max_iter, Error_type, Error_tolerance);
/*kew*/

// for special ns_supg solver

  //TODO: try to move ns supg init to another file
  if(Solver_type == -100){

	lsv_mkb_core_solver[Solver_id].subsystem[0].level[0].Solver=Solver_type;
	lsv_mkb_core_solver[Solver_id].subsystem[0].level[0].Precon=MULTI_GRID_AMG;



#ifdef NS_SUPG_MKB_CORE_EXTENSIONS

	printf("Entering NS_SUPG_MKB_CORE_EXTENSIONS for solver %d:\n",lsv_mkb_core_cur_solver_id);

	// read control data for ns_supg solver
	int info = lsr_ns_supg_ext_init(lsv_mkb_core_cur_solver_id, Parallel, Max_num_levels_p,
					Filename, Max_iter, Error_type, Error_tolerance,
					Monitoring_level);
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[0].Storage_type = LAC_STORAGE_PETSC;
	if(info < 0){
	  // error handling
	}
#endif

  } // end ns_supg solver
  else{ // classic mkb_core


	if(0 != strlen(Filename)){
	  /* read parameters from input file */
	  fp=fopen(Filename,"r");
	  if(fp==NULL) {
	printf("ERROR in opening solver input file %s !\n", Filename);
	// version 1 - restrictive
	exit(-1); // may be commented out
	// version 2 - permissive
	printf("Using default parameter values.\n");
	  }
	}

	if (0 == strlen(Filename) || fp==NULL){

	  if(nr_level <= 1){

	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].solver_id = lsv_mkb_core_cur_solver_id;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].parallel = Parallel;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].nr_subsystems = 1;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].cur_subsystem = 0;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].nr_level = 1;
	mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[0]);

	//mkb_core_level->Solver = Solver_type; - not yet implemented
	mkb_core_level->Solver = GMRES;
	mkb_core_level->SM_and_LV_id = -1;
	mkb_core_level->Storage_type = LSC_STORAGE_UNDEFINED; // selected later by lad_intf
														  // based on block size
	if(Max_iter<0) mkb_core_level->Max_iter = 100;
	else mkb_core_level->Max_iter = Max_iter;
	if(Error_type<0) mkb_core_level->Conv_type = REL_RES_INI;
	else  mkb_core_level->Conv_type = Error_type;
	if(Error_tolerance<0.0) mkb_core_level->Conv_meas = 1e-6;
	else mkb_core_level->Conv_meas = Error_tolerance;
	if(Monitoring_level < 0) mkb_core_level->Monitor=LSC_SILENT;
	else mkb_core_level->Monitor = Monitoring_level;

	/* default solver is single level GMRES with 50 Krylov vectors */
	mkb_core_level->GMRES_type = STANDARD_GMRES;
	mkb_core_level->Krylbas = 50;
	mkb_core_level->Pdeg = -1;
	mkb_core_level->Precon = BLOCK_GS; /* preconditioning by a single iteration */
	mkb_core_level->Nr_prec = 1;       /* of block Gauss-Seidel */
	mkb_core_level->ILU_k = 0;       /* not used for block Gauss-Seidel */
	mkb_core_level->Block_type = 1; /* single DOF structure subdomains */
	mkb_core_level->Nr_pre_smooth = 1;
	mkb_core_level->Nr_post_smooth = 0;

	printf("MKB_CORE iterative solver not received correct control file\n");
	printf("Using default single level GMRES solver with simple block-GS preconditioning\n");
	if(Max_iter<0){
	  printf("Using default value %d for max_iter\n", mkb_core_level->Max_iter);
	}
	if(Error_type<0){
	  printf("Using default REL_RES_INI (0) - relative wrt residual convergence measure\n");
	}
	if(Error_tolerance<0.0){
	  printf("Using default convergence measure %lf\n", mkb_core_level->Conv_meas);
	}
	if(Monitoring_level < 0){
	  printf("Using default monitoring level (0-silent)\n");
	}


	  } // end if single level requested
	  else {

	printf("No default values for multigrid solver %d:\n",lsv_mkb_core_cur_solver_id);
	printf("Specify linear solver parameters in proper input file (format: mkb_core.dat) \n");
	exit(0);

	  }

	} /* end if no file with solver data */
	else {

	  fscanf(fp,"%s",keyword);  lsr_mkb_core_util_skip_rest_of_line(fp);

	  if(!strcmp(keyword,"AMG"))
	  {
		solver = MULTI_GRID_AMG;
		  fscanf(fp,"%d", &nr_level);
		  lsr_mkb_core_util_skip_rest_of_line(fp);

		lsv_mkb_core_solver[lsv_mkb_cur_solver_id].solver_id = lsv_mkb_cur_solver_id;
		lsv_mkb_core_solver[lsv_mkb_cur_solver_id].parallel = Parallel;
		lsv_mkb_core_solver[lsv_mkb_cur_solver_id].nr_subsystems = 1;
		lsv_mkb_core_solver[lsv_mkb_cur_solver_id].cur_subsystem = 0;
		lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].nr_level = nr_level;
		mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].level[nr_level - 1]);
		mkb_core_level->Storage_type = LAC_STORAGE_PETSC;

		mkb_core_level->Solver = solver;
		mkb_core_level->Pdeg = -1;
		mkb_core_level->SM_and_LV_id = nr_level - 1;
		int coarsening_algorithm;
		int interpolation_algorithm;
		double strength_threshold;
		fscanf(fp,"%lf",&strength_threshold);
		lsr_mkb_core_util_skip_rest_of_line(fp);
		fscanf(fp,"%d %d",&(mkb_core_level->Nr_pre_smooth), &(mkb_core_level->Nr_post_smooth));
		lsr_mkb_core_util_skip_rest_of_line(fp);
		fscanf(fp,"%d", &coarsening_algorithm);
		lsr_mkb_core_util_skip_rest_of_line(fp);
		fscanf(fp,"%d", &interpolation_algorithm);
		lsr_mkb_core_util_skip_rest_of_line(fp);
		fscanf(fp,"%d %d %lg", &(mkb_core_level->Max_iter),
		  &(mkb_core_level->Conv_type),&(mkb_core_level->Conv_meas) );
		lsr_mkb_core_util_skip_rest_of_line(fp);
		fscanf(fp,"%d",&(mkb_core_level->Monitor));

		if(Max_iter>0) mkb_core_level->Max_iter = Max_iter;
		if(Error_type>0) mkb_core_level->Conv_type = Error_type;
		if(Error_tolerance>0) mkb_core_level->Conv_meas = Error_tolerance;
		if(Monitoring_level > 0) mkb_core_level->Monitor = Monitoring_level;

		int level_index;
		for(level_index=0; level_index<lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].nr_level-1; level_index++)
		{
		  mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].level[level_index]);

		  mkb_core_level->Solver = MULTI_GRID_AMG;
		  mkb_core_level->Max_iter = 10000 ;
		  mkb_core_level->Conv_type = REL_RES_INI;
		  mkb_core_level->Conv_meas = 1e-15 ;
		  mkb_core_level->Monitor = 0;
		  mkb_core_level->Nr_pre_smooth = lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].level[nr_level - 1].Nr_pre_smooth;
		  mkb_core_level->Nr_post_smooth = lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].level[nr_level - 1].Nr_post_smooth;
		  mkb_core_level->Storage_type = LAC_STORAGE_PETSC;
		  mkb_core_level->SM_and_LV_id = level_index;
		}
		#ifdef AMG_MKB_CORE_EXTENSIONS
		init_amg(coarsening_algorithm, interpolation_algorithm, strength_threshold, nr_level);
		#endif
		fclose(fp);
	  } else if(!strcmp(keyword,"NS_SUPG_SOLVER_DATA")) {
		  fscanf(fp,"%d %d %d", &solver, &nr_level, &krylbas );
		  lsr_mkb_core_util_skip_rest_of_line(fp);
		  //fclose(fp); // close input file (it can be reopended by ns_supg_ext_init)

	  //*Max_num_levels_p = 1;
	  //*Storage_type = lsv_mkb_solver[lsv_mkb_cur_solver_id].level[0].storage_type;
		  printf("Entering NS_SUPG_MKB_EXTENSIONS for solver %d:\n",lsv_mkb_cur_solver_id);

		gmres_type = STANDARD_GMRES;

		if(nr_level==1) {
		  if(solver==MULTI_GMRES) solver = GMRES;
		  if(solver==MULTI_GRID) solver = STANDARD_IT;
		}

		if(solver==GMRES||solver==STANDARD_IT){
		  /* single level GMRES or standard iterative method */

		  lsv_mkb_core_solver[lsv_mkb_cur_solver_id].solver_id = lsv_mkb_cur_solver_id;
		  lsv_mkb_core_solver[lsv_mkb_cur_solver_id].parallel = Parallel;
		  lsv_mkb_core_solver[lsv_mkb_cur_solver_id].nr_subsystems = 5;
		  lsv_mkb_core_solver[lsv_mkb_cur_solver_id].cur_subsystem = 0;
		  lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].nr_level = 1;
		  mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_cur_solver_id].subsystem[0].level[0]);

		  mkb_core_level->Solver = solver;
		  mkb_core_level->GMRES_type = gmres_type;
		  mkb_core_level->Krylbas = krylbas;
		  mkb_core_level->Pdeg = -1;

		  fscanf(fp,"%d %d %d",&(mkb_core_level->Precon),
			 &(mkb_core_level->Nr_prec),&(mkb_core_level->Block_type));
		  lsr_mkb_core_util_skip_rest_of_line(fp);
		  fscanf(fp,"%d %d %lg", &(mkb_core_level->Max_iter),
			 &(mkb_core_level->Conv_type),&(mkb_core_level->Conv_meas) );
		  lsr_mkb_core_util_skip_rest_of_line(fp);
		  fscanf(fp,"%d %d %d",
			 &(mkb_core_level->Monitor),
			 &(mkb_core_level->Nr_pre_smooth), &(mkb_core_level->Nr_post_smooth));
		  lsr_mkb_core_util_skip_rest_of_line(fp);

		  if(Max_iter>0) mkb_core_level->Max_iter = Max_iter;
		  if(Error_type>0) mkb_core_level->Conv_type = Error_type;
		  if(Error_tolerance>0) mkb_core_level->Conv_meas = Error_tolerance;
		  if(Monitoring_level > 0) mkb_core_level->Monitor = Monitoring_level;

		  #ifdef AMG_MKB_CORE_EXTENSIONS
		  //TODO: discuss and implement better initialization in ns_supg
		  init_amg(1, 1, 0.5, nr_level);
		  mkb_core_level->Storage_type = LSC_STORAGE_PETSC;
		  #endif
		  fclose(fp);
		}
	  }

	  else {

	  if(strcmp(keyword,"FINE_LEVEL") != 0) {
	printf("ERROR in reading solver input file for FINE_LEVEL !\n");
	exit(-1);
	  }

	  fscanf(fp,"%d %d %d", &solver, &nr_level, &krylbas );
	  lsr_mkb_core_util_skip_rest_of_line(fp);


	  if(nr_level > *Max_num_levels_p){
#ifdef DEBUG_SIM
	printf("Solver called with max_num_levels = %d.\n", *Max_num_levels_p);
	printf("Read from input file nr_levels = %d.\n", nr_level);
	//printf("Press any key (twice) to continue with enforced nr_levels = %d\n",
	//	     *Max_num_levels_p);
	//printf("or press CTRL C to stop\n");
	//getchar();getchar();
#endif
	nr_level = *Max_num_levels_p;
	  }

	  gmres_type = STANDARD_GMRES;

	  if(nr_level==1) {
	if(solver==MULTI_GMRES) solver = GMRES;
	if(solver==MULTI_GRID) solver = STANDARD_IT;
	  }

	  if(solver==GMRES||solver==STANDARD_IT){
	/* single level GMRES or standard iterative method */

	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].solver_id = lsv_mkb_core_cur_solver_id;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].parallel = Parallel;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].nr_subsystems = 1;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].cur_subsystem = 0;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].nr_level = 1;
	mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[0]);

	mkb_core_level->Solver = solver;
	mkb_core_level->SM_and_LV_id = -1;
	mkb_core_level->GMRES_type = gmres_type;
	mkb_core_level->Krylbas = krylbas;
	mkb_core_level->Pdeg = -1;

	fscanf(fp,"%d %d %d",&(mkb_core_level->Precon),
		   &(mkb_core_level->Nr_prec),&(mkb_core_level->Block_type));
	if(mkb_core_level->Block_type<0){
	  mkb_core_level->Storage_type = abs(mkb_core_level->Block_type);
	  mkb_core_level->Block_type = -1;
	}
	else mkb_core_level->Storage_type = LSC_STORAGE_BLOCK;

	if( (mkb_core_level->Precon>=10*MULTI_ILU && mkb_core_level->Precon<10*MULTI_ILU+9) ||
		(mkb_core_level->Precon>=10*BLOCK_ILU && mkb_core_level->Precon<10*BLOCK_ILU+9) ){
	  mkb_core_level->ILU_k = mkb_core_level->Precon%10;
	  mkb_core_level->Precon = mkb_core_level->Precon/10;
	}
	else  mkb_core_level->ILU_k = 0;

	lsr_mkb_core_util_skip_rest_of_line(fp);
	fscanf(fp,"%d %d %lg", &(mkb_core_level->Max_iter),
		   &(mkb_core_level->Conv_type),&(mkb_core_level->Conv_meas) );
	lsr_mkb_core_util_skip_rest_of_line(fp);
	fscanf(fp,"%d %d %d",
		   &(mkb_core_level->Monitor),
		   &(mkb_core_level->Nr_pre_smooth), &(mkb_core_level->Nr_post_smooth));
	lsr_mkb_core_util_skip_rest_of_line(fp);

	if(Max_iter>0) mkb_core_level->Max_iter = Max_iter;
	if(Error_type>0) mkb_core_level->Conv_type = Error_type;
	if(Error_tolerance>0) mkb_core_level->Conv_meas = Error_tolerance;
	if(Monitoring_level > 0) mkb_core_level->Monitor = Monitoring_level;

	fclose(fp);

	  } /* end if single level solver */
	  else if(solver==MULTI_GMRES||solver==MULTI_GRID){

	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].solver_id = lsv_mkb_core_cur_solver_id;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].parallel = Parallel;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].nr_subsystems = 1;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].cur_subsystem = 0;
	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].nr_level = nr_level;

	/* fine level data */
	mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[nr_level-1]);

	mkb_core_level->Solver = solver;
	mkb_core_level->SM_and_LV_id = -1;
	mkb_core_level->GMRES_type = gmres_type;
	mkb_core_level->Krylbas = krylbas;
	mkb_core_level->Pdeg = -1;

	ilev = nr_level-1;

	fscanf(fp,"%d %d %d",&(mkb_core_level->Precon),
		   &(mkb_core_level->Nr_prec),&(mkb_core_level->Block_type));
	lsr_mkb_core_util_skip_rest_of_line(fp);
	if(mkb_core_level->Block_type<0){
	  mkb_core_level->Storage_type = abs(mkb_core_level->Block_type);
	  mkb_core_level->Block_type = -1;
	}
	else mkb_core_level->Storage_type = LSC_STORAGE_BLOCK;
	int storage_type = mkb_core_level->Storage_type;
	if( (mkb_core_level->Precon>=10*MULTI_ILU && mkb_core_level->Precon<10*MULTI_ILU+9) ||
		(mkb_core_level->Precon>=10*BLOCK_ILU && mkb_core_level->Precon<10*BLOCK_ILU+9) ){
	  mkb_core_level->ILU_k = mkb_core_level->Precon%10;
	  mkb_core_level->Precon = mkb_core_level->Precon/10;
	}
	else  mkb_core_level->ILU_k = 0;

	fscanf(fp,"%d %d %lg", &(mkb_core_level->Max_iter),
		   &(mkb_core_level->Conv_type),&(mkb_core_level->Conv_meas) );
	lsr_mkb_core_util_skip_rest_of_line(fp);
	fscanf(fp,"%d %d %d",
		   &(mkb_core_level->Monitor),
		   &(mkb_core_level->Nr_pre_smooth), &(mkb_core_level->Nr_post_smooth));
	lsr_mkb_core_util_skip_rest_of_line(fp);

	if(Max_iter>0) mkb_core_level->Max_iter = Max_iter;
	if(Error_type>0) mkb_core_level->Conv_type = Error_type;
	if(Error_tolerance>0) mkb_core_level->Conv_meas = Error_tolerance;
	if(Monitoring_level > 0) mkb_core_level->Monitor = Monitoring_level;

	fscanf(fp,"%s",keyword); lsr_mkb_core_util_skip_rest_of_line(fp);
#ifdef DEBUG_LSM
	if(strcmp(keyword,"COARSE_LEVEL") != 0) {
	  printf("ERROR in reading solver input file for COARSE_LEVEL !\n");
	  exit(1);
	}
#endif

	if(nr_level>1){

	  /*coarse level data */
	  mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[0]);

	  mkb_core_level->SM_and_LV_id = -1;
	  mkb_core_level->Storage_type = storage_type; // the same for all levels

	  /* read data from a file for the coarse solver */
	  fscanf(fp,"%d %d %d", &(mkb_core_level->Solver),
		 &(mkb_core_level->Pdeg), &(mkb_core_level->Krylbas) );
	  lsr_mkb_core_util_skip_rest_of_line(fp);
	  fscanf(fp,"%d %d %d", &(mkb_core_level->Precon),
		 &(mkb_core_level->Nr_prec),&(mkb_core_level->Block_type) );
	  lsr_mkb_core_util_skip_rest_of_line(fp);
	  if( (mkb_core_level->Precon>=10*MULTI_ILU && mkb_core_level->Precon<10*MULTI_ILU+9) ||
		  (mkb_core_level->Precon>=10*BLOCK_ILU && mkb_core_level->Precon<10*BLOCK_ILU+9) ){
		mkb_core_level->ILU_k = mkb_core_level->Precon%10;
		mkb_core_level->Precon = mkb_core_level->Precon/10;
	  }
	  else  mkb_core_level->ILU_k = 0;


	  fscanf(fp,"%d %d %lg", &(mkb_core_level->Max_iter),
		 &(mkb_core_level->Conv_type),&(mkb_core_level->Conv_meas) );
	  lsr_mkb_core_util_skip_rest_of_line(fp);
	  fscanf(fp,"%d",
		 &(mkb_core_level->Monitor));
	  lsr_mkb_core_util_skip_rest_of_line(fp);
	  mkb_core_level->Nr_pre_smooth = 1;
	  mkb_core_level->Nr_post_smooth = 0;

	  fscanf(fp,"%s",keyword); lsr_mkb_core_util_skip_rest_of_line(fp);
#ifdef DEBUG_LSM
	  if(strcmp(keyword,"INTER_LEVELS") != 0) {
		printf("ERROR in reading solver input file for INTER_LEVELS !\n");
		exit(1);
	  }
#endif
	} /* end if nr_level>1 */

	if(nr_level>2){

	  /* intermediate levels data */
	  mkb_core_level = &temp_level;
	  mkb_core_level->Storage_type = storage_type; // the same for all levels
	  /* read data from a file for the intermediate solvers */
	  fscanf(fp,"%d %d %d",
		 &(mkb_core_level->Nr_pre_smooth),&(mkb_core_level->Nr_post_smooth),
		 &(mkb_core_level->Pdeg));
	  lsr_mkb_core_util_skip_rest_of_line(fp);
	  fscanf(fp,"%d %d %d",&(mkb_core_level->Precon),
		 &(mkb_core_level->Nr_prec),&(mkb_core_level->Block_type) );
	  lsr_mkb_core_util_skip_rest_of_line(fp);
	  if( (mkb_core_level->Precon>=10*MULTI_ILU && mkb_core_level->Precon<10*MULTI_ILU+9) ||
		  (mkb_core_level->Precon>=10*BLOCK_ILU && mkb_core_level->Precon<10*BLOCK_ILU+9) ){
		mkb_core_level->ILU_k = mkb_core_level->Precon%10;
		mkb_core_level->Precon = mkb_core_level->Precon/10;
	  }
	  else  mkb_core_level->ILU_k = 0;

	  fscanf(fp,"%d",
		 &(mkb_core_level->Monitor) );
	  lsr_mkb_core_util_skip_rest_of_line(fp);

	}

	fclose(fp);

	for(ilev=1;
		ilev<lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].nr_level-1;
		ilev++){

	  /* for the time being solver works in a multigrid fashion, creating */
	  /* intermediate solver structures at each mesh level */
	  /* (contrary to a domain decomposition fashion were solver structures */
	  /* are usually only for extreme coarse and fine mesh levels) */
	  lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].cur_level=ilev;
	  mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[ilev]);

	  mkb_core_level->Solver = 10 ;
	  mkb_core_level->SM_and_LV_id = -1;
	  mkb_core_level->GMRES_type = -1 ;
	  mkb_core_level->Krylbas = -1 ;
	  mkb_core_level->Pdeg = temp_level.Pdeg ;
	  mkb_core_level->Precon = temp_level.Precon ;
	  mkb_core_level->ILU_k = temp_level.ILU_k ;
	  mkb_core_level->Nr_prec = temp_level.Nr_prec ;
	  mkb_core_level->Storage_type = temp_level.Storage_type ;
	  mkb_core_level->Block_type = temp_level.Block_type ;
	  mkb_core_level->Max_iter = 10000 ;
	  mkb_core_level->Conv_type = REL_RES_INI ;
	  mkb_core_level->Conv_meas = 1e-15 ;
	  mkb_core_level->Monitor = temp_level.Monitor ;
	  mkb_core_level->Nr_pre_smooth = temp_level.Nr_pre_smooth;
	  mkb_core_level->Nr_post_smooth = temp_level.Nr_post_smooth;

	} /* end loop over levels: ilev */

	if(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[nr_level-1].Monitor>LSC_ERRORS)
	  printf("Linear solver parameters read from file for solver %d:\n",
		 lsv_mkb_core_cur_solver_id);

	  } /* end if multi-level solver */
	  } /* end if FINE_LEVEL */
	} /* end if found file with solver data */


	for(ilev=0;ilev<nr_level;ilev++){

	  /* check data for single subsystem solver */
	  lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].cur_level=ilev;
	  mkb_core_level = &lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[ilev];

	  if(nr_level>1 && mkb_core_level->Block_type==0)
	mkb_core_level->Block_type=1;
	  if(mkb_core_level->Solver==GMRES && mkb_core_level->GMRES_type==MATRIX_FREE
	 && mkb_core_level->Block_type==0 ) mkb_core_level->Block_type=1;
	  if(mkb_core_level->Precon == ADD_SCHWARZ && mkb_core_level->Block_type==0 )
	mkb_core_level->Block_type=1;
	  if(mkb_core_level->Precon == MULTI_ILU ||
	 mkb_core_level->Precon == BLOCK_ILU) mkb_core_level->Block_type=1;
	  if(mkb_core_level->Precon == NO_PRECON) mkb_core_level->Block_type=0;

	}

	/*ok_kbw*/
	if(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[nr_level-1].Monitor>LSC_ERRORS){

	  printf("\nInitiated data for solver (id = %d): %d mesh levels with parameters:\n",
		 lsv_mkb_core_cur_solver_id,
		 lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].nr_level);

	  for(ilev=0;ilev<lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].nr_level;ilev++){

	lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].cur_level=ilev;
	mkb_core_level = &(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[ilev]);

	printf("\nLevel: %d\n", ilev);
	printf("Solver type %d, GMRES type  %d, Krylbas %d, Pdeg %d\n",
		   mkb_core_level->Solver, mkb_core_level->GMRES_type, mkb_core_level->Krylbas,
		   mkb_core_level->Pdeg);
	printf("Preconditioner %d, Number of sweeps %d, Block type/size %d, Storage_type %d\n",
		   mkb_core_level->Precon, mkb_core_level->Nr_prec, mkb_core_level->Block_type, mkb_core_level->Storage_type);
	printf("Max_iter %d, Conv_type %d, Conv_meas %20.15lf\n",
		   mkb_core_level->Max_iter, mkb_core_level->Conv_type, mkb_core_level->Conv_meas);
	printf("Monitor %d, Nr_pre_smmoth %d, Nr_post_smooth %d, ILU_k %d \n",
		   mkb_core_level->Monitor,
		   mkb_core_level->Nr_pre_smooth, mkb_core_level->Nr_post_smooth,
		   mkb_core_level->ILU_k );
	  }
	}
	/*kew*/
	*Max_num_levels_p = lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].nr_level;

	} // end if classic mkb_core (not ns_supg solver)


  return(lsv_mkb_core_solver[lsv_mkb_core_cur_solver_id].subsystem[0].level[0].Storage_type);
}


/*---------------------------------------------------------
  lsr_mkb_core_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
---------------------------------------------------------*/
extern int lsr_mkb_core_solve( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
						/* <0 - error code */
	int Solver_id,  /* in: solver ID */
	int Comp_type,  /* in: indicator for the scope of computations: */
					/*   LSC_SOLVE - solve the system */
					/*   LSC_RESOLVE - resolve for the new right hand side */
	int* L_matrix_id, /* in: list of identifiers of SM and LV in lad_... module */
	int* L_ndof, 	/* in: 	the number of degrees of freedom */
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

  int nr_levels, ilev, nrdofgl, nrdofbl, max_nrdof;
  double storage, total_storage; /* storage in MBytes */
  int i,j,k, iaux, ient;

  lst_mkb_core_levels *mkb_core_level;

  // initiate SM and LV IDs for each level
  nr_levels = lsv_mkb_core_solver[Solver_id].subsystem[0].nr_level;

  for(ilev=0;ilev<nr_levels;ilev++){

	mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[0].level[ilev];

	mkb_core_level->Nrdofgl = L_ndof[ilev];
	if(mkb_core_level->SM_and_LV_id < 0){ // SM and LV IDs not yet initialized
	  mkb_core_level->SM_and_LV_id = L_matrix_id[ilev];
	}
	else{
	  if(mkb_core_level->SM_and_LV_id != L_matrix_id[ilev]){
	printf("Wrong SM and LV ID (%d!=%d) for level %d passed to lsr_mkb_core_solve. Exiting.\n",
		   mkb_core_level->SM_and_LV_id, L_matrix_id[ilev], ilev );
	exit(-1);
	  }
	}
  }

/*ok_kbw*/
  printf("\nEntering lsr_mkb_core_solve: Solver id %d, Comp_type %d, nr_levels %d\n",
	 Solver_id, Comp_type, nr_levels);
  printf("SM_and_LV_id\tNdof\n");
  for(ilev=0;ilev<nr_levels;ilev++){

	mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[0].level[ilev];

	printf("%d (%d)\t\t%d (%d)\n",
	   mkb_core_level->SM_and_LV_id, L_matrix_id[ilev],
	   mkb_core_level->Nrdofgl, L_ndof[ilev]);
  }
/*kew*/


/*------------------------------------------------------------*/
/* when monitoring the run compute the storage */
/*------------------------------------------------------------*/
  if(Monitor>LSC_ERRORS){

	int nroffbl, jaux, kaux;
	double daux;

	nr_levels = lsv_mkb_core_solver[Solver_id].subsystem[0].nr_level;

	total_storage=0.0;

	for(ilev=0;ilev<nr_levels;ilev++){

	  lsv_mkb_core_solver[Solver_id].subsystem[0].cur_level=ilev;
	  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[0].level[ilev];

	  storage = 0.0; daux = 0.0;

	  printf("\nLEVEL %d\n",ilev);

	  // get storage from LA package
	  daux = lar_get_storage( mkb_core_level->SM_and_LV_id );
	  storage += daux;
	  printf("Storage for SM, LV and preconditioner matrices: %f MBytes\n",daux);

	  if(mkb_core_level->Solver==GMRES||mkb_core_level->Solver==MULTI_GMRES)
	storage+=(mkb_core_level->Krylbas+1.0)*mkb_core_level->Nrdofgl
	  *sizeof(double)/1024.0/1024.0;
	  printf("Storage for Krylov vectors %f MBytes\n",storage-daux);
	  daux=storage;

	  storage+= 2.0*mkb_core_level->Nrdofgl
	*sizeof(double)/1024.0/1024.0;
	  printf("Storage for initial guess and solution vectors %f MBytes\n",
		 storage-daux);

	  printf("Total storage for level: over %f MBytes\n",storage);
	  total_storage += storage;
	  printf("Total storage for solver: over %f MBytes\n",total_storage);

	} /* end loop over levels */
  } /* end if computing storage */

/*------------------------------------------------------------*/
/* SOLVE THE PROBLEM FOR A GIVEN LEVEL (WHEN NEEDED)          */
/*------------------------------------------------------------*/

  nr_levels = lsv_mkb_core_solver[Solver_id].subsystem[0].nr_level;
  nrdofgl = lsv_mkb_core_solver[Solver_id].subsystem[0].level[nr_levels-1].Nrdofgl;
  //v_sol = (double *)malloc(nrdofgl*sizeof(double));

  for(ilev=0;ilev<nr_levels;ilev++){

	lsv_mkb_core_solver[Solver_id].subsystem[0].cur_level=ilev;
	mkb_core_level = &(lsv_mkb_core_solver[Solver_id].subsystem[0].level[ilev]);

	/* if it is the last (or the only) mesh or if we solve to get initial guess */
	//if( ilev==nr_levels-1 || FULL_MULTIGRID!=0 ){
	if( ilev==nr_levels-1  ){

	  int iaux,jaux,kaux;
	  double caux,daux,faux,eaux;

	  /* take into account input parameters */
	  if(*Nr_iter>0) mkb_core_level->Max_iter = *Nr_iter;
	  if(*Toler>0) mkb_core_level->Conv_meas = *Toler;
	  if(Monitor>0) mkb_core_level->Monitor = Monitor;

	  /* specify maximal number of iterations*/
	  kaux=mkb_core_level->Max_iter;

	  /* rewrite initial guess to enable computing of error in update */
	  //for (i=0;i<mkb_core_level->Nrdofgl;i++) v_sol[i] = X[i];

	  /*------- AUXILIARY STANDARD ITERATIONS SOLVER  --------*/
	  if(ilev==nr_levels-1 &&
	 ( mkb_core_level->Solver==STANDARD_IT || mkb_core_level->Solver==MULTI_GRID || mkb_core_level->Solver == MULTI_GRID_AMG) ){

	/* specify kind and value of tolerance */
	daux = mkb_core_level->Conv_meas;

	iaux = lsr_mkb_core_standard(Solver_id, 0, ilev,
					 mkb_core_level->Nrdofgl, Ini_zero, X,
					 NULL, /* RHS vector taken from block data */
					 &kaux, &daux, mkb_core_level->Monitor,
					 &caux);

	/* prepare parameters for return */
	if(Toler!=NULL) *Toler = daux;
	if(Conv_rate!=NULL) *Conv_rate = caux;
	/* number of performed iterations */
	if(Nr_iter!=NULL) *Nr_iter = kaux;

	if(mkb_core_level->Monitor>=LSC_INFO){
	  if(iaux)printf("Convergence in GS after %d iterations!!!\n",kaux);
	  else printf("No convergence in GS after %d iterations!!!\n",kaux);
	  printf("Total convergence rate (average decrease in update per iteration) %f\n",
		 caux);
	  printf("Norm of update %15.12f\n",daux);
	}

	  }
	  /*------- DEFAULT GMRES SOLVER  --------*/
	  else{

	/* specify kind and value of tolerance */
	daux = 0.0;     /* relative to initial residual tolerance */
	eaux = 0.0;	/* residual tolerance */
	faux = 0.0;	/* relative to rhs residual tolerance */
	if(mkb_core_level->Conv_type==REL_RES_INI) daux = mkb_core_level->Conv_meas;
	else if(mkb_core_level->Conv_type==ABS_RES) eaux = mkb_core_level->Conv_meas;
	else if(mkb_core_level->Conv_type==REL_RES_RHS) faux = mkb_core_level->Conv_meas;

	/*------- M A I N  S O L U T I O N  P R O C E D U R E --------*/

	//printf("calling gmres with parameters: \n:");
	//printf("Solver_id %d, mkb_core_level->Nrdofgl %d, Ini_zero %d\n",
	//       Solver_id, mkb_core_level->Nrdofgl, Ini_zero);
	//printf("mkb_core_level->Krylbas %d, kaux %d, eaux %15.12lf, daux %15.12lf, faux %15.12lf\n",
	//       mkb_core_level->Krylbas, kaux, eaux, daux, faux);

	iaux = lsr_mkb_core_gmres(Solver_id, 0, ilev,
				  mkb_core_level->Nrdofgl, Ini_zero, X,
				  NULL, /* RHS vector taken from block data */
				  mkb_core_level->Krylbas, &kaux, &eaux, &daux, &faux,
				  mkb_core_level->Monitor, &caux);




	if(mkb_core_level->Monitor>=LSC_INFO){
	  if(iaux) printf("Convergence in GMRES after %d iterations!!!\n",kaux);
	  else  printf("No convergence in GMRES after %d iterations!!!\n",kaux);
	  printf("Total convergence rate (average decrease in residual per iteration) %f\n",
		 caux);
	  printf("Norm of preconditioned residual %15.12f, relative decrease %15.12f\n",
		 eaux,daux);
	  printf("Current ratio (norm of preconditioned residual)/(norm of rhs) = %15.12f\n",faux);
	}

	/* prepare parameters for return */
	if(Toler!=NULL){
	  if(mkb_core_level->Conv_type==REL_RES_INI) *Toler = daux;
	  else if(mkb_core_level->Conv_type==ABS_RES) *Toler = eaux;
	  else if(mkb_core_level->Conv_type==REL_RES_RHS) *Toler = faux;
	}
	if(Conv_rate!=NULL) *Conv_rate = caux;
	/* number of performed iterations */
	if(Nr_iter!=NULL) *Nr_iter = kaux;

	  }
 /*kbw
#ifdef DEBUG_LSM
	printf("level %d, nrdof_glob %d, Solution:\n", ilev, mkb_core_level->Nrdofgl );
	for(i=0;i<mkb_core_level->Nrdofgl;i++)printf("%20.15lf",X[i]);
	printf("\n");
	getchar();
#endif
/*kew*/

	} /* end if solving the problem for a given level */

  } /* the end of the main loop over levels */

}

/*---------------------------------------------------------
  lsr_mkb_core_destroy - to destroy a particular instance of the solver
  (should be done in reverse order wrt creation)
---------------------------------------------------------*/
int lsr_mkb_core_destroy(
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
  )
{

  // currently nothing to do

  return(1);
}


/*---------------------------------------------------------
lsr_mkb_core_create_precon - to create preconditioner
---------------------------------------------------------*/
int lsr_mkb_core_create_precon( /* returns:   <=0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,           /* in: level ID */
  int SM_and_LV_id
  )
{

  int info=1;
  int nr_blocks_dia=0;

  int solver_type = lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Solver;
  int storage_type = lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Storage_type;

  // branch if ns_suupg extensions
  if(solver_type == MKB_CORE_NS_SUPG_SOLVER){

#ifdef NS_SUPG_MKB_EXTENSIONS
	// entry to ns_supg extensions
	info = lsr_ns_supg_ext_create_precon(Solver_id, Level_id);
#endif

  }
  // Solver_types 1..100 are reserved for mkb_core (including ns_supg_extensions)
  else if(solver_type>0 && solver_type<100){

	lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].SM_and_LV_id = SM_and_LV_id;

	int precon = lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Precon;

	if(precon != NO_PRECON){


	  int second_arg = lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Block_type;
	  if(storage_type!=LAC_STORAGE_BLOCK && (precon==MULTI_ILU || precon==BLOCK_ILU)){
	second_arg = lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].ILU_k;
	  }

	  lar_allocate_preconditioner(SM_and_LV_id, precon, second_arg);

	}

  } // end if classic mkb_core (NOT one of direct solvers through MKB interface)
  else if(solver_type == MULTI_GRID_AMG){
	  lar_allocate_preconditioner(SM_and_LV_id,
			  lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Precon,
			  lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Block_type);
  }
  else{

	// place anything else....

	printf("Wrong solver type %d in lsr_mkb_core_create_precon. Exiting.\n", solver_type);
	exit(-1);

  }

  return(info);
}


/*---------------------------------------------------------
  lsr_mkb_core_fill_precon - to prepare preconditioner by factorizing the stiffness
					matrix, either only diagonal blocks or block ILU(0)
---------------------------------------------------------*/
int lsr_mkb_core_fill_precon(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,           /* in: level ID */
  int SM_and_LV_id
  )
{

  int info = 0;

/*kbw
  printf("before calling lar_fill_preconditioner: solver %d, level %d\n",
	 Solver_id, Level_id);
/*kew*/

  int solver_type = lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Solver;

  // Solver_types 1..100 are reserved for mkb_core (including ns_supg_extensions)
  if(solver_type>0 && solver_type<100){

	if(lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Precon != NO_PRECON){
	  info = lar_fill_preconditioner( SM_and_LV_id );
	}

  } // end if classic mkb_core not ns_supg extensions
  else{

	// place for amg or anything else....

	printf("Wrong solver type %d in lsr_mkb_core_fill_precon. Exiting.\n", solver_type);
	exit(-1);

  }

  return(info);
}

/*---------------------------------------------------------
lsr_mkb_core_destroy_precon - to free preconditioner data structure
---------------------------------------------------------*/
extern int lsr_mkb_core_destroy_precon(
						 /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,           /* in: level ID */
  int SM_and_LV_id
  )
{

  if(lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id].Precon != NO_PRECON){
	lar_free_preconditioner(SM_and_LV_id);
  }
  return 0;
}


///////////////////////////// Internal  solver algorithms ///////////////

/*----------------------------------------------------------------
lsr_mkb_core_comp_norm_rhs - to compute the norm of the preconditioned rhs
				   v = || M^-1 *  b ||
----------------------------------------------------------------*/
double lsr_mkb_core_comp_norm_rhs ( /* return: the norm of the rhs vector */
	int Solver_id,      /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
		int Level_id       /* in: level data structure to be used */
	)
{

/* constants */
  int ione = 1;

/* local variables */
  lst_mkb_core_levels* mkb_core_level;  /* pointer to current level data structure */

/* auxiliary variables */
  int i, iaux, jaux, ndof;
  double daux, *vtemp, *vtemp1;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* take proper solver data structure */
#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_comp_norm_rhs! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
#endif

  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];
  ndof = mkb_core_level->Nrdofgl;
  vtemp = lsr_mkb_core_util_dvector(ndof,"temporary vector in comp_norm_rhs");
  vtemp1 = lsr_mkb_core_util_dvector(ndof,"temporary vector in comp_norm_rhs");
  lsr_mkb_core_util_d_zero(vtemp,ndof);

/* just compute initial preconditioned(!) residual with zero initial guess */
  iaux=1; /* take into account right hand side vector */
  jaux=1; /* set initial guess to zero */
  lsr_mkb_core_compreres(Solver_id, Subsystem_id, Level_id,
			iaux, jaux, ndof, vtemp, NULL, vtemp1);
  daux = dnrm2_(&ndof, vtemp1, &ione);


  free(vtemp);
  free(vtemp1);

  return(daux);
}

/*---------------------------------------------------------
lsr_mkb_core_standard - to solve the problem by standard iterations
			  x_out ~= A^-1 * b, || x^k - x^k-1 ||_max < Toler
---------------------------------------------------------*/
int lsr_mkb_core_standard( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
	int Solver_id,      /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Ndof, 	/* in: 	the number of degrees of freedom */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* B,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int*  Iter, 	/* in:	the maximum iterations to be performed */
			/* out:	actual number of iterations performed */
	double* Toler, 	/* in:	tolerance level for max norm of update */
			/* out:	the final value of max norm of update */
	int Monitor,	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data, 3 - 2+iteration data */
	double* Pconvr	/* out: convergence rate */
	)
{

/* constants */
  int ione = 1;

/* local variables */
  lst_mkb_core_levels* mkb_core_level;  /* pointer to current level data structure */
  double *uit;	        /* auxiliary vector */
  int it;		/* iteration counter */
  double diff;	        /* max norm of the difference between two iterates */
  double normini;	/* norm of initial correction */
  int iconv = 0;	/* convergence indicator */

/* auxiliary variables */
  int i, iaux;
  double *vtemp, daux, *exact;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* take proper solver data structure */
#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_standard! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
#endif

  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];
  if(Ndof != mkb_core_level->Nrdofgl) {
	printf("Not proper number of unknowns (dofs) in lsr_mkb_core_standard: %d!=%d\n",
	   Ndof,mkb_core_level->Nrdofgl);
  }

/* write parameters */
  if(Monitor>=LSC_INFO){
	printf("\nStandard iterative solver - entered level %d with tolerance: %14.12f\n",
	   Level_id, *Toler);
  }

/*kbw
  printf("Level %d - norm of X on entrance %lf\n",Level_id,dnrm2_(&Ndof, X, &ione));
  if(B!=NULL) printf("norm B on entrance %lf\n",dnrm2_(&Ndof, B, &ione));
  if(Level_id<lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].nr_level-1){
	printf("X on entrance\n");
	for(i=0;i<Ndof;i++) printf("%10.6lf",X[i]);
	printf("\n");
	printf("B on entrance\n");
	for(i=0;i<Ndof;i++) printf("%10.6lf",B[i]);
	printf("\n");
  }
/*kew*/

  uit = lsr_mkb_core_util_dvector(Ndof,"usol in GS");

  it=0; diff = *Toler+1;
  double global_diff = 0;
  while((diff>*Toler || global_diff>*Toler) &&it<*Iter){

	it++;

/* rewrite x to uit - to be able to compute the norm of update*/
	dcopy_(&Ndof, X, &ione, uit, &ione);

	iaux=1; /* take into account rhs (possibly from data structure B==NULL */
	if(it>1) Ini_zero=0; /* X is non-zero */

	if(mkb_core_level->Solver==STANDARD_IT){

	  lsr_mkb_core_smooth(Solver_id,Subsystem_id,Level_id,iaux,Ini_zero,X,B);

	}
	else if(mkb_core_level->Solver==MULTI_GRID){

	  /* v-cycle is used as multigrid solver */
	  lsr_mkb_core_vcycle(Solver_id,Subsystem_id,Level_id,iaux,Ini_zero,X,B);

	} /* for multi-level GS  */
	#ifdef AMG_MKB_CORE_EXTENSIONS
	else if(mkb_core_level->Solver==MULTI_GRID_AMG){
		int nr_pre = lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id].Nr_pre_smooth;
		int nr_post = lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id].Nr_post_smooth;
		lsr_mkb_vcycle_amg(Solver_id,Subsystem_id,Level_id,iaux,Ini_zero,X,B,nr_pre,nr_post);
	}
	#endif
	else{

	  printf("Standard iterations called with wrong Solver parameter\n");
	  exit(1);

	}

	diff=0;
	for(i=0;i<Ndof;i++) if(fabs(uit[i]-X[i])>diff) diff=fabs(uit[i]-X[i]);

	if(it==1) normini = diff;

	if(Monitor>=LSC_INFO){
	  printf("GS iteration %d, max norm of difference %.12f\n",
		 it,diff);
	}

#ifdef AMG_MKB_CORE_EXTENSIONS
	  if(mkb_core_level->Solver==MULTI_GRID_AMG){
	   global_diff = amg_get_global_diff(diff);
	  }
#endif

  }

  free(uit);

  if(it<*Iter) iconv = 1;
  else iconv = 0;

  *Iter=it;
  *Toler=diff;
  *Pconvr = exp(log(diff/normini)/(it));

  return(iconv);
}

/*---------------------------------------------------------
lsr_mkb_core_vcycle - to perform one V-cycle of multigrid (as multi-level
			smoother): x_out = x_in + M^-1 * ( b - A * x_in )
---------------------------------------------------------*/
void lsr_mkb_core_vcycle (
	int Solver_id,      /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Use_rhs,	/* in: flag for considering RHS */
	int Ini_zero,	/* in: flag for zero initial guess */
	double* X, 	/* in/out: initial guess vector and solution */
		double* B	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	)
{

/* constants */
  int ione = 1;

/* auxiliary variables */
  lst_mkb_core_levels* mkb_core_level; /* in: pointer to current level data structure */
  int i, iaux, jaux, ndof, nr_pre, nr_post;
  double daux, *vtemp, *vtemp1;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kbw
  double *vtest;
/*kew*/
#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_vcycle! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
#endif

/* take proper solver data structure */
  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];
  ndof = mkb_core_level->Nrdofgl;

/*kbw
  vtest = lsr_mkb_core_util_dvector(ndof,"temporary vector in vcycle");
  dcopy_(&ndof, X, &ione, vtest, &ione);
/*kew*/

/*kbw
  if(Level_id>=0){
	printf("level %d X on entrance\n", Level_id);
	for(i=0;i<40;i++) printf("%10.6lf",X[i]);
	for(i=ndof-40;i<ndof;i++) printf("%10.6lf",X[i]);
	printf("\n");
	if(B!=NULL){
	  printf("level %d B on entrance\n", Level_id);
	  for(i=0;i<40;i++) printf("%10.6lf",B[i]);
	  for(i=ndof-40;i<ndof;i++) printf("%10.6lf",B[i]);
	  printf("\n");
	}
	getchar();getchar();
  }
/*kew*/

/* if not on the coarsest level */
  if(Level_id>0){

	vtemp = lsr_mkb_core_util_dvector(ndof,"temporary vector in vcycle");
	vtemp1 = lsr_mkb_core_util_dvector(ndof,"temporary vector in vcycle");

	nr_pre=mkb_core_level->Nr_pre_smooth;
	nr_post=mkb_core_level->Nr_post_smooth;

/*kbw
	printf("nr_pre %d, nr_post %d\n", nr_pre, nr_post);
	getchar();getchar();
/*kew*/

/* perform presmoothing using proper algorithms */
	for(i=0;i<nr_pre;i++) {
	  lsr_mkb_core_smooth(Solver_id, Subsystem_id, Level_id, Use_rhs, Ini_zero, X, B);
	}

/* indicate initial guess is no longer zero */
	Ini_zero=0;

/*kbw
  if(Level_id>=0){
	printf("M^-1*(b-Ax) after %d smoothing steps: \n",nr_pre);
	//for(i=0;i<ndof;i++) printf("%10.6lf",X[i]-vtest[i]);
	for(i=0;i<40;i++) printf("%10.6lf",X[i]);
	for(i=ndof-40;i<ndof;i++) printf("%10.6lf",X[i]);
	printf("\n");
	getchar();
  }
/*kew*/

	if(mkb_core_level->Monitor>LSC_INFO)
	  printf("\nProjecting error from level %d to level %d\n",Level_id,Level_id-1);

/* transpose of L2 projection from coarse to fine */
/* compute the residual of the not-preconditioned system vtemp = b - Ax */
	//lsr_mkb_core_compres(mkb_core_level, Use_rhs, Ini_zero, ndof, X, B, vtemp);
	lar_compute_residual(mkb_core_level->SM_and_LV_id, Use_rhs, Ini_zero, ndof, X, B, vtemp);

/*kbw
  if(Level_id>=0){
	printf("residual on fine grid: \n");
	for(i=0;i<40;i++) printf("%10.6lf",vtemp[i]);
	for(i=ndof-40;i<ndof;i++) printf("%10.6lf",vtemp[i]);
	printf("\n");
	getchar();
  }
/*kew*/

/* project residual on coarse grid - stored in vtemp1 */
	lsr_mkb_core_fem_proj_sol_lev(Solver_id,Subsystem_id,Level_id,
				  vtemp,Level_id-1,vtemp1);

/*kbw
  if(Level_id>=0){
	printf("residual projected on coarse grid (10 entries only): \n");
	for(i=0;i<10;i++) printf("%10.6lf",vtemp1[i]);
	printf("\n");
	getchar();
  }
/*kew*/

/* solve problem on coarser grid */
/* set current solver data structure number */
	lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level=Level_id-1;

/* vtemp1 is the rhs vector for the coarse grid correction problem */
/* vtemp is the initial guess and is set to zero */
	iaux=1; lsr_mkb_core_util_d_zero(vtemp,ndof);

/* v-cycle is used for coarse grid correction (if one wants w-cycle */
/* lsr_mkb_core_standard should be used here) - vtemp is corrected vtemp1*/
	lsr_mkb_core_vcycle(Solver_id,Subsystem_id,Level_id-1,iaux,iaux,vtemp,vtemp1);

/*kbw
  if(Level_id>=0){
	printf("after coarse solve: \n");
	for(i=0;i<ndof;i++) printf("%10.6lf",vtemp[i]);
	printf("\n");
	getchar();
  }
/*kew*/

/* reset current solver data structure number */
	lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level=Level_id;

/* project back to finer grid */
	if(mkb_core_level->Monitor>LSC_INFO)
	  printf("\nProjecting solution from level %d to level %d\n",Level_id-1,Level_id);

	lsr_mkb_core_fem_proj_sol_lev(Solver_id,Subsystem_id,Level_id-1,
				  vtemp,Level_id,vtemp1);


/*kbw
  printf("coarse grid correction projected on fine grid\n");
  for(i=0;i<ndof;i++) printf("%10.6lf",vtemp1[i]);
  printf("\n");
  getchar();
/*kew*/


/* sum up the smoothed initial guess and coarse grid correction */
/* daxpy: X:= daux*vtemp + X */
	daux=1;
	daxpy_(&ndof,&daux,vtemp1,&ione,X,&ione);
	Ini_zero=0; /* indicate initial guess is no longer zero */

	free(vtemp);
	free(vtemp1);

/* perform post-smoothing using proper algorithms */
	for(i=0;i<nr_post;i++) {
	  lsr_mkb_core_smooth(Solver_id, Subsystem_id, Level_id, Use_rhs, Ini_zero, X, B);
	}

/*kbw
  if(Level_id>=0){
	printf("V after smoothing %d\n", nr_post);
	for(i=0;i<ndof;i++) printf("%10.6lf",X[i]-vtest[i]);
	printf("\n");
	getchar();
  }
  free(vtest);
/*kew*/

  } /* end if not on the coarsest level */
  else {

	lsr_mkb_core_solve_coarse(Solver_id, Subsystem_id, Level_id, Ini_zero, X, B);

/*kbw
  if(Level_id>=0){
	printf("X after solving coarse problem\n");
	for(i=0;i<ndof;i++) printf("%10.6lf",X[i]);
	printf("\n");
	getchar();
  }
/*kew*/

  }

  return;
}

#ifdef AMG_MKB_CORE_EXTENSIONS
void lsr_mkb_vcycle_amg (
	int Solver_id,
	int Subsystem_id,
	int Level_id,
	int Use_rhs,
	int Ini_zero,
	double* X,
	double* B,
	int nr_pre,
	int nr_post
	)
{

/* auxiliary variables */
  int i;
/*++++++++++++++++ executable statements ++++++++++++++++*/
  if(Level_id == 0)
  {
	if(lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].nr_level > 1)
	{
		amg_project_solution_from_level(Level_id + 1);
	}
	lsr_mkb_core_smooth(Solver_id, Subsystem_id, Level_id, Use_rhs, Ini_zero, X, B);
	if(lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].nr_level > 1)
	{
		amg_project_solution_to_level(Level_id + 1);;
	}
  }
  else
  {
	  if(Level_id + 1 < lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].nr_level)
		  amg_project_solution_from_level(Level_id + 1);
	  for(i=0;i<nr_pre;i++)
	  {
		lsr_mkb_core_smooth(Solver_id, Subsystem_id, Level_id, Use_rhs, Ini_zero, X, B);
		Ini_zero=0;
	  }
	  lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level=Level_id-1;
	  lsr_mkb_vcycle_amg(Solver_id,Subsystem_id,Level_id-1,Use_rhs,Ini_zero,X,B,nr_pre,nr_post);
	  lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level=Level_id;

	  for(i=0;i<nr_post;i++)
		lsr_mkb_core_smooth(Solver_id, Subsystem_id, Level_id, Use_rhs, Ini_zero, X, B);
	  if(Level_id + 1 < lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].nr_level)
		  amg_project_solution_to_level(Level_id + 1);

  }
  return;
}
#endif

/*---------------------------------------------------------
lsr_mkb_core_smooth - to perform one smoothing step using different algorithms
			x_out = x_in + M^-1 * ( b - A * x_in )
---------------------------------------------------------*/
void lsr_mkb_core_smooth(
	int Solver_id,      /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Use_rhs,	/* in: flag for considering RHS */
	int Ini_zero,	/* in: flag for zero initial guess */
	double* X, 	/* in/out: initial guess vector and solution */
		double* B	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	)
{

/* constants */
  int ione = 1;

/* auxiliary variables */
  lst_mkb_core_levels* mkb_core_level; /* in: pointer to current level data structure */
  int i, iaux, jaux, ndof;
  double daux, *vtemp, *vtemp1;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_smooth! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
#endif
/* take proper solver data structure */
  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];
  ndof = mkb_core_level->Nrdofgl;

/*kbw
  if(Level_id>=0){
	printf("level %d X on entrance to lsr_mkb_core_smooth\n", Level_id);
	for(i=0;i<40;i++) printf("%20.15lf",X[i]);
	printf("\n");
	//for(i=0;i<ndof;i++) printf("%5d%15.10lf",i/4+1,X[i]);
	for(i=ndof-40;i<ndof;i++) printf("%20.15lf",X[i]);
	printf("\n");
	printf("after level %d X on entrance\n", Level_id);
  }
/*kew*/
/*kbw
  if(Level_id>=0){
	if(B!=NULL){
	  printf("level %d B on entrance to lsr_mkb_core_smooth\n", Level_id);
	  for(i=0;i<40;i++) printf("%20.15lf",B[i]);
	  printf("\n");
	  //for(i=0;i<ndof;i++) printf("%5d%15.10lf",i/4+1,B[i]);
	  for(i=ndof-40;i<ndof;i++) printf("%20.15lf",B[i]);
	  printf("\n");
	}
	printf("after level %d B on entrance\n", Level_id);
	getchar();
  }
/*kew*/

/* different preconditioner options used for smoothing */
  if(mkb_core_level->Solver==GMRES&&mkb_core_level->GMRES_type==MATRIX_FREE){

	// if(mkb_core_level->Precon==BLOCK_GS){
	  // lsr_mkb_core_mfmiter(mkb_core_level,Use_rhs,Ini_zero,ndof,X,B); // not implemented for ModFEM
	// }
	// else if(mkb_core_level->Precon==BLOCK_JACOBI){
	  // lsr_mkb_core_mfaiter(mkb_core_level,Use_rhs,Ini_zero,ndof,X,B); // not implemented for ModFEM
	// }

  }
  else{

	if(mkb_core_level->Precon==MULTI_ILU || mkb_core_level->Precon==BLOCK_ILU){

/* incomplete LU preconditioning : */
	  vtemp = lsr_mkb_core_util_dvector(ndof,"temporary vector in smooth");
	  vtemp1 = lsr_mkb_core_util_dvector(ndof,"temporary vector in smooth");

	  iaux = 0; /* initial guess not zero */
	  //lsr_mkb_core_compres(mkb_core_level, Use_rhs, iaux, ndof, X, B, vtemp);
	  lar_compute_residual ( mkb_core_level->SM_and_LV_id,  Use_rhs, iaux, ndof, X, B, vtemp);
/* routine returned: vtemp = B - A * X */

	  //lsr_mkb_core_rhsub(mkb_core_level,ndof,vtemp1,vtemp);
	  lar_perform_rhsub( mkb_core_level->SM_and_LV_id, ndof, vtemp1, vtemp);
/* routine returned: vtemp1 = M^-1 * vtemp = M^-1 * ( B - A * X ) */

/* add X to get  X:= X + M^-1 ( B - A * X )*/
	  daux = 1.;
	  daxpy_(&ndof, &daux, vtemp1, &ione, X, &ione);

	  free(vtemp1);
	  free(vtemp);

/*kbw
  if(Level_id>=0){
	printf("X before exchange dofs\n");
	//for(i=0;i<ndof;i++) printf("%5d%15.10lf",i/4+1,X[i]);
	for(i=0;i<40;i++) printf("%20.15lf",X[i]);
	printf("\n");
	getchar();
  }
/*kew*/

/*||begin||*/
	  lsr_mkb_core_fem_exchange_dofs(Solver_id,Subsystem_id,Level_id,X);
/*||end||*/

/*kbw
  if(Level_id>=0){
	printf("X after exchange dofs\n");
	//for(i=0;i<ndof;i++) printf("%5d%15.10lf",i/4+1,X[i]);
	for(i=0;i<ndof;i++) printf("%20.15lf",X[i]);
	printf("\n");
	getchar();
  }
/*kew*/

	}
	else if(mkb_core_level->Precon==BLOCK_JACOBI || mkb_core_level->Precon==BLOCK_GS){

	int nr_iter = mkb_core_level->Nr_prec;
	int nr_outer_iter = 1;

	if(nr_iter<0){
	  nr_outer_iter = -nr_iter;
	}

	int i_iter;
	for(i_iter=0; i_iter < nr_outer_iter; i_iter++){

	  if(mkb_core_level->Nr_prec<0 && i_iter%2==0) nr_iter=-1;
	  else if(mkb_core_level->Nr_prec<0 && i_iter%2==1) nr_iter=-2;

/*kbw
  if(Level_id>=0){
	printf("smooth (level %d): X before BJ or GS\n", Level_id);
	//for(i=0;i<ndof;i++) printf("%5d%15.10lf",i/4+1,X[i]);
	printf("from %d to %d\n", 0, 40);
	for(i=0;i<40;i++) printf("%20.15lf",X[i]);
	printf("\n");
	printf("from %d to %d\n", 1160, 1200);
	for(i=1160;i<1200;i++) printf("%20.15lf",X[i]);
	printf("\n");
	printf("from %d to %d\n", ndof-40, ndof);
	for(i=ndof-40;i<ndof;i++) printf("%20.15lf",X[i]);
	printf("\n");
	//getchar();
  }
/*kew*/

	  lar_perform_BJ_or_GS_iterations( mkb_core_level->SM_and_LV_id,
					   Use_rhs, Ini_zero, nr_iter,
					   ndof, X, B);

	  /* if(abs(mkb_core_level->Block_type)==0) { */
	  /* 	lsr_mkb_core_blsiter(mkb_core_level,Use_rhs,Ini_zero,ndof,X,B); */
	  /* } */
	  /* else { */
	  /* 	lsr_mkb_core_blliter(mkb_core_level,Use_rhs,Ini_zero,ndof,X,B); */
	  /* } */

/*kbw
  if(Level_id==0){
	printf("smooth (level %d): X before exchange dofs\n", Level_id);
	for(i=0;i<ndof;i++) printf("%15.10lf",X[i]);
	//for(i=0;i<ndof;i++) printf("%5d%15.10lf",i/4+1,X[i]);
	printf("from %d to %d\n", 0, 40);
	for(i=0;i<40;i++) printf("%20.15lf",X[i]);
	printf("\n");
	printf("from %d to %d\n", 1160, 1200);
	for(i=1160;i<1200;i++) printf("%20.15lf",X[i]);
	printf("\n");
	printf("from %d to %d\n", ndof-40, ndof);
	for(i=ndof-40;i<ndof;i++) printf("%20.15lf",X[i]);
	printf("\n");
	//getchar();
  }
/*kew*/

/*||begin||*/
	  lsr_mkb_core_fem_exchange_dofs(Solver_id,Subsystem_id,Level_id,X);
/*||end||*/

/*kbw
  if(Level_id==0){
	printf("X after exchange dofs\n");
	for(i=0;i<ndof;i++) printf("%15.10lf",X[i]);
	printf("from %d to %d\n", 0, 40);
	for(i=0;i<40;i++) printf("%20.15lf",X[i]);
	printf("\n");
	printf("from %d to %d\n", 1160, 1200);
	for(i=1160;i<1200;i++) printf("%20.15lf",X[i]);
	printf("\n");
	printf("from %d to %d\n", ndof-40, ndof);
	for(i=ndof-40;i<ndof;i++) printf("%20.15lf",X[i]);
	printf("\n");
	//getchar();
  }
/*kew*/
	}

	}
	else if(mkb_core_level->Solver == MULTI_GRID_AMG)
	{
		lar_perform_BJ_or_GS_iterations( Level_id,
						 Use_rhs, Ini_zero, mkb_core_level->Nr_prec,
						 ndof, X, B);
	}
	else{

	  printf("Unknown option for smoothing (not preconditioning!) %d\n",
		 mkb_core_level->Precon);
	  exit(1);

	}

  }




  return;

}


/*---------------------------------------------------------
lsr_mkb_core_precon - to perform preconditioning using different algorithms
			x_out = M^-1 * b
---------------------------------------------------------*/
void lsr_mkb_core_precon(
	int Solver_id,      /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	double* X, 	/* in/out: initial guess vector and solution */
		double* B	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	)
{

/* constants */
  int ione = 1;

/* auxiliary variables */
  lst_mkb_core_levels* mkb_core_level; /* in: pointer to current level data structure */
  int i, iaux, jaux, ndof;
  double daux, *vtemp, *vtemp1;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_precon! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
#endif
/* take proper solver data structure */
  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];
  ndof = mkb_core_level->Nrdofgl;

/* different preconditioner options */
  if(mkb_core_level->Precon==ADD_SCHWARZ){

	iaux = 1; /* initial guess is zero for preconditioning */
	jaux = 1; /* use supplied B as rhs */
	lsr_mkb_core_util_d_zero(X,ndof);

/*kbw
	printf("lsr_mkb_core_precon - V on entrance (ndof == %d)\n", ndof);
  int i;
  for(i=0;i<20;i++) printf("%20.15lf",X[i]);
  //for(i=0;i<Ndof;i++) printf("%10.6lf",V[i]);
  printf("\n");
  getchar();
/*kew*/

	//lsr_mkb_core_blliter(mkb_core_level,jaux,iaux,ndof,X,B);
	lar_perform_BJ_or_GS_iterations( mkb_core_level->SM_and_LV_id,
					 jaux, iaux,  mkb_core_level->Nr_prec, ndof, X, B);

  }
  else if(mkb_core_level->Precon==BLOCK_ILU){

/* incomplete LU preconditioning : */
	//lsr_mkb_core_rhsub(mkb_core_level,ndof,X,B);
	lar_perform_rhsub( mkb_core_level->SM_and_LV_id, ndof, X, B);
/* routine returned: X = M^-1 * B */

  }
  else{

	printf("Unknown option for preconditioning (not smoothing!) %d\n",
	   mkb_core_level->Precon);
	exit(1);

  }

/*||begin||*/
  lsr_mkb_core_fem_exchange_dofs(Solver_id,Subsystem_id,Level_id,X);
/*||end||*/

  return;

}

/*---------------------------------------------------------
lsr_mkb_core_compreres - to compute the residual of the left preconditioned
	system of equations, v = M^-1 * ( b - Ax )
		( used also to compute the product v = M^-1 * Ax)
---------------------------------------------------------*/
void lsr_mkb_core_compreres (
	int Solver_id,   /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Control, /* in: indicator whether to compute residual (1-use RHS) */
			 /*	or matrix-vector product (0-do not use RHS) */
	int Ini_zero,/* in: flag for zero initial guess */
	int Ndof,    /* in: number of unknowns (components of x) */
	double* X,   /* in: initial guess vector */
		double* B,   /* in:  the rhs vector, if NULL take rhs */
			 /*      from block data structure */
	double* V    /* out: initial residual, v = M^-1*(b-Ax) */
	)
{

/* constants */
  int ione = 1;

/* auxiliary variables */
  lst_mkb_core_levels* mkb_core_level; /* in: pointer to current level data structure */
  int i, iaux, jaux;
  double daux, *vtemp;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_compreres! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
#endif
/* take proper solver data structure */
  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];

  if(Ndof != mkb_core_level->Nrdofgl) {
	printf("Not proper number of unknowns (dofs) in lsr_mkb_core_compreres: %d!=%d\n",
	   Ndof,mkb_core_level->Nrdofgl);
  }

/*kbw
  if(Level_id>=0){
	printf("it_compreres: Solver_id %d, Subsystem_id %d, Level_id %d, Control %d, Ini_zero %d, Ndof %d\n",
	   Solver_id, Subsystem_id, Level_id, Control, Ini_zero, Ndof);

	printf("X on entrance\n");
	for(i=0;i<40;i++) printf("%10.6lf",X[i]);
	for(i=Ndof-40;i<Ndof;i++) printf("%10.6lf",X[i]);
	printf("\n");
	if(B!=NULL){
	  printf("B on entrance\n");
	  for(i=0;i<40;i++) printf("%10.6lf",B[i]);
	  for(i=Ndof-40;i<Ndof;i++) printf("%10.6lf",B[i]);
	  printf("\n");
	}
	getchar();
  }
/*kew*/

/*kbw
  if(Level_id>=0){
	if(Ndof>=426822){
	  printf("X on entrance in compreres 77124 %lf  426822 %lf\n", X[77124], X[426822]);
	} else {
	  printf("X on entrance in compreres (Ndof %d) 77124 %lf\n", Ndof, X[77124]);

	}

	//printf("X on entrance in compreres\n");
	//for(i=0;i<20;i++) printf("%10.6lf",X[i]);
	//for(i=0;i<Ndof;i++) printf("%10.6lf",X[i]);
	//printf("\nafter X on entrance in compreres\n");
	//getchar();getchar();
  }
/*kew*/

  // first classic mkb_core
  if(lsv_mkb_core_solver[Solver_id].nr_subsystems==1){

	/* if there is no preconditioning */
	if(mkb_core_level->Precon == NO_PRECON){

	  /* just compute the residual of the not-preconditioned system */
	  //lsr_mkb_core_compres(mkb_core_level, Control, Ini_zero, Ndof, X, B, V);
	  lar_compute_residual ( mkb_core_level->SM_and_LV_id, Control, Ini_zero, Ndof, X, B, V);

	}
	/* single level preconditioning */
	else if(mkb_core_level->Precon==ADD_SCHWARZ||mkb_core_level->Precon==BLOCK_ILU){

	  /* incomplete LU preconditioning or */
	  /* textbook additive Schwarz preconditioning with */
	  /* 1. multiplying X by A (or computing not preconditioned residual) */
	  /* 2. standard iterations with residual as rhs and zero initial guess */

	  vtemp = lsr_mkb_core_util_dvector(Ndof,"temporary vector in compreres");

	  //lsr_mkb_core_compres(mkb_core_level, Control, Ini_zero, Ndof, X, B, vtemp);
	  lar_compute_residual ( mkb_core_level->SM_and_LV_id, Control, Ini_zero, Ndof, X, B, vtemp);
	  /* the result is vtemp = B - A * X */

/*kbw
	  printf("entering precon on the last level!!!!!!!!!\n");
  if(Level_id>=0){
	printf("V before precon\n");
	for(i=0;i<40;i++) printf("%10.6lf",vtemp[i]);
	printf("\n");
  }
	  getchar();getchar();
/*kew*/

	  /* if on the coarsest or the only one level */
	  if(Level_id==0){

	/* perform preconditioning */
	lsr_mkb_core_precon(Solver_id, Subsystem_id, Level_id, V, vtemp);
	/* the result is V = M^-1 * ( B - A * X ) */

	  } /* end if on the coarsest level */
	  else {

	/* perform multi-level preconditioning */
	printf("Multi-level preconditioning only based on smoothing. Exiting\n");
	exit(0);

	  }

/*kbw
	  printf("after precon on the last level!!!!!!!!!\n");
  if(Level_id>=0){
	printf("V after precon\n");
	for(i=0;i<40;i++) printf("%10.6lf",V[i]);
	printf("\n");
  }
	  getchar();getchar();
/*kew*/

	  free(vtemp);

	}
	/* single or multi level preconditioning using smoothing */
	else{

	  /* copy X to V */
	  dcopy_(&Ndof, X, &ione, V, &ione);

	  /* if on the coarsest or the only one level */
	  if(Level_id==0){

/*kbw
	  printf("entering smooth on the last level!!!!!!!!!\n");
  if(Level_id>=0){
	printf("V before smoothing\n");
	for(i=0;i<40;i++) printf("%10.6lf",V[i]);
	printf("\n");
  }
	  getchar();getchar();
/*kew*/


	/* perform smoothing as preconditioning */
	lsr_mkb_core_smooth(Solver_id, Subsystem_id, Level_id, Control, Ini_zero, V, B);
	/* the result is V = X + M^-1 * ( B - A * X ) */

/*kbw
  if(Level_id>=0){
	printf("V after smoothing (%d times)\n", mkb_core_level->Nr_pre_smooth);
	for(i=0;i<Ndof;i++) printf("%10.6lf",V[i]);
	//for(i=0;i<20;i++) printf("%10.6lf",V[i]);
	printf("\n");
  }
	  getchar();getchar();
/*kew*/

	  } /* end if on the coarsest level */
	  else {

	/* one v-cycle multigrid is used as approximate solver-preconditioner */
	lsr_mkb_core_vcycle(Solver_id, Subsystem_id, Level_id, Control, Ini_zero, V, B);
	/* the result is V = X + M^-1 * ( B - A * X ) */

	  }

	  /* subtract X from V to get V = M^-1 * ( B - A * X ) */
	  daux = -1.;
	  daxpy_(&Ndof, &daux, X, &ione, V, &ione);

	} /* end if there is preconditioning */

	if(Control==0){
	  /* multiply by -1 to get v = M^-1 *  Ax  */
	  daux= -1;
	  dscal_(&Ndof, &daux, V, &ione);
	}

  } // end if classic mkb (nr_subsystems==1)
  else if(lsv_mkb_core_solver[Solver_id].nr_subsystems==5 || lsv_mkb_core_solver[Solver_id].nr_subsystems==6){

#ifdef NS_SUPG_MKB_CORE_EXTENSIONS
	// entry to ns_supg extensions
	lsr_ns_supg_ext_compreres(
	Solver_id,   /* in: solver data structure to be used */
		Subsystem_id,  /* in: subsystem data structure to be used */
	Level_id,	/* in: index of the mesh (level) */
	Control, /* in: indicator whether to compute residual (1-use RHS) */
			 /*	or matrix-vector product (0-do not use RHS) */
	Ini_zero,/* in: flag for zero initial guess */
	Ndof,    /* in: number of unknowns (components of x) */
	X,   /* in: initial guess vector */
		B,   /* in:  the rhs vector, if NULL take rhs */
			 /*      from block data structure */
	V    /* out: initial residual, v = M^-1*(b-Ax) */
				  );
#endif


  } // end if ns_supg extensions


/*kbw
  if(Level_id==0){
	//if(Ndof>=426822){
	//  printf("X on entrance in compreres 77124 %lf  426822 %lf\n", V[77124], V[426822]);
	//} else {
	//  printf("X on entrance in compreres (Ndof %d) 77124 %lf\n", Ndof, V[77124]);
	//
	//}
	printf("V after compreres\n");
	for(i=0;i<Ndof;i++) printf("%10.6lf",V[i]);
	//for(i=0;i<20;i++) printf("%10.6lf",V[i]);
	printf("\n");
	//getchar();getchar();
  }
/*kew*/

  return;
}


void lsr_mkb_core_compres (
	int Solver_id,   /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Control, /* in: indicator whether to compute residual (1-use RHS) */
			 /*	or matrix-vector product (0-do not use RHS) */
	int Ini_zero,/* in: flag for zero initial guess */
	int Ndof,    /* in: number of unknowns (components of x) */
	double* X,   /* in: initial guess vector */
		double* B,   /* in:  the rhs vector, if NULL take rhs */
			 /*      from block data structure */
	double* V    /* out: initial residual, v = M^-1*(b-Ax) */
	)
{

/* constants */
  int ione = 1;

/* auxiliary variables */
  lst_mkb_core_levels* mkb_core_level; /* in: pointer to current level data structure */
  int i, iaux, jaux;
  double daux, *vtemp;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_compreres! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
#endif
/* take proper solver data structure */
  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];

  if(Ndof != mkb_core_level->Nrdofgl) {
	printf("Not proper number of unknowns (dofs) in lsr_mkb_core_compreres: %d!=%d\n",
	   Ndof,mkb_core_level->Nrdofgl);
  }

  // first classic mkb_core
  if(lsv_mkb_core_solver[Solver_id].nr_subsystems==1){

	  /* just compute the residual of the not-preconditioned system */
	  //lsr_mkb_core_compres(mkb_core_level, Control, Ini_zero, Ndof, X, B, V);
	  lar_compute_residual ( mkb_core_level->SM_and_LV_id, Control, Ini_zero, Ndof, X, B, V);
	  //rpgwars

  } // end if classic mkb (nr_subsystems==1)
  else if(lsv_mkb_core_solver[Solver_id].nr_subsystems==5 || lsv_mkb_core_solver[Solver_id].nr_subsystems==6){

#ifdef NS_SUPG_MKB_CORE_EXTENSIONS
	  lsr_ns_supg_ext_compres_vector(mkb_core_level->SM_and_LV_id, Ndof, X, B, V);
#endif

  } // end if ns_supg extensions

  return;
}


/*---------------------------------------------------------
lsr_mkb_core_solve_coarse - to launch solver for the coarse grid problem
---------------------------------------------------------*/
int lsr_mkb_core_solve_coarse(  /* returns: 1 - success; <=0 - error code*/
	int Solver_id,        /* in: solver data structure to be used */
		int Subsystem_id,  /* in: subsystem data structure to be used */
	int Level_id,	/* in: index of the mesh (level) */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X,	/* in/out: initial guess and solution vector */
	double* B	/* in: rhs vector */
	)
{

/* auxiliary variables */
  lst_mkb_core_levels *mkb_core_level;
  int i, iaux, jaux, kaux, iblock, ierr;
  double daux, eaux, faux, gaux, normb;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG
  if(Level_id!=lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level){
	printf("Wrong cur_level ( %d != %d ) in lsr_mkb_core_solve_coarse! Exiting!",
	   Level_id,lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].cur_level);
	exit(0);
  }
  if(Level_id!=0){
	printf("Wrong level ( %d != 0 ) in lsr_mkb_core_solve_coarse! Exiting!",
	   Level_id);
	exit(0);
  }
#endif
/* take proper solver data structure */
  mkb_core_level = &lsv_mkb_core_solver[Solver_id].subsystem[Subsystem_id].level[Level_id];

/* specify maximal number of iterations */
  kaux=mkb_core_level->Max_iter;

  if(mkb_core_level->Solver==STANDARD_IT||mkb_core_level->Solver==MULTI_GRID){

/* specify kind and value of tolerance */
	daux = mkb_core_level->Conv_meas;

/*------- STANDARD ITERATIONS SOLVER  --------*/
	iaux=lsr_mkb_core_standard(Solver_id, Subsystem_id, Level_id, mkb_core_level->Nrdofgl,
			   Ini_zero, X, B, &kaux, &daux, mkb_core_level->Monitor, &gaux);

	if(mkb_core_level->Monitor>=LSC_INFO){
	  if(iaux)printf("Convergence in GS after %d iterations!!!\n",kaux);
	  else printf("No convergence in GS after %d iterations!!!\n",kaux);
	  printf("Total convergence rate (average decrease in update per iteration) %f\n",
		 gaux);
	  printf("Norm of update %15.12f\n",daux);
	}

  }
  else{

/* specify kind and value of tolerance */
	daux = 0.0;       /* relative to initial residual tolerance */
	eaux = 0.0;	/* residual tolerance */
	faux = 0.0;	/* relative to rhs residual tolerance */
	if(mkb_core_level->Conv_type==REL_RES_INI) daux = mkb_core_level->Conv_meas;
	else if(mkb_core_level->Conv_type==ABS_RES) eaux = mkb_core_level->Conv_meas;
	else if(mkb_core_level->Conv_type==REL_RES_RHS) faux = mkb_core_level->Conv_meas;

/* compute the norm of rhs (excluded since it is not parallel) */
	iaux=1;
	normb=1.0;

/*kbw
	printf("Before entering GMRES: Solver_id %d, Subsystem_id %d, Nrdofgl %d\n",
	   Solver_id, Subsystem_id, mkb_core_level->Nrdofgl);
	printf("  , Ini_zero %d, Krylbas %d, Monitor %d\n",
	   Ini_zero, mkb_core_level->Krylbas, mkb_core_level->Monitor);
printf("rhs:\n");
for(i=0;i<40;i++) printf("%10.6lf",B[i]);
//for(i=0;i<mkb_core_level->Nrdofgl;i++) printf("%10.6lf",B[i]);
printf("\nini_guess:\n");
for(i=0;i<40;i++) printf("%10.6lf",X[i]);
//for(i=0;i<mkb_core_level->Nrdofgl;i++) printf("%10.6lf",X[i]);
printf("\n");
/*kew*/

/*--- GMRES SOLVER ---*/
	iaux = lsr_mkb_core_gmres(Solver_id, Subsystem_id, Level_id, mkb_core_level->Nrdofgl,
			Ini_zero, X, B, mkb_core_level->Krylbas,
			&kaux, &eaux, &daux, &faux, mkb_core_level->Monitor, &gaux);

	if(mkb_core_level->Monitor>=LSC_INFO){
	  if(iaux) printf("Convergence in GMRES after %d iterations!!!\n",kaux);
	  else  printf("No convergence in GMRES after %d iterations!!!\n",kaux);
	  printf("Total convergence rate (average decrease in residual per iteration) %f\n",
		 gaux);
	  printf("Norm of preconditioned residual %15.12f, relative decrease %15.12f\n",eaux,daux);
	  printf("Current ratio (norm of preconditioned residual)/(norm of rhs) = %15.12f\n",faux);
	}
  }

  if(mkb_core_level->Monitor>=LSC_INFO){
	printf("\nCoarse grid correction problem Solved for level %d\n", Level_id);
  }


/*kbw
printf("mesh level %d, nrdof %d, precon %d\n",
Level_id,mkb_core_level->Nrdofgl,mkb_core_level->Precon);
printf("rhs:\n");
for(i=0;i<mkb_core_level->Nrdofgl;i++) printf("%10.6lf",B[i]);
printf("\nsolution:\n");
for(i=0;i<mkb_core_level->Nrdofgl;i++) printf("%10.6lf",X[i]);
printf("\n");
getchar();
/*kew*/
/* check whether B = A * X */
#ifdef DEBUG
  if(Level_id==0){

#ifdef PARALLEL
	/*||begin||*/
	lsr_mkb_core_fem_exchange_dofs(Solver_id,Subsystem_id,Level_id,X);
	/*||end||*/
#endif

	// just compute the residual
	double *vtemp;
	vtemp = (double *) malloc(mkb_core_level->Nrdofgl*sizeof(double));

	int matrix_id = mkb_core_level->SM_and_LV_id;
	lar_compute_residual ( matrix_id, 1, 0, mkb_core_level->Nrdofgl, X, B, vtemp);
	ierr=0;
	for(i=0; i<mkb_core_level->Nrdofgl; i++){
	  if(fabs(vtemp[i])>TOL){
	ierr=1;
	printf("Error in coarse solve: entry %d - residual %lf \n",
		   i,vtemp[i]);
	  }
	}
	if(!ierr) printf("Coarse problem solved correctly !\n");
	free(vtemp);
  }
  else{
	printf("Level_id != 0 in coarse_solve ! Exiting !\n");
	exit(-10);
  }


#endif

  return(1);
}


/*---------------------------------------------------------
lsr_mkb_core_gmres - basic restarted GMRES solver procedure
---------------------------------------------------------*/
int lsr_mkb_core_gmres(	/* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - no convergence */
	int Solver_id,  /* in: solver data structure indicator to be passed
						   to data structure dependent routines */
		int Subsystem_id,  /* in: subsystem data structure ID (passed) */
	int Level_id,	/* in: index of the mesh (level) (passed) */
	int ndof, 	/* in: 	the number of degrees of freedom */
	int ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* x, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* b,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int krylbas, 	/* in: 	number of Krylov space vectors */
	int*  iter, 	/* in:	the maximum GMRES iterations to be performed */
			/* out:	actual number of iterations performed */
	double* resid, 	/* in:	tolerance level for residual */
			/* out:	the final value of residual */
	double* res_rel, /* in:	tolerance level for relative residual */
			/* out:	the final value of relative residual */
	double* res_rhs, /* in:	tolerance level for ratio rhs/residual */
			/* out:	the final value of the ratio rhs/residual */
	int monitor,	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data 3 - iteration data */
	double* pconvr	/* out: pointer to convergence rate */
	)
{


/* local variables (for meaning of symbols see "Templates..." */
  double* v; 	/* workspace for Krylov space basis vectors v */
  double* h; 	/* upper Hessenberg matrix h */
		/* (first column for vectors s and y) */
  double* giv;	/* stored Givens rotations */
  double convrate;/* convergence rate */
  double rnrm2;	/* L2 norm of the residual vector */
  double bnrm2;	/* L2 norm of the right hand side vector (or just 1.0) */
  double beta0;	/* L2 norm of the initial residual vector */
  double beta1; 	/* to compute convergence in last restart */
  double tol_res,tol_rel,tol_rhs;	/* tolerances for error measures */
  int restrt; 	/* number of restarts */
  int restit = 0;	/* counter for restarts */
  int gmiter = 0;	/* counter for GMRES iterations */
  int ikryl = 0;	/* counter for Krylov space base vectors */

/* auxiliary variables */
  int i,k,iaux,jaux;
  double daux,eaux,faux;

/* constants */
  int IONE = 1;
  int IZERO = 0;
  double DONE = 1.0;
  const double MAX_CONVRATE = 0.999;


/*++++++++++++++++ executable statements ++++++++++++++++*/

/* rewrite input parameters */
  tol_res = *resid; if(tol_res<1.e-12) tol_res=1.e-12;
  tol_rel = *res_rel; if(tol_rel<1.e-12) tol_rel=1.e-12;
  tol_rhs = *res_rhs;
  if(tol_rhs > 0) {
	printf("RHS based convergence indicator not yet supported in parallel versionof solver ! \nResetting to 0.0\n");
	//exit(-1);
	tol_rhs = 0.0;
  }

/* write parameters */
  if(monitor>=2){
	printf("\nGMRES: entered with tolerances for: absolute norm of residual %14.12f\n",
	   tol_res);
	printf("       (relative to RHS %14.12f, relative to initial %14.12f)\n",
	   tol_rhs,tol_rel);
  }

/* allocate and initialize arrays */
  v = lsr_mkb_core_util_dvector(ndof*(krylbas+1),"v in gmres");
  h = lsr_mkb_core_util_dvector((krylbas+1)*(krylbas+1),"h in gmres");
  giv = lsr_mkb_core_util_dvector(2*(krylbas+1),"giv in gmres");

/* compute number of restarts */
  restrt = (*iter) / krylbas + 1 ;

  double current_residual = 0;
/* restarts' loop */
  do{

/* compute initial preconditioned(!) (and not preconditioned) residual (iaux=1) */
	iaux=1; /* take into account right hand side vector */
	if(gmiter==0) jaux=ini_zero; /* for first restart take input vector x */
	else jaux=0; /* for subsequent restarts assume initial guess is not zero */

	//get real residual at the beginning of each restart
	lsr_mkb_core_compres(Solver_id,Subsystem_id,Level_id,iaux,jaux,ndof,x,b,v);
	double previous_residual = current_residual;
	current_residual = lsr_mkb_core_fem_vec_norm(Solver_id,Subsystem_id,Level_id,ndof, v);


	lsr_mkb_core_compreres(Solver_id,Subsystem_id,Level_id,iaux,jaux,ndof,x,b,v);
	/*||begin||*/
	rnrm2 = lsr_mkb_core_fem_vec_norm(Solver_id,Subsystem_id,Level_id,ndof, v);


/*kbw
	printf("rnrm2 = %lf\n",rnrm2);
/*kew*/
	//    rnrm2 = dnrm2_(&ndof, v, &ione);

	bnrm2=0.0;
	if(ini_zero&&gmiter==0) {
	  bnrm2 = rnrm2;
	}
	if (bnrm2 < SMALL*SMALL) {
	  bnrm2 = 1.;
	}

/* compute and print convergence rate for the restart */
	if(restit==0) {
	  convrate = 0.0;
	  faux = 0.0;
	  beta0 = rnrm2;
	  beta1 = beta0;
	}
	else  if(rnrm2 < SMALL*SMALL || gmiter==0 ||
			 beta0 < SMALL*SMALL || beta1 < SMALL*SMALL  ) {
	  convrate = 0.0;
	  faux = 0.0;
	  beta1 = 0.0;
	}
	else{
	  convrate = exp(log(rnrm2/beta0)/(gmiter));
	  faux = exp(log(rnrm2/beta1)/(krylbas));
	  beta1 = rnrm2;
	}
	if(monitor>=2){
	  if(restit==0){
	printf("\nGMRES: norm of initial residual %10.8f (relative to RHS %10.8f)\n",
		   rnrm2,rnrm2/bnrm2);
	printf("       norm of RHS vector %14.12f\n",
		   bnrm2);
	  }
	  else{
	printf("\nGMRES: after restart %d (%d vectors), norm of preconditioned residual %14.12f\n",
		   restit,krylbas,rnrm2);
	printf("       (relative to RHS %14.12f, relative to initial %14.12f)\n",
		   rnrm2/bnrm2,rnrm2/beta0);
	printf("       total convergence rate %14.12f, in last restart %14.12f\n",
		   convrate,faux);
	  }
	}

/* check for convergence*/
	*resid = rnrm2;
	if(beta0!=0.0) *res_rel = rnrm2/beta0; else *res_rel = 0;
	if(bnrm2!=0.0) *res_rhs = rnrm2/bnrm2; else *res_rhs = 0;
	if ( rnrm2 < tol_res  ||  *res_rel < tol_rel
	 ||  *res_rhs < tol_rhs || rnrm2 < SMALL*SMALL )  goto convergence;
	if( convrate > MAX_CONVRATE || faux > MAX_CONVRATE){
	  printf("\nConvergence rate above %f - abandoning GMRES real residual decrease in this restart %lf \n",MAX_CONVRATE,
			  exp(log(current_residual/previous_residual)/(gmiter-(restit-1)*krylbas)));
	  goto noconvergence;
	}

/* construct the first Krylov space basis vector (scaled initial residual) */
	daux = 1. / rnrm2;
	dscal_(&ndof, &daux, &v[0], &IONE);

/* initialize the elementary vector e1 scaled by norm of preconditioned residual */
/* and store it in the first column (index 0) of h */
	h[0] = rnrm2;
	for(i=1;i<=krylbas+1;i++) h[i] = 0;

/* Krylov space basis vectors' loop - GMRES iterations - Arnoldi process */
	for(ikryl=1;ikryl<=krylbas;ikryl++){

	  gmiter++;

/* compute preconditioned matrix-vector product M^-1 * A * v_i */
	  iaux=0; /* do not take into account right hand side vector */
	  jaux=0; /* assume initial guess is not zero */
	  lsr_mkb_core_compreres(Solver_id, Subsystem_id, Level_id, iaux, jaux, ndof,
		   &v[(ikryl-1)*ndof], b, &v[ikryl*ndof]);

/* create a column of matrix h */
	  k=(krylbas+1)*ikryl;
	  for(i=1;i<=ikryl;i++){
	/*||begin||*/
	h[k] = lsr_mkb_core_fem_sc_prod(Solver_id,Subsystem_id,Level_id,
				   ndof, &v[(i-1)*ndof], &v[ikryl*ndof]);
/*kbw
	printf("h[%d] = %lf\n",k,h[k]);
/*kew*/

	/*||end||*/
	// h[k] = ddot_(&ndof, &v[(i-1)*ndof], &IONE, &v[ikryl*ndof], &IONE);
	k++;
	  }

/* subtract the last column of v from previous columns */
	  k=(krylbas+1)*ikryl;
	  for(i=1;i<=ikryl;i++){
	daux = -h[k];
	daxpy_(&ndof, &daux, &v[(i-1)*ndof], &IONE, &v[ikryl*ndof], &IONE);
	k++;
	  }

	  k=(krylbas+1)*ikryl+ikryl;
	  /*||begin||*/
	  h[k] = lsr_mkb_core_fem_vec_norm(Solver_id,Subsystem_id,Level_id,
				  ndof, &v[ikryl*ndof]);
/*kbw
	printf("h[%d] = %lf\n",k,h[k]);
/*kew*/
	  /*||end||*/
	  // h[k] = dnrm2_(&ndof, &v[ikryl*ndof], &IONE);

/* normalize kolumn of v to get next Krylov space basis vector */
	  daux = 1. / h[k];
	  dscal_(&ndof, &daux, &v[ikryl*ndof], &IONE);

/* apply stored Givens rotations to the i-th column of h */
	  for(i=1;i<ikryl;i++){
	drot_(&IONE, &h[i-1+ikryl*(krylbas+1)], &IONE,
		  &h[i+ikryl*(krylbas+1)], &IONE,
		  &giv[2*(i-1)], &giv[2*(i-1)+1]);
	  }

/* construct the i-th rotation matrix and apply it to h */
/* so that h(i+1,i) = 0. */
	  daux=h[ikryl-1+ikryl*(krylbas+1)];
	  eaux=h[ikryl+ikryl*(krylbas+1)];
	  drotg_(&daux, &eaux,
		 &giv[2*(ikryl-1)], &giv[2*(ikryl-1)+1]);

	  drot_(&IONE,&h[ikryl-1+ikryl*(krylbas+1)], &IONE,
		&h[ikryl+ikryl*(krylbas+1)], &IONE,
		&giv[2*(ikryl-1)], &giv[2*(ikryl-1)+1]);

/* apply the i-th rotation matrix to [ s(i), s(i+1) ]' */
/* this gives an approximation of the residual norm */
	  drot_(&IONE, &h[ikryl-1], &IONE, &h[ikryl],
		&IONE, &giv[2*(ikryl-1)], &giv[2*(ikryl-1)+1]);

/* compute convergence indicator and convergence rate */
	  daux = *resid;
	  *resid = fabs(h[ikryl]);
	  if(beta0!=0.0) *res_rel = *resid/beta0; else *res_rel = 0;
	  if(bnrm2!=0.0) *res_rhs = *resid/bnrm2; else *res_rhs = 0;
	  if(*res_rel<SMALL*SMALL){
	convrate=0.0;
	  }
	  else{
	convrate=exp(log(*res_rel)/(gmiter));
	  }

/* monitor results */
	  if(monitor>=1&&(*resid)>daux){
	printf("\nGMRES warning: diverging residuals\n");
	  }
	  if(monitor>=3){
	printf("\nGMRES: iteration %d, norm of preconditioned residual %14.12f\n",
		   gmiter,*resid);
	printf("       (relative to RHS %14.12f, relative to initial %14.12f)\n",
		   *res_rhs,*res_rel);
	printf("       convergence rate: total %14.12f, in last iteration %14.12f\n",
		   convrate, *resid/daux);
	  }

/* if conditions are met to stop Arnoldi process */
	  if ( ikryl==krylbas || gmiter== *iter || *resid < tol_res  ||
	   *res_rel < tol_rel ||  *res_rhs < tol_rhs ) {

/* compute the vector of coefficients y (stored in first column of h) */
	iaux=krylbas+1;
	dtrsv_("U", "N", "N", &ikryl, &h[krylbas+1], &iaux, &h[0], &IONE);

/* monitor results */
	if(monitor>=4){
	  printf("GMRES: computed vector y:\n");
	  for(i=0;i<krylbas+1;i++) printf("%e  ",h[i]);
	  printf("\n");
	}

/* compute current solution vector x */
	dgemv_("N", &ndof, &ikryl, &DONE, &v[0], &ndof, &h[0], &IONE,
		   &DONE, &x[0], &IONE);

/* monitor results */
	if ( gmiter== *iter || *resid < tol_res  ||
		 *res_rel < tol_rel ||  *res_rhs < tol_rhs ) {
	  if(monitor>=2){
		faux = exp(log(*resid/beta1)/(gmiter-restit*krylbas));
		printf("\nGMRES: after restart %d (%d vectors), norm of preconditioned residual %14.12f\n",
		   restit+1,gmiter-restit*krylbas,*resid);
		printf("       (relative to RHS %14.12f, relative to initial %14.12f)\n",
		   *res_rhs,*res_rel);
		printf("       convergence rate: total %14.12f, in last restart %14.12f\n",
		   convrate,faux);
	  }
	}

	if ( *resid < tol_res  || *res_rel < tol_rel
		 || *res_rhs < tol_rhs ) goto convergence;
	else if(gmiter== *iter) goto noconvergence;
	else break;

	  }
/* ^the end of breaking procedures */

	}
/* ^the end of Arnoldi process */

	restit++;

  }while(restit<=restrt);
/* ^the end of restarts' loop */

  printf("Hi, you should not be here at all...\n");

convergence:{
	if(monitor>=2){
	  printf("\nGMRES convergence achieved in %d iterations\n\n",gmiter);
	}

	*iter=gmiter;
	*pconvr=convrate;

	free(v);
	free(h);
	free(giv);

	return(IONE);
  }

noconvergence:{
	if(monitor>=2){
	  printf("\nNo GMRES convergence after %d iterations\n\n",gmiter);
	}

	*iter=gmiter;
	*pconvr=convrate;

	free(v);
	free(h);
	free(giv);

	return(IZERO);
  }

}


/*---------------------------------------------------------
lsr_mkb_core_util_dvector - to allocate a double vector: name[0..ncom-1]:
				  name=lsr_mkb_core_util_dvector(ncom,error_text)
---------------------------------------------------------*/
double *lsr_mkb_core_util_dvector( /* return: pointer to allocated vector */
	int ncom,  	/* in: number of components */
	char error_text[]/* in: error text to be printed */
	)
{

  double *v;

  v = (double *) malloc (ncom*sizeof(double));
  if(!v){
	printf("Not enough space for allocating vector: %s ! Exiting\n", error_text);
	exit(1);
  }
  return v;
}


/*---------------------------------------------------------
lsr_mkb_core_util_ivector - to allocate an integer vector: name[0..ncom-1]:
				  name=lsr_mkb_core_util_ivector(ncom,error_text)
---------------------------------------------------------*/
int *lsr_mkb_core_util_ivector(    /* return: pointer to allocated vector */
	int ncom, 	/* in: number of components */
	char error_text[]/* in: error text to be printed */
	)
{

  int *v;

  v = (int *) malloc (ncom*sizeof(int));
  if(!v){
	printf("Not enough space for allocating vector: %s ! Exiting\n", error_text);
	exit(1);
  }
  return v;
}

/*---------------------------------------------------------
lsr_mkb_core_util_imatrix - to allocate an integer matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
---------------------------------------------------------*/
int **lsr_mkb_core_util_imatrix( /* returns: pointer to array of pointers to integers */
	int Nrow, 	/* in: number of rows */
	int Ncol, 	/* in: number of columns */
	char Error_text[]/* in: text to print in case of error */
	)
{

  int i;
  int **m;

  m = (int **) malloc (Nrow*sizeof(int *));
  if(!m){
	printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
	exit(1);
  }
  for(i=0;i<Nrow;i++){
	m[i] = (int *) malloc (Ncol*sizeof(int));
	if(!m[i]){
	  printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
	  exit(1);
	}
  }
  return m;
}

/*---------------------------------------------------------
lsr_mkb_core_util_dmatrix - to allocate a double matrix name[0..nrow-1][0..ncol-1]:
				  name=imatrix(nrow,ncol,error_text)
---------------------------------------------------------*/
double **lsr_mkb_core_util_dmatrix( /* returns: pointer to array of pointers to doubles */
	int Nrow, 	/* in: number of rows */
	int Ncol, 	/* in: number of columns */
	char Error_text[]/* in: text to print in case of error */
	)
{

  int i;
  double **m;

  m = (double **) malloc (Nrow*sizeof(double *));
  if(!m){
	printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
	exit(1);
  }
  for(i=0;i<Nrow;i++){
	m[i] = (double *) malloc (Ncol*sizeof(double));
	if(!m[i]){
	  printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
	  exit(1);
	}
  }
  return m;
}

/*---------------------------------------------------------
lsr_mkb_core_util_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int lsr_mkb_core_util_chk_list(	/* returns: */
			/* >0 - position on the list */
					/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	)
{

  int i, il;

  for(i=0;i<Ll;i++){
	if((il=List[i])==0) break; // 0 indicates the end of list
	/* found on the list on (i+1) position */
	if(Num==il) return(i+1);
  }
  /* not found on the list */
  return(0);
}

/*---------------------------------------------------------
lsr_mkb_core_util_put_list - to put Num on the list List with length Ll
	(filled with numbers and zeros at the end)
---------------------------------------------------------*/
int lsr_mkb_core_util_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
				/*   0 - put on the list */
				/*  -1 - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
	)
{

  int i, il;

  for(i=0;i<Ll;i++){
	if((il=List[i])==0) break;
	/* found on the list on (i+1) position */
	if(Num==il) return(i+1);
  }
  /* if list is full return error message */
  if(i==Ll) return(-1);
  /* update the list and return*/
  List[i]=Num;
  return(0);
}

/*---------------------------------------------------------
lsr_mkb_core_util_d_zero - to zero a double vector
---------------------------------------------------------*/
void lsr_mkb_core_util_d_zero(double *Vec, int Num)
{

  int i;

  for(i=0;i<Num;i++){
	Vec[i]=0;
  }

}

/*---------------------------------------------------------
lsr_mkb_core_util_i_zero - to zero an integer vector
---------------------------------------------------------*/
void lsr_mkb_core_util_i_zero(int *Vec, int Num)
{

  int i;

  for(i=0;i<Num;i++){
	Vec[i]=0;
  }

}


/*---------------------------------------------------------
lsr_mkb_core_util_sort - to heap-sort an array (code taken from fortran...)
---------------------------------------------------------*/
void lsr_mkb_core_util_sort(
   int    *Ind_array,    /* in/out: index array for sorting */
   double *Val_array     /* in: array of values used for sorting */
   )
{

  int i,j,l,ir,index;
  double q;

  l = Ind_array[0]/2+1;
  ir =  Ind_array[0];

  s20: {}

  if(l>1){
	l--;
	index=Ind_array[l];
	q=Val_array[index];
  }
  else{
	index=Ind_array[ir];
	q=Val_array[index];
	Ind_array[ir]=Ind_array[1];
	ir--;
	if(ir==1){
	  Ind_array[1]=index;
	  return;
	}
  }

  i=l;
  j=2*l;

  s30: {}

  if(j>ir) goto s40;
  if((j<ir)&&(Val_array[Ind_array[j]]<Val_array[Ind_array[j+1]])) j++;
  if(q<Val_array[Ind_array[j]]){
	Ind_array[i]=Ind_array[j];
	i=j;
	j*=2;;
  }
  else{
	j=ir+1;
  }

  goto s30;

  s40: {}

  Ind_array[i]=index;

  goto s20;

}


/************************************************************
lsr_mkb_core_util_dgetrf - quasi-LU decomposition of a matrix (only for small
  matrices that fit into cache)
*************************************************************/
void lsr_mkb_core_util_dgetrf(double* a, int m, int* ips)
/*
in:
	a - matrix to decompose
	m - number of rows
out:
	a - decomposed matrix
	ips - partial pivoting information storage
*/
{

  int i,j,k,mm1,kp1,ipiv,kpiv,ipivot;
  double temp,pivot,big;

  for(i=0;i<m;i++){ips[i]=i;}

  mm1=m-1;
  for(k=0;k<mm1;k++){
	big=0.0;
	for(i=k;i<m;i++){
	  ipiv=ips[i];
	  if(big<=fabs(a[ipiv+m*k])){
	big=fabs(a[ipiv+m*k]);
	ipivot=i;
	  }
	}
	if(big<SMALL*SMALL){
	  printf("Zero pivot 1 in LU decomposition (%.15lf)\n",big);
	  exit(1);
	}
	kpiv=ips[ipivot];
	if(ipivot!=k){
	  ips[ipivot]=ips[k];
	  ips[k]=kpiv;
	}
	pivot=a[kpiv+m*k];
	kp1=k+1;
	for(i=kp1;i<m;i++){
	  ipiv=ips[i];
	  temp=a[ipiv+m*k]/pivot;
	  a[ipiv+m*k]=temp;
	  for(j=kp1;j<m;j++) a[ipiv+m*j] -= temp*a[kpiv+m*j];
	}
  }
  if(fabs(a[ips[mm1]+m*mm1])<SMALL*SMALL){
	printf("Zero pivot 2 in LU decomposition (%.15lf)\n",
	   a[ips[mm1]+m*mm1]);
	exit(1);
  }
}

/************************************************************
lsr_mkb_core_util_dgetrs - to perform forward reduction and back substitution
	of the RHS vector for solving a system of linear equations
************************************************************/
void lsr_mkb_core_util_dgetrs(double* a, int m, double* b, double* x, int* ips)
/*
in:
	a - LU decomposed matrix
	m - number of rows
	b - right hand side vector
out:
	x - solution of a system of linear equations
	ips - partial pivoting information storage
*/

{
  int i,j,ipiv;
  double sum;

  x[0]=b[ips[0]];
  for(i=1;i<m;i++){
	ipiv=ips[i];
	sum=b[ipiv];
	for(j=0;j<i;j++){ sum -= a[ipiv+m*j]*x[j]; }
	x[i]=sum;
  }
  x[m-1]=x[m-1]/a[ipiv+m*(m-1)];
  for(i=m-2;i>=0;i--){
	ipiv=ips[i];
	sum=x[i];
	for(j=i+1;j<m;j++){ sum -= a[ipiv+m*j]*x[j]; }
	x[i]= sum / a[ipiv+m*i];
  }
}


/*---------------------------------------------------------
lsr_mkb_core_util_skip_rest_of_line - to allow for comments in input files
---------------------------------------------------------*/
void lsr_mkb_core_util_skip_rest_of_line(
			  FILE *Fp  /* in: input file */
			  )
{
  while(fgetc(Fp)!='\n'){}
  return;
}

///////////////////////////// GEOMETRIC MULTIGRID UTILITY //////////////////////

/*---------------------------------------------------------
  lsr_mkb_core_get_pdeg_coarse - to get enforced pdeg for the coarse mesh
---------------------------------------------------------*/
int lsr_mkb_core_get_pdeg_coarse( // returns: enforced pdeg for the coarse mesh
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  int Level_id /* in: level number */
  )
{

  lst_mkb_core_levels *mkb_core_level =
	&lsv_mkb_core_solver[Solver_id].subsystem[0].level[Level_id];

  return(mkb_core_level->Pdeg);

}


