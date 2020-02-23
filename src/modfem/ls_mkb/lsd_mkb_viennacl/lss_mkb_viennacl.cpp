/*******viennacl_ext******************************************************

History:
  10.2016 - Kazimierz Chlon
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<assert.h>

// lsd_mkb interface - including interface for direct solvers
#include <modfem/ls_mkb/lsh_mkb_intf.h>

/* internal information for the solver module */
#include "./lsh_mkb_viennacl.h"

// API for linear algebra routines supporting MKB solver
#include <modfem/ls_mkb/lah_intf.h>
// SUPERLU interface needs access to internal data structures of SM
//#include "../lad_crs_generic/lah_crs_generic.h"
//#include "../lad_crs/lah_crs.h"


//#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/gmres.hpp"
//#include "viennacl/linalg/bicgstab.hpp"
//#include "viennacl/linalg/prod.hpp"
//#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
//#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/direct_solve.hpp"
//#include "viennacl/matrix.hpp"
//#include "viennacl/linalg/lu.hpp"
//#include "viennacl/misc/bandwidth_reduction.hpp"


#include "viennacl/linalg/direct_solve.hpp"
#include "viennacl/linalg/lu.hpp"

#include "uth_system.h"
#include<vector>

#ifdef VIENNACL_WITH_OPENCL

  #include "tmh_intf.h"

  // interface of opencl implementation
  #include "../../tmd_opencl/tmh_ocl.h"
  #include "../../tmd_opencl/tmh_ocl_num_int.h"

#endif


#ifdef __cplusplus
extern "C" {
#endif


/* GLOBAL VARIABLES */
//extern int lsv_mkb_superlu_cur_solver_id;   /* ID of the current solver */
lst_mkb_viennacl_solvers lsv_mkb_viennacl_solver[LSC_MAX_NUM_SOLV];  /* array of solvers */

/**-----------------------------------------------------------
   lsr_mkb_direct_init - to create a new solver instance, read its control
   parameters and initialize its data structure
------------------------------------------------------------*/
int lsr_mkb_direct_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  char* Filename,  /* in: name of the file with control parameters */
  int Monitoring_level /* Level of output, -1 for Filename */
  ){

/*kbw/
  printf("\nIn lsr_mkb_direct_init: Solver id %d, Filename %s, Monitor %d\n",
  Solver_id,  Filename, Monitoring_level);
/*kew*/

#ifdef PARALLEL
  //mf_log_warn("Using direct solver in distributed memory environment!");
  printf("MKB direct solvers not adapted to parallel (distributed memory) execution. Exiting!\n");
  exit(-1);
#endif


  printf("\nViennacl\n");

  //lsv_mkb_viennacl_cur_solver_id = Solver_id;
  lsv_mkb_viennacl_solver[Solver_id].solver_id = Solver_id;
  lsv_mkb_viennacl_solver[Solver_id].monitor = Monitoring_level;



 // ViennaCL with Opencl
#ifdef VIENNACL_WITH_OPENCL

  // ----------- Get platform configuration ----------------

  // choose the platform
  int platform_index = tmr_ocl_get_current_platform_index();

  int device_tmc_type = tmr_ocl_get_current_device_type();

  // choose device_index
  int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);

  // choose device_index_cl
  cl_device_id device_index_cl = tmr_ocl_select_device_id(platform_index, device_index);

  // choose the context
  cl_context context = tmr_ocl_select_context(platform_index, device_index);

  // choose the command queue
  cl_command_queue command_queue =
	tmr_ocl_select_command_queue(platform_index, device_index);

	// -----------End Get platform configuration ----------------


  //supply existing context 'context'
  //with one device and one queue to ViennaCL using id '0':
  viennacl::ocl::setup_context(0, context, device_index_cl, command_queue);

  //std::cout << viennacl::ocl::current_device().info() << std::endl;
  //printf("\n\n VCL context = %d \n\n",viennacl::ocl::current_context().handle().get());

#endif


  return(Solver_id);
}


  /**-----------------------------------------------------------
   lsr_mkb_direct_create - to create a new solver instance, read its control
   parameters and initialize its data structure
------------------------------------------------------------*/
int lsr_mkb_direct_create( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id,   /* in: solver ID (used to identify the subproblem) */
  char* Filename,  /* in: name of the file with control parameters */
  int Monitoring_level /* Level of output, -1 for Filename */
  ){

/*kbw/
  printf("\nIn lsr_mkb_direct_init: Solver id %d, Filename %s, Monitor %d\n",
  Solver_id,  Filename, Monitoring_level);
/*kew*/

#ifdef PARALLEL
  //mf_log_warn("Using direct solver in distributed memory environment!");
  printf("MKB direct solvers not adapted to parallel (distributed memory) execution. Exiting!\n");
  exit(-1);
#endif


  /* initialize data structure */
  lsv_mkb_viennacl_solver[Solver_id].crs_row=NULL;
  lsv_mkb_viennacl_solver[Solver_id].crs_col=NULL;
  lsv_mkb_viennacl_solver[Solver_id].crs_val=NULL;
  lsv_mkb_viennacl_solver[Solver_id].rhs=NULL;
  lsv_mkb_viennacl_solver[Solver_id].offset=0;
  lsv_mkb_viennacl_solver[Solver_id].free_flag=0;

  return(Solver_id);
}



typedef viennacl::compressed_matrix<double> ViennaMatrixType;


/**--------------------------------------------------------
  lsr_mkb_direct_solve - to solve a system of equations, given previously constructed
			 system matrix, preconditioner
---------------------------------------------------------*/
 int lsr_mkb_direct_solve( /* returns: convergence indicator: */
			/* 1 - convergence */
			/* 0 - noconvergence */
						/* <0 - error code */
	int Solver_id,  /* in: solver ID */
	int Comp_type,  /* in: indicator for the scope of computations: */
					/*   LSC_SOLVE - solve the system */
					/*   LSC_RESOLVE - resolve for the new right hand side */
	int Matrix_id,  /* in: matrix ID */
	int Ndof, 	/* in: 	the number of degrees of freedom */
	double* X, 	/* in: 	the initial guess */
			/* out:	the iterated solution */
		double* B,	/* in:  the rhs vector, if NULL take rhs */
			/*      from block data structure */
	int Monitor	/* in:	flag to determine monitoring level */
			/*	0 - silent run, 1 - warning messages */
			/*	2 - 1+restart data, 3 - 2+iteration data */
	)
{


  int i;
  int j,k;
  int nnz, offset;
  int *crs_row=NULL, *crs_col=NULL;
  double *crs_val=NULL, *rhs=NULL;


	//temporary vector for rhs and sol_vect
  std::vector<double> stl_vec(Ndof);

  viennacl::vector<double> vcl_sol_vect(Ndof);

  double time_solve=time_clock();


// ViennaCL with Opencl
#ifdef VIENNACL_WITH_OPENCL


  cl_mem ocl_crs_col_ind;
  cl_mem ocl_crs_row_ptr;
  cl_mem ocl_crs_val;
  cl_mem ocl_rhs;


  // get pointers to crs structures on ACCELERATOR
  tmr_ocl_get_solver_structures_crs(&nnz, &ocl_crs_col_ind,
					  &ocl_crs_row_ptr, &ocl_crs_val, &ocl_rhs);

  // Wrapp Opencl memory to ViennaCL //
  viennacl::vector<double> vcl_rhs_vect(ocl_rhs, Ndof);
  viennacl::vcl_size_t crs_nnz = nnz;
  viennacl::vcl_size_t crs_size = Ndof;

  ViennaMatrixType vcl_dov_matrix(ocl_crs_row_ptr,
								  ocl_crs_col_ind,
								  ocl_crs_val,
								  crs_size,
								  crs_size,
								  crs_nnz,
								  viennacl::ocl::current_context());



// ViennaCL with Openmp and sequential
#else


  // get data from lar module
  offset=lsv_mkb_viennacl_solver[Solver_id].offset;

  lsv_mkb_viennacl_solver[Solver_id].free_flag = lar_get_SM_and_LV_crs( Matrix_id, offset,
	  &(lsv_mkb_viennacl_solver[Solver_id].crs_row),
	  &(lsv_mkb_viennacl_solver[Solver_id].crs_col),
	  &(lsv_mkb_viennacl_solver[Solver_id].crs_val),
	  &(lsv_mkb_viennacl_solver[Solver_id].rhs));

  crs_row=lsv_mkb_viennacl_solver[Solver_id].crs_row;
  crs_col=lsv_mkb_viennacl_solver[Solver_id].crs_col;
  crs_val=lsv_mkb_viennacl_solver[Solver_id].crs_val;
  rhs=lsv_mkb_viennacl_solver[Solver_id].rhs;

  nnz=crs_row[Ndof];


  viennacl::compressed_matrix<double> vcl_dov_matrix(Ndof,Ndof,nnz);
  viennacl::vector<double> vcl_rhs_vect(Ndof);


   //set ViennaCL structures from CRS
  vcl_dov_matrix.set((void*) crs_row,(void*) crs_col, crs_val,Ndof,Ndof,nnz);


  //set ViennaCL rhs vector
  for (i = 0;i < Ndof; ++i){
	stl_vec[i]=rhs[i];
  }
  copy(stl_vec.begin(), stl_vec.end(), vcl_rhs_vect.begin());

#endif

  time_solve=time_clock()-time_solve;
  printf("\n ViennaCL: prepare: %lfs\n",time_solve);

  time_solve=time_clock();


// ilu preconditioner
  viennacl::linalg::block_ilu_precond< ViennaMatrixType, viennacl::linalg::ilu0_tag> vcl_precon(vcl_dov_matrix,viennacl::linalg::ilu0_tag());


  //viennacl::linalg::cg_tag custom(1e-10, 100);
  viennacl::linalg::gmres_tag custom(1e-11, 800, 50);
  vcl_sol_vect = viennacl::linalg::solve(vcl_dov_matrix, vcl_rhs_vect, custom,vcl_precon);

  time_solve=time_clock()-time_solve;

  printf("\n ViennaCL: No. of iters: %d. Est. error: %.12lf in %lfs\n",custom.iters(),custom.error(),time_solve);


  time_solve=time_clock();

// get solution from ViennaCL structures to std:vector
  copy(vcl_sol_vect.begin(), vcl_sol_vect.end(), stl_vec.begin());

// rewrite the solution
  for (i = 0;i < Ndof; ++i) {
	X[i]=stl_vec[i];
	//printf("wynik=%20.10lf \n",std_rhs_vect[i]);
  }


  time_solve=time_clock()-time_solve;
  printf("\n ViennaCL: rewrite solution: %lfs\n",time_solve);



 return(0);

}

  /**--------------------------------------------------------
  lsr_mkb_direct_free - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
int lsr_mkb_direct_free( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
			 )
{


  if(lsv_mkb_viennacl_solver[Solver_id].free_flag != 0){
	free(lsv_mkb_viennacl_solver[Solver_id].crs_row);
	free(lsv_mkb_viennacl_solver[Solver_id].crs_col);
	free(lsv_mkb_viennacl_solver[Solver_id].crs_val);
  }


  return(0);
}



/**--------------------------------------------------------
  lsr_mkb_direct_destroy - to destroy a particular instance of the direct solver
------------------------------------------------------------*/
int lsr_mkb_direct_destroy( /* returns: >0 - solver ID, <0 - error code */
  int Solver_id   /* in: solver ID (used to identify the subproblem) */
			 )
{

	// ViennaCL with Opencl
#ifdef VIENNACL_WITH_OPENCL

  //tmr_cleanup_crs();
  //tmr_ocl_cleanup();

#endif

  // currently nothing to do

  return(0);
}



#ifdef __cplusplus
}
#endif
