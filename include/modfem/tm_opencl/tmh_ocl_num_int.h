/************************************************************************
File tmh_ocl_num_int.h - OpenCL interface to numerical integration and solution procedures

Contains declarations of constants, types and interface routines:

  tmr_ocl_create_num_int_kernel - for a selected platform, device and kernel index
  tmr_ocl_create_assembly_structures - create all necessary for numerical integration 
      and assembling data structures in host and accelerator memory
  tmr_ocl_free_assembly_structures - free numerical integration 
      and assembling data structures created by tmr_ocl_create_assembly_structures
  tmr_ocl_create_assemble_stiff_mat_elem - create local element stiffness matrices
        and possibly assemble them to the global stiffness matrix
  tmr_ocl_create_solver_structures_crs - to create CRS structure on ACCELERATOR
  tmr_ocl_free_solver_structures_crs - to release CRS structure on ACCELERATOR

History:        
	02.2013 - Krzysztof Banas, initial version
        10.2016 - Jan Bielanski, Kazimierz Chlon - assembling into crs on gpu

------------------------------  			
History:        
	02.2013 - Krzysztof Banas, initial version

*************************************************************************/

#ifndef TMH_OPENCL_NUM_INT_H
#define TMH_OPENCL_NUM_INT_H

#include<stdlib.h>
#include<stdio.h>

#include <CL/cl.h>

/* Constants */
// float versus double (MUST BE COMPATIBLE WITH KERNEL SWITCH!!!!)
// data type for integration and assembly
//#define SCALAR float
#define SCALAR double


#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  // Classic tables
  double *crs_val_cpu;
  double *crs_val_gpu;
  unsigned int *crs_col_ind;
  unsigned int *crs_row_ptr;
  double *rhs_val_cpu;
  double *rhs_val_gpu;

  int Nnz; //Non zero values in crs
  int Nrdof; //Global number of DOFs
  
  // OpenCL CRS memory objects
  cl_mem ocl_crs_val;
  cl_mem ocl_crs_external_val;
  cl_mem ocl_crs_col_ind;
  cl_mem ocl_crs_row_ptr;
  
  // CRS size of data
  int ocl_crs_val_bytes;
  int ocl_crs_col_ind_bytes;
  int ocl_crs_row_ptr_bytes;

  // OpenCL RHS vector
  cl_mem ocl_rhs_val;
  cl_mem ocl_rhs_external_val;
    
  // RHS size of data
  int ocl_rhs_bytes;
  
  
} tmt_ocl_crs_struct;

typedef struct {

  // OpenCL assembly table memory objects
  cl_mem ocl_asse_pos_first_dof_int_ent; 
  cl_mem ocl_assembly_table;
  
  int nr_asse_blocks_all_int_ent;
  int ocl_asse_pos_first_dof_int_ent_bytes;
  int ocl_assembly_table_bytes;
  
  cl_mem ocl_local_to_global;
  cl_mem ocl_global_to_posglob;
  cl_mem ocl_pos_first_dof_int_ent;
  int ocl_local_to_global_bytes;

} tmt_ocl_assembly_struct;

typedef struct {

  // Geometric data - stored by elements
  int nr_geo; int ocl_el_geo_data_in_size;
  int ocl_el_geo_data_in_bytes;
  cl_mem ocl_el_geo_data_in; 

  // Geometric data - strored by dofs id
  int ocl_el_geo_data_by_dofs_in_size;
  int ocl_el_geo_data_by_dofs_in_bytes;
  cl_mem ocl_el_geo_data_by_dofs_in; 
  
  // Gauss data
  int ngauss; int ocl_gauss_data_in_size;
  int ocl_gauss_data_in_bytes;
  cl_mem ocl_gauss_data_in;
  
  // Shape fun data
  int nshape; int ocl_shape_fun_data_in_size;
  int ocl_shape_fun_data_in_bytes;
  cl_mem ocl_shape_fun_data_in;
  
  //PDE data
  int all_el_pde_coeff_size;
  int one_el_pde_coeff_size;
  int one_int_p_pde_coeff_size;
  int ocl_el_pde_dat_in_bytes;
  cl_mem ocl_el_pde_data_in;
  
    // CPU only:
  int* ngeo_color_pos;
  
} tmt_ocl_num_int_struct;

typedef struct {

  // OpenCL solution table memory objects
  cl_mem ocl_dofs_vector_current; 
  cl_mem ocl_dofs_vector_prev_iter;
  cl_mem ocl_dofs_vector_prev_step;
  
  int nrdofs;

} tmt_ocl_solution_struct;



typedef struct {

  tmt_ocl_crs_struct crs_struct;
  tmt_ocl_assembly_struct assembly_struct;
  tmt_ocl_num_int_struct num_int_struct;
  tmt_ocl_solution_struct solution_struct;

} tmt_ocl_problem_struct;

/* A single (for the time being) global structure for solving a single problem */
extern tmt_ocl_problem_struct tmv_ocl_problem_struct;

/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
/* -------------- specific procedures for numerical integration  ------------------ */

/*---------------------------------------------------------
  tmr_ocl_create_num_int_kernel - for a selected platform, device and kernel index
---------------------------------------------------------*/
extern int tmr_ocl_create_num_int_kernel(
  int Platform_index,
  int Device_index,
  int Kernel_index,
  char* Kernel_name,
  const char* Kernel_file,
  int Monitor
    );


/*---------------------------------------------------------
  tmr_ocl_create_assembly_structures - create all necessary for numerical integration 
      and assembling data structures in host and accelerator memory
----------------------------------------------------------*/
extern int tmr_ocl_create_assembly_structures(
				       );
      
/*---------------------------------------------------------
  tmr_ocl_free_assembly_structures - free numerical integration 
      and assembling data structures created by tmr_ocl_create_assembly_structures
----------------------------------------------------------*/
extern int tmr_ocl_free_assembly_structures(
				     );

/*---------------------------------------------------------
  tmr_ocl_create_assemble_stiff_mat_elem - create local element stiffness matrices
        and possibly assemble them to the global stiffness matrix
----------------------------------------------------------*/
extern int tmr_ocl_create_assemble_stiff_mat_elem(
  int Problem_id, 
  //int Level_id, 
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  
  //int Nr_int_ent,
  //int* L_int_ent_type,
  //int* L_int_ent_id,
  
  //int Nr_colors_elems,
  //int* L_color_index_elems,
  
  //int Nr_colors_accel,
  //int* Asse_pos_first_dof_int_ent,
  //int* Assembly_table,
  //int* Pos_first_dof_int_ent,
  
  //int* Local_to_global,
  //int Max_dofs_int_ent

  // Number of elements per color
  int Nr_elems,

  // Index of first element in assembly table for color
  int First_asse_elem_index,

  // Index of first geo for color
  int First_geo_index,

  // Index of first pde for color
  int First_pde_index
					   );
/*---------------------------------------------------------
  tmr_ocl_create_solver_structures_crs - to create CRS structure on ACCELERATOR
----------------------------------------------------------*/
extern int tmr_ocl_create_solver_structures_crs(
					 );

  
/*---------------------------------------------------------
  tmr_ocl_free_solver_structures_crs - to release CRS structure on ACCELERATOR
----------------------------------------------------------*/
extern int tmr_ocl_free_solver_structures_crs(
				       );

/*---------------------------------------------------------
tmr_ocl_get_solver_structures_crs - to get CRS structure on ACCELERATOR
/*---------------------------------------------------------*/
extern int tmr_ocl_get_solver_structures_crs(
  int* nnz,
  cl_mem* crs_col_ind,
  cl_mem* crs_row_ptr,
  cl_mem* crs_val,
  cl_mem* rhs );


//tmr_ocl_send_data_to_gpu_crs ??????????

//tmr_ocl_prepare_final_crs ????????????? tmr_ocl_merge_solver_structures

//tmr_ocl_solve_viennacl_gpu ????????????

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!! EVERYTHING BELOW ARE OLD OBSOLETE VERSIONS !!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  
/* /\**-------------------------------------------------------- */
/*   tmr_ocl_num_int_cleanup - discards created OpenCL resources */
/* ---------------------------------------------------------*\/ */
/* extern void tmr_ocl_num_int_cleanup(); */



/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
/* -------------- TEMPORARY jb_kc_opencl version ------------------ */

/* /\*--------------------------------------------------------- */
/*   tmr_perform_creation_crs - create crs structure on gpu */
/* ---------------------------------------------------------*\/ */
/* int tmr_perform_creation_crs( */
/*  int Problem_id, */
/*  int Nr_int_ent, */
/*  int* Asse_pos_first_dof_int_ent, */
/*  int* Assembly_table, */
/*  int Max_dofs_int_ent */
/*  //,tmt_ocl_crs_struct *tmv_ocl_crs_struct - currently structure is global defined */
/* ); */


/* /\*--------------------------------------------------------- */
/*   tmr_send_data_to_gpu_crs - send  */
/* ---------------------------------------------------------*\/ */
/* int tmr_prepare_final_crs(); */

/* /\*--------------------------------------------------------- */
/*   tmr_cleanup_crs - free tmr crs structure */
/* ---------------------------------------------------------*\/ */
/* int tmr_cleanup_crs(); */

/* /\*--------------------------------------------------------- */
/*   tmr_get_crs_structure - get crs structure from lsd */
/* ---------------------------------------------------------*\/ */
/* int tmr_get_crs_structure( */
/*  int Nrdof_glob, */
/*  int nnz, */
/*  int* crs_col_ind, */
/*  int* crs_row_ptr, */
/*  double* crs_val, */
/*  double* rhs */
/* ); */

/* /\*--------------------------------------------------------- */
/*   tmr_solve_viennacl_gpu - solve on gpu */
/* ---------------------------------------------------------*\/ */
/* int tmr_solve_viennacl_gpu(double* X); */





/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
/* -------------- TEMPORARY PORTED FROM include/tmh_accel_intf.h ------------------ */


/* /\*--------------------------------------------------------- */
/*   tmr_send_geo_data_accel_test - send geo data to accelerator */
/* ---------------------------------------------------------*\/ */
/* int tmr_send_geo_data_accel_test( */
/*                      int Problem_id, */
/*                      int Nr_geo, */
/*                      SCALAR* El_geo_dat_host_p // switch to scalar in new code */
/* ); */


/* /\*--------------------------------------------------------- */
/*   tmr_cleanup_geo_data_accel_test - remove geo data from accelerator */
/* ---------------------------------------------------------*\/ */
/* int tmr_cleanup_geo_data_accel_test( */
/*                      int Problem_id, */
/*                      int Nr_geo, */
/*                      SCALAR* El_geo_dat_host_p */
/* ); */

/* /\*--------------------------------------------------------- */
/*   tmr_rewrite_geo_by_color - rewrite geo from dof by dof organization */
/*                     to elem by elem organization */
/* ---------------------------------------------------------*\/ */
/* int tmr_rewrite_geo_by_color_test( */
/* 			 int Problem_id, */
/* 			 int Nr_elems, */
/* 			 int First_elem_in_color_index, */
/* 			 int First_geo_in_color_index			  */
/* ); */
		    


/* -------------- TEMPORARY PORTED FROM include/tmh_accel_intf.h ------------------ */
/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

#ifdef __cplusplus
}
#endif


#endif
