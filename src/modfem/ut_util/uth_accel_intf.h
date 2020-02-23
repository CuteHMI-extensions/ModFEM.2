/***********************************************************************
File uth_accel_intf - utility routines supporting streaming processing
					  on accelerators (common to all problem modules)

Contains declarations of routines:

  utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
						   and assemble them to the global SM using ACCELERATOR
  utr_create_assembly_structures_accel - to create data structurs in ACCELERATOR memory
	 that are further used in numerical integration and assembly (no solver structures)

  ///////////&&&&&&&&&&&&$$$$$$$$$$$$ OBSOLETE $$$$$$$$$$$$$&&&&&&&&&&&&&&&//////////
  utr_prepare_create_assemble_accel
  utr_perform_create_assemble_accel_for_color
  utr_create_solver_crs_accel - create crs structure on accelerator
  utr_get_solver_crs_structure_accel - get crs structure from lsd
  utr_prepare_final_solver_crs_accel
  utr_cleanup_solver_crs_accel - free crs structure in accelerator memory
  utr_send_geo_data_accel - send geo data to accelerator
  utr_cleanup_geo_data_accel - remove geo data from accelerator
  utr_send_dofs_data_accel - send dofs data to accelerator
  utr_cleanup_dofs_data_accel - remove dofs data from accelerator


------------------------------
History:
	08.2015 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
*************************************************************************/

#ifndef _uth_accel_intf_
#define _uth_accel_intf_

#include <modfem/tm_opencl/tmh_ocl.h>
#include <modfem/tm_opencl/tmh_ocl_num_int.h>

#ifdef __cplusplus
extern "C"{
#endif

/** Constants */



/*** FUNCTION DECLARATIONS - headers for external functions ***/



/*------------------------------------------------------------
 utr_create_assembly_structures_accel - to create data structurs in ACCELERATOR memory
	that are further used in numerical integration and assembly (no solver structures)
------------------------------------------------------------*/
extern int utr_create_assembly_structures_accel(
                        );

/*------------------------------------------------------------
  utr_free_assembly_structures_accel - to free data structurs in ACCELERATOR memory
------------------------------------------------------------*/
extern int utr_free_assembly_structures_accel(

                     );

/*------------------------------------------------------------
  utr_get_solver_structures_crs_accel - to free data structurs in ACCELERATOR memory
------------------------------------------------------------*/
extern int utr_get_solver_structures_crs_accel(
  int* nnz,
  cl_mem* crs_col_ind,
  cl_mem* crs_row_ptr,
  cl_mem* crs_val,
  cl_mem* rhs );

/*------------------------------------------------------------
 utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
								 and assemble them to the global SM using ACCELERATOR
------------------------------------------------------------*/
extern int utr_create_assemble_stiff_mat_elem_accel(
  int Problem_id,
  int Level_id,
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  int Nr_int_ent,
  int* L_int_ent_type,
  int* L_int_ent_id,
  int Nr_colors_elems,
  int* L_color_index_elems,
  int Nr_colors_accel,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int* Pos_first_dof_int_ent,
  int* Local_to_global,
  int Max_dofs_int_ent
                     );

/*------------------------------------------------------------
utr_merge_solver_structures_accel to merge parts of SM and LV created on CPU and GPU
------------------------------------------------------------*/
extern int utr_merge_solver_structures_accel(

                       );

 ///////////&&&&&&&&&&&&$$$$$$$$$$$$ OBSOLETE $$$$$$$$$$$$$&&&&&&&&&&&&&&&//////////


/*------------------------------------------------------------
 utr_prepare_geo_dat_for_all_color_streaming - to create table with geo data for all colors
------------------------------------------------------------*/
void  utr_prepare_geo_dat_for_all_color_streaming (
  int Problem_id,
  int Nr_colors_accel,
  int* L_color_index_elems,
  int* L_int_ent_id,

  int* Ngauss_p,
  SCALAR** Gauss_dat_host_p,

  int* Nshape_p,
  SCALAR** Shape_fun_dat_host_p,

  int* Ngeo_p,
  int** Ngeo_color_p,
  SCALAR** El_geo_dat_host_p,

  int* All_el_pde_coeff_size,
  int* One_el_pde_coeff_size,
  int* One_int_p_pde_coeff_size,
  SCALAR** El_pde_dat_host_p
                           );

/*---------------------------------------------------------
  utr_perform_creation_crs - create crs structure on gpu
---------------------------------------------------------*/
int utr_perform_creation_crs(
 int Problem_id,
 int Nr_int_ent,
 int* Asse_pos_first_dof_int_ent,
 int* Assembly_table,
 int* Pos_first_dof_int_ent,
 int* Local_to_global,
 int Max_dofs_int_ent
 //,tmt_ocl_crs_struct *tmv_ocl_crs_struct - currently structure is global defined
);


/*---------------------------------------------------------
  utr_send_data_to_gpu_crs - send
---------------------------------------------------------*/
int utr_prepare_final_crs();

/*---------------------------------------------------------
  utr_cleanup_crs - free tmr crs structure
---------------------------------------------------------*/
int utr_cleanup_crs();


/*---------------------------------------------------------
  utr_solve_viennacl_gpu - solve on gpu
---------------------------------------------------------*/
int utr_solve_viennacl_gpu(double* X);



/*---------------------------------------------------------
  utr_send_geo_data_accel_test - send geo data to accelerator
---------------------------------------------------------*/
int utr_send_geo_data_accel_test(
                     int Problem_id,
                     int Nr_geo,
                     SCALAR* El_geo_dat_host_p // switch to scalar in new code
);


/*---------------------------------------------------------
  utr_cleanup_geo_data_accel_test - remove geo data from accelerator
---------------------------------------------------------*/
int utr_cleanup_geo_data_accel_test(
                     int Problem_id,
                     int Nr_geo,
                     SCALAR* El_geo_dat_host_p
);

/*---------------------------------------------------------
  utr_rewrite_geo_by_color - rewrite geo from dof by dof organization
					to elem by elem organization
---------------------------------------------------------*/
int utr_rewrite_geo_by_color_test(
             int Problem_id,
             int Nr_elems,
             int First_elem_in_color_index,
             int First_geo_in_color_index
);


/*---------------------------------------------------------
  utr_prepare_create_assemble_accel
----------------------------------------------------------*/
int utr_prepare_create_assemble_accel(
  int Problem_id,
  int Monitor
                      );

/*---------------------------------------------------------
  utr_perform_create_assemble_accel_for_color
----------------------------------------------------------*/
int utr_perform_create_assemble_accel_for_color(
  int Problem_id,
  int Level_id,
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_elems_this_color,
  int First_elem_in_color_index,
  int* L_elem_id_this_color,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int* Pos_first_dof_int_ent,
  int* Local_to_global,
  int Ngauss,
  SCALAR* Gauss_dat_host,
  int Num_shap,
  SCALAR* Shape_fun_dat_host,
  int Num_geo_dofs,
  int El_geo_dat_position_for_color,
  SCALAR* Geo_dofs_global,
  int Num_dofs,
  SCALAR* U_n_dofs_global,
  SCALAR* U_k_dofs_global,
  int Size_global_pde_data_int,  // size for all elements
  SCALAR* Global_pde_dat_host_int,
  int Size_global_pde_data_scalar,  // size for all elements
  SCALAR* Global_pde_dat_host_scalar,
  int Size_per_element_pde_data_int, // size for one element and all integration point
  int* El_pde_dat_host_int,
  int Size_per_element_pde_data_scalar,  // size for one element and one integration point
  SCALAR* El_pde_dat_host_scalar,
  int Max_dofs_int_ent
                      );


/*---------------------------------------------------------
  utr_create_solver_crs_accel - create crs structure on accelerator
---------------------------------------------------------*/
int utr_create_solver_crs_accel(
  int Problem_id,
  int Nr_int_ent,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int* Pos_first_dof_int_ent,
  int* Local_to_global,
  int Max_dofs_int_ent
  //,tmt_ocl_crs_struct *tmv_ocl_crs_struct - currently structure is global defined
);

/*---------------------------------------------------------
  utr_get_solver_crs_structure_accel - get crs structure from lsd
---------------------------------------------------------*/
int utr_get_solver_crs_structure_accel(
 int Nrdof_glob,
 int nnz,
 int* crs_col_ind,
 int* crs_row_ptr,
 double* crs_val,
 double* rhs
);

/*---------------------------------------------------------
  utr_prepare_final_solver_crs_accel
---------------------------------------------------------*/
int utr_prepare_final_solver_crs_accel();

/*---------------------------------------------------------
  utr_cleanup_solver_crs_accel - free crs structure in accelerator memory
---------------------------------------------------------*/
int utr_cleanup_solver_crs_accel();


/*---------------------------------------------------------
  utr_send_geo_data_accel - send geo data to accelerator
---------------------------------------------------------*/
int utr_send_geo_data_accel(
                     int Problem_id,
                     int Nr_elements,
                     int Nr_geo,
                     double* El_geo_dat_host_p
);

/*---------------------------------------------------------
  utr_cleanup_geo_data_accel - remove geo data from accelerator
---------------------------------------------------------*/
int utr_cleanup_geo_data_accel(
                     int Problem_id,
                     int Nr_elements,
                     int Nr_geo,
                     double* El_geo_dat_host_p
);


/*---------------------------------------------------------
  utr_send_dofs_data_accel - send dofs data to accelerator
---------------------------------------------------------*/
int utr_send_dofs_data_accel(
                     int Problem_id,
                     int Nr_dofs,
                     double* El_dofs_dat_host_p
);


/*---------------------------------------------------------
  utr_cleanup_dofs_data_accel - remove dofs data from accelerator
---------------------------------------------------------*/
int utr_cleanup_dofs_data_accel(
                     int Problem_id,
                     int Nr_dofs,
                     double* El_dofs_dat_host_p
);



#ifdef __cplusplus
}
#endif

#endif

