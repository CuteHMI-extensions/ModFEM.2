/***********************************************************************
File uts_accel_intf - utility routines supporting streaming processing
					  on accelerators (common to all problem modules)

Contains definitions of routines:

local:
  utr_create_assembly_structures_accel - to create data structurs in ACCELERATOR memory
  utr_free_assembly_structures_accel - to free data structurs in ACCELERATOR memory
  utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
						   and assemble them to the global SM using ACCELERATOR
  utr_merge_solver_structures_accel to merge parts of SM and LV created on CPU and GPU

------------------------------
History:
	08.2015 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
		10.2016 - Jan Bielanski, Kazimierz Chlon - assembling into crs on gpu
*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <signal.h>

#ifdef _OPENMP
#include <omp.h>
#endif


/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>

/* interface for all approximation modules */
#include <modfem/aph_intf.h>

/* interface for general purpose utilities - for all problem dependent modules*/
#include <modfem/uth_intf.h>
// interface for accelerator utilities - for all problem dependent modules*/
#include "uth_accel_intf.h"

/* interface for linear algebra packages */
#include <modfem/lin_alg_intf.h>

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>

#include <modfem/pdh_intf.h>

#include <modfem/tmh_intf.h>


// general OpenCL API
#include <CL/cl.h>

// interface of opencl implementation
#include <modfem/tm_opencl/tmh_ocl.h>
#include <modfem/tm_opencl/tmh_ocl_num_int.h>

#define TIME_TEST_ACCEL_INT_ASSE
#ifdef TIME_TEST_ACCEL_INT_ASSE
double t_begin;
double t_end;
double total_time;
#endif


/*------------------------------------------------------------
  utr_create_assembly_structures_accel - to create data structurs in ACCELERATOR memory
	 that are further used in numerical integration and assembly (no solver structures)
------------------------------------------------------------*/
int utr_create_assembly_structures_accel(
		const int Problem_id,
		const int nr_int_ent,
		const int nr_dof_ent,
		const int max_dofs_int_ent,
		const int nrdofs_glob,
		const int block_size,
		double * geo_dofs_vector,
		double * dofs_vector_current,
		double * dofs_vector_prev_iter,
		double * dofs_vector_prev_step,
		int nr_dof_blocks_all_int_ent,
		int * pos_first_dof_int_ent,
		int * local_to_global,
		int * global_to_posglob,
		const int nr_asse_blocks_all_int_ent,
		int * asse_pos_first_dof_int_ent,
		int * assembly_table,
		const int nr_colors_elems,
		int * l_color_index_elems,
		int * l_int_ent_id
)
{


	// ----------- PREPARE GAUSS/SHAPE_FUN/GEO DATA FOR ALL COLORS -----------

	int ngauss;
	SCALAR * gauss_dat_all_color_host;

	int nshape;
	SCALAR * shape_fun_dat_all_color_host;

	int ngeo;
	int * ngeo_color_p;
	SCALAR * el_geo_dat_all_color_host;

	int all_el_pde_coeff_size;
	int one_el_pde_coeff_size;
	int one_int_p_pde_coeff_size; // for debugging ONLY
	SCALAR * el_pde_dat_host;


	// Read GAUSS/SHAPE_FUN/GEO data

	utr_prepare_geo_dat_for_all_color_streaming (
			Problem_id,
			nr_colors_elems,
			l_color_index_elems,
			l_int_ent_id,
			&ngauss,
			&gauss_dat_all_color_host,
			&nshape,
			&shape_fun_dat_all_color_host,
			&ngeo,
			&ngeo_color_p,
			&el_geo_dat_all_color_host,
			&all_el_pde_coeff_size,
			&one_el_pde_coeff_size,
			&one_int_p_pde_coeff_size,
			&el_pde_dat_host
	);

	int nr_elems_all_gpu_colors = l_color_index_elems[nr_colors_elems] -
			l_color_index_elems[0];


	// call OpenCL routine
#ifdef OPENCL_GPU
	tmr_ocl_create_assembly_structures(
			Problem_id,
			nr_int_ent,
			nr_dof_ent,
			max_dofs_int_ent,
			nrdofs_glob,
			block_size,
			geo_dofs_vector,
			dofs_vector_current,
			dofs_vector_prev_iter,
			dofs_vector_prev_step,
			nr_dof_blocks_all_int_ent,
			pos_first_dof_int_ent,
			local_to_global,
			global_to_posglob,
			nr_asse_blocks_all_int_ent,
			asse_pos_first_dof_int_ent,
			assembly_table,
			ngauss,
			gauss_dat_all_color_host,
			nshape,
			shape_fun_dat_all_color_host,
			ngeo,
			ngeo_color_p,
			el_geo_dat_all_color_host,
			nr_elems_all_gpu_colors,
			all_el_pde_coeff_size,
			one_el_pde_coeff_size,
			one_int_p_pde_coeff_size,
			el_pde_dat_host
	);
#endif

};


/*------------------------------------------------------------
  utr_free_assembly_structures_accel - to free data structurs in ACCELERATOR memory
------------------------------------------------------------*/
int utr_free_assembly_structures_accel(

)
{
	// call OpenCL routine
#ifdef OPENCL_GPU
	tmr_ocl_free_assembly_structures();
#endif

}


/*------------------------------------------------------------
  utr_get_solver_structures_crs_accel - to get CRS structure on ACCELERATOR
------------------------------------------------------------*/
int utr_get_solver_structures_crs_accel(
		int * nnz,
		cl_mem * crs_col_ind,
		cl_mem * crs_row_ptr,
		cl_mem * crs_val,
		cl_mem * rhs )
{
	// call OpenCL routine
#ifdef OPENCL_GPU
	tmr_ocl_get_solver_structures_crs(nnz, crs_col_ind, crs_row_ptr, crs_val, rhs);
#endif

}



/*------------------------------------------------------------
 utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
								 and assemble them to the global SM using ACCELERATOR
------------------------------------------------------------*/
int utr_create_assemble_stiff_mat_elem_accel(
		int Problem_id,
		int Level_id,
		int Comp_type,         /* in: indicator for the scope of computations: */
		//extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
		//extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
		//extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
		//extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
		int * Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
		int Nr_int_ent,
		int * L_int_ent_type,
		int * L_int_ent_id,
		int Nr_colors_elems,
		int * L_color_index_elems,
		int Nr_colors_accel,
		int * Asse_pos_first_dof_int_ent,
		int * Assembly_table,
		int * Pos_first_dof_int_ent,
		int * Local_to_global,
		int Max_dofs_int_ent
) {


	int icolor;
	int first_pde_index = 0;
	int first_asse_elem_index = 0;
	int nr_elems_this_color = 0;


	for (icolor = 0; icolor < Nr_colors_accel; icolor++) {

		nr_elems_this_color = L_color_index_elems[icolor + 1] - L_color_index_elems[icolor];
		first_asse_elem_index = L_color_index_elems[icolor];

		printf("Assembling for color %d\n", icolor);

		tmr_ocl_create_assemble_stiff_mat_elem(Problem_id, Comp_type, Pdeg_coarse_p,
				nr_elems_this_color, first_asse_elem_index,
				icolor, first_pde_index);

		first_pde_index += nr_elems_this_color;

	}



	return (0);

}

/*------------------------------------------------------------
utr_merge_solver_structures_accel to merge parts of SM and LV created on CPU and GPU
------------------------------------------------------------*/
int utr_merge_solver_structures_accel(

)
{

	tmr_ocl_merge_solver_structures();

	/*   tmr_prepare_final_crs(); */
	/* #ifndef VIENNACL_WITH_OPENCL */
	/*     tmr_cleanup_crs(); */
	/* #endif */

	return 0;
}

/*------------------------------------------------------------
 utr_prepare_geo_dat_for_all_color_streaming - to create table with geo data for all colors
------------------------------------------------------------*/
void  utr_prepare_geo_dat_for_all_color_streaming (
		int Problem_id,
		int Nr_colors_accel,
		int * L_color_index_elems,
		int * L_int_ent_id,

		int * Ngauss_p,
		SCALAR ** Gauss_dat_host_p,

		int * Nshape_p,
		SCALAR ** Shape_fun_dat_host_p,

		int * Ngeo_p,
		int ** Ngeo_color_p,
		SCALAR ** El_geo_dat_host_p,

		int * All_el_pde_coeff_size,
		int * One_el_pde_coeff_size,
		int * One_int_p_pde_coeff_size,
		SCALAR ** El_pde_dat_host_p
)
{
	int i, j;
	int icolor, ielem;

	int mesh_id = pdr_ctrl_i_params(Problem_id, 2);
	int field_id = pdr_ctrl_i_params(Problem_id, 3);

	int nr_elems_all_gpu_colors; // Number of elements for all colors
	int nr_elems_this_color; // Number of elements per color
	int * L_elem_id; // Element id

	SCALAR * gauss_dat_host;
	int ref_el_quadr_dat_size;

	SCALAR * shape_fun_host;
	int ref_el_shape_fun_size;
	int shape_fun_host_bytes;

	int num_geo_dofs;
	int one_el_geo_dat_size;
	int position_geo_dofs;
	SCALAR * el_geo_dat;
	int * el_geo_pos_color;


	// Number of elements for all colors
	nr_elems_all_gpu_colors = L_color_index_elems[Nr_colors_accel] - L_color_index_elems[0];

	// 1. ------------- QUADRATURE DATA AND JACOBIAN TERMS -------------

	L_elem_id = &L_int_ent_id[L_color_index_elems[0]];

	// get the size of quadrature data (we assume all elements are of the same type)
	int base = apr_get_base_type(field_id, L_elem_id[0]);
	int ngauss;            /* number of gauss points */
	double xg[300];   	 /* coordinates of gauss points in 3D */
	double wg[100];       /* gauss weights */

	// choose an example element for a given pdeg and color
	int el_id = L_elem_id[0];
	int pdeg = apr_get_el_pdeg(field_id, el_id, NULL);
	int num_shap = apr_get_el_pdeg_numshap(field_id, el_id, &pdeg);
	apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);

	// we may need quadrature data for the reference element
	// for each gauss point its coordinates and weight are sent for reference element
	ref_el_quadr_dat_size = ngauss * 4;

	gauss_dat_host = (SCALAR *) malloc(ngauss * 4 * sizeof(SCALAR));

	/*kbw*/
	printf("Allocated %d bytes for gauss_dat_host %lu\n",
			ngauss * 4 * sizeof(SCALAR), gauss_dat_host);
	/*kew*/


	// This is the version when we send only necessary variables but we do not care about sending
	for (i = 0; i < ngauss; i++) {
		gauss_dat_host[4 * i] = xg[3 * i];
		gauss_dat_host[4 * i + 1] = xg[3 * i + 1];
		gauss_dat_host[4 * i + 2] = xg[3 * i + 2];
		gauss_dat_host[4 * i + 3] = wg[i];
	}

	/*kbw
	  for(i=0;i<5;i++){

	  printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
	  //printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);

	  }
	  for(i=ngauss-5;i<ngauss;i++){

	  printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
	  //printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);

	  }
	//kew*/

	// 2. ------------- SHAPE FUNCTIONS ---------------

	// space for element shape functions' values and derivatives in global memory
	// we need all shape functions and their derivatives
	// at all integration points for the reference element
	ref_el_shape_fun_size = 4 * num_shap * ngauss; // for the reference element

	// !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
	// fill and send buffers with shape function values for the reference element
	shape_fun_host_bytes = ref_el_shape_fun_size * sizeof(SCALAR);
	shape_fun_host = (SCALAR *)malloc(shape_fun_host_bytes);

	/*kbw*/
	printf("Allocated %d bytes for shape_fun_host %lu\n",
			shape_fun_host_bytes, shape_fun_host);
	/*kew*/

	// shape functions and derivatives are computed here but used also later on
	double base_phi_ref[APC_MAXELVD];    /* basis functions */
	double base_dphix_ref[APC_MAXELVD];  /* x-derivatives of basis function */
	double base_dphiy_ref[APC_MAXELVD];  /* y-derivatives of basis function */
	double base_dphiz_ref[APC_MAXELVD];  /* z-derivatives of basis function */

	// we need all shape functions and their derivatives
	// at all integration points for the reference element
	int ki;
	for (ki = 0; ki < ngauss; ki++) {

		int temp = apr_shape_fun_3D(base, pdeg, &xg[3 * ki],
						base_phi_ref, base_dphix_ref, base_dphiy_ref, base_dphiz_ref);
		assert(temp == num_shap);

		for (i = 0; i < num_shap; i++) {

			shape_fun_host[ki * 4 * num_shap + 4 * i] = base_phi_ref[i];
			shape_fun_host[ki * 4 * num_shap + 4 * i + 1] = base_dphix_ref[i];
			shape_fun_host[ki * 4 * num_shap + 4 * i + 2] = base_dphiy_ref[i];
			shape_fun_host[ki * 4 * num_shap + 4 * i + 3] = base_dphiz_ref[i];
		}
	}

	/*kbw
		for(i=0;i<5;i++){

		printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
		//printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);

		}
		for(i=num_shap-5;i<num_shap;i++){

		printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
		//printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);

		}
	//kew*/

	// 3. ------------- GEO_DOFS ---------------

	// Create data structure for geo data
	L_elem_id = &L_int_ent_id[L_color_index_elems[0]];
	//num_geo_dofs = mmr_el_node_coor(mesh_id, L_elem_id[0], NULL, NULL);
	num_geo_dofs = apr_get_el_geo_dofs(field_id, L_elem_id[0], NULL, NULL, NULL);


	one_el_geo_dat_size = 3 * num_geo_dofs;

	el_geo_dat = (SCALAR *)malloc(3 * num_geo_dofs * nr_elems_all_gpu_colors * sizeof(SCALAR));
	el_geo_pos_color = (int *)malloc((Nr_colors_accel + 1) * sizeof(int));

	/*kbw*/
	printf("Allocated %d bytes for el_geo_dat_host %lu\n",
			3 * num_geo_dofs * nr_elems_all_gpu_colors * sizeof(SCALAR), el_geo_dat);
	/*kew*/

	position_geo_dofs = 0;

	// Loop over colors for accelerator
	for (icolor = 0; icolor < Nr_colors_accel; icolor++) {

		// Elements per color
		nr_elems_this_color = L_color_index_elems[icolor + 1] - L_color_index_elems[icolor];
		L_elem_id = &L_int_ent_id[L_color_index_elems[icolor]]; // Element id

		el_geo_pos_color[icolor] = position_geo_dofs;

		for (ielem = 0; ielem < nr_elems_this_color; ielem++) {

			// element ID
			int el_id = L_elem_id[ielem];

			// checking whether this element has the same data as assumed for this color
			assert( pdeg == apr_get_el_pdeg(field_id, el_id, NULL) );
			assert( num_shap == apr_get_el_pdeg_numshap(field_id, el_id, &pdeg) );
			assert( base == apr_get_base_type(field_id, el_id) );

			// element geo_dofs for computing Jacobian terms on GPU

			// IDs of element vertices (nodes) and their coordinates as geo_dofs
			double geo_dofs[3 * APC_MAX_GEO_DOFS]; /* coord of nodes of El */
			//mmr_el_node_coor(mesh_id, el_id, NULL, geo_dofs);
			num_geo_dofs = apr_get_el_geo_dofs(field_id, el_id, NULL, NULL, geo_dofs);


			for (i = 0; i < num_geo_dofs; i++) {
				el_geo_dat[position_geo_dofs] = geo_dofs[3 * i];
				el_geo_dat[position_geo_dofs + 1] = geo_dofs[3 * i + 1];
				el_geo_dat[position_geo_dofs + 2] = geo_dofs[3 * i + 2];
				position_geo_dofs += 3;
			}
		}

	}

	el_geo_pos_color[Nr_colors_accel] = position_geo_dofs;

	// prepare PDE dependent data for streaming assembly (if necessary?)
	// (depending on problem name we decide whether new coefficients are necessary)
	pdr_prepare_pde_coeff_streaming(Problem_id,
			nr_elems_all_gpu_colors,
			&L_int_ent_id[L_color_index_elems[0]],
			All_el_pde_coeff_size,
			One_el_pde_coeff_size,
			One_int_p_pde_coeff_size,
			El_pde_dat_host_p);

	/*jbw*/
	printf("Number of elements processed by GPU %d\n", nr_elems_all_gpu_colors);
	printf("Coefficient common for all elements: %d\n", (*All_el_pde_coeff_size));
	printf("Coefficient for one elements: %d\n", (*One_el_pde_coeff_size));
	printf("Coefficient for gauss point: %d\n", (*One_int_p_pde_coeff_size));
	/*jbw*/

	/*kbw*/
	printf("\nData structure sizes:\n");
	printf("\tQuadrature data for reference element %d\n",
			ref_el_quadr_dat_size);
	printf("\tShape functions and derivatives for reference element %d\n",
			ref_el_shape_fun_size);
	printf("\tGeo dofs for each element: %d\n", one_el_geo_dat_size);
	/*kew*/

	*Ngauss_p = ngauss;
	*Gauss_dat_host_p = gauss_dat_host;
	*Nshape_p = num_shap;
	*Shape_fun_dat_host_p = shape_fun_host;
	*Ngeo_p = num_geo_dofs;
	*Ngeo_color_p = el_geo_pos_color;
	*El_geo_dat_host_p = el_geo_dat;

}

