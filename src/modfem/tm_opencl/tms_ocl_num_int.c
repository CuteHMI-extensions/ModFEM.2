/***********************************************************************
File tms_ocl_num_int - OpenCl routines supporting streaming processing
					   on accelerators (common to all problem modules)

Contains definitions of routines:
  tmr_ocl_create_assembly_structures - create all necessary for numerical integration
	  and assembling data structures in host and accelerator memory
  tmr_ocl_create_assemble_stiff_mat_elem - create local element stiffness matrices
		and possibly assemble them to the global stiffness matrix
  tmr_ocl_free_assembly_structures - free numerical integration
	  and assembling data structures created by tmr_ocl_create_assembly_structures
  tmr_ocl_create_solver_structures_crs - to create CRS structure on ACCELERATOR
  tmr_ocl_free_solver_structures_crs - to release CRS structure on ACCELERATOR
  tmr_ocl_create_num_int_kernel - for a selected platform, device and kernel index

------------------------------
History:
	08.2016 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
		10.2016 - Jan Bielanski, Kazimierz Chlon - assembling into crs on gpu
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<signal.h>
#include<limits.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include <CL/cl.h>

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>

/* interface for all approximation modules */
#include <modfem/aph_intf.h>

/* interface for general purpose utilities - for all problem dependent modules*/
#include <modfem/uth_intf.h>
#include <modfem/uth_system.h>

/* debugging/logging/checking macros */
#include <modfem/uth_log.h>

/* interface for linear algebra packages */
#include <modfem/lin_alg_intf.h>

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>

#include <modfem/pdh_intf.h>

// generic thread management interface
#include <modfem/tmh_intf.h>

// interface of opencl implementation
#include <modfem/tm_opencl/tmh_ocl.h>
#include <modfem/tm_opencl/tmh_ocl_num_int.h>

/* A single global problem structure */
tmt_ocl_problem_struct tmv_ocl_problem_struct;

#ifdef TUNING
FILE * optf, *resf;
FILE * headf; //header file only for result titles
//#define COUNT_OPER
#endif

// SWITCH 0:
// Master switch: GPU versus PHI versus CPU - controlled by compilation options !!!!!

// SWITCH 1: Opencl_HSA is for Heterogenous System Architecture with Shared Virtual Memory eg. APU
#ifdef OPENCL_HSA
#define OPENCL_GPU
#endif

// SWITCH 2: one_el_one_thread strategy versus one_el_one_workgroup strategy
#define ONE_EL_ONE_THREAD
//#define ONE_EL_ONE_WORKGROUP
//#define ONE_EL_TWO_THREADS

// SWITCH 3: size for  work-group
#ifdef OPENCL_CPU
#define WORK_GROUP_SIZE 8
#elif defined OPENCL_GPU
#define WORK_GROUP_SIZE 64
#elif defined OPENCL_PHI
#define WORK_GROUP_SIZE 16
#else
#define WORK_GROUP_SIZE 0
#endif

// SWITCH 4: How many workgroups will be enough to make GPU cores (comp_units) busy?
#ifdef OPENCL_CPU
#define NR_WORK_GROUPS_PER_COMP_UNIT  1 // ????
#elif defined OPENCL_GPU
#define NR_WORK_GROUPS_PER_COMP_UNIT  4 // ????
#elif defined OPENCL_PHI
#define NR_WORK_GROUPS_PER_COMP_UNIT  1 // ????
#else
#define NR_WORK_GROUPS_PER_COMP_UNIT 0
#endif

// this should always be defined for GPUs
#ifdef OPENCL_GPU
#define COAL_WRITE
#endif

//#define TIME_TEST_OCL
#ifdef TIME_TEST_OCL
double ocl_time_total_kernel_execution = 0.0;
double ocl_time_auxilliary_operation = 0.0;
double ocl_time_total_transfer = 0.0;
double ocl_total_time = 0.0;
double t_begin, t_end;
#endif

// !!!!!!!!! For debug only !!!!!!!!!! //
//#define DEBUG_TMD
// !!!!!!!!! For debug only !!!!!!!!!! //

/*---------------------------------------------------------
  tmr_ocl_create_assembly_structures - create all necessary for numerical integration
	  and assembling data structures in host and accelerator memory
----------------------------------------------------------*/
int tmr_ocl_create_assembly_structures(
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
		const int ngauss,
		SCALAR * gauss_dat_all_color_host,
		const int nshape,
		SCALAR * shape_fun_dat_all_color_host,
		const int ngeo,
		int * ngeo_color_p,
		SCALAR * el_geo_dat_all_color_host,
		const int nr_elems,
		const int all_el_pde_coeff_size,
		const int one_el_pde_coeff_size,
		const int one_int_p_pde_coeff_size,
		SCALAR * el_pde_dat_host
)
{

	// ----------- Get platform configuration ----------------

	// choose the platform
	int platform_index = tmr_ocl_get_current_platform_index();

	int device_tmc_type = tmr_ocl_get_current_device_type();

	// choose device_index
	int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);

	// choose the context
	cl_context context = tmr_ocl_select_context(platform_index, device_index);

	// choose the command queue
	cl_command_queue command_queue =
			tmr_ocl_select_command_queue(platform_index, device_index);

	// -----------End Get platform configuration ----------------

	// ----------ALLOCATE MEMEORY FOR ASSEMBLY TABLE ------------

	// Calculate the block size
	int block_size_int_ent = (max_dofs_int_ent + 1) * block_size * max_dofs_int_ent * block_size;
	int size_assemble_table_col = block_size_int_ent * nr_int_ent; // compute size of assembly table (size of block for int ent * number of int ent)
	int global_to_posblob_bytes = (nr_dof_ent + 1) * sizeof(int);
	tmv_ocl_problem_struct.assembly_struct.nr_asse_blocks_all_int_ent = nr_asse_blocks_all_int_ent;
	tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table_bytes = size_assemble_table_col * sizeof(int);
	tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent_bytes = (nr_int_ent + 1) * sizeof(int);

	tmv_ocl_problem_struct.assembly_struct.ocl_local_to_global_bytes = nr_dof_blocks_all_int_ent * sizeof(int);

	// Create assembly table buffer on GPU
	tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent_bytes, NULL, NULL);
	tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table_bytes, NULL, NULL);

	tmv_ocl_problem_struct.assembly_struct.ocl_local_to_global =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.assembly_struct.ocl_local_to_global_bytes, NULL, NULL);
	tmv_ocl_problem_struct.assembly_struct.ocl_global_to_posglob =
			clCreateBuffer(context, CL_MEM_READ_ONLY, global_to_posblob_bytes, NULL, NULL);
	tmv_ocl_problem_struct.assembly_struct.ocl_pos_first_dof_int_ent =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent_bytes, NULL, NULL);


	// Write data on GPU
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent, CL_TRUE, 0,
			tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent_bytes, asse_pos_first_dof_int_ent, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table, CL_TRUE, 0,
			tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table_bytes, assembly_table, 0, NULL, NULL);

	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.assembly_struct.ocl_local_to_global, CL_TRUE, 0,
			tmv_ocl_problem_struct.assembly_struct.ocl_local_to_global_bytes, local_to_global, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.assembly_struct.ocl_global_to_posglob, CL_TRUE, 0,
			global_to_posblob_bytes, global_to_posglob, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.assembly_struct.ocl_pos_first_dof_int_ent, CL_TRUE, 0,
			tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent_bytes, pos_first_dof_int_ent, 0, NULL, NULL);


	//----------ALLOCATE MEMEMORY FOR GEO/GAUSS/SHAPE DATA --------

	// Geometric data - stored by elements
	tmv_ocl_problem_struct.num_int_struct.nr_geo = ngeo;
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in_size = nr_elems * ngeo * 3;
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in_bytes = nr_elems * ngeo * 3 * sizeof(SCALAR);

	// Create buffer
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in_bytes, NULL, NULL);

	// Write data
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in, CL_TRUE, 0,
			tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in_bytes, el_geo_dat_all_color_host, 0, NULL, NULL);


	// Geometric data - strored by dofs id
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in_size = nr_dof_ent * 3;
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in_bytes = nr_dof_ent * 3 * sizeof(SCALAR);

	// Create buffer
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in_bytes, NULL, NULL);

	// Write data
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in, CL_TRUE, 0,
			tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in_bytes, geo_dofs_vector, 0, NULL, NULL);


	// Gauss data
	tmv_ocl_problem_struct.num_int_struct.ngauss = ngauss;
	tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in_size = ngauss * 4;
	tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in_bytes = ngauss * 4 * sizeof(SCALAR);

	// Create buffer
	tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in_bytes, NULL, NULL);

	// Write data
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in, CL_TRUE, 0,
			tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in_bytes, gauss_dat_all_color_host, 0, NULL, NULL);

	// Shape fun data
	tmv_ocl_problem_struct.num_int_struct.nshape = nshape;
	tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in_size = nshape * ngauss * 4;
	tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in_bytes = nshape * ngauss * 4 * sizeof(SCALAR);

	// Create buffer
	tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in_bytes, NULL, NULL);

	// Write data
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in, CL_TRUE, 0,
			tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in_bytes, shape_fun_dat_all_color_host, 0, NULL, NULL);

	//PDE data
	tmv_ocl_problem_struct.num_int_struct.all_el_pde_coeff_size = all_el_pde_coeff_size;
	tmv_ocl_problem_struct.num_int_struct.one_el_pde_coeff_size = one_el_pde_coeff_size;
	tmv_ocl_problem_struct.num_int_struct.one_int_p_pde_coeff_size = one_int_p_pde_coeff_size;

	tmv_ocl_problem_struct.num_int_struct.ocl_el_pde_dat_in_bytes = (all_el_pde_coeff_size +
					nr_elems * (one_el_pde_coeff_size + ngauss * one_int_p_pde_coeff_size)) * sizeof(SCALAR);

	// Create buffer
	tmv_ocl_problem_struct.num_int_struct.ocl_el_pde_data_in =
			clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_problem_struct.num_int_struct.ocl_el_pde_dat_in_bytes, NULL, NULL);

	// Write data
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.num_int_struct.ocl_el_pde_data_in, CL_TRUE, 0,
			tmv_ocl_problem_struct.num_int_struct.ocl_el_pde_dat_in_bytes, el_pde_dat_host, 0, NULL, NULL);


	// CPU only, geo position for color
	tmv_ocl_problem_struct.num_int_struct.ngeo_color_pos = ngeo_color_p ;

	//----------ALLOCATE MEMEMORY solve vec --------
	tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_current =
			clCreateBuffer(context, CL_MEM_READ_ONLY, (nrdofs_glob * sizeof(double)), NULL, NULL);

	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_current, CL_TRUE, 0,
			(nrdofs_glob * sizeof(double)), dofs_vector_current, 0, NULL, NULL);

	tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_prev_iter =
			clCreateBuffer(context, CL_MEM_READ_ONLY, (nrdofs_glob * sizeof(double)), NULL, NULL);

	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_prev_iter, CL_TRUE, 0,
			(nrdofs_glob * sizeof(double)), dofs_vector_prev_iter, 0, NULL, NULL);

	tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_prev_step =
			clCreateBuffer(context, CL_MEM_READ_ONLY, (nrdofs_glob * sizeof(double)), NULL, NULL);

	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_prev_step, CL_TRUE, 0,
			(nrdofs_glob * sizeof(double)), dofs_vector_prev_step, 0, NULL, NULL);

	return (1);
}

/*---------------------------------------------------------
  tmr_ocl_free_assembly_structures - free numerical integration
	  and assembling data structures created by tmr_ocl_create_assembly_structures
----------------------------------------------------------*/
int tmr_ocl_free_assembly_structures(
)
{
	// Release assembly struct
	tmv_ocl_problem_struct.assembly_struct.nr_asse_blocks_all_int_ent = 0;
	tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table_bytes = 0;
	tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent_bytes = 0;

	clReleaseMemObject(tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent);
	clReleaseMemObject(tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table);


	// Release GEO DATA
	tmv_ocl_problem_struct.num_int_struct.nr_geo = 0;
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in_size = 0;
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in_bytes = 0;

	clReleaseMemObject(tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in);

	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in_size = 0;
	tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in_bytes = 0;

	clReleaseMemObject(tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in);

	// Release GAUSS DATA
	tmv_ocl_problem_struct.num_int_struct.ngauss = 0;
	tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in_size = 0;
	tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in_bytes = 0;

	clReleaseMemObject(tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in);

	// Release SHAPE FUN DATA
	tmv_ocl_problem_struct.num_int_struct.nshape = 0;
	tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in_size = 0;
	tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in_bytes = 0;

	clReleaseMemObject(tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in);

	// Release PDE DATA
	tmv_ocl_problem_struct.num_int_struct.all_el_pde_coeff_size = 0;
	tmv_ocl_problem_struct.num_int_struct.one_el_pde_coeff_size = 0;
	tmv_ocl_problem_struct.num_int_struct.one_int_p_pde_coeff_size = 0;

	clReleaseMemObject(tmv_ocl_problem_struct.num_int_struct.ocl_el_pde_data_in);

	return (1);
}

/*---------------------------------------------------------
  tmr_ocl_create_assemble_stiff_mat_elem - create local element stiffness matrices
		and possibly assemble them to the global stiffness matrix
----------------------------------------------------------*/
int tmr_ocl_create_assemble_stiff_mat_elem(
		int Problem_id,
		//int Level_id,
		int Comp_type,         /* in: indicator for the scope of computations: */

		//extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
		//extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
		//extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
		//extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
		int * Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes

		//
		//int Nr_int_ent,
		//int* L_int_ent_type,
		//int* L_int_ent_id,

		//
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
)
{

	// OpenCL profiling
#ifdef TIME_TEST_OCL
	cl_ulong startTime;
	cl_ulong endTime;
#endif

	// OpenCL operation counter
#ifdef COUNT_OPER
	SCALAR count_oper[3] = {0.0, 0.0, 0.0};
	cl_mem ocl_count_oper;
#endif

	// ----------- [BEGIN] Get platform configuration [BEGIN] ----------------

	// choose the platform
	int platform_index = tmr_ocl_get_current_platform_index();

	int device_tmc_type = tmr_ocl_get_current_device_type();

	// choose device_index
	int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);

	// choose the context
	cl_context context = tmr_ocl_select_context(platform_index, device_index);

	// choose the command queue
	cl_command_queue command_queue =
			tmr_ocl_select_command_queue(platform_index, device_index);

	// THERE MAY BE SEVERAL KERNELS FOR THE CODE, THE INDEX IN OPENCL DATA
	// STRUCTURE FOR THE NUMERICAL INTEGRATION KERNEL IS ASSIGNED IN (tmd_ocl/tmh_ocl.h)
	int kernel_index = TMC_OCL_KERNEL_NUM_INT_INDEX;

	// Get kernel
	cl_kernel kernel = tmr_ocl_select_kernel(platform_index, device_index, kernel_index);

	// OpenCL device characteristics stored in data structure
	tmt_ocl_device_struct device_struct =
			tmv_ocl_struct.list_of_platforms[platform_index].list_of_devices[device_index];
	int max_num_comp_units = device_struct.max_num_comp_units;
	int max_work_group_size = device_struct.max_work_group_size;

	// ----------- [END] Get platform configuration [END] ----------------


	// ----------- [BEGIN] Select kernel algorithm [BEGIN] ----------------

	int kernel_version_alg = 0 ;
	First_geo_index = tmv_ocl_problem_struct.num_int_struct.ngeo_color_pos[First_geo_index];

	// KERNEL ALGORITHM LIST
#define ONE_EL_ONE_THREAD_KERNEL 0
#define ONE_EL_TWO_THREADS_KERNEL 1
#define ONE_EL_ONE_WORKGROUP_KERNEL 2

#ifdef ONE_EL_ONE_THREAD
	kernel_version_alg = ONE_EL_ONE_THREAD_KERNEL;
#ifdef DEBUG_TMD
	mf_log_info("Selected: kernel_index %d, version %d ONE_EL_ONE_THREAD)\n", kernel_index, kernel_version_alg);
#endif
#elif defined(ONE_EL_TWO_THREADS)
	kernel_version_alg = ONE_EL_TWO_THREADS_KERNEL;
#ifdef DEBUG_TMD
	mf_log_info("Selected: kernel_index %d, version %d ONE_EL_TWO_THREAD)\n", kernel_index, kernel_version_alg);
#endif
#elif defined ONE_EL_ONE_WORKGROUP
	kernel_version_alg = ONE_EL_ONE_WORKGROUP_KERNEL;
#ifdef DEBUG_TMD
	mf_log_info("Selected: kernel_index %d, version %d ONE_EL_ONE_ONE_EL_ONE_WORKGROUP)\n", kernel_index, kernel_version_alg);
#endif
#else
	mf_fatal_err("Wrong kernel version specified!!!). Exiting.\n");
	exit(-1);
#endif

	// ----------- [END] Select kernel algorithm [END] ---------------


	// ----------- [BEGIN] OPENCL TEST (CONTEXT/COMMAND_QUEUE/KERNEL) [BEGIN] ----------------

	if (context == NULL || command_queue == NULL || kernel == NULL) {

		mf_fatal_err("Failed to restore kernel for platform %d, device %d, kernel %d\n context %lu, command queue %lu, kernel %lu\n",
				platform_index, device_index, kernel_index, context, command_queue, kernel);
		exit(-1);
	}

	// ----------- [END] OPENCL TEST (CONTEXT/COMMAND_QUEUE/KERNEL) [END] ----------------


	// ----------- [BEGIN] Prepare task for execution on DEVICE  [BEGIN] ----------------

#define NR_EXEC_PARAMS 16  // size of array with execution parameters

	// In the current version the procedure is called for a single color
	int nr_elems = Nr_elems;

	// ASSUMED WORK_GROUP_SIZE
	int work_group_size = WORK_GROUP_SIZE;

	// Usually for GPUs it is good to maximize the number of work_groups and threads
	int nr_work_groups = NR_WORK_GROUPS_PER_COMP_UNIT * max_num_comp_units;
	int nr_elems_per_work_group;

	// Multiplicator
	double arbitrary_multiplication_constant = 1.2;

	// SELECT STRATEGY OF DIVIDE TASK BETWEEN THREADS ON DEVICE
	// Calculate number of elements per work group for differents strategies
	if (kernel_version_alg == ONE_EL_ONE_THREAD_KERNEL)
	{

		nr_elems_per_work_group =
				ceil(((double)nr_elems) / nr_work_groups / work_group_size) * work_group_size;

	}
	else if (kernel_version_alg == ONE_EL_TWO_THREADS_KERNEL)
	{

		nr_elems_per_work_group =
				ceil(((double)nr_elems) / nr_work_groups / work_group_size) * work_group_size;
	}
	else if (kernel_version_alg == ONE_EL_ONE_WORKGROUP_KERNEL)
	{

		nr_elems_per_work_group =
				ceil(((double)nr_elems) / nr_work_groups);

	}
	else
	{
		mf_fatal_err("Wrong kernel version specified!!!). Exiting.\n");
		exit(-1);
	}

	// DEFINE MINIMAL NUMBER OF ELEMENTS FOR DIFFRENTS STARTEGIES
#ifdef ONE_EL_ONE_THREAD
#define MIN_NR_ELEMS_PER_WORK_GROUP WORK_GROUP_SIZE
#endif

#ifdef ONE_EL_TWO_THREADS
#define MIN_NR_ELEMS_PER_WORK_GROUP (WORK_GROUP_SIZE/2)
#endif

#ifdef ONE_EL_ONE_WORKGROUP
#define MIN_NR_ELEMS_PER_WORK_GROUP 1
#endif

	if (nr_elems_per_work_group < MIN_NR_ELEMS_PER_WORK_GROUP)
	{
		mf_fatal_err("Too low number of elements per workgroup %d\n" \
				"For this strategies minimal number of elements is %d\n" \
				"Check parameters of execution. Exiting!!!\n",
				nr_elems_per_work_group, MIN_NR_ELEMS_PER_WORK_GROUP);
		exit - 1;
	}

	// CALCULATE NUMBER OF THREADS

	// Finally we decide how many threads will perform calculations
	int nr_threads = nr_work_groups * work_group_size;
	int nr_elems_per_thread = 0;

	// Set number of elements per thread
	if (kernel_version_alg == ONE_EL_ONE_THREAD_KERNEL)
	{
		nr_elems_per_thread = nr_elems_per_work_group / work_group_size;
	}
	else if (kernel_version_alg == ONE_EL_TWO_THREADS_KERNEL)
	{
		nr_elems_per_thread = nr_elems_per_work_group;
	}
	else if (kernel_version_alg == ONE_EL_ONE_WORKGROUP_KERNEL)
	{
		nr_elems_per_thread = nr_elems_per_work_group;
	}
	else
	{
		mf_fatal_err("Wrong kernel version specified!!!). Exiting.\n");
		exit(-1);
	}


	// Extended log info
#ifdef DEBUG_TMD
	if (kernel_version_alg == ONE_EL_ONE_THREAD_KERNEL)
	{
		mf_log_info("Kernel algorithm version: ONE_EL_ONE_THREAD_KERNEL\n");
	}
	else if (kernel_version_alg == ONE_EL_TWO_THREADS_KERNEL)
	{
		mf_log_info("Kernel algorithm version: ONE_EL_TWO_THREADS_KERNEL\n");
	}
	else if (kernel_version_alg == ONE_EL_ONE_WORKGROUP_KERNEL)
	{
		mf_log_info("Kernel algorithm version: ONE_EL_ONE_WORKGROUP_KERNEL\n");
	}
	mf_log_info(" --> Threads mesh on DEVICE:\n"	\
			"Number of threads: %d\n"	\
			"Workgroup size: %d\n" \
			"Number of workgroups: %d\n",
			nr_threads, work_group_size, nr_work_groups);

	mf_log_info(" --> Divide of elements on DEVICE:\n"	\
			"Number of elements per KERNEL: %d\n"	\
			"Number of elements per WORK GROUP: %d\n" \
			"Number of elements per THREAD: %d\n",
			nr_elems, nr_elems_per_work_group, nr_elems_per_thread);
#endif


	//
	// CREATE BUFFER WITH EXECUTION PARAMETERS
	//


	// Size of execution parameters table
	const size_t execution_parameters_dev_bytes = NR_EXEC_PARAMS * sizeof(int);

	// Create memory object for store execution parameters table
	cl_mem execution_parameters_dev;

	// Create buffer on device for store execution parameters table
	execution_parameters_dev = clCreateBuffer(context, CL_MEM_READ_ONLY,
					execution_parameters_dev_bytes, NULL, NULL);


	//
	// FILL BUFFER WITH EXECUTION PARAMETERS
	//
	int execution_parameters_host[NR_EXEC_PARAMS];

	// Number of elements per thread for this kernel
	execution_parameters_host[0] = nr_elems_per_thread;

	// Number of elements for this kernel
	execution_parameters_host[1] = nr_elems;

	// Index of first geo for color
	execution_parameters_host[2] = First_geo_index;

	// Index of first pde for color
	execution_parameters_host[3] = First_pde_index;

	// Index of first element in assembly table for color
	execution_parameters_host[4] = First_asse_elem_index;

#ifdef DEBUG_TMD
	printf("\n[DEBUG] Execution_parameters_host %d", NR_EXEC_PARAMS);
	printf("\n[DEBUG]  - execution_parameters_host[0] (nr_elems_per_thread) = %d", nr_elems_per_thread);
	printf("\n[DEBUG]  - execution_parameters_host[1] (nr_elems) = %d", nr_elems);
	printf("\n[DEBUG]  - execution_parameters_host[2] (First_geo_index) = %d", First_geo_index);
	printf("\n[DEBUG]  - execution_parameters_host[3] (First_pde_index) = %d", First_pde_index);
	printf("\n[DEBUG]  - execution_parameters_host[4] (First_asse_elem_index) = %d", First_asse_elem_index);
	printf("\n");
#endif

	// send data to device - non-blocking call (CL_FALSE)
	clEnqueueWriteBuffer(command_queue, execution_parameters_dev, CL_TRUE, 0,
			execution_parameters_dev_bytes, execution_parameters_host,
			0, NULL, NULL);

	// ----------- [END] Prepare task for execution on DEVICE  [END] ----------------


	// ----------- [BEGIN] Set KERNEL arguments  [BEGIN] ----------------

	// Data will be sent to ACCELERATOR in following order:
	// 1) EXECUTION PARAMETERS
	// 2) INTEGRATION POINTS DATA (GAUSS DATA)
	// 3) REFERENCE ELEMENT SHAPE FUNCTION
	// 4) GEOMETRIC COORDINATES OF ELEMENTS
	// 5) PDE DATA (old el_data_in) ......
	// 6/7) ASSEMBY TABLE
	// 8/9) RHS and CRS
	// 10/11/12) SOLUTION VECTORS
	//

	// OpenCL test
	cl_int retval = CL_SUCCESS;

	// Kernel parameters
	retval |= clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &execution_parameters_dev);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 0, (__LINE__ - 3), __func__);
#endif

	// Intergration points data
	retval |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.num_int_struct.ocl_gauss_data_in);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 1, (__LINE__ - 3), __func__);
#endif

	// Reference element shape function data
	retval |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.num_int_struct.ocl_shape_fun_data_in);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 2, (__LINE__ - 3), __func__);
#endif

	// Geometry of elements data
#ifdef GEO_DATA_STORED_BY_DOFS
	retval |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_by_dofs_in);
#else
	retval |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.num_int_struct.ocl_el_geo_data_in);
#endif
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 3, (__LINE__ - 3), __func__);
#endif

	// PDE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	retval |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.num_int_struct.ocl_el_pde_data_in);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 4, (__LINE__ - 3), __func__);
#endif
	// PDE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// Assembly tables
	retval |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 5, (__LINE__ - 3), __func__);
#endif
	retval |= clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 6, (__LINE__ - 3), __func__);
#endif

	// Local to global table
	retval |= clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.assembly_struct.ocl_pos_first_dof_int_ent);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 7, (__LINE__ - 3), __func__);
#endif
	retval |= clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.assembly_struct.ocl_local_to_global);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 8, (__LINE__ - 3), __func__);
#endif

	retval |= clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.assembly_struct.ocl_global_to_posglob);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 9, (__LINE__ - 3), __func__);
#endif

	// RHS and CRS
	retval |= clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.crs_struct.ocl_rhs_val);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 10, (__LINE__ - 3), __func__);
#endif
	retval |= clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.crs_struct.ocl_crs_val);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 11, (__LINE__ - 3), __func__);
#endif

	// Solution vectors
	retval |= clSetKernelArg(kernel, 12, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_current);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 12, (__LINE__ - 3), __func__);
#endif
	retval |= clSetKernelArg(kernel, 13, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_prev_iter);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 13, (__LINE__ - 3), __func__);
#endif
	retval |= clSetKernelArg(kernel, 14, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.solution_struct.ocl_dofs_vector_prev_step);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 14, (__LINE__ - 3), __func__);
#endif

#ifdef COUNT_OPER
	ocl_count_oper = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 3 * sizeof(SCALAR), NULL, NULL);
	retval |= clSetKernelArg(kernel, 15, sizeof(cl_mem), (void *) &ocl_count_oper);
#ifdef DEBUG_TMD
	tmr_set_kernel_arguments_log(retval, 15, (__LINE__ - 3), __func__);
#endif
#endif

	// Correctness test
	if (retval != CL_SUCCESS) {
		mf_fatal_err("Failed to Set the kernel argument %d.\n", 0);
		exit(-1);
	}


	// ----------- [END] Set KERNEL arguments  [END] ----------------



	// ----------- [BEGIN] EXECUTE KERNEL  [BEGIN] ----------------


	// Set kernel invocation parameters
	size_t globalWorkSize[1] = { 0 };
	size_t localWorkSize[1] = { 0 };

	globalWorkSize[0] = nr_threads ;
	localWorkSize[0] = work_group_size ;

	// Events
	cl_event ndrEvt;

	//Queue the kernel up for execution
	retval = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
					globalWorkSize, localWorkSize,
					0, NULL,  &ndrEvt);
#ifdef DEBUG_TMD
	tmr_kernel_execution_log(retval, (__LINE__ - 4), __func__);
#endif

#ifdef TIME_TEST_OCL

	// Wait for events and finish command queue
	clWaitForEvents(1, &ndrEvt);
	clFinish(command_queue);

	// Get kernel profiling info
	clGetEventProfilingInfo(ndrEvt,
			CL_PROFILING_COMMAND_START,
			sizeof(cl_ulong),
			&startTime,
			0);
	clGetEventProfilingInfo(ndrEvt,
			CL_PROFILING_COMMAND_END,
			sizeof(cl_ulong),
			&endTime,
			0);
	double time_internal = ((double)endTime - (double)startTime) * 1.0e-9;
	utv_time.GPU_integration_and_assembly += time_internal;

	mf_log_info("KERNEL EXECUTION TIME (internal): %lf\n", time_internal);


#endif

#ifdef COUNT_OPER
	clEnqueueReadBuffer(command_queue, ocl_count_oper, CL_TRUE, 0, 3 * sizeof(SCALAR), count_oper, 0, NULL, NULL);
	mf_log_info("CRS Arthmetic operations: %.0lf, Local mem access: %.0lf, Global mem access: %.0lf\n", count_oper[0], count_oper[1], count_oper[2]);
	clReleaseMemObject(ocl_count_oper);
#endif


	// ----------- [END] EXECUTE KERNEL  [END] ----------------


	return (1);
}

/*---------------------------------------------------------
tmr_ocl_create_solver_structures_crs - to create CRS structure on ACCELERATOR
/*---------------------------------------------------------*/
int tmr_ocl_create_solver_structures_crs(
		int Nrdof_glob,
		int nnz,
		int * crs_col_ind,
		int * crs_row_ptr,
		double * crs_val,
		double * rhs )
{

	// ----------- Get platform configuration ----------------

	// choose the platform
	int platform_index = tmr_ocl_get_current_platform_index();

	int device_tmc_type = tmr_ocl_get_current_device_type();

	// choose device_index
	int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);

	// choose the context
	cl_context context = tmr_ocl_select_context(platform_index, device_index);

	// choose the command queue
	cl_command_queue command_queue =
			tmr_ocl_select_command_queue(platform_index, device_index);

	// -----------End Get platform configuration ----------------


	// ----------- Set crs data ----------------

	tmv_ocl_problem_struct.crs_struct.crs_col_ind = crs_col_ind;
	tmv_ocl_problem_struct.crs_struct.crs_row_ptr = crs_row_ptr;
	tmv_ocl_problem_struct.crs_struct.crs_val_cpu = crs_val;
	tmv_ocl_problem_struct.crs_struct.rhs_val_cpu = rhs;
	tmv_ocl_problem_struct.crs_struct.crs_val_gpu = (SCALAR *) malloc (nnz * sizeof(SCALAR));
	tmv_ocl_problem_struct.crs_struct.rhs_val_gpu = (SCALAR *) malloc (Nrdof_glob * sizeof(SCALAR));

	// For opencl 1.1
	memset(tmv_ocl_problem_struct.crs_struct.crs_val_gpu, 0.0, nnz * sizeof(SCALAR));
	memset(tmv_ocl_problem_struct.crs_struct.rhs_val_gpu, 0.0, Nrdof_glob * sizeof(SCALAR));

	tmv_ocl_problem_struct.crs_struct.Nnz = nnz;
	tmv_ocl_problem_struct.crs_struct.Nrdof = Nrdof_glob;
	tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes = nnz * sizeof(SCALAR);
	tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind_bytes = (nnz + 1) * sizeof(int);
	tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr_bytes = (Nrdof_glob + 1) * sizeof(int);
	tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes = Nrdof_glob * sizeof(SCALAR);



	/* ------------- ALLOCATE MEMORY CRS STRUCTURE -------------------- */

	tmv_ocl_problem_struct.crs_struct.ocl_crs_val = clCreateBuffer(context, CL_MEM_READ_WRITE,
					tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes, NULL, NULL);

	tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind = clCreateBuffer(context, CL_MEM_READ_ONLY,
					tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind_bytes, NULL, NULL);


	tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr = clCreateBuffer(context, CL_MEM_READ_ONLY,
					tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr_bytes, NULL, NULL);

	tmv_ocl_problem_struct.crs_struct.ocl_rhs_val = clCreateBuffer(context, CL_MEM_READ_WRITE,
					tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes, NULL, NULL);


	// Write CRS data column indices and rows pointer

	/* structure must be created before write */

	// Write column index table
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind, CL_TRUE, 0,
			tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind_bytes, tmv_ocl_problem_struct.crs_struct.crs_col_ind,
			0, NULL, NULL);

	// Write row pointer
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr, CL_TRUE, 0,
			tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr_bytes, tmv_ocl_problem_struct.crs_struct.crs_row_ptr,
			0, NULL, NULL);

	// Write val buffer on gpu for opencl 1.1
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_crs_val, CL_TRUE, 0,
			tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes, tmv_ocl_problem_struct.crs_struct.crs_val_gpu,
			0, NULL, NULL);
	// Write val buffer on gpu for opencl 1.1
	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_rhs_val, CL_TRUE, 0,
			tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes, tmv_ocl_problem_struct.crs_struct.rhs_val_gpu,
			0, NULL, NULL);


	// Write val buffer on gpu for new opencl > 1.1
	// double fill_data = 0.0;
	// clEnqueueFillBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_crs_val,&fill_data, sizeof(double), 0,
	// tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes, 0, NULL, NULL);
	// Write val buffer on gpu for new opencl > 1.1
	// double fill_data = 0.0;
	// clEnqueueFillBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_rhs_val,&fill_data, sizeof(double), 0,
	// tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes, 0, NULL, NULL);



	return (1);
}

/*---------------------------------------------------------
tmr_ocl_get_solver_structures_crs - to get CRS structure on ACCELERATOR
/*---------------------------------------------------------*/
int tmr_ocl_get_solver_structures_crs(
		int * nnz,
		cl_mem * crs_col_ind,
		cl_mem * crs_row_ptr,
		cl_mem * crs_val,
		cl_mem * rhs )
{

	*crs_col_ind = tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind;
	*crs_row_ptr = tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr;
	*crs_val = tmv_ocl_problem_struct.crs_struct.ocl_crs_val;
	*rhs = tmv_ocl_problem_struct.crs_struct.ocl_rhs_val;
	*nnz = tmv_ocl_problem_struct.crs_struct.Nnz;

	return (1);
}



/*---------------------------------------------------------
tmr_ocl_free_solver_structures_crs - to release CRS structure on ACCELERATOR
/*---------------------------------------------------------*/
int tmr_ocl_free_solver_structures_crs(				       )
{


	// Free CRS data
	if (tmv_ocl_problem_struct.crs_struct.crs_val_gpu != NULL) {
		free(tmv_ocl_problem_struct.crs_struct.crs_val_gpu);
		tmv_ocl_problem_struct.crs_struct.crs_val_gpu = NULL;
	}
	if (tmv_ocl_problem_struct.crs_struct.rhs_val_gpu != NULL) {
		free(tmv_ocl_problem_struct.crs_struct.rhs_val_gpu);
		tmv_ocl_problem_struct.crs_struct.rhs_val_gpu = NULL;
	}
	tmv_ocl_problem_struct.crs_struct.crs_col_ind = NULL;
	tmv_ocl_problem_struct.crs_struct.crs_row_ptr = NULL;

	// Release GPU memory objects
	//clReleaseMemObject(tmv_ocl_problem_struct.assembly_struct.ocl_asse_pos_first_dof_int_ent);
	//clReleaseMemObject(tmv_ocl_problem_struct.assembly_struct.ocl_assembly_table);
	clReleaseMemObject(tmv_ocl_problem_struct.crs_struct.ocl_crs_val);
	//clReleaseMemObject(tmv_ocl_problem_struct.crs_struct.ocl_crs_external_val);
	clReleaseMemObject(tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind);
	clReleaseMemObject(tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr);
	clReleaseMemObject(tmv_ocl_problem_struct.crs_struct.ocl_rhs_val);
	//clReleaseMemObject(tmv_ocl_problem_struct.crs_struct.ocl_rhs_external_val);

	// Cleaning
	tmv_ocl_problem_struct.crs_struct.Nnz = 0;
	tmv_ocl_problem_struct.crs_struct.Nrdof = 0;
	tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes = 0;
	tmv_ocl_problem_struct.crs_struct.ocl_crs_col_ind_bytes = 0;
	tmv_ocl_problem_struct.crs_struct.ocl_crs_row_ptr_bytes = 0;
	tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes = 0;


	return (1);
}

//tmr_ocl_send_data_to_gpu_crs ??????????

//tmr_ocl_prepare_final_crs ????????????? tmr_ocl_merge_solver_structures

//tmr_ocl_solve_viennacl_gpu ????????????


/*---------------------------------------------------------
  tmr_ocl_merge_solver_structures - to merge parts of SM and LV created on CPU and GPU
---------------------------------------------------------*/
int tmr_ocl_merge_solver_structures() {

#ifdef TIME_TEST_OCL
	double kernel_execution_time = 0.0;
#endif

#ifdef COUNT_OPER
	SCALAR count_oper[3];
	cl_mem ocl_count_oper;
#endif

	// ----------- Get platform configuration ----------------

	int i;

	// choose the platform
	int platform_index = tmr_ocl_get_current_platform_index();

	int device_tmc_type = tmr_ocl_get_current_device_type();

	// choose device_index
	int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);

	// OpenCL device characteristics stored in data structure
	tmt_ocl_device_struct device_struct =
			tmv_ocl_struct.list_of_platforms[platform_index].list_of_devices[device_index];
	int max_num_comp_units = device_struct.max_num_comp_units;
	int max_work_group_size = device_struct.max_work_group_size;

	// ASSUMED WORK_GROUP_SIZE
	int work_group_size = WORK_GROUP_SIZE;

	// usually for GPUs it is good to maximize the number of work_groups and threads
	int nr_work_groups = NR_WORK_GROUPS_PER_COMP_UNIT * max_num_comp_units;
	int nr_threads = nr_work_groups * work_group_size;

	// choose the context
	cl_context context = tmr_ocl_select_context(platform_index, device_index);

	// choose the command queue
	cl_command_queue command_queue =
			tmr_ocl_select_command_queue(platform_index, device_index);

	// THERE MAY BE SEVERAL KERNELS FOR THE CODE, THE INDEX IN OPENCL DATA
	// STRUCTURE FOR THE NUMERICAL INTEGRATION KERNEL IS ASSIGNED IN (tmd_ocl/tmh_ocl.h)
	int kernel_index = TMC_OCL_KERNEL_CRS_FINALIZE;

	// !!!!! FOR THE TIME BEING WE ALWAYS READ AND COMPILE tmr_ocl_num_int_el_and_assembling.cl
	cl_kernel kernel = tmr_ocl_select_kernel(platform_index, device_index, kernel_index);

#ifndef NR_EXEC_PARAMS
#define NR_EXEC_PARAMS 16  // size of array with execution parameters
#endif

	// WE ALLOCATE OPENCL BUFFERS FOR nr_elems_per_kernel ELEMENTS
	// this should be good for all kernel invocations and all colors !!!
	// (actual kernel calculations will be for nr_elems_this_kercall)

	cl_int retval = CL_SUCCESS;

#ifdef TIME_TEST_OCL
	t_begin = time_clock();
	//t_begin = time_clock();
#endif


	// there are NR_EXEC_PARAMS execution parameters (ALLWAYS CHECK!!!)

	const size_t execution_parameters_dev_bytes = NR_EXEC_PARAMS * sizeof(int);

	cl_mem execution_parameters_dev;

	execution_parameters_dev = clCreateBuffer(context, CL_MEM_READ_ONLY,
					execution_parameters_dev_bytes, NULL, NULL);
	retval |= clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &execution_parameters_dev);

	// Create OpenCL memory objects (only create and copy necessary data without kernel assigning

	/* ------------- ALLOCATE MEMORY CRS STRUCTURE -------------------- */

	tmv_ocl_problem_struct.crs_struct.ocl_crs_external_val = clCreateBuffer(context, CL_MEM_READ_ONLY,
					tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes, NULL, NULL);

	tmv_ocl_problem_struct.crs_struct.ocl_rhs_external_val = clCreateBuffer(context, CL_MEM_READ_ONLY,
					tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes, NULL, NULL);

	/* structure must be created before write */

	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_crs_external_val, CL_TRUE, 0,
			tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes, tmv_ocl_problem_struct.crs_struct.crs_val_cpu,
			0, NULL, NULL);

	clEnqueueWriteBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_rhs_external_val, CL_TRUE, 0,
			tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes, tmv_ocl_problem_struct.crs_struct.rhs_val_cpu,
			0, NULL, NULL);

#ifdef COUNT_OPER
	ocl_count_oper = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 3 * sizeof(SCALAR), NULL, NULL);
#endif

	/* ------------- BIND DATA -------------------- */

	// Bind externel RHS vector
	retval |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.crs_struct.ocl_rhs_external_val);

	// Bind externel CRS structure
	retval |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.crs_struct.ocl_crs_external_val);


	// Bind RHS vector
	retval |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.crs_struct.ocl_rhs_val);

	// Bind CRS structure
	retval |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &tmv_ocl_problem_struct.crs_struct.ocl_crs_val);

#ifdef COUNT_OPER
	retval |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &ocl_count_oper);
#endif

	/* ------------- PERFORM KERNEL -------------------- */

	int execution_parameters_host[NR_EXEC_PARAMS];
	execution_parameters_host[0] = tmv_ocl_problem_struct.crs_struct.Nrdof;
	execution_parameters_host[1] = tmv_ocl_problem_struct.crs_struct.Nnz;
	execution_parameters_host[2] = tmv_ocl_problem_struct.crs_struct.Nrdof / (nr_work_groups * work_group_size);
	execution_parameters_host[3] = tmv_ocl_problem_struct.crs_struct.Nrdof % (nr_work_groups * work_group_size);
	execution_parameters_host[4] = tmv_ocl_problem_struct.crs_struct.Nnz / (nr_work_groups * work_group_size);
	execution_parameters_host[5] = tmv_ocl_problem_struct.crs_struct.Nnz % (nr_work_groups * work_group_size);

	// send data to device - non-blocking call (CL_FALSE)
	clEnqueueWriteBuffer(command_queue, execution_parameters_dev, CL_TRUE, 0,
			execution_parameters_dev_bytes, execution_parameters_host,
			0, NULL, NULL);


#ifdef TIME_TEST_OCL
	clFinish(command_queue);
	t_end = time_clock();
	//t_end = time_clock();

	ocl_total_time += t_end - t_begin;
	mf_log_info("CRS: copying input crs buffer: %.12lf, total time: %.12lf\n",
			t_end - t_begin, ocl_total_time);
#endif

	// FINALLY EXECUTE THE KERNEL
	// set kernel invocation parameters
	size_t globalWorkSize[1] = { 0 };
	size_t localWorkSize[1] = { 0 };
	globalWorkSize[0] = nr_threads ;
	localWorkSize[0] = work_group_size ;

	cl_event ndrEvt;

#ifdef TIME_TEST_OCL
	t_begin = time_clock();
	//t_begin = time_clock();
#endif

	// Queue the kernel up for execution
	retval = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
					globalWorkSize, localWorkSize,
					0, NULL,  &ndrEvt);

#ifdef DEBUG_TMD
	tmr_kernel_execution_log(retval, (__LINE__ - 5), "tmr_ocl_merge_solver_structures");
#endif

#ifdef TIME_TEST_OCL
	clWaitForEvents(1, &ndrEvt);
	clFinish(command_queue);

	t_end = time_clock();

	// Get kernel profiling info
	cl_ulong startTime;
	cl_ulong endTime;
	clGetEventProfilingInfo(ndrEvt,
			CL_PROFILING_COMMAND_START,
			sizeof(cl_ulong),
			&startTime,
			0);
	clGetEventProfilingInfo(ndrEvt,
			CL_PROFILING_COMMAND_END,
			sizeof(cl_ulong),
			&endTime,
			0);
	double time_internal = ((double)endTime - (double)startTime) * 1.0e-9;

	// double kernel_execution_time=0.0;
	kernel_execution_time = t_end - t_begin;
	#pragma omp atomic
	ocl_time_total_kernel_execution += (t_end - t_begin);

	mf_log_info("CRS: executing addition kernel: %.12lf (internal %.12lf)\n",
			kernel_execution_time, time_internal);

	ocl_total_time += t_end - t_begin;
	t_begin = time_clock();
#endif


#ifdef GPU_DEBUG
	memset(tmv_ocl_problem_struct.crs_struct.crs_val_cpu, 0, tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes);
	memset(tmv_ocl_problem_struct.crs_struct.rhs_val_cpu, 0, tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes);
#endif

	clEnqueueReadBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_crs_val, CL_TRUE, 0, tmv_ocl_problem_struct.crs_struct.ocl_crs_val_bytes, tmv_ocl_problem_struct.crs_struct.crs_val_cpu, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, tmv_ocl_problem_struct.crs_struct.ocl_rhs_val, CL_TRUE, 0, tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes, tmv_ocl_problem_struct.crs_struct.rhs_val_cpu, 0, NULL, NULL);

	/*
	  for(i=0; i<tmv_ocl_problem_struct.crs_struct.Nnz; i++) {
		tmv_ocl_problem_struct.crs_struct.crs_val_cpu[i]= tmv_ocl_problem_struct.crs_struct.crs_val_gpu[i];
	  }

	  for(i=0; i<(tmv_ocl_problem_struct.crs_struct.ocl_rhs_bytes/sizeof(double)); i++) {
		tmv_ocl_problem_struct.crs_struct.rhs_val_cpu[i]=tmv_ocl_problem_struct.crs_struct.rhs_val_gpu[i];
	  }
	*/

#ifdef COUNT_OPER
	clEnqueueReadBuffer(command_queue, ocl_count_oper, CL_TRUE, 0, 3 * sizeof(SCALAR), count_oper, 0, NULL, NULL);
#endif

#ifdef TIME_TEST_OCL
	clFinish(command_queue);
	t_end = time_clock();
	ocl_total_time += t_end - t_begin;
	#pragma omp atomic
	ocl_time_total_transfer += ocl_total_time;

	mf_log_info("CRS: copying output crs buffer: %.12lf, total time: %.12lf\n",
			t_end - t_begin, ocl_total_time);
#endif

#ifdef COUNT_OPER
	mf_log_info("CRS: %.0lf, Local mem access: %.0lf, Global mem access: %.0lf\n", count_oper[0], count_oper[1], count_oper[2]);
	clReleaseMemObject(ocl_count_oper);
#endif

	return (0);

}

/*---------------------------------------------------------
  tmr_ocl_create_num_int_kernel - for a selected platform, device and kernel index
---------------------------------------------------------*/
int tmr_ocl_create_num_int_kernel(
		int Platform_index,
		int Device_index,
		int Kernel_index,
		char * Kernel_name,
		const char * Kernel_file,
		int Monitor
)
{

	FILE * Interactive_output = stdout;
	cl_int retval;

	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];

	// check the device !!!!!!!!!!!!!!!!! (or at least its index)
	if (Device_index < 0) {
		printf("Wrong device_index %d passed to tmr_ocl_create_kernel! Exiting.\n",
				Device_index);
		exit(-1);
	}

	cl_device_id device = platform_struct.list_of_devices[Device_index].id;
	cl_context context = platform_struct.list_of_contexts[
				 platform_struct.list_of_devices[Device_index].context_index
		 ];

	if (Monitor > TMC_PRINT_INFO) {
		printf("Program file is: %s\n", Kernel_file);
	}

#ifdef TUNING
	//for assembler
	system("rm -rf ~/.nv");
#endif

	// read source file into data structure
	const char * source = tmr_ocl_readSource(Interactive_output, Kernel_file);

	cl_program program = clCreateProgramWithSource(context, 1,
					&source,
					NULL, NULL);
	if (program == NULL)
	{
		printf("Failed to create CL program from source.\n");
		exit(-1);
	}

#ifdef TUNING

	optf = fopen("options.txt", "r");
	if (!optf) {
		printf("Could not open options file!\n");
		exit(-1);
	}
	resf = fopen("result.csv", "a+");
	if (!resf) {
		printf("Could not open results file!\n");
		exit(-1);
	}

#define NOPT 9

	char opt[NOPT][40];
	char buffer[40 * (NOPT + 10)];
	int options[NOPT];
	int tmp;
	int i;
	for (i = 0; i < 40 * (NOPT + 10); i++)
		buffer[i] = 0;



	sprintf(&opt[0], " -D COAL_READ");
	sprintf(&opt[1], " -D COAL_WRITE");
	sprintf(&opt[2], " -D COMPUTE_ALL_SHAPE_FUN_DER");
	sprintf(&opt[3], " -D USE_WORKSPACE_FOR_PDE_COEFF");
	sprintf(&opt[4], " -D USE_WORKSPACE_FOR_GEO_DATA");
	sprintf(&opt[5], " -D USE_WORKSPACE_FOR_SHAPE_FUN");
	sprintf(&opt[6], " -D USE_WORKSPACE_FOR_STIFF_MAT");
	sprintf(&opt[7], " -D WORKSPACE_PADDING=0");
	sprintf(&opt[8], " -D WORKSPACE_PADDING=1");

//	  sprintf(&opt[0]," -D USE_WORKSPACE");
//	  sprintf(&opt[1]," -D USE_REGISTERS_FOR_COEFF");
//	  sprintf(&opt[2]," -D USE_SHAPE_FUN_WORKSPACE");
//	  sprintf(&opt[3]," -D USE_REGISTERS_FOR_SHAPE_FUN");
//	  sprintf(&opt[4]," -D USE_SHAPE_FUN_REF_DIRECTLY");
//	  sprintf(&opt[5]," -D STIFF_MAT_IN_SHARED");
//	  sprintf(&opt[6]," -D PADDING=0");
//	  sprintf(&opt[7]," -D PADDING=1");

	//sprintf(&opt[7]," -D COAL_READ");
	//sprintf(&opt[8]," -D COAL_WRITE");
	//sprintf(&opt[9]," -D CONSTANT_COEFF");
	//sprintf(&opt[9]," -D COUNT_OPER");

	unsigned long line_count = 0;
	int c;

	while ( (c = fgetc(resf)) != EOF ) {
		if ( c == '\n' )
			line_count++;
	}

	printf("Result file has %u lines\n", line_count);
	if (line_count == 0)
	{
		headf = fopen("header.csv", "w");
		if (!headf) {
			printf("Could not open header file!\n");
			exit(-1);
		}
		system("touch num");
	}

	for (i = 0; i < NOPT * line_count; i++)
		fscanf(optf, "%d", &tmp);

	for (i = 0; i < NOPT; i++)
		fscanf(optf, "%d", &options[i]);

	printf("Options indicators:\n");
	for (i = 0; i < NOPT; i++)
	{
		if (line_count == 0)
		{
			fprintf(headf, "%s,", opt[i]);
		}
		printf("%d ,", options[i]);
		fprintf(resf, "%d ,", options[i]);
	}
	//fprintf(resf,";");
	printf("Options:\n");
	strcpy(buffer, " ");
	for (i = 0; i < NOPT; i++)
	{
		if (options[i] == 1)
		{
			strcat(buffer, opt[i]);
			if (i == 1)
			{
				system("touch coal_write");
			}
		}
		else
			strcat(buffer, " ");
	}

#ifdef OPENCL_CPU
	strcat(buffer, " -D WORK_GROUP_SIZE=8");
#elif defined OPENCL_GPU
	strcat(buffer, " -D WORK_GROUP_SIZE=64");
#elif defined OPENCL_PHI
	strcat(buffer, " -D WORK_GROUP_SIZE=16");
#endif

#ifdef LAPLACE
	strcat(buffer, " -D LAPLACE");
#elif TEST_SCALAR
	strcat(buffer, " -D TEST_NUMINT");
#elif HEAT
	strcat(buffer, " -D HEAT");
#endif

#ifdef COUNT_OPER
	strcat(buffer, " -D COUNT_OPER");
#endif

#ifdef OPENCL_HSA
	strcat(buffer, " -cl-std=CL2.0");
#endif

	if (tmr_ocl_check_vendor() == 1)
		strcat(buffer, " -cl-nv-verbose");

	if (tmr_ocl_check_vendor() == 3)
	{
		//strcat(buffer," -fno-bin-source -fno-bin-llvmir -fno-bin-amdil -save-temps=");
		strcat(buffer, " -save-temps=");
		system("mkdir kernele");
		char name[50];
		int word;
		FILE * fp;
		for (i = 0; i < 50; i++)
			name[i] = 0;

		strcat(name, "kernele/");
		for (i = 0; i < NOPT; i++)
		{
			if (options[i] == 1)
				strcat(name, "1");
			else
				strcat(name, "0");
		}


#ifdef OPENCL_GPU
#ifdef OPENCL_HSA
		strcat(name, "_HSA");
#else
		strcat(name, "_GPU");
#endif
#endif
#ifdef OPENCL_CPU
		strcat(name, "_CPU");
#endif

		word = wc("tmr_ocl_num_int_el.cl", "//QSS");
		if (word == 1)
			strcat(name, "_QSS");
		word = wc("tmr_ocl_num_int_el.cl", "//SQS");
		if (word == 1)
			strcat(name, "_SQS");
		word = wc("tmr_ocl_num_int_el.cl", "//SSQ");
		if (word == 1)
			strcat(name, "_SSQ");

		strcat(buffer, name);

	}

	printf("%s\n", buffer);

	retval = clBuildProgram(program, 0, NULL, buffer, NULL, NULL);

#else

#ifdef OPENCL_HSA
	retval = clBuildProgram(program, 0, NULL, "-cl-std=CL2.0", NULL, NULL);
#else

	if (tmr_ocl_check_vendor() == 1)
	{
		// TO GET INFO FROM NVIDIA COMPILER
		retval = clBuildProgram(program, 0, NULL, "-cl-nv-verbose", NULL, NULL);
		// TO FORCE NVIDIA COMPILER TO LIMIT THE NUMBER OF USED REGISTERS
		// retval = clBuildProgram(program, 0, NULL, "-cl-nv-maxrregcount=32", NULL, NULL);
	}
	else
	{
		// generic call - no compiler options
		retval = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
		//retval = clBuildProgram(program, 0, NULL, "-cl-nv-verbose", NULL, NULL);
	}

#endif

#endif //TUNING


	char * buildLog;
	size_t size_of_buildLog;
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
			0, NULL, &size_of_buildLog);
	buildLog = (char *)malloc(size_of_buildLog + 1);
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
			size_of_buildLog, buildLog, NULL);
	buildLog[size_of_buildLog] = '\0';
	if (Monitor > TMC_PRINT_INFO) {
		printf("Kernel buildLog: %s\n", buildLog);
	}
	if (retval != CL_SUCCESS)
	{
		printf("Error building program in tmr_ocl_create_kernel\n");
#ifdef TUNING
		//fprintf(resf,"Kernel buildLog: %s", buildLog);
		fprintf(resf, "Kernel error - too much shared data\n");
		fclose(resf);
#endif
		exit(-1);
	}


#ifdef TUNING

	if (tmr_ocl_check_vendor() == 1)
	{
		char * s;
		int t = 0;
		int stack;
		int stores;
		int loads;
		int registers;
		int smem;
		int cmem;
		int cmem2;
		int nr;

		s = strstr(buildLog, "stack");
		if (s != NULL)
		{
			t = (int)(s - buildLog);
			//printf("Found string at index = %d\n", t );
			sscanf(&buildLog[t - 10], "%d", &stack);
			//printf("stack=%d\n",stack);
		}
		else
		{
			//printf("String not found\n");
			stack = 0;
		}
		if (line_count == 0)
		{
			fprintf(headf, "bytes stack frame,");
		}
		fprintf(resf, "%d,", stack);
		s = strstr(buildLog, "bytes spill stores");
		if (s != NULL)
		{
			t = (int)(s - buildLog);
			//printf("Found string at index = %d\n", t );
			sscanf(&buildLog[t - 10], "%*s %d", &stores);
			//printf("stores=%d\n",stores);
		}
		else
		{
			//printf("String not found\n");
			stores = 0;
		}
		if (line_count == 0)
		{
			fprintf(headf, "bytes spill stores,");
		}
		fprintf(resf, "%d,", stores);
		s = strstr(buildLog, "bytes spill loads");
		if (s != NULL)
		{
			t = (int)(s - buildLog);
			//printf("Found string at index = %d\n", t );
			sscanf(&buildLog[t - 10], "%*s %d", &loads);
			//printf("loads=%d\n",loads);
		}
		else
		{
			//printf("String not found\n");
			loads = 0;
		}
		if (line_count == 0)
		{
			fprintf(headf, "bytes spill loads,");
		}
		fprintf(resf, "%d,", loads);
		s = strstr(buildLog, "registers");
		if (s != NULL)
		{
			t = (int)(s - buildLog);
			//printf("Found string at index = %d\n", t );
			sscanf(&buildLog[t - 7], "%*s %d", &registers);
			//printf("registers=%d\n",registers);
		}
		else
		{
			//printf("String not found\n");
			registers = 0;
		}
		if (line_count == 0)
		{
			fprintf(headf, "registers,");
		}
		fprintf(resf, "%d,", registers);
		s = strstr(buildLog, "smem");
		if (s != NULL)
		{
			t = (int)(s - buildLog);
			//printf("Found string at index = %d\n", t );
			sscanf(&buildLog[t - 20], "%*s %d", &smem);
			//printf("smem=%d\n",smem);
		}
		else
		{
			//printf("String not found\n");
			smem = 0;
		}
		if (line_count == 0)
		{
			fprintf(headf, "bytes smem,");
		}
		fprintf(resf, "%d,", smem);
		s = strstr(buildLog, "cmem[0]");
		if (s != NULL)
		{
			t = (int)(s - buildLog);
			//printf("Found string at index = %d\n", t );
			sscanf(&buildLog[t - 10], "%d", &cmem);
			//printf("cmem[0]=%d\n",cmem);
		}
		else
		{
			//printf("String not found\n");
			cmem = 0;
		}
		if (line_count == 0)
		{
			fprintf(headf, "bytes cmem[0],");
		}
		fprintf(resf, "%d,", cmem);
		s = strstr(&buildLog[t + 9], "cmem[");
		if (s != NULL)
		{
			t = (int)(s - buildLog);
			//printf("Found string at index = %d\n", t );
			sscanf(&buildLog[t - 10], "%d", &cmem2);
			sscanf(&buildLog[t + 5], "%d", &nr);
			//printf("cmem[%d]=%d\n",nr,cmem2);
		}
		else
		{
			//printf("String not found\n");
			nr = -1;
			cmem2 = 0;
		}
		if (line_count == 0)
		{
			fprintf(headf, "bytes cmem[%d],", nr);
		}
		fprintf(resf, "%d,", cmem2);

	}//end nvidia

#endif

	free(buildLog);

	// Create OpenCL kernel
	cl_kernel kernel = clCreateKernel(program, Kernel_name, NULL);
	if (kernel == NULL)
	{
		printf("Failed to create kernel.\n");
		exit(-1);
	}

	if (Monitor > TMC_PRINT_INFO) {
		printf("Created kernel for platform %d, device %d, kernel index %d\n",
				Platform_index, Device_index, Kernel_index);
	}

#ifdef TUNING

	if (tmr_ocl_check_vendor() == 1)
	{
		system("mkdir kernele");
		char name[50];
		FILE * fp;
		size_t bin_sz;
		int word;
		for (i = 0; i < 50; i++)
			name[i] = 0;

		strcat(name, "kernele/");
		for (i = 0; i < NOPT; i++)
		{
			if (options[i] == 1)
				strcat(name, "1");
			else
				strcat(name, "0");
		}

		word = wc("tmr_ocl_num_int_el.cl", "//QSS");
		if (word == 1)
			strcat(name, "_QSS");
		word = wc("tmr_ocl_num_int_el.cl", "//SQS");
		if (word == 1)
			strcat(name, "_SQS");
		word = wc("tmr_ocl_num_int_el.cl", "//SSQ");
		if (word == 1)
			strcat(name, "_SSQ");

		strcat(name, "_kernel.ptx");

		printf("Kernel: %s\n", name);

		clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin_sz, NULL);

		// Read binary (PTX file) to memory buffer
		unsigned char * bin = (unsigned char *)malloc(bin_sz);
		clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &bin, NULL);

		printf("TU!!!\n");


		// Save PTX to add_vectors_ocl.ptx
		fp = fopen(name, "wb");
		fwrite(bin, sizeof(char), bin_sz, fp);
		fclose(fp);
		free(bin);


		word = wc(name, "add.f64");
		if (line_count == 0)
		{
			fprintf(headf, "add.f64,");
		}
		fprintf(resf, "%d,", word);
		//printf("add.f64=%d\n",word);
		word = wc(name, "sub.f64");
		if (line_count == 0)
		{
			fprintf(headf, "sub.f64,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "mul.f64");
		if (line_count == 0)
		{
			fprintf(headf, "mul.f64,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "fma.rn.f64");
		if (line_count == 0)
		{
			fprintf(headf, "fma.rn.f64,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "neg.f64");
		if (line_count == 0)
		{
			fprintf(headf, "neg.f64,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "rcp.rn.f64");
		if (line_count == 0)
		{
			fprintf(headf, "rcp.rn.f64,");
		}
		fprintf(resf, "%d,", word);

		word = wc(name, "global.f64");
		if (line_count == 0)
		{
			fprintf(headf, "global.f64,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "shared.f64");
		if (line_count == 0)
		{
			fprintf(headf, "shared.f64,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "const.f64");
		if (line_count == 0)
		{
			fprintf(headf, "const.f64,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "local.f64");
		if (line_count == 0)
		{
			fprintf(headf, "local.f64,");
		}
		fprintf(resf, "%d,", word);


	}
	else if (tmr_ocl_check_vendor() == 2)
	{
		system("mkdir kernele");
		char name[50];
		int word;
		for (i = 0; i < 50; i++)
			name[i] = 0;

		strcat(name, "kernele/");
		for (i = 0; i < NOPT; i++)
		{
			if (options[i] == 1)
				strcat(name, "1");
			else
				strcat(name, "0");
		}

		word = wc("tmr_ocl_num_int_el.cl", "//QSS");
		if (word == 1)
			strcat(name, "_QSS");
		word = wc("tmr_ocl_num_int_el.cl", "//SQS");
		if (word == 1)
			strcat(name, "_SQS");
		word = wc("tmr_ocl_num_int_el.cl", "//SSQ");
		if (word == 1)
			strcat(name, "_SSQ");

#ifdef OPENCL_CPU
		strcat(name, "_kernel_CPU.ptx");
#endif
#ifdef OPENCL_PHI
		strcat(name, "_kernel_PHI.ptx");
#endif


		printf("Kernel: %s\n", name);

		char compile[500];

		for (i = 0; i < 500; i++)
			compile[i] = 0;
		strcat(compile, "ioc64 -cmd=build -input=tmr_ocl_num_int_el.cl -device=");
#ifdef OPENCL_CPU
		strcat(compile, "cpu -asm=");
#endif
#ifdef OPENCL_PHI
		strcat(compile, "co -asm=");
#endif
		strcat(compile, name);
		strcat(compile, " -bo=\"");
		strcat(compile, buffer);
		strcat(compile, "\"");

		printf("Compile: %s\n", compile);

		system(compile);
//PHI
		word = wc(name, "vsubrpd");
		if (line_count == 0)
		{
			fprintf(headf, "vsubrpd,");
		}
		fprintf(resf, "%d,", word);
		//printf("add.f64=%d\n",word);
		word = wc(name, "vfmadd213pd");
		if (line_count == 0)
		{
			fprintf(headf, "vfmadd213pd,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "vfnmadd231pd");
		if (line_count == 0)
		{
			fprintf(headf, "vfnmadd231pd,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "vfmadd231pd");
		if (line_count == 0)
		{
			fprintf(headf, "vfmadd231pd,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "vfmadd132pd");
		if (line_count == 0)
		{
			fprintf(headf, "vfmadd132pd,");
		}
		fprintf(resf, "%d,", word);
//common
		word = wc(name, "vmulpd");
		if (line_count == 0)
		{
			fprintf(headf, "vmulpd,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "vaddpd");
		if (line_count == 0)
		{
			fprintf(headf, "vaddpd,");
		}
		fprintf(resf, "%d,", word);

//cpu
		word = wc(name, "vsubsd");
		if (line_count == 0)
		{
			fprintf(headf, "vsubsd,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "vmulsd");
		if (line_count == 0)
		{
			fprintf(headf, "vmulsd,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "vaddsd");
		if (line_count == 0)
		{
			fprintf(headf, "vaddsd,");
		}
		fprintf(resf, "%d,", word);
		word = wc(name, "vsubpd");
		if (line_count == 0)
		{
			fprintf(headf, "vsubpd,");
		}
		fprintf(resf, "%d,", word);

	}
	else if (tmr_ocl_check_vendor() == 3)
	{

		//spectre
		/*
		v_add_f64
		v_mul_f64
		v_div_scale_f64
		v_rcp_f64
		v_fma_f64
		v_div_fmas_f64
		v_div_fixup_f64
		*/

	}



#endif


	platform_struct.list_of_devices[Device_index].program[Kernel_index] = program;
	platform_struct.list_of_devices[Device_index].kernel[Kernel_index] = kernel;

	return (1);
}



#ifdef __cplusplus
}
#endif
