/************************************************************************
File uts_util - utility routines (common to all problem modules,
				possibly used also by other modules)

Contains definitions of routines:
  utr_set_interactive - to set up names for interactive input and output
	that control the run (possibly stdin and/or stdout for on-line control)
  utr_initialize_mesh - to initialize mesh of a specified type
  utr_initialize_field - to initialize field of a specified type
  utr_export_mesh - to export mesh of a specified type in a specified format
  utr_export_field - to export field of a specified type in a specified format

  utr_get_ent_dofs - to read geometry or solution DOFs for a given DOF entity
				 from suitable (approximation or mesh) data structures

  utr_average_sol_el - to compute the average of solution over element

  utr_dvector - to allocate a double vector: name[0..ncom-1]:
				  name=ut_dvector(ncom,error_text)
  utr_ivector - to allocate an integer vector: name[0..ncom-1]:
				  name=ut_ivector(ncom,error_text)
  utr_d_zero - to zero a double vector of length Num
  utr_chk_list - to check whether Num is on the list List
	with length Ll
  utr_put_list - to put Num on the list List with length Ll

	(filled with numbers and zeros at the end)
  utr_mat3_inv - to invert a 3x3 matrix (stored as a vector!)
  utr_vec3_prod - to compute vector product of 3D vectors
  utr_vec3_mxpr - to compute mixed vector product of 3D vectors
  utr_vec3_length - to compute length of a 3D vector
  utr_mat3vec - to compute matrix vector product in 3D space
  utr_mat3mat - to compute matrix matrix product in 3D space
	(all matrices are stored by rows as vectors!)
  utr_mat_det - to compute determinant of a matrix
  utr_skip_rest_of_line - to allow for comments in input files

------------------------------
History:
	08.2008 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<signal.h>

#ifdef MULTITHREADED
#include <omp.h>
#endif

#include <modfem/svnversion.h>
#include <modfem/gitversion.h>
#include <modfem/mf_version.h>
//#include "mf_version.h"

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>

/* interface for all approximation modules */
#include <modfem/aph_intf.h>

#ifdef PARALLEL
/* interface for parallel mesh manipulation modules */
#include "mmph_intf.h"

/* interface with parallel communication library */
#include "pch_intf.h"
#endif

#include <modfem/uth_log.h>
#include <modfem/uth_mesh.h>
#include <modfem/uth_io_results.h>
#include <modfem/uth_io_compression.h>
#include <modfem/uth_system.h>
#include <modfem/uth_io_files.h>

/* interface with multithreading management library */
#include <modfem/tmh_intf.h>

/* interface for general purpose utilities - for all problem dependent modules*/
#include <modfem/uth_intf.h>

/* interface for linear algebra packages */
#include <modfem/lin_alg_intf.h>

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>

#include <modfem/pdh_intf.h>

utt_time_measurement utv_time;


// maximal number of solution components
#define UTC_MAXEQ PDC_MAXEQ

/****************************************/
/* Ctrl+C handling functtion            */
/****************************************/
// flag indicating if SIGINT (ctrl-C) is caught
int utv_SIGINT_not_caught = 1;
void utr_ctrl_c_signal_handler(int param)
{
	char c = 0;

	signal(SIGINT, SIG_DFL);

	mf_log_warn("\nCought Ctrl+C signal!\nPress key and confirm it by pressing [Enter]:\n[c] or [q] to quit now\n[m] to finalize current step and display menu\n[r] or any other key to resume(default) : ");
	c = getchar();
	switch (c) {
		case 3: // Ctrl-C code
		case '^':
		case 'q':
		case 'c': {
			mf_log_info("\nExiting application!\n");
			exit(-1);
		}
		break;
		case 'm': {
			utv_SIGINT_not_caught = 0;
			mf_log_info("\nFinalizing current progress and entering menu...");
		}
		break;
		case 27: // esc
		default: {
			mf_log_info("\nResuming execution...");
			utv_SIGINT_not_caught = 1;
		}
		break;
	} //!switch

	signal(SIGINT, utr_ctrl_c_signal_handler);

	return;
}


/**--------------------------------------------------------
  utr_reset_time_counters
----------------------------------------------------------*/
void utr_reset_time_counters( int Control)
{

	utv_time.initialize_mesh = 0.0;
	utv_time.initialize_field = 0.0;
	utv_time.decompress_files = 0.0;

	utv_time.export_mesh = 0.0;
	utv_time.export_field = 0.0;
	utv_time.compress_files = 0.0;

	utv_time.create_int_ent_list = 0.0;
	utv_time.create_dof_ent_list = 0.0;

	utv_time.renumber_dofs = 0.0;
	utv_time.total_coloring = 0.0;
	utv_time.color_elements = 0.0;
	utv_time.color_faces = 0.0;

	utv_time.total_integration_and_assembly = 0.0;
	utv_time.handle_solver_interface_data_structures = 0.0;
	utv_time.handle_solver_interface_GPU_data_structures = 0.0;
	utv_time.create_solver_data_structures = 0.0;

	utv_time.total_sir_create = 0.0;

	utv_time.read_geo_dofs_vector = 0.0;
	utv_time.read_single_sol_dofs_vector = 0.0;
	utv_time.read_all_sol_dofs_vectors = 0.0;
	utv_time.CPU_integration_elems = 0.0;
	utv_time.CPU_assembly_faces = 0.0;
	utv_time.CPU_integration_faces = 0.0;
	utv_time.CPU_assembly_faces = 0.0;
	utv_time.GPU_integration_and_assembly = 0.0;
	utv_time.CPU_GPU_solver_ds_merging = 0.0;
	utv_time.create_sm_uncon = 0.0;
	utv_time.adding_constraints_info = 0.0;

	utv_time.create_patches = 0.0;
	utv_time.recover_derivatives = 0.0;
	utv_time.calculate_error = 0.0;
	utv_time.refine_mesh = 0.0;
	utv_time.project_field = 0.0;
	utv_time.manual_refinemts = 0.0;
	utv_time.total_outside_solver = 0.0;

	utv_time.read_initial_guess = 0.0;
	utv_time.solve = 0.0;
	utv_time.write_solution = 0.0;
	utv_time.total_sir_solve = 0.0;

	return;
}

/**--------------------------------------------------------
  utr_display_time_counters
----------------------------------------------------------*/
void utr_display_time_counters( int Control) {

	printf("\n###<<<---***---!!!--- Final performance data report begin ---!!!---***--->>>###\n");
	printf(">>> Total time from beginning of simulations: \t%lf <<<<<<<<<<<<<<<<<<<<<<\n",
			utv_time.total);


	printf("\n>>>>>>>>>>>>>>> Before main menu loop:\t\t\t\t\t<<<<<<<<\n");

	printf("Time for mesh initialization: \t\t\t\t\t%lf\n",
			utv_time.initialize_mesh);
	printf("Time for field initialization: \t\t\t\t\t%lf\n",
			utv_time.initialize_field);
	printf("\tincluding time for decompressing files: \t\t   %lf\n",
			utv_time.decompress_files);
	printf(">>>>>>>>>>>>>>> Total time for problem initialization: \t%lf <<<<<<<<<<<<<<<\n",
			utv_time.problem_initialization);

	printf("\n>>>>>>>>>>>>>>> Inside sir_create:\t\t\t\t\t<<<<<<<<\n");

	//printf("Time for creating int_ent list: \t\t\t\t%lf\n", utv_time.create_int_ent_list);
	//printf("Time for creating dof_ent list: \t\t\t\t%lf\n", utv_time.create_dof_ent_list);
	printf("Time for renumbering DOFs: \t\t\t\t\t%lf\n",
			utv_time.renumber_dofs);

	printf("\nTotal time for coloring: \t\t\t\t\t%lf\n",
			utv_time.total_coloring);
	printf("\tincluding time for coloring elements: \t\t\t   %lf\n",
			utv_time.color_elements);
	printf("\tincluding time for coloring faces: \t\t\t   %lf\n",
			utv_time.color_faces);

	printf("\nTotal time for handling solver interface data structures:\t%lf\n",
			utv_time.handle_solver_interface_data_structures);
	printf("\t(dof_struct, local_to_global, assembly_table)\n");
	printf("Total time for handling solver interface GPU data structures:\t%lf\n",
			utv_time.handle_solver_interface_GPU_data_structures);
	printf("\tincluding time for reading geo_dofs vector:\t\t   %lf\n",
			utv_time.read_geo_dofs_vector);
	printf("\tincluding time for reading single sol_dofs vector:\t   %lf\n",
			utv_time.read_single_sol_dofs_vector);
	printf("\tincluding time for reading all sol_dofs vectors:\t   %lf\n",
			utv_time.read_all_sol_dofs_vectors);

	printf("Total time for creating solver data structures:\t\t\t%lf\n",
			utv_time.create_solver_data_structures);

	printf(">>>>>>>>>>>>>>> Total time in sir_create: \t\t%lf <<<<<<<<<<<<<<<\n",
			utv_time.total_sir_create);

	printf("\n>>>>>>>>>>>>>>> Inside sir_solve:\t\t\t\t\t<<<<<<<<\n");

	printf("Time for reading initial guess:\t\t\t\t\t%lf\n",
			utv_time.read_initial_guess);
	printf("Total time for integration and assembly: \t\t\t%lf\n",
			utv_time.total_integration_and_assembly);
	printf("\tincluding concurrent time for GPU integration and assembly:\t   %lf\n",
			utv_time.GPU_integration_and_assembly);
	printf("\tincluding time for merging CPU and GPU linear system arrays:\t   %lf\n",
			utv_time.CPU_GPU_solver_ds_merging);
	printf("\tincluding time (sum over all threads) for:\n");
	printf("\t\tCPU elements integration: \t\t\t   %lf\n",
			utv_time.CPU_integration_elems);
	printf("\t\tCPU elements assembly: \t\t\t\t   %lf\n",
			utv_time.CPU_assembly_elems);
	printf("\t\tCPU faces integration: \t\t\t\t   %lf\n",
			utv_time.CPU_integration_faces);
	printf("\t\tCPU faces assembly: \t\t\t\t   %lf\n",
			utv_time.CPU_assembly_faces);
	printf("\t\t\tcreating SMs and LVs unconstrained:\t      %lf\n",
			utv_time.create_sm_uncon);
	printf("\t\t\tadding lists of real dof_ent:\t\t      %lf\n",
			utv_time.adding_constraints_info);



	printf("Time for solving the system of linear equations:\t\t%lf\n",
			utv_time.solve);
	printf("Time for writing solution to data structure:\t\t\t%lf\n",
			utv_time.write_solution);
	printf(">>>>>>>>>>>>>>> Total time in sir_solve: \t\t%lf <<<<<<<<<<<<<<<\n",
			utv_time.total_sir_solve);


	printf("\n>>>>>>>>>>>>>>> In main menu loop (outside solver)\t\t\t<<<<<<<<\n");

	printf("Time for calculating error: \t\t\t\t\t%lf\n",
			utv_time.calculate_error);
	printf("Time for exporting mesh: \t\t\t\t\t%lf\n",
			utv_time.export_mesh);
	printf("Time for exporting field: \t\t\t\t\t%lf\n",
			utv_time.export_field);
	printf("\tincluding time for compressing files: \t\t\t   %lf\n",
			utv_time.compress_files);

	printf(">>>>>>>>>>>>>>> Total time in main menu actions:\t%lf <<<<<<<<<<<<<<<\n",
			utv_time.total_outside_solver);

	printf("\n###<<<---***---!!!---- Final performance data report end ----!!!---***--->>>###\n");


	return;
}

/*---------------------------------------------------------
  utr_set_interactive - to set up names for interactive input and output
	that control the run (possibly stdin and/or stdout for on-line control)
----------------------------------------------------------*/
int utr_set_interactive(
		char * Work_dir,
		int argc,
		char ** argv,
		FILE ** interactive_input_p,
		FILE ** interactive_output_p
)
{

	time_init();
	char interactive_input_name[300];
	char interactive_output_name[300], tmp[100];
	FILE * interactive_input = NULL;
	FILE * interactive_output = NULL;
#ifdef PARALLEL
	int nr_proc, my_proc_id;
	int info;
#endif

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	utr_reset_time_counters(0);

	/* initial setting for interactive_input and interactive_output */
	sprintf(interactive_input_name, "%s/input_interactive.txt", Work_dir);
	/*kbw
	  printf("input file: %s\n", interactive_input_name);
	//kew*/
	interactive_input = fopen(interactive_input_name, "r");
	if (interactive_input == NULL) {

		interactive_input = stdin;
		interactive_output = stdout;
		/****************************************/
		/* Ctrl+C signal handling (only stdin)  */
		//signal(SIGINT,utr_ctrl_c_signal_handler);
		/****************************************/

	} else {
		fscanf(interactive_input, "%s\n", tmp);
		if (strcmp(tmp, "STDOUT") == 0 ||
				strcmp(tmp, "stdout") == 0) interactive_output = stdout;
		else {
			/*kbw
			  printf("output file: %s\n", interactive_output_name);
			//kew*/
			sprintf(interactive_output_name, "%s/%s", Work_dir, tmp);
			//interactive_output = fopen(interactive_output_name,"w");
		}
	}
	/*kbw
	  printf("output file: %s\n", interactive_output_name);
	//kew*/

#ifdef PARALLEL
	/* first initiate parallel execution */
	info = pcr_init_parallel(&argc, argv, Work_dir,
					interactive_output_name, &interactive_output,
					&nr_proc, &my_proc_id);
	/*kbw
	  printf("output file: %s\n", interactive_output_name);
	//kew*/

	if (info) {
		fprintf(interactive_output,
				"Unable to initiate parallel communication module! Exiting.\n");
		exit(0);
	}
	/*ok_kbw*/
	fprintf(interactive_output,
			"\nInitiated parallel communication: nr_proc %d, my_proc_id %d\n",
			nr_proc, my_proc_id);
//kew*/
	if (my_proc_id == pcr_print_master()) {
		if (interactive_input == stdin) {
			printf("\nWARNING! Input from STDIN in PARALLEL mode!\n");
			printf("May lead to deadlock in batch execution!\n");
			printf("Press any key to continue.\n");
			getchar();
			getchar();
		}
	}
#endif

#ifndef PARALLEL
	if (interactive_output != stdout) {
		interactive_output = fopen(interactive_output_name, "w");
		/*kbw
		  printf("output file: %s\n", interactive_output_name);
		//kew*/
	}
#endif

	if (interactive_input == NULL || interactive_output == NULL) {
		fprintf(stderr,
				"Cannot establish interactive input and/or interactive output. Exiting\n");
		exit(0);
	}
	// Setting global pointer to output file for logging purposes.
	utv_log_out = interactive_output;
	//

	if (strcmp(Work_dir, ".") == 0) {
		mf_log_info("No special working directory specified in command line");
		mf_log_info("Current directory: %s , assumed as working directory", Work_dir);
	} else {
		mf_log_info("Working directory specified in command line: %s", Work_dir);
	}

	if (interactive_input == stdin) {
		mf_log_info("Interactive input from terminal (stdin)\nInteractive output to stdout");
	} else {
		if (interactive_output == stdout) {
			mf_log_info("Interactive input from file: %s\nInteractive output to stdout",
					interactive_input_name);
		} else {
			mf_log_info("Interactive input from file: %s",
					interactive_input_name);
			mf_log_info("Interactive output to file: %s",
					interactive_output_name);
		}
	}

	*interactive_input_p = interactive_input;
	*interactive_output_p = interactive_output;


#ifdef MULTITHREADED

	omp_set_nested(0);

	#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
		{
			mf_log_info("OpenMP multithreading enabled with %d threads!", omp_get_num_threads());
			if (omp_get_num_threads() == 1) {
				// 1 thread is absolutely correct and used for testing !!!!!!!
			}
		}
	}


#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)

	int monitor = TMC_PRINT_INFO + 1;
// control not used - decided at compile time !!!!!
	int control = -1;
	tmr_init_multithreading_opencl(Work_dir, argc, argv,
			interactive_input, interactive_output,
			control, monitor);
	mf_log_info("\nOpenCL multithreading enabled!");

#endif // if OpenMP + OpenCL

#endif // if OpenMP alone

	//mf_log_info("ModFEM ver. %s(%d)", MF_GITVERSION , SVNVERSION);
	mf_log_info("ModFEM ver. %s.%s (%s)", MF_VERSION, GITVERSION, __DATE__);
	//mf_log_info("ModFEM ver. %s.%s (%s)", MF_VERSION, GITVERSION, __DATE__);
	//mf_log_info("ModFEM ver. %s (%s)", MF_GITVERSION,__DATE__);

	return (1);
}


/*---------------------------------------------------------
utr_initialize_mesh - to initialize mesh of a specified type
---------------------------------------------------------*/
int utr_initialize_mesh( /* returns mesh_id */
		FILE * Interactive_output, /* file or stdout to write messages */
		const char * Work_dir, // path to working directory
		char Mesh_type, /* letter symbol denoting mesh_type (j, p, t or h) */
		const char * Mesh_file /* mesh file name - conforming to naming convention */
)
{


	return utr_io_initialize_mesh(Interactive_output, Work_dir, Mesh_type, Mesh_file);
}


/*---------------------------------------------------------
utr_initialize_field - to initialize field of a specified type
---------------------------------------------------------*/
int utr_initialize_field(  /* returns: >0 - field ID, <0 - error code */
		FILE * interactive_output, /* file or stdout to write messages */
		char Field_type, /* options: s-standard, c-discontinuous with complete basis */
		/*          t - discontinuous with tensor product basis */
		char Control,	 /* in: control variable: */
		/*      z - to initialize the field to zero */
		/*      r - to read field values from the file */
		/*      i - to initialize the field using function */
		/*                 provided by the problem dependent module */
		int Mesh_id,	 /* in: ID of the corresponding mesh */
		int Nreq,	 /* in: number of equations - solution vector components */
		int Nr_sol,	 /* in: number of solution vectors for each dof entity */
		int Pdeg,      /* in: degree of approximating polynomial */
		char * Filename, /* in: name of the file to read approximation data */
		double (*Fun_p)(int, double *, int) /* pointer to function that provides */
		/* problem dependent initial condition data */
)
{

	int field_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	double time_begin = time_clock();

	/*kbw
	  printf("Initializing field: type %c, control %c\n", Field_type, Control);
	  printf("Mesh_id %d, Nreq %d, Nr_sol %d, Pdeg %d\n",
		 Mesh_id, Nreq, Nr_sol, Pdeg);
	//kew*/

	if ( (Filename != NULL)  && (Control == 'r') ) {
		char decompressed_file[255] = {0};
		utr_io_decompress_file(NULL, Filename, decompressed_file);
		strcpy(Filename, decompressed_file);
	}

	/* reading field values saved by previous runs */
	if (Control == 'r') {
		field_id = apr_init_field(Field_type, APC_READ, Mesh_id, Nreq, Nr_sol, Pdeg, Filename, NULL);
	}
	/* initializing field values according to the specified function */
	else if (Control == 'i') {
		field_id = apr_init_field(Field_type, APC_INIT, Mesh_id, Nreq, Nr_sol, Pdeg, NULL, Fun_p);
	}
	/* initializing field values to zero */
	else {
		field_id = apr_init_field(Field_type, APC_ZERO, Mesh_id, Nreq, Nr_sol, Pdeg, NULL, NULL);
	}

	#pragma omp atomic
	utv_time.initialize_field += time_clock() - time_begin;

	return (field_id);
}

/**--------------------------------------------------------
utr_export_field - to export mesh of a specified type
---------------------------------------------------------*/
int utr_export_field(/** returns number of files exported */
		int Field_id,    /* in: field ID */
		int Nreq,        /* in: number of equations (scalar dofs) */
		int Select,      /* in: parameter to select written vectors */
		int Accuracy,    /* in: parameter specyfying accuracy - significant digits */
		/* (put 0 for full accuracy in "%g" format) */
		char * Filename  /* in: name of the file to write field data */
)
{

	double time_begin = time_clock();

	mf_check_mem(Filename);

	char filename_correct[1000] = {0};

#ifdef PARALLEL
	sprintf(filename_correct, "%s_%d", Filename, pcv_my_proc_id);
#else
	sprintf(filename_correct, "%s", Filename);
#endif // PARALLEL

	int nf =  apr_write_field(Field_id, Nreq, Select, Accuracy, filename_correct);

#ifndef PARALLEL
	utr_io_compress_file(NULL, filename_correct);
#endif

	#pragma omp atomic
	utv_time.export_field += time_clock() - time_begin;

	return nf;
}

/*---------------------------------------------------------
utr_export_mesh - to export mesh of a specified type in a specified format
---------------------------------------------------------
int utr_export_mesh(
  FILE* Interactive_output,
  int Mesh_id,
  char * Filename,
  char Format)
{
  int result=1;

  switch(Format)
	{
	case 'p':
	  {
	char arg[255]={0};
	sprintf(arg,"%s.vtu",Filename);
	result=utr_write_paraview_mesh(Mesh_id,arg);
	fprintf(Interactive_output,"Writing file %s\n",arg);
	  }
	  break;
	}
  return result;
}*/

/*---------------------------------------------------------
utr_export_field - to export field of a specified type in a specified format
---------------------------------------------------------
int utr_export_field(
  FILE* Interactive_output,
  int Mesh_id,
  int Field_id,
  const char * Filename,
  char Format,
  const char ** Desc)
{
  int result=1;
  static const char * defaultDesc[2]={"SCALARS pressure double 1\n","VECTORS velocity double\n"};
  if(Desc == NULL) {
	Desc=defaultDesc;
  }

  switch(Format)
	{
	case 'p':
	  {
	char arg[255]={0};
	sprintf(arg,"%s.vtu",Filename);
	result=utr_write_paraview_field(Mesh_id,Field_id,arg,Desc);
	fprintf(Interactive_output,"Writing file %s\n",arg);
	  }
	  break;
	}
  return result;
}*/


/*---------------------------------------------------------
  utr_get_ent_dofs - to read geometry or solution DOFs for a given DOF entity
				 from suitable (approximation or mesh) data structures
----------------------------------------------------------*/
int utr_get_ent_dofs( /* returns: 0 - success, -1 - no such data structure */
		int Problem_id,     /** in: solver ID (used to identify the subproblem) */
		int Dof_ent_type,
		int Dof_ent_id,
		int Nrdofs,
		int Vect_id,       /** in: vector ID to select one of solution vectors */
		/*      or geometrical DOFs */
		double * Vect_dofs /** in: dofs to be read */
)
{

	int field_id = pdr_ctrl_i_params(Problem_id, 3);
	int mesh_id = apr_get_mesh_id(field_id);
	//int mesh_id = pdr_ctrl_i_params(Problem_id, 2);
	int nr_sol =  pdr_ctrl_i_params(Problem_id, 4);

	if (Vect_id == Geometry_dofs_ID) {

		// current version for geometrically multilinear elements only
		// (geo dofs only for vertices == nodes)
		mmr_node_coor(mesh_id, Dof_ent_id, Vect_dofs);

	}
	else {

		if (Vect_id > nr_sol) return (-1);

		apr_read_ent_dofs(field_id, Dof_ent_type, Dof_ent_id,
				Nrdofs, Vect_id, Vect_dofs);
	}

	return (0);
}

/*---------------------------------------------------------
utr_average_sol_el - to compute the average of solution over element
----------------------------------------------------------*/
double utr_average_sol_el( /* returns: the average of solution over element */
		int Field_id,      /* in: data structure to be used  */
		int El          /* in: element number */
)
{

	/* local variables */
	int base_q;           /* type of basis functions */
	int el_type;          /* element type */
	int nreq;             /* number of equations */
	int pdeg;             /* degree of polynomial */
	int num_shap;         /* number of element shape functions */
	int ndofs;            /* local dimension of the problem */
	int ngauss;           /* number of gauss points */
	double xg[3000];      /* coordinates of gauss points in 3D */
	double wg[1000];      /* gauss weights */
	double determ;        /* determinant of jacobi matrix */
	double vol;           /* volume for integration rule */
	double xcoor[3];      /* global coord of gauss point */
	double u_val[PDC_MAXEQ]; /* computed solution */
	double base_phi[APC_MAXELVD];    /* basis functions */
	int el_nodes[APC_MAX_GEO_DOFS + 1];      /* list of nodes of El */
	int el_nodes_type[APC_MAX_GEO_DOFS + 1];      /* list of nodes type of El */
	double node_coor[3 * APC_MAX_GEO_DOFS]; /* coord of nodes of El */
	double dofs_loc[APC_MAXELSD]; /* element solution dofs */

	/* auxiliary variables */
	int mesh_id;
	int i, ki, iaux;
	double volume, average;

	/*++++++++++++++++ executable statements ++++++++++++++++*/


	mesh_id = apr_get_mesh_id(Field_id);
	nreq = apr_get_nreq(Field_id);
	base_q = apr_get_base_type(Field_id, El);

	/* find element type */
	el_type = mmr_el_type(mesh_id, El);

	/* find degree of polynomial and number of element scalar dofs */
	apr_get_el_pdeg(Field_id, El, &pdeg);
	num_shap = apr_get_el_pdeg_numshap(Field_id, El, &pdeg);
	ndofs = num_shap * nreq;

	/* get the coordinates of element nodes in the right order */
	//mmr_el_node_coor(mesh_id,El,el_nodes,node_coor);
	apr_get_el_geo_dofs(Field_id, El, el_nodes, el_nodes_type, node_coor);

	/* get the current solution degrees of freedom */
	i = 1; // sol_vec_id
	apr_get_el_dofs(Field_id, El, i, dofs_loc);

	/* prepare data for gaussian integration */
	/*kb!!!*/
#ifdef DEBUG
	printf("Average computed for pdeg=303\n");
#endif

	if (base_q == APC_BASE_COMPLETE_DG) pdeg = 3;
	else pdeg = 303;
	apr_set_quadr_3D(base_q, &pdeg, &ngauss, xg, wg);
	/*kb!!!*/

	/*kbw
	  printf("In err-comp for element %d, type %d\n",El,el_type);
	  printf("pdeg %d, ngauss %d\n",pdeg,ngauss);
	  printf("nreq %d, ndof %d, local_dim %d\n",nreq,num_shap,ndofs);
	  printf("%d nodes with coordinates:\n",el_nodes[0]);
	  for(i=0;i<el_nodes[0];i++){
	  printf("node %d: x - %f, y - %f, z - %f\n",
	  el_nodes[i+1],node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
	  }
	  printf("solution dofs:\n");
	  for(i=0;i<ndofs;i++){
	  printf("%20.12lf",dofs_loc[i]);
	  }
	  printf("\n");
	  getchar();
	  kew*/

	volume = 0.0;
	average = 0.0;
	for (ki = 0; ki < ngauss; ki++) {

		/* at the gauss point, compute basis functions, determinant etc*/
		iaux = 2; /* calculations with jacobian but not on the boundary */

		determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
						&xg[3 * ki], node_coor, dofs_loc,
						base_phi, NULL, NULL, NULL,
						xcoor, u_val, NULL, NULL, NULL, NULL);
		vol = determ * wg[ki];

		/*kbw
		  printf("at gauss point %d, local coor %lf, %lf, %lf\n",
		  ki,xg[3*ki],xg[3*ki+1],xg[3*ki+2]);
		  printf("weight %lf, determ %lf, coeff %lf\n",
		  wg[ki],determ,vol);
		  printf("u_val = %f\n", u_val[0]);
		  kew*/

		volume += vol;
		average += u_val[0] * vol;

	} /* ki */

	/*kbw
	  printf("Average in element %d = %lf (at center = %lf)\n",
	  EL,average/volume,(dofs_loc[0]+dofs_loc[1]+dofs_loc[2])/3.0);
	  kew*/

	return (average / volume);

}


/*---------------------------------------------------------
utr_dvector - to allocate a double vector: name[0..ncom-1]:
				  name=ut_dvector(ncom,error_text)
---------------------------------------------------------*/
double * utr_dvector( /* return: pointer to allocated vector */
		int ncom,  	/* in: number of components */
		char error_text[]/* in: error text to be printed */
)
{

	double * v;

	v = (double *) malloc (ncom * sizeof(double));
	if (!v) {
		printf("Not enough space for allocating vector: %s ! Exiting\n", error_text);
		exit(1);
	}
	return v;
}


/*---------------------------------------------------------
  utr_ivector - to allocate an integer vector: name[0..ncom-1]:
  name=ut_ivector(ncom,error_text)
---------------------------------------------------------*/
int * utr_ivector(   /* return: pointer to allocated vector */
		int ncom, 	/* in: number of components */
		char error_text[]/* in: error text to be printed */
)
{

	int * v;

	v = (int *) malloc (ncom * sizeof(int));
	if (!v) {
		printf("Not enough space for allocating vector: %s ! Exiting\n", error_text);
		exit(1);
	}
	return v;
}

/*---------------------------------------------------------
utr_d_zero - to zero a double vector of length Num
---------------------------------------------------------*/
void utr_d_zero(
		double * Vec, 	/* in, out: vector to initialize */
		int Num		/* in: vector length */
)
{
	int i;

	for (i = 0; i < Num; i++) {
		Vec[i] = 0;
	}

	return;
}

/*---------------------------------------------------------
utr_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int utr_chk_list(	/* returns: */
		/* >0 - position on the list */
		/* 0 - not found on the list */
		int Num, 	/* number to be checked */
		int * List, 	/* list of numbers */
		int Ll		/* length of the list */
)
{

	int i, il;

	for (i = 0; i < Ll; i++) {
		if ((il = List[i]) == 0) break;
		/* found on the list on (i+1) position */
		if (Num == il) return (i + 1);
	}
	/* not found on the list */
	return (0);
}

/*---------------------------------------------------------
utr_put_list - to put Num on the list List with length Ll
	(filled with numbers and zeros at the end)
---------------------------------------------------------*/
int utr_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
		/*  <0 - position at which put on the list */
		/*   0 - list full, not found on the list */
		int Num, 	/* in: number to put on the list */
		int * List, 	/* in: list */
		int Ll		/* in: total list's lengths */
)
{

	int i, il;

	for (i = 0; i < Ll; i++) {
		if ((il = List[i]) == UTC_LIST_END_MARK) break;
		/* found on the list on (i+1) position */
		if (Num == il) return (i + 1);
	}
	/* if list is full return error message */
	if (i == Ll) return (0);
	/* update the list and return*/
	List[i] = Num;
	return (-(i + 1));
}


/*---------------------------------------------------------
utr_mat3_inv - to invert a 3x3 matrix (stored as a vector!)
---------------------------------------------------------*/
double utr_mat3_inv(	/* returns: determinant of matrix to invert */
		const double * mat,	/* matrix to invert */
		double * mat_inv	/* inverted matrix */
)
{
	double s0, s1, s2, rjac, rjac_inv;

	rjac = mat[0] * mat[4] * mat[8] + mat[3] * mat[7] * mat[2]
			+ mat[6] * mat[1] * mat[5] - mat[6] * mat[4] * mat[2]
			- mat[0] * mat[7] * mat[5] - mat[3] * mat[1] * mat[8];
	rjac_inv = 1 / rjac;

	s0 =   mat[4] * mat[8] - mat[7] * mat[5];
	s1 = - mat[3] * mat[8] + mat[6] * mat[5];
	s2 =   mat[3] * mat[7] - mat[6] * mat[4];

	mat_inv[0] = s0 * rjac_inv;
	mat_inv[3] = s1 * rjac_inv;
	mat_inv[6] = s2 * rjac_inv;

	s0 =   mat[7] * mat[2] - mat[1] * mat[8];
	s1 =   mat[0] * mat[8] - mat[6] * mat[2];
	s2 = - mat[0] * mat[7] + mat[6] * mat[1];

	mat_inv[1] = s0 * rjac_inv;
	mat_inv[4] = s1 * rjac_inv;
	mat_inv[7] = s2 * rjac_inv;

	s0 =   mat[1] * mat[5] - mat[4] * mat[2];
	s1 = - mat[0] * mat[5] + mat[3] * mat[2];
	s2 =   mat[0] * mat[4] - mat[3] * mat[1];

	mat_inv[2] = s0 * rjac_inv;
	mat_inv[5] = s1 * rjac_inv;
	mat_inv[8] = s2 * rjac_inv;

	return (rjac);
}

/*---------------------------------------------------------
utr_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
void utr_vec3_prod(
		const double * vec_a, 	/* in: vector a */
		const double * vec_b, 	/* in: vector b */
		double * vec_c	/* out: vector product axb */
)
{

	vec_c[0] = vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
	vec_c[1] = vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2];
	vec_c[2] = vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

	return;
}

/*---------------------------------------------------------
utr_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
double utr_vec3_mxpr( /* returns: mixed product [a,b,c] */
		const double * vec_a, 	/* in: vector a */
		const double * vec_b, 	/* in: vector b */
		const double * vec_c	/* in: vector c */
)
{
	double daux;

	daux  = vec_c[0] * (vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1]);
	daux += vec_c[1] * (vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2]);
	daux += vec_c[2] * (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

	return (daux);
}

/*---------------------------------------------------------
utr_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
double utr_vec3_length(	/* returns: vector length */
		const double * vec	/* in: vector */
)
{
	return (sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
}

//----------------------------------------------------------
//utr_vec3_add - to compute sum of 3D vectors
//----------------------------------------------------------

void utr_vec3_add(
		const double vec_a[3], //IN
		const double vec_b[3], //IN
		double vec_c[3])	//OUT
{
	vec_c[0] = vec_a[0] + vec_b[0];
	vec_c[1] = vec_a[1] + vec_b[1];
	vec_c[2] = vec_a[2] + vec_b[2];
}

//----------------------------------------------------------
//utr_vec3_subst - to compute the difference of 3D vectors
//----------------------------------------------------------

void utr_vec3_subst(
		const double vec_a[3],	//IN
		const double vec_b[3],	// IN
		double vec_c[3])	// OUT
{
	vec_c[0] = vec_a[0] - vec_b[0];
	vec_c[1] = vec_a[1] - vec_b[1];
	vec_c[2] = vec_a[2] - vec_b[2];
}


//----------------------------------------------------------
//utr_vec3_dot - to compute vector dot product of 3D vectors
//----------------------------------------------------------
double utr_vec3_dot(
		const double vec_a[3],
		const double vec_b[3])
{
	return vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1] + vec_a[2] * vec_b[2];
}

///----------------------------------------------------------
///utr_vec3_angle - to compute vector dot product of 3D vectors
///                      scalar_prod(a,b)
/// angle(a,b) = arc cos ---------------------
///                      length(a) * length(b)
///----------------------------------------------------------
double utr_vec3_angle( ///< returns angle in radians
		const double vec_a[3],
		const double vec_b[3])
{
	return acos(  utr_vec3_dot(vec_a, vec_b) /
					(utr_vec3_length(vec_a) * utr_vec3_length(vec_b)));
}

/*---------------------------------------------------------
utr_mat3vec - to compute matrix vector product in 3D space
---------------------------------------------------------*/
void utr_mat3vec(
		const double * m1, 	/* in: matrix (stored by rows as a vector!) */
		const double * v1, 	/* in: vector */
		double * v2	/* out: resulting vector */
)
{

	v2[0] = m1[0] * v1[0] + m1[1] * v1[1] + m1[2] * v1[2] ;
	v2[1] = m1[3] * v1[0] + m1[4] * v1[1] + m1[5] * v1[2] ;
	v2[2] = m1[6] * v1[0] + m1[7] * v1[1] + m1[8] * v1[2] ;
	return;
}

/*---------------------------------------------------------
utr_mat3mat - to compute matrix matrix product in 3D space
	(all matrices are stored by rows as vectors!)
---------------------------------------------------------*/
void utr_mat3mat(
		const double * m1,	/* in: matrix */
		const double * m2,	/* in: matrix */
		double * m3	/* out: matrix m1*m2 */
)
{

	m3[0] = m1[0] * m2[0] + m1[1] * m2[3] + m1[2] * m2[6] ;
	m3[1] = m1[0] * m2[1] + m1[1] * m2[4] + m1[2] * m2[7] ;
	m3[2] = m1[0] * m2[2] + m1[1] * m2[5] + m1[2] * m2[8] ;

	m3[3] = m1[3] * m2[0] + m1[4] * m2[3] + m1[5] * m2[6] ;
	m3[4] = m1[3] * m2[1] + m1[4] * m2[4] + m1[5] * m2[7] ;
	m3[5] = m1[3] * m2[2] + m1[4] * m2[5] + m1[5] * m2[8] ;

	m3[6] = m1[6] * m2[0] + m1[7] * m2[3] + m1[8] * m2[6] ;
	m3[7] = m1[6] * m2[1] + m1[7] * m2[4] + m1[8] * m2[7] ;
	m3[8] = m1[6] * m2[2] + m1[7] * m2[5] + m1[8] * m2[8] ;

}

/*---------------------------------------------------------
  utr_mat_det - to compute determinant of a matrix
---------------------------------------------------------*/
/*  returns determinant
	Matrix m wil NOT change!
*/
double utr_mat_det(const double * m, int n, char store, double * det)
{
	int aux = 0, i = 0;
	int  * pivots = (int *) malloc(n * sizeof(int));
	double * M = (double *) malloc(n * n * sizeof(double));

	*det = 1.0;

	if (M != NULL && pivots != NULL) {
		if (store == 'R' || store == 'r') {
			int j = 0;
			for (i = 0; i < n; ++i) {
				for (j = 0; j < n; ++j) {
					M[n * j + i] = m[n * i + j];
				}
			}
		}
		else {	// 'c' or 'C' column store schema
			memcpy(M, m, n * n * sizeof(double));
		}
		dgetrf_(&n, &n, M, &n, pivots, &aux);
		mf_check(aux >= 0, "DGETRF exited with error %d", aux);
		for (i = 0; i < n; ++i) {
			*det *= M[n * i + i] * (pivots[i] != (i + 1) ? -1.0 : 1.0);
		}
		free(M);
		free(pivots);
	}
	return (*det);
}


/*---------------------------------------------------------
utr_skip_rest_of_line - to allow for comments in input files
---------------------------------------------------------*/
void utr_skip_rest_of_line(
		FILE * Fp /* in: input file */
)
{
	while (fgetc(Fp) != '\n') {}
	return;
}


/*-------------------------------------------------------------
utr_fprintf_double -	prints double to stream in "%f" fashion
			with specified significant digits number
-------------------------------------------------------------*/
int utr_fprintf_double(
		FILE * stream,	// output stream
		int sdigits,	// number of significant digits to write - put 0 for
		// full accuracy (.12) written in "%g" fashion
		double value	// double to write
)
{
	int dignum;
	double fpart;
	double ipart;
	double res;
	double pw;

	if (sdigits == 0)
	{
		fprintf(stream, "%.12g", value);
		return 0;
	}

	fpart = modf(value, &ipart);
	if (ipart != 0.0)
	{
		dignum = (int)(log10 (fabs(ipart))) + 1;
		if (dignum <= sdigits)
			fprintf(stream, "%.*f", sdigits - dignum, value);
		else
		{
			pw = pow(10.0, (double)(dignum - sdigits));
			//res = ipart/temp;
			//	  res = round(ipart/pw)*pw;
			res = floor(0.5 + ipart / pw) * pw;
			fprintf(stream, "%.0f", res);
		}
	}
	else
	{
		dignum = -((int)(log10 (fabs(fpart))));
		fprintf(stream, "%.*f", sdigits + dignum, value);
	}

	return 0;
}

/*---------------------------------------------------------
utr_heap_sort - to heap-sort an array (code taken from fortran...)
---------------------------------------------------------*/
void utr_heap_sort(
		int  *  Ind_array,    /* in/out: index array for sorting */
		/* Ind_array[0] - the length of both arrays */
		double * Val_array    /* in: array of values used for sorting */
)
{

	int i, j, l, ir, index;
	double q;

	l = Ind_array[0] / 2 + 1;
	ir =  Ind_array[0];

s20:
	{}

	if (l > 1) {
		l--;
		index = Ind_array[l];
		q = Val_array[index];
	}
	else {
		index = Ind_array[ir];
		q = Val_array[index];
		Ind_array[ir] = Ind_array[1];
		ir--;
		if (ir == 1) {
			Ind_array[1] = index;
			return;
		}
	}

	i = l;
	j = 2 * l;

s30:
	{}

	if (j > ir) goto s40;
	if ((j < ir) && (Val_array[Ind_array[j]] < Val_array[Ind_array[j + 1]])) j++;
	if (q < Val_array[Ind_array[j]]) {
		Ind_array[i] = Ind_array[j];
		i = j;
		j *= 2;;
	}
	else {
		j = ir + 1;
	}

	goto s30;

s40:
	{}

	Ind_array[i] = index;

	goto s20;

}

/*---------------------------------------------------------
utr_insert_sort - to sort an array by insertion (for small integer arrays)
---------------------------------------------------------*/
int utr_insert_sort(
		int * A, // array
		int p,  //first index
		int k   // last index
)
{
	int i, j;
	int t;

	for (i = p + 1; i <= k; i++) {
		t = A[i];
		j = i - 1;
		while ( j >= p && A[j] > t ) {
			A[j + 1] = A[j];
			j--;
		}
		A[j + 1] = t;
	}
	return 0;
}

int utr_ctrl_pts_init(const char * Work_dir, const char * Filename, const int N_solutions,
		const char * Target_ctrl_pts_filename)
{
	utr_io_result_add_column(RESULT_CONTROL_POINTS);
	utr_io_result_ctrl_pts_read_file(Work_dir, Filename);
	utr_io_result_ctrl_pts_set_sol_len(N_solutions);
	utr_io_result_ctrl_pts_set_filename(Work_dir, Target_ctrl_pts_filename);
	return 1;
}

///
/// \brief utr_ctrl_pts_write - to get compute values at control points and send to the
/// utr_io_result submodule.
/// \param Field_id
/// \return
///
int utr_ctrl_pts_add_values(const int Field_id)
{
	const int mesh_id = apr_get_mesh_id(Field_id);

	int elems[10] = {0};
	static int was_initialized = 0;
	double glob_coords[3] = {0.0};
	double loc_coords[100] = {0.0}, solution[100] = {0.0};
	const int n_val = apr_get_nreq(Field_id);
	int i;

	if (was_initialized == 0) {
		was_initialized = 1;
		for (i = 0; i < utr_io_result_ctrl_pts_nr(); ++i) {
			utr_io_result_ctrl_pt_get_coords(i, glob_coords);

			//  apr_sol_xglob(Field_id,glob_coords,1,elems,loc_coords,solution,
			//              NULL,NULL,NULL,APC_CLOSE,APE_SOL_XGLOB_DEFAULT);
			apr_sol_xglob(Field_id, glob_coords, 1, elems, loc_coords, solution,
					NULL, NULL, NULL, APC_CLOSE, APE_SOL_XGLOB_DEFAULT);
			if (elems[0] > 0) {
				utr_io_result_ctrl_pt_set_elemID(i, mesh_id, elems[1]);
			}
			else {
				mf_fatal_err("Control point %d [%lf,%lf,%lf] was not found inside given mesh(%d)",
						i, glob_coords[0], glob_coords[1], glob_coords[2], mesh_id);
			}
		}
	}

	for (i = 0; i < utr_io_result_ctrl_pts_nr(); ++i) {

		utr_io_result_ctrl_pt_get_coords(i, glob_coords);

		// regular check with known elem ID
		elems[0] = utr_io_result_ctrl_pt_get_elemID(i, mesh_id);
		apr_sol_xglob(Field_id, glob_coords, 1, elems, loc_coords, solution,
				NULL, NULL, NULL, APC_CLOSE, APE_SOL_XGLOB_DEFAULT
				| APE_SOL_XGLOB_CHECK_ONLY_GIVEN_ELEMENT);
		//              NULL,NULL,NULL,0);
		//            NULL,NULL,NULL,APC_CLOSE, APE_SOL_XGLOB_DEFAULT
		//          | APE_SOL_XGLOB_CHECK_ONLY_GIVEN_ELEMENT);

		utr_io_result_add_ctrl_pt_values(i, n_val, solution);

	}

	return 1;
}
