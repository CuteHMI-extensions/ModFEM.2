/****************************************************************
File tms_ocl_intf.c - implementation of OpenMP/OpenCL multithreading
					  version for thread management - general purpose routines
					  (specific routines for particular stages of computations
					   are in separate files, e.g. )

Contains definitions of interface routines:
  tmr_init_multithreading_opencl - to initialize (multi)thread management
  tmr_ocl_create_command_queues - for a selected platform and device
  tmr_ocl_create_kernel - for a selected platform, device and kernel index

  tmr_ocl_get_current_platform_index() - to return the index of current OpenCL platform
  tmr_ocl_get_current_device_type
//   returns the current OpenCL device type specified based on compiler switches
  tmr_ocl_select_device
//   returns OpenCL device index (for local data structures) or -1 if device
//   is not available (not existing or not serviced) for the specified platform
  tmr_ocl_device_type - returns LOCAL device type (INTEGER)
  tmr_ocl_select_context - selects context for a given device
  tmr_ocl_select_command_queue - selects command queue for a given device
  tmr_ocl_select_kernel - selects kernel for a given device and kernel index
  tmr_ocl_cleanup - discards created OpenCL resources

Function for printing debug info:
  tmr_kernel_execution_log - print info about error returned by clEnqueueNDRangeKernel
  tmr_set_kernel_arguments_log - print info about error returned by clSetKernelArg



------------------------------
History:
	02.2013 - Krzysztof Banas, initial version
****************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

// generic thread management interface
#include <modfem/tmh_intf.h>

#include <CL/cl.h>

// interface of opencl implementation
#include <modfem/tm_opencl/tmh_ocl.h>

/* debugging/logging/checking macros */
#include <modfem/uth_log.h>


/* A single global variable */
tmt_ocl_struct tmv_ocl_struct;


#ifdef __cplusplus
extern "C" {
#endif

// DECLARATIONS OF LOCAL FUNCTIONS
void DisplayPlatformInfo(
		FILE * Interactive_output, /* file or stdout to write messages */
		cl_platform_id id,
		cl_platform_info name,
		char * str);

void DisplayDeviceInfo(
		FILE * Interactive_output, /* file or stdout to write messages */
		cl_device_id id,
		cl_device_info name,
		char * str);


/*---------------------------------------------------------
  tmr_init_multithreading_opencl - to initialize (multi)thread management
---------------------------------------------------------*/
extern int tmr_init_multithreading_opencl(
		char * Work_dir,
		int Argc,
		char ** Argv,
		FILE * Interactive_input,
		FILE * Interactive_output,
		int Control,   // not used !!!
		int Monitor
)
{

	//Monitor = TMC_PRINT_INFO + 1;

	int ocl_device = -1;
#ifdef OPENCL_PHI
	ocl_device = TMC_OCL_DEVICE_ACCELERATOR; // should be substituted in some way
#endif
#ifdef OPENCL_GPU
	ocl_device = TMC_OCL_DEVICE_GPU;
#endif
#ifdef OPENCL_CPU
	ocl_device = TMC_OCL_DEVICE_CPU; // should be substituted in some way
#endif
#ifdef OPENCL_HSA
	ocl_device = TMC_OCL_DEVICE_GPU;
#endif


	tmr_ocl_init(Work_dir, Interactive_input, Interactive_output, ocl_device, Monitor);


	return (1);
}

/*---------------------------------------------------------
  tmr_ocl_init - to initialize OpenCL thread management
---------------------------------------------------------*/
int tmr_ocl_init(
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output,
		int Control,  // OpenCL device type
		int Monitor
)
{

	// MONITOR SHOULD BE PASSED AS PARAMETER !!!!!!!!!!!!!!!!!!!!!!
	// for all operations indicate explicit info messages
	//int Monitor = TMC_PRINT_INFO + 1;

#ifdef TIME_TEST_OCL
	double t_begin, t_end;
	t_begin = time_clock();
#endif

	// Create OpenCL contexts on all available platforms
	// contexts are stored in table with indices: 0-CPU, 1-GPU, 2-ACCELERATOR
	// table entry is NULL if there is no such context for a given platform
	int number_of_platforms = tmr_ocl_create_contexts(Interactive_output,
					TMC_OCL_ALL_PLATFORMS, Monitor);

#ifdef DEBUG_TMM
	if (number_of_platforms > 1) {
		fprintf(Interactive_output,
				"\nMore than one platform in a system! Check whether it's OK with the code.\n");
	}
#endif

	if (number_of_platforms > 1) {
		if (Interactive_input == stdin) {
#ifdef OPENCL_PLATFORM_ID
			tmv_ocl_struct.current_platform_index = OPENCL_PLATFORM_ID;
#else
			fprintf(Interactive_output,
					"\nSelect platform ID (index): ");
			fscanf(Interactive_input, "%d", &tmv_ocl_struct.current_platform_index);
			//tmv_ocl_struct.current_platform_index = 0;
#endif


		}
		else { // if input from file - first platform selected

#ifdef OPENCL_PLATFORM_ID
			tmv_ocl_struct.current_platform_index = OPENCL_PLATFORM_ID;
#else
			tmv_ocl_struct.current_platform_index = 0;
#endif
		}

	}
	else { // if only one platform - first platform selected

		tmv_ocl_struct.current_platform_index = 0;

	}

	fprintf(Interactive_output,
			"\nSelected platform ID (index): %d",
			tmv_ocl_struct.current_platform_index);


#ifdef TIME_TEST_OCL
	t_end = time_clock();
	fprintf(Interactive_output, "EXECUTION TIME: creating contexts: %.12lf\n", t_end - t_begin);
#endif

#ifdef TIME_TEST_OCL
	t_begin = time_clock();
#endif


	int platform_index = tmv_ocl_struct.current_platform_index;
	// or better
	// int platform_index = tmr_ocl_get_current_platform_index();

	// create command queues on all devices
	tmr_ocl_create_command_queues(Interactive_output, platform_index,
			TMC_OCL_ALL_DEVICES, Monitor);

#ifdef TIME_TEST_OCL
	t_end = time_clock();
	fprintf(Interactive_output, "EXECUTION TIME: creating command queues: %.12lf\n",
			t_end - t_begin);
#endif



	//cl_platform_id platform = tmv_ocl_struct.list_of_platforms[platform_index].id;
	tmt_ocl_platform_struct platform_struct =
			tmv_ocl_struct.list_of_platforms[platform_index];

	/*-------------- KERNEL CREATION PHASE--------------------*/

	int device_tmc_type = Control;
	tmv_ocl_struct.current_device_type = Control;
	int device_index;
	int kernel_index;

#ifdef TIME_TEST_OCL
	t_begin = time_clock();
#endif


	if (Monitor > TMC_PRINT_INFO) {
		if (device_tmc_type == 0) {
			fprintf(Interactive_output,
					"\nSelecting OpenCL device type %d (CPU) for platform %d\n",
					device_tmc_type, platform_index);
		}
		else if (device_tmc_type == 1) {
			fprintf(Interactive_output,
					"\nSelecting OpenCL device type %d (GPU) for platform %d\n",
					device_tmc_type, platform_index);
		}
		else if (device_tmc_type == 2) {
			fprintf(Interactive_output,
					"\nSelecting OpenCL device type %d (ACCELERATOR) for platform %d\n",
					device_tmc_type, platform_index);
		}
		else {
			fprintf(Interactive_output,
					"\nUnknown OpenCL device type %d for platform %d. Exiting!\n",
					device_tmc_type, platform_index);
			exit(1);
		}
	}

	// choose device_index
	device_index = tmr_ocl_select_device(Interactive_output, platform_index, device_tmc_type);

	// if required context exists
	if (device_index >= 0) {

		// save selected device index
		tmv_ocl_struct.current_device = device_index;

		// create the kernel for numerical integration
		kernel_index = 0;

		if (Monitor > TMC_PRINT_INFO) {
			fprintf(Interactive_output,
					"\nCreating OpenCL kernel %d for device type %d and device index %d\n",
					kernel_index, device_tmc_type, device_index);
		}

		// HARDCODED DEFAULT
		int kernel_version_hw = 0;
		int kernel_version_alg = 0;

		char arg[255] = {0};

		sprintf(arg, "%s/tmr_ocl_num_int_el_and_assembling.cl", Work_dir);

		// BASED ON INPUT PARAMETERS SELECT OPENCL KERNEL FOR NUMERICAL INTEGRATION !!!!!!!
		if (kernel_index == 0) {

			// specific for numerical integration
			tmr_ocl_create_num_int_kernel(platform_index, device_index, kernel_index,
					// kernel name:         , file:
					//  "tmr_ocl_num_int_el", "tmr_ocl_num_int_el.cl", Monitor);
					"tmr_ocl_num_int_el", arg, Monitor);

			// generic procedure
			/* tmr_ocl_create_kernel(Interactive_output,platform_index,device_index,kernel_index, */
			/*  // kernel name:         , file: */
			/*  //  "tmr_ocl_num_int_el", "tmr_ocl_num_int_el.cl", Monitor); */
			/*   "tmr_ocl_num_int_el", arg, Monitor); */
		}

		// create kernel for creating the final CRS
		kernel_index = 1;

		if (Monitor > TMC_PRINT_INFO) {
			fprintf(Interactive_output,
					"\nCreating OpenCL kernel %d for device type %d and device index %d\n",
					kernel_index, device_tmc_type, device_index);
		}

		// specific for finalize the crs creation
		tmr_ocl_create_num_int_kernel(platform_index, device_index, kernel_index,
				// kernel name:         , file:
				//  "tmr_ocl_num_int_el", "tmr_ocl_num_int_el.cl", Monitor);
				"tmr_ocl_prepare_final_crs", arg, Monitor);



		/* #ifdef GPU_SEPARATED_GEO_DATA */

		/* // create kernel for rewrite geo dofs */
		/* kernel_index = 3; */

		/* if(Monitor>TMC_PRINT_INFO){ */
		/*   fprintf(Interactive_output,  */
		/*          "\nCreating OpenCL kernel %d for device type %d and device index %d\n",  */
		/*          kernel_index, device_tmc_type, device_index); */
		/* } */

		/* // BASED ON INPUT PARAMETERS SELECT OPENCL KERNEL FOR REWRITING GEO DOFS !!!!!!! */
		/* if(kernel_index == 3){ */
		/*   // specific for rewrite geo dofs */
		/*   tmr_ocl_create_num_int_kernel(platform_index,device_index,kernel_index, */
		/* 				    // kernel name:         , file: */
		/* 				    //  "tmr_ocl_num_int_el", "tmr_ocl_num_int_el.cl", Monitor); */
		/* 				    "tmr_ocl_rewrite_geo_dofs", arg, Monitor); */
		/* } */

		/* #endif */

#ifdef TIME_TEST_OCL
		t_end = time_clock();
		fprintf(Interactive_output, "EXECUTION TIME: creating OpenCL kernel for integration: %.12lf\n", t_end - t_begin);
#endif

	}
	else {
		fprintf(Interactive_output, "Requested device not supported for this platform.\n");
		exit(-1);
	}

	return (1);
}


/*---------------------------------------------------------
  DisplayPlatformInfo - utility local procedure
---------------------------------------------------------*/
void DisplayPlatformInfo(
		FILE * Interactive_output, /* file or stdout to write messages */
		cl_platform_id id,
		cl_platform_info name,
		char * str)
{
	cl_int retval;
	size_t paramValueSize;

	retval = clGetPlatformInfo(
					id,
					name,
					0,
					NULL,
					&paramValueSize);
	if (retval != CL_SUCCESS) {
		fprintf(Interactive_output, "Failed to find OpenCL platform %s.\n", str);
		return;
	}

	char * info = (char *)malloc(sizeof(char) * paramValueSize);
	retval = clGetPlatformInfo(
					id,
					name,
					paramValueSize,
					info,
					NULL);
	if (retval != CL_SUCCESS)  {
		fprintf(Interactive_output, "Failed to find OpenCL platform %s.\n", str);
		return;
	}

	fprintf(Interactive_output, "\t%s:\t%s\n", str, info );
	free(info);
}

/*---------------------------------------------------------
  DisplayDeviceInfo - utility local procedure
---------------------------------------------------------*/
void DisplayDeviceInfo(
		FILE * Interactive_output, /* file or stdout to write messages */
		cl_device_id id,
		cl_device_info name,
		char * str)
{
	cl_int retval;
	size_t paramValueSize;

	retval = clGetDeviceInfo(
					id,
					name,
					0,
					NULL,
					&paramValueSize);
	if (retval != CL_SUCCESS) {
		fprintf(Interactive_output, "Failed to find OpenCL device info %s.\n", str);
		return;
	}

	char * info = (char *)malloc(sizeof(char) * paramValueSize);
	retval = clGetDeviceInfo(
					id,
					name,
					paramValueSize,
					info,
					NULL);

	if (retval != CL_SUCCESS) {
		fprintf(Interactive_output, "Failed to find OpenCL device info %s.\n", str);
		return;
	}

	fprintf(Interactive_output, "\t\t%s:\t%s\n", str, info );
	free(info);
};

/*---------------------------------------------------------
//  tmr_ocl_create_contexts - to create OpenCL contexts on a selected platform
---------------------------------------------------------*/
int tmr_ocl_create_contexts(
		FILE * Interactive_output, /* file or stdout to write messages */
		int Platform_id_in,
		int Monitor
)
{
	cl_int retval;
	cl_uint numPlatforms;
	cl_platform_id * platformIds;
	cl_context context = NULL;
	cl_uint iplat, jdev, k;

	// First, query the total number of platforms
	retval = clGetPlatformIDs(0, (cl_platform_id *) NULL, &numPlatforms);

	// allocate memory for local platform structures
	tmv_ocl_struct.number_of_platforms = numPlatforms;
	tmv_ocl_struct.list_of_platforms =
			(tmt_ocl_platform_struct *) malloc( sizeof(tmt_ocl_platform_struct)
					* numPlatforms);

	// Next, allocate memory for the installed platforms, and qeury
	// to get the list.
	platformIds = (cl_platform_id *)malloc(sizeof(cl_platform_id) * numPlatforms);
	retval = clGetPlatformIDs(numPlatforms, platformIds, NULL);

	if (Monitor >= TMC_PRINT_INFO) {
		fprintf(Interactive_output, "\nNumber of OpenCL platforms: \t%d\n", numPlatforms);
	}

	// Iterate through the list of platforms displaying associated information
	for (iplat = 0; iplat < numPlatforms; iplat++) {

		if (Monitor > TMC_PRINT_INFO) {
			fprintf(Interactive_output, "\n");
			fprintf(Interactive_output, "Platform %d:\n", iplat);
		}

		tmv_ocl_struct.list_of_platforms[iplat].id = platformIds[iplat];
		//clGetPlatformInfo(platformIds[iplat], CL_PLATFORM_NAME, size_of_name???,
		//		      tmv_ocl_struct.list_of_platforms[iplat].name, (size_t *) NULL);

		if (Monitor > TMC_PRINT_INFO) {

			// First we display information associated with the platform
			DisplayPlatformInfo(Interactive_output,
					platformIds[iplat],
					CL_PLATFORM_NAME,
					"CL_PLATFORM_NAME");
			DisplayPlatformInfo(Interactive_output,
					platformIds[iplat],
					CL_PLATFORM_PROFILE,
					"CL_PLATFORM_PROFILE");
			DisplayPlatformInfo(Interactive_output,
					platformIds[iplat],
					CL_PLATFORM_VERSION,
					"CL_PLATFORM_VERSION");
			DisplayPlatformInfo(Interactive_output,
					platformIds[iplat],
					CL_PLATFORM_VENDOR,
					"CL_PLATFORM_VENDOR");
		}

		// Now query the set of devices associated with the platform
		cl_uint numDevices;
		retval = clGetDeviceIDs(
						platformIds[iplat],
						CL_DEVICE_TYPE_ALL,
						0,
						NULL,
						&numDevices);


		tmv_ocl_struct.list_of_platforms[iplat].number_of_devices = numDevices;
		tmv_ocl_struct.list_of_platforms[iplat].list_of_devices =
				(tmt_ocl_device_struct *) malloc( sizeof(tmt_ocl_device_struct)
						* numDevices);

		cl_device_id * devices =
				(cl_device_id *) malloc (sizeof(cl_device_id) * numDevices);

		retval = clGetDeviceIDs(
						platformIds[iplat],
						CL_DEVICE_TYPE_ALL,
						numDevices,
						devices,
						NULL);

		if (Monitor >= TMC_PRINT_INFO) {
			fprintf(Interactive_output, "\n\tNumber of devices: \t%d\n", numDevices);
		}
		// Iterate through each device, displaying associated information
		for (jdev = 0; jdev < numDevices; jdev++)
		{

			if (Monitor > TMC_PRINT_INFO) {
				fprintf(Interactive_output, "\tDevice %d:\n", jdev);
			}
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].id =
					devices[jdev];
			clGetDeviceInfo(devices[jdev], CL_DEVICE_TYPE, sizeof(cl_device_type),
					&tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].type, NULL);

			if (tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].type == CL_DEVICE_TYPE_CPU) {
				tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].tmc_type = TMC_OCL_DEVICE_CPU;
			}
			if (tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].type == CL_DEVICE_TYPE_GPU) {
				tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].tmc_type = TMC_OCL_DEVICE_GPU;
			}
			if (tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].type == CL_DEVICE_TYPE_ACCELERATOR) {
				tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].tmc_type = TMC_OCL_DEVICE_ACCELERATOR;
			}

			cl_ulong mem_size_ulong = 0;
			int err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_GLOBAL_MEM_SIZE,
							sizeof(cl_ulong), &mem_size_ulong, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].global_mem_bytes =
					(double)mem_size_ulong;

			err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_MAX_MEM_ALLOC_SIZE,
							sizeof(cl_ulong), &mem_size_ulong, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].global_max_alloc =
					(double)mem_size_ulong;

			err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_LOCAL_MEM_SIZE,
							sizeof(cl_ulong), &mem_size_ulong, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].shared_mem_bytes =
					(double)mem_size_ulong;

			err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,
							sizeof(cl_ulong), &mem_size_ulong, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].cache_bytes =
					(double)mem_size_ulong;

			err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
							sizeof(cl_ulong), &mem_size_ulong, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].constant_mem_bytes =
					(double)mem_size_ulong;

			cl_uint cache_line_size = 0;
			err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,
							sizeof(cl_uint), &cache_line_size, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].cache_line_bytes =
					(int) cache_line_size;

			cl_uint max_num_comp_units = 0;
			err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_MAX_COMPUTE_UNITS,
							sizeof(cl_uint), &max_num_comp_units, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].max_num_comp_units =
					(int) max_num_comp_units;

			size_t max_work_group_size = 0;
			err_num = clGetDeviceInfo(devices[jdev], CL_DEVICE_MAX_WORK_GROUP_SIZE,
							sizeof(size_t), &max_work_group_size, NULL);
			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].max_work_group_size =
					(int) max_work_group_size;

			// possible further inquires:
			//CL_DEVICE_MAX_WORK_GROUP_SIZE,
			//CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, CL_DEVICE_MAX_WORK_ITEM_SIZES
			//CL_DEVICE_MAX_CONSTANT_ARGS
			//CL_DEVICE_MAX_PARAMETER_SIZE
			//CL_DEVICE_PREFERRED_VECTOR_WIDTH_ - char, int, float, double etc.
			//CL_DEVICE_MEM_BASE_ADDR_ALIGN, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE

			tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].command_queue = 0;

			for (k = 0; k < TMC_OCL_MAX_NUM_KERNELS; k++) {
				tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].program[k] = 0;
				tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].kernel[k] = 0;
			}

			//clGetDeviceInfo(devices[jdev], CL_DEVICE_NAME, sizeof(device_name?),
			//&tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].name, NULL);

			if (Monitor > TMC_PRINT_INFO) {

				DisplayDeviceInfo(Interactive_output,
						devices[jdev],
						CL_DEVICE_NAME,
						"CL_DEVICE_NAME");

				DisplayDeviceInfo(Interactive_output,
						devices[jdev],
						CL_DEVICE_VENDOR,
						"CL_DEVICE_VENDOR");

				DisplayDeviceInfo(Interactive_output,
						devices[jdev],
						CL_DEVICE_VERSION,
						"CL_DEVICE_VERSION");
				fprintf(Interactive_output, "\t\tdevice global memory size (MB) = %lf\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].global_mem_bytes / 1024 / 1024);
				fprintf(Interactive_output, "\t\tdevice global max alloc size (MB) = %lf\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].global_max_alloc / 1024 / 1024);
				fprintf(Interactive_output, "\t\tdevice local memory size (kB) = %lf\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].shared_mem_bytes / 1024);
				fprintf(Interactive_output, "\t\tdevice constant memory size (kB) = %lf\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].constant_mem_bytes / 1024);
				fprintf(Interactive_output, "\t\tdevice cache memory size (kB) = %lf\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].cache_bytes / 1024);
				fprintf(Interactive_output, "\t\tdevice cache line size (B) = %d\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].cache_line_bytes);
				fprintf(Interactive_output, "\t\tdevice maximal number of comptme units = %d\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].max_num_comp_units);
				fprintf(Interactive_output, "\t\tdevice maximal number of work units in work group = %d\n",
						tmv_ocl_struct.list_of_platforms[iplat].list_of_devices[jdev].max_work_group_size);

				fprintf(Interactive_output, "\n");
			}
		}

		free(devices);

		// Next, create OpenCL contexts on platforms
		cl_context_properties contextProperties[] = {
			CL_CONTEXT_PLATFORM,
			(cl_context_properties)platformIds[iplat],
			0
		};

		if (Platform_id_in == TMC_OCL_ALL_PLATFORMS || Platform_id_in == iplat) {

			if (Monitor > TMC_PRINT_INFO) {
				fprintf(Interactive_output, "\tCreating CPU context (index=0) on platform %d\n", iplat);
			}

			tmv_ocl_struct.list_of_platforms[iplat].list_of_contexts[0] =
					clCreateContextFromType(contextProperties,
							CL_DEVICE_TYPE_CPU, NULL, NULL, &retval);

			if (Monitor >= TMC_PRINT_INFO && retval != CL_SUCCESS) {
				fprintf(Interactive_output, "\tCould not create CPU context on platform %d\n", iplat);
			}

			if (Monitor > TMC_PRINT_INFO) {
				fprintf(Interactive_output, "\tCreating GPU context (index=1) on platform %d\n", iplat);
			}

			tmv_ocl_struct.list_of_platforms[iplat].list_of_contexts[1] =
					clCreateContextFromType(contextProperties,
							CL_DEVICE_TYPE_GPU, NULL, NULL, &retval);

			if (Monitor >= TMC_PRINT_INFO && retval != CL_SUCCESS) {
				fprintf(Interactive_output, "\tCould not create GPU context on platform %d\n", iplat);
			}

			if (Monitor > TMC_PRINT_INFO) {
				fprintf(Interactive_output, "\tCreating ACCELERATOR context (index=2) on platform %d\n", iplat);
			}

			tmv_ocl_struct.list_of_platforms[iplat].list_of_contexts[2] =
					clCreateContextFromType(contextProperties,
							CL_DEVICE_TYPE_ACCELERATOR, NULL, NULL, &retval);
			if (Monitor >= TMC_PRINT_INFO && retval != CL_SUCCESS) {
				fprintf(Interactive_output, "\tCould not create ACCELERATOR context on platform %d\n", iplat);
			}

		}
	}

	free(platformIds);
	return numPlatforms;
}


/*---------------------------------------------------------
  tmr_ocl_create_command_queues - for a selected platform and device
---------------------------------------------------------*/
int tmr_ocl_create_command_queues(
		FILE * Interactive_output, /* file or stdout to write messages */
		int Platform_index,
		int Device_type,
		int Monitor
)
{

	// in a loop over all platforms
	int platform_index;
	for (platform_index = 0;
			platform_index < tmv_ocl_struct.number_of_platforms;
			platform_index++) {

		// shortctm for global platform structure
		tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[platform_index];

		// if creating contexts for all platforms or just this one
		if (Platform_index == TMC_OCL_ALL_PLATFORMS ||
				Platform_index == platform_index) {

			// in a loop over all devices
			int idev;
			for (idev = 0;
					idev < platform_struct.number_of_devices;
					idev++) {

				// variable for storing device_id
				cl_device_id device = 0;

				// select context for the device (CPU context for CPU device, etc.)
				// (contexts are already created!,
				// icon is just the index in the platform structure)
				int icon;

				// check whether this is a CPU device - then context is no 0
				if (platform_struct.list_of_devices[idev].type ==
						CL_DEVICE_TYPE_CPU) {

					if (Device_type == TMC_OCL_ALL_DEVICES ||
							Device_type == TMC_OCL_DEVICE_CPU) {

						device = platform_struct.list_of_devices[idev].id;
						platform_struct.list_of_devices[idev].tmc_type = TMC_OCL_DEVICE_CPU;
						icon = 0;

					}
					else {

						device = NULL;

					}

				}
				// check whether this is a GPU device - then context is no 1
				else if (platform_struct.list_of_devices[idev].type ==
						CL_DEVICE_TYPE_GPU) {

					if (Device_type == TMC_OCL_ALL_DEVICES ||
							Device_type == TMC_OCL_DEVICE_GPU) {

						device = platform_struct.list_of_devices[idev].id;
						platform_struct.list_of_devices[idev].tmc_type = TMC_OCL_DEVICE_GPU;
						icon = 1;

					}
					else {

						device = NULL;

					}

				}
				// check whether this is an ACCELERATOR device - then context is no 2
				else if (platform_struct.list_of_devices[idev].type ==
						CL_DEVICE_TYPE_ACCELERATOR) {

					if (Device_type == TMC_OCL_ALL_DEVICES ||
							Device_type == TMC_OCL_DEVICE_ACCELERATOR) {

						device = platform_struct.list_of_devices[idev].id;
						platform_struct.list_of_devices[idev].tmc_type = TMC_OCL_DEVICE_ACCELERATOR;
						icon = 2;

					}
					else {

						device = NULL;

					}

				}

				if (device != NULL) {

					// choose OpenCL context selected for a device
					cl_context context = platform_struct.list_of_contexts[icon];
					platform_struct.list_of_devices[idev].context_index = icon;

					// if context exist
					if (context != NULL) {

						if (Monitor > TMC_PRINT_INFO) {
							if (platform_struct.list_of_devices[idev].tmc_type == TMC_OCL_DEVICE_CPU) {
								fprintf(Interactive_output, "\nCreating command queue for CPU context %d, device index %d, platform %d\n",
										icon, idev, platform_index);
							}
							if (platform_struct.list_of_devices[idev].tmc_type == TMC_OCL_DEVICE_GPU) {
								fprintf(Interactive_output, "\nCreating command queue for GPU context %d, device index %d, platform %d\n",
										icon, idev, platform_index);
							}
							if (platform_struct.list_of_devices[idev].tmc_type == TMC_OCL_DEVICE_ACCELERATOR) {
								fprintf(Interactive_output, "\nCreating command queue for ACCELERATOR context %d, device index %d, platform %d\n",
										icon, idev, platform_index);
							}
						}

						// Create a command-queue on the device for the context
						cl_command_queue_properties prop = 0;
						prop |= CL_QUEUE_PROFILING_ENABLE;
						platform_struct.list_of_devices[idev].command_queue =
								clCreateCommandQueue(context, device, prop, NULL);
						if (platform_struct.list_of_devices[idev].command_queue == NULL)
						{
							fprintf(Interactive_output, "Failed to create command queue for context %d, device %d, platform %d\n",
									icon, idev, platform_index);
							exit(-1);
						}

					} // end if context exist for a given device

				} // end if device is of specified type

			} // end loop over devices

		} // end if platform is of specified type

	} // end loop over platforms

	return (1);
}

/*---------------------------------------------------------
// auxiliary procedure for reading source files
----------------------------------------------------------*/
char * tmr_ocl_readSource(
		FILE * Interactive_output, /* file or stdout to write messages */
		const char * kernelPath
) {

	cl_int status;
	FILE * fp;
	char * source;
	long int size;

	fp = fopen(kernelPath, "rb");
	if (!fp) {
		fprintf(Interactive_output, "Could not open kernel file\n");
		exit(-1);
	}
	status = fseek(fp, 0, SEEK_END);
	if (status != 0) {
		fprintf(Interactive_output, "Error seeking to end of file\n");
		exit(-1);
	}
	size = ftell(fp);
	if (size < 0) {
		fprintf(Interactive_output, "Error getting file position\n");
		exit(-1);
	}

	rewind(fp);

	source = (char *)malloc(size + 1);

	int i;
	for (i = 0; i < size + 1; i++) {
		source[i] = '\0';
	}

	if (source == NULL) {
		fprintf(Interactive_output, "Error allocating space for the kernel source\n");
		exit(-1);
	}

	fread(source, 1, size, fp);
	source[size] = '\0';

	return source;
}

/*---------------------------------------------------------
// procedure for vendor checking needed for kernel creation
 * 1 - NVIDIA
 * 2 - INTEL
 * 3 - AMD
----------------------------------------------------------*/

int tmr_ocl_check_vendor()
{
	cl_platform_id * platformIds;
	size_t paramValueSize;

	platformIds = (cl_platform_id *)malloc(sizeof(cl_platform_id));
	clGetPlatformIDs(tmv_ocl_struct.number_of_platforms, platformIds, NULL);



	clGetPlatformInfo(
			platformIds[tmv_ocl_struct.current_platform_index],
			CL_PLATFORM_VENDOR,
			0,
			NULL,
			&paramValueSize);

	char * info = (char *)malloc(sizeof(char) * paramValueSize);
	clGetPlatformInfo(
			platformIds[tmv_ocl_struct.current_platform_index],
			CL_PLATFORM_VENDOR,
			paramValueSize,
			info,
			NULL);
	//fprintf(Interactive_output,"%s\n", info );

	int test;
	test = strcmp("NVIDIA Corporation", info);
	if (test == 0)
		return 1;
	test = strcmp("Intel(R) Corporation", info);
	if (test == 0)
		return 2;
	test = strcmp("Advanced Micro Devices, Inc.", info);
	if (test == 0)
		return 3;

	return 0;
}

int wc(char * file_path, char * word) {
	FILE * fp;
	int count = 0;
	int ch, len;

	if (NULL == (fp = fopen(file_path, "r")))
		return -1;
	len = strlen(word);
	for (;;) {
		int i;
		if (EOF == (ch = fgetc(fp))) break;
		if ((char)ch != *word) continue;
		for (i = 1; i < len; ++i) {
			if (EOF == (ch = fgetc(fp))) goto end;
			if ((char)ch != word[i]) {
				fseek(fp, 1 - i, SEEK_CUR);
				goto next;
			}
		}
		++count;
next:
		;
	}
end:
	fclose(fp);
	return count;
}

/*---------------------------------------------------------
  tmr_ocl_create_kernel - for a selected platform, device and kernel index
						  (generic procedure - for particular kernels may exist
						   specific procedures, e.g. tmr_ocl_create_num_int_kernel
						   for numerical integration)
---------------------------------------------------------*/
int tmr_ocl_create_kernel(
		FILE * Interactive_output, /* file or stdout to write messages */
		int Platform_index,
		int Device_index,
		int Kernel_index,
		char * Kernel_name,
		const char * Kernel_file,
		int Monitor
)
{

	cl_int retval;

	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];

	// check the device !!!!!!!!!!!!!!!!! (or at least its index)
	if (Device_index < 0) {
		fprintf(Interactive_output, "Wrong device_index %d passed to tmr_ocl_create_kernel! Exiting.\n",
				Device_index);
		exit(-1);
	}

	cl_device_id device = platform_struct.list_of_devices[Device_index].id;
	cl_context context = platform_struct.list_of_contexts[
				 platform_struct.list_of_devices[Device_index].context_index
		 ];

	if (Monitor > TMC_PRINT_INFO) {
		fprintf(Interactive_output, "Program file is: %s\n", Kernel_file);
	}

	// read source file into data structure
	const char * source = tmr_ocl_readSource(Interactive_output, Kernel_file);

	cl_program program = clCreateProgramWithSource(context, 1,
					&source,
					NULL, NULL);
	if (program == NULL)
	{
		fprintf(Interactive_output, "Failed to create CL program from source.\n");
		exit(-1);
	}


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


	char * buildLog;
	size_t size_of_buildLog;
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
			0, NULL, &size_of_buildLog);
	buildLog = (char *)malloc(size_of_buildLog + 1);
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
			size_of_buildLog, buildLog, NULL);
	buildLog[size_of_buildLog] = '\0';
	if (Monitor > TMC_PRINT_INFO) {
		fprintf(Interactive_output, "Kernel buildLog: %s\n", buildLog);
	}
	if (retval != CL_SUCCESS) {
		fprintf(Interactive_output, "Error building program in tmr_ocl_create_kernel\n");
		exit(-1);
	}

	// Create OpenCL kernel
	cl_kernel kernel = clCreateKernel(program, Kernel_name, NULL);
	if (kernel == NULL)
	{
		fprintf(Interactive_output, "Failed to create kernel.\n");
		exit(-1);
	}

	if (Monitor > TMC_PRINT_INFO) {
		fprintf(Interactive_output, "Created kernel for platform %d, device %d, kernel index %d\n",
				Platform_index, Device_index, Kernel_index);
	}


	platform_struct.list_of_devices[Device_index].program[Kernel_index] = program;
	platform_struct.list_of_devices[Device_index].kernel[Kernel_index] = kernel;

	return (1);
}

/*---------------------------------------------------------
tmr_ocl_get_current_platform_index() - to return the index of current OpenCL platform
---------------------------------------------------------*/
int tmr_ocl_get_current_platform_index()
{
	return (tmv_ocl_struct.current_platform_index);
}

/*---------------------------------------------------------
  tmr_ocl_get_current_device_type
//   returns the current OpenCL device type specified based on compiler switches
---------------------------------------------------------*/
int tmr_ocl_get_current_device_type()
{
	return (tmv_ocl_struct.current_device_type);
}

/**--------------------------------------------------------
tmr_ocl_get_current_device() - to return the index of current OpenCL device
(selected based on compiler switches)
---------------------------------------------------------*/
int tmr_ocl_get_current_device()
{
	return (tmv_ocl_struct.current_device);
}

/*---------------------------------------------------------
  tmr_ocl_select_device
//   returns OpenCL device index (for local data structures) or -1 if device
//   is not available (not existing or not serviced) for the specified platform
---------------------------------------------------------*/
int tmr_ocl_select_device(
		FILE * Interactive_output, /* file or stdout to write messages */
		int Platform_index,
		int Device_tmc_type
)
{
	int device_index = -1;

	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];

	// in a loop over all devices
	int idev;
	int idev_type;
	int num_devices_of_the_same_type = 0;

	// count devices
	for (idev = 0; idev < platform_struct.number_of_devices; idev++) {

		fprintf(Interactive_output, "--> platform %d, idev %d, type %d, inptm_type %d\n",
				Platform_index, idev, platform_struct.list_of_devices[idev].tmc_type, Device_tmc_type);

		// check device type
		if (platform_struct.list_of_devices[idev].tmc_type == Device_tmc_type) {
			num_devices_of_the_same_type++;
			device_index = idev;
			//break;
		}
	}

#ifdef PARALLEL

	int my_id = pcr_my_proc_id();
	int my_device = my_id % num_devices_of_the_same_type;

	for (idev = 0, idev_type = 0; idev < platform_struct.number_of_devices; idev++) {
		if (platform_struct.list_of_devices[idev].tmc_type == Device_tmc_type) {
			if (my_device == idev_type) {
				device_index = idev;
				break;
			}
			idev_type++;
		}
	}

	fprintf(Interactive_output, "--> selected device for process %d is %d\n",
			my_id, device_index);

#else

	for (idev = 0; idev < platform_struct.number_of_devices; idev++) {

		// check device type
		if (platform_struct.list_of_devices[idev].tmc_type == Device_tmc_type) {
			num_devices_of_the_same_type++;
			device_index = idev;
			break;
		}
	}

	fprintf(Interactive_output, "--> selected device is %d\n",
			device_index);

#endif


	return (device_index);
}


/*---------------------------------------------------------
  tmr_ocl_device_type - returns LOCAL device type (INTEGER)
---------------------------------------------------------*/
int tmr_ocl_device_type(
		int Platform_index,
		int Device_index
)
{
	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];
	return (platform_struct.list_of_devices[Device_index].tmc_type);
}


/*---------------------------------------------------------
  tmr_ocl_select_device_id - selects command queue for a given device
---------------------------------------------------------*/
cl_device_id tmr_ocl_select_device_id(
		int Platform_index,
		int Device_index
)
{
	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];
	return (platform_struct.list_of_devices[Device_index].id);
}


/*---------------------------------------------------------
  tmr_ocl_select_context - selects context for a given device
---------------------------------------------------------*/
cl_context tmr_ocl_select_context(
		int Platform_index,
		int Device_index
)
{
	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];
	int context_index = platform_struct.list_of_devices[Device_index].context_index;
	return (platform_struct.list_of_contexts[context_index]);
}

/*---------------------------------------------------------
  tmr_ocl_select_command_queue - selects command queue for a given device
---------------------------------------------------------*/
cl_command_queue tmr_ocl_select_command_queue(
		int Platform_index,
		int Device_index
)
{
	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];
	return (platform_struct.list_of_devices[Device_index].command_queue);
}

/*---------------------------------------------------------
  tmr_ocl_select_kernel - selects kernel for a given device and kernel index
---------------------------------------------------------*/
cl_kernel tmr_ocl_select_kernel(
		int Platform_index,
		int Device_index,
		int Kernel_index
)
{
	// choose the platform
	tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];
	return (platform_struct.list_of_devices[Device_index].kernel[Kernel_index]);
}

/*---------------------------------------------------------
  tmr_ocl_cleanup - discards created OpenCL resources
---------------------------------------------------------*/
void tmr_ocl_cleanup()
{
	int i, j, k;

	for (i = 0; i < tmv_ocl_struct.number_of_platforms; i++) {

		for (j = 0; j < tmv_ocl_struct.list_of_platforms[i].number_of_devices; j++) {

			if (tmv_ocl_struct.list_of_platforms[i].list_of_devices[j].command_queue != 0) {
				clReleaseCommandQueue(
						tmv_ocl_struct.list_of_platforms[i].list_of_devices[j].command_queue);
			}

			for (k = 0; k < TMC_OCL_MAX_NUM_KERNELS; k++) {
				if (tmv_ocl_struct.list_of_platforms[i].list_of_devices[j].kernel[k] != 0) {
					clReleaseKernel(tmv_ocl_struct.list_of_platforms[i].list_of_devices[j].kernel[k]);
				}
			}

			for (k = 0; k < TMC_OCL_MAX_NUM_KERNELS; k++) {
				if (tmv_ocl_struct.list_of_platforms[i].list_of_devices[j].program[k] != 0) {
					clReleaseProgram(tmv_ocl_struct.list_of_platforms[i].list_of_devices[j].program[k]);
				}
			}

		}

		free(tmv_ocl_struct.list_of_platforms[i].list_of_devices);

		if (tmv_ocl_struct.list_of_platforms[i].list_of_contexts[0] != 0)
			clReleaseContext(tmv_ocl_struct.list_of_platforms[i].list_of_contexts[0]);

		if (tmv_ocl_struct.list_of_platforms[i].list_of_contexts[1] != 0)
			clReleaseContext(tmv_ocl_struct.list_of_platforms[i].list_of_contexts[1]);

		if (tmv_ocl_struct.list_of_platforms[i].list_of_contexts[2] != 0)
			clReleaseContext(tmv_ocl_struct.list_of_platforms[i].list_of_contexts[2]);

	}

	free(tmv_ocl_struct.list_of_platforms);

}

/*---------------------------------------------------------
tmr_kernel_execution_log - print info about error returned by clEnqueueNDRangeKernel
---------------------------------------------------------*/
void tmr_kernel_execution_log(cl_int error_code, int line_code, const char * function_name)
{
	switch (error_code)
	{

		case CL_INVALID_PROGRAM_EXECUTABLE:
			mf_fatal_err("CL_INVALID_PROGRAM_EXECUTABLE function: %s code line: %d \n Description: if there is no successfully built program executable available for device associated with command_queue.", function_name, line_code);
			break;

		case CL_INVALID_COMMAND_QUEUE:
			mf_fatal_err("CL_INVALID_COMMAND_QUEUE function: %s code line: %d \n Description: if command_queue is not a valid command-queue.", function_name, line_code);
			break;

		case CL_INVALID_KERNEL:
			mf_fatal_err("CL_INVALID_KERNEL function: %s code line: %d \n Description: if kernel is not a valid kernel object.", function_name, line_code);
			break;
		case CL_INVALID_CONTEXT:
			mf_fatal_err("CL_INVALID_CONTEXT function: %s code line: %d Description: if context associated with command_queue and kernel is not the same or if the context associated with command_queue and events in event_wait_list are not the same.\n ", function_name, line_code);
			break;

		case CL_INVALID_KERNEL_ARGS:
			mf_fatal_err("CL_INVALID_KERNEL_ARGS function: %s code line: %d \n Description: if the kernel argument values have not been specified.", function_name, line_code);
			break;

		case CL_INVALID_WORK_DIMENSION:
			mf_fatal_err("CL_INVALID_WORK_DIMENSION function: %s code line: %d \n Description: if work_dim is not a valid value (i.e. a value between 1 and 3).", function_name, line_code);
			break;

		case CL_INVALID_WORK_GROUP_SIZE:
			mf_fatal_err("CL_INVALID_WORK_GROUP_SIZE function: %s code line: %d \n Description 1: if local_work_size is specified and number of work-items specified by global_work_size is not evenly divisable by size of work-group given by local_work_size or does not match the work-group size specified for kernel using the __attribute__((reqd_work_group_size(X, Y, Z))) qualifier in program source. \n Description 2: if local_work_size is specified and the total number of work-items in the work-group computed as local_work_size[0] *... local_work_size[work_dim - 1] is greater than the value specified by CL_DEVICE_MAX_WORK_GROUP_SIZE in the table of OpenCL Device Queries for clGetDeviceInfo. \n Description 3: if local_work_size is NULL and the __attribute__((reqd_work_group_size(X, Y, Z))) qualifier is used to declare the work-group size for kernel in the program source.", function_name, line_code);
			break;

		case CL_INVALID_WORK_ITEM_SIZE:
			mf_fatal_err("CL_INVALID_WORK_ITEM_SIZE function: %s code line: %d \n Description: if the number of work-items specified in any of local_work_size[0], ... local_work_size[work_dim - 1] is greater than the corresponding values specified by CL_DEVICE_MAX_WORK_ITEM_SIZES[0], .... CL_DEVICE_MAX_WORK_ITEM_SIZES[work_dim - 1].", function_name, line_code);
			break;

		case CL_INVALID_GLOBAL_OFFSET:
			mf_fatal_err("CL_INVALID_GLOBAL_OFFSET function: %s code line: %d \n Description: if global_work_offset is not NULL.", function_name, line_code);
			break;

		case CL_OUT_OF_RESOURCES:
			mf_fatal_err("CL_OUT_OF_RESOURCES function: %s code line: %d \n Description: if there is a failure to queue the execution instance of kernel on the command-queue because of insufficient resources needed to execute the kernel. For example, the explicitly specified local_work_size causes a failure to execute the kernel because of insufficient resources such as registers or local memory. Another example would be the number of read-only image args used in kernel exceed the CL_DEVICE_MAX_READ_IMAGE_ARGS value for device or the number of write-only image args used in kernel exceed the CL_DEVICE_MAX_WRITE_IMAGE_ARGS value for device or the number of samplers used in kernel exceed CL_DEVICE_MAX_SAMPLERS for device.", function_name, line_code);
			break;

		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			mf_fatal_err("CL_MEM_OBJECT_ALLOCATION_FAILURE function: %s code line: %d \n Description: if there is a failure to allocate memory for data store associated with image or buffer objects specified as arguments to kernel.", function_name, line_code);
			break;

		case CL_INVALID_EVENT_WAIT_LIST:
			mf_fatal_err("CL_INVALID_EVENT_WAIT_LIST function: %s code line: %d \n Description: if event_wait_list is NULL and num_events_in_wait_list > 0, or event_wait_list is not NULL and num_events_in_wait_list is 0, or if event objects in event_wait_list are not valid events.", function_name, line_code);
			break;

		case CL_OUT_OF_HOST_MEMORY:
			mf_fatal_err("CL_OUT_OF_HOST_MEMORY function: %s code line: %d \n Description: if there is a failure to allocate resources required by the OpenCL implementation on the host.", function_name, line_code);
			break;

		case CL_SUCCESS:
			mf_log_info("SUCCESS function: %s code line: %d ", function_name, line_code);
			break;

		default:
			mf_log_warn("UNKNOWN FLAG");
			break;
	};
}

/*---------------------------------------------------------
tmr_set_kernel_arguments_log - print info about error returned by clSetKernelArg
---------------------------------------------------------*/
void tmr_set_kernel_arguments_log(cl_int error_code, int argument_index, int line_code, const char * function_name)
{
	switch (error_code)
	{

		case CL_INVALID_KERNEL:
			mf_fatal_err("CL_INVALID_KERNEL argument: %d function: %s code line: %d \n Description: if kernel is not a valid kernel object.", argument_index, function_name, line_code);
			break;

		case CL_INVALID_ARG_INDEX:
			mf_fatal_err("CL_INVALID_ARG_INDEX argument: %d function: %s code line: %d \n Description: if arg_index is not a valid argument index.", argument_index, function_name, line_code);
			break;

		case CL_INVALID_ARG_VALUE:
			mf_fatal_err("CL_INVALID_ARG_VALUE  argument: %d function: %s code line: %d \n Description: if arg_value specified is NULL for an argument that is not declared with the __local qualifier or vice-versa.", argument_index, function_name, line_code);
			break;

		case CL_INVALID_MEM_OBJECT:
			mf_fatal_err("CL_INVALID_MEM_OBJECT argument: %d function: %s code line: %d \n Description: for an argument declared to be a memory object when the specified arg_value is not a valid memory object.", argument_index, function_name, line_code);
			break;

		case CL_INVALID_SAMPLER :
			mf_fatal_err("CL_INVALID_SAMPLER  argument: %d function: %s code line: %d \n Description: for an argument declared to be of type sampler_t when the specified arg_value is not a valid sampler object.", argument_index, function_name, line_code);
			break;

		case CL_INVALID_ARG_SIZE :
			mf_fatal_err("CL_INVALID_ARG_SIZE  argument: %d function: %s code line: %d \n Description: if arg_size does not match the size of the data type for an argument that is not a memory object or if the argument is a memory object and arg_size != sizeof(cl_mem) or if arg_size is zero and the argument is declared with the __local qualifier or if the argument is a sampler and arg_size != sizeof(cl_sampler).", argument_index, function_name, line_code);
			break;

		case CL_SUCCESS:
			mf_log_info("SUCCESS  argument: %d function: %s code line: %d ", argument_index, function_name, line_code);
			break;

		default:
			mf_log_warn("UNKNOWN FLAG");
			break;
	};
}

#ifdef __cplusplus
}
#endif
