/************************************************************************
File tmh_ocl.h - OpenCL interface to thread management routines

Contains declarations of constants, types and interface routines:
  tmr_ocl_create_command_queues - for a selected platform and device
  tmr_ocl_create_kernel - for a selected platform, device and kernel index
  tmr_ocl_get_current_device_type 
//   returns the current OpenCL device type specified based on compiler switches
  tmr_ocl_select_device 
//   returns OpenCL device index (for local data structures) or -1 if device
//   is not available (not existing or not serviced) for the specified platform
  tmr_ocl_get_current_device() - to return the index of current OpenCL device
  tmr_ocl_device_type - returns LOCAL device type (INTEGER)
  tmr_ocl_select_context - selects context for a given device
  tmr_ocl_select_command_queue - selects command queue for a given device
  tmr_ocl_select_kernel - selects kernel for a given device and kernel index
  tmr_ocl_cleanup - discards created OpenCL resources

Function for printing debug info:
  tmr_kernel_execution_log - print info about error returned by clEnqueueNDRangeKernel 

------------------------------  			
History:        
	02.2013 - Krzysztof Banas, initial version

*************************************************************************/

#ifndef TMH_MTH_OPENCL_H
#define TMH_MTH_OPENCL_H

#include<stdlib.h>
#include<stdio.h>

#include <CL/cl.h>

/* Constants */

#define TMC_OCL_ALL_PLATFORMS -1 // for future use

#define TMC_OCL_ALL_DEVICES -1 // for future use
#define TMC_OCL_DEVICE_CPU 0 // (for OPENCL_CPU compile time switch)
#define TMC_OCL_DEVICE_GPU 1 // (for OPENCL_GPU compile time switch)
#define TMC_OCL_DEVICE_ACCELERATOR 2 // (for OPENCL_PHI compile time switch)

#define TMC_OCL_MAX_NUM_KERNELS 10

// index for storing numerical integration kernel in data structure
#define TMC_OCL_KERNEL_NUM_INT_INDEX  0
#define TMC_OCL_KERNEL_CRS_FINALIZE 1
#define TMC_OCL_KERNEL_SOLVE 2
#define TMC_OCL_KERNEL_REWRITE_GEO 3
// place for further kernels

#ifdef __cplusplus
extern "C" {
#endif


/* Datatypes */
typedef struct {
  //  char name[128];
  cl_device_id id;
  int context_index;
  int tmc_type;
  cl_device_type type;
  double global_mem_bytes; // in B
  double global_max_alloc; // in B
  double shared_mem_bytes; // in B
  double constant_mem_bytes; // in B
  double cache_bytes; // in B
  int cache_line_bytes; // in B
  int max_num_comp_units;
  int max_work_group_size;
  cl_command_queue command_queue;
  int number_of_kernels;
  cl_program program[TMC_OCL_MAX_NUM_KERNELS];
  cl_kernel kernel[TMC_OCL_MAX_NUM_KERNELS];
} tmt_ocl_device_struct;

typedef struct {
  //  char name[128];
  cl_platform_id id;
  // cl_uint number_of_devices;
  int number_of_devices;
  tmt_ocl_device_struct *list_of_devices;
  cl_context list_of_contexts[3]; // always: [0]-CPU, [1]-GPU, [2]-accel
} tmt_ocl_platform_struct;

typedef struct {
  //cl_uint preferred_alignment = 16;    
  //cl_uint number_of_platforms;
  int number_of_platforms;
  tmt_ocl_platform_struct* list_of_platforms;
  int current_platform_index; // we always choose the first platform?
  int current_device_type;
  int current_device; // selected device for the platform
  char kernel_source_directory;
} tmt_ocl_struct;


/* A single global variable for storing data for OpenCL management*/
extern tmt_ocl_struct tmv_ocl_struct;


/* Declarations of interface routines: */
/**--------------------------------------------------------
  tmr_ocl_init - to initialize OpenCL thread management
---------------------------------------------------------*/
extern int tmr_ocl_init(
  char* Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output,
  int Control,   // not used
  int Monitor
);


/**--------------------------------------------------------
//  tmr_ocl_create_contexts - to create OpenCL contexts on a selected platform
---------------------------------------------------------*/
extern int tmr_ocl_create_contexts(
  FILE *Interactive_output, /* file or stdout to write messages */
  int Chosen_platform_id,
  int Monitor
  );

/**--------------------------------------------------------
  tmr_ocl_create_command_queues - for a selected platform and device
---------------------------------------------------------*/
extern int tmr_ocl_create_command_queues(
    FILE *Interactive_output, /* file or stdout to write messages */
    int Chosen_platform_index,
    int Chosen_device_type,
    int Monitor
  );

/**--------------------------------------------------------
  tmr_ocl_create_kernel - for a selected platform, device and kernel index
---------------------------------------------------------*/
extern int tmr_ocl_create_kernel(
  FILE *Interactive_output, /* file or stdout to write messages */
  int Platform_index,
  int Device_index,
  int Kernel_index,
  char* Kernel_name,
  const char* FileName,
  int Monitor
);

/*---------------------------------------------------------
// auxiliary procedure for reading source files
----------------------------------------------------------*/
extern char* tmr_ocl_readSource(
  FILE *Interactive_output, /* file or stdout to write messages */
  const char* kernelPath
			 );

/**--------------------------------------------------------
tmr_ocl_get_current_platform_index() - to return the index of current OpenCL platform
---------------------------------------------------------*/
extern int tmr_ocl_get_current_platform_index();

/**--------------------------------------------------------
tmr_ocl_get_current_device_type() - to return the type of current OpenCL device
(selected based on compiler switches)
---------------------------------------------------------*/
extern int tmr_ocl_get_current_device_type();

/**--------------------------------------------------------
tmr_ocl_get_current_device() - to return the index of current OpenCL device
(selected based on compiler switches)
---------------------------------------------------------*/
extern int tmr_ocl_get_current_device();


/**--------------------------------------------------------
  tmr_ocl_select_device 
//   returns OpenCL device index (for local data structures) or -1 if device
//   is not available (not existing or not serviced) for the specified platform
---------------------------------------------------------*/
extern int tmr_ocl_select_device(
  FILE* Interactive_output, /* file or stdout to write messages */
  int Platform_index,
  int Device_tmc_type
			  );

/**--------------------------------------------------------
  tmr_ocl_device_type - returns LOCAL device type (INTEGER)
---------------------------------------------------------*/
extern int tmr_ocl_device_type(
  int Platform_index,
  int Device_index
);

/**---------------------------------------------------------
  tmr_ocl_select_device_id - selects command queue for a given device
---------------------------------------------------------*/
extern cl_device_id tmr_ocl_select_device_id(
  int Platform_index,
  int Device_index
);

/**--------------------------------------------------------
  tmr_ocl_select_context - selects context for a given device
---------------------------------------------------------*/
extern cl_context tmr_ocl_select_context(
  int Platform_index,
  int Device_index
);

/**--------------------------------------------------------
  tmr_ocl_select_command_queue - selects command queue for a given device
---------------------------------------------------------*/
extern cl_command_queue tmr_ocl_select_command_queue(
  int Platform_index,
  int Device_index
);

/**--------------------------------------------------------
  tmr_ocl_select_kernel - selects kernel for a given device and kernel index
---------------------------------------------------------*/
extern cl_kernel tmr_ocl_select_kernel(
  int Platform_index,
  int Device_index,
  int Kernel_index
);


/**--------------------------------------------------------
  tmr_ocl_cleanup - discards created OpenCL resources
---------------------------------------------------------*/
extern void tmr_ocl_cleanup();

  
/*---------------------------------------------------------
tmr_kernel_execution_log - print info about error returned by clEnqueueNDRangeKernel 
---------------------------------------------------------*/
extern void tmr_kernel_execution_log(cl_int error_code, int line_code, const char *function_name);

/*---------------------------------------------------------
tmr_set_kernel_arguments_log - print info about error returned by clSetKernelArg
---------------------------------------------------------*/
extern void tmr_set_kernel_arguments_log(cl_int error_code, int argument_index, int line_code, const char *function_name);


#ifdef __cplusplus
}
#endif


#endif
