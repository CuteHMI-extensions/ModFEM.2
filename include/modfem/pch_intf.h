/************************************************************************
File pch_intf.h - interface of the parallel communication library

  pcr_init_parallel - to initialize parallel communication library
  pcr_print_master - to return the id of read/write master processor
  pcr_send_buffer_open - to open a buffer for sending data
  pcr_buffer_pack_int - to pack an array of integers to a buffer
  pcr_buffer_pack_double - to pack an array of doubles to a buffer
  pcr_buffer_pack_char - to pack an array of characters to a buffer
  pcr_buffer_send - to send a buffer to destination processor
  !pcr_buffer_start_send - to start sending a buffer to destination processor
  !pcr_buffer_finish_send - to finish sending a buffer to destination processor
  !pcr_pack_send_int - to pack and send an array of integers
  !pcr_receive_int - to receive an array of integers
  !pcr_pack_send_double - to pack and send an array of doubles 
  !pcr_receive_double - to receive an array of doubles 
  pcr_buffer_bcast - to broadcast the content of a buffer 
  pcr_buffer_receive - to receive a buffer 
  !pcr_buffer_start_receive - to start receiving a buffer
  !pcr_buffer_finish_receive - to finish receiving a buffer 
  pcr_buffer_unpack_int - to unpack an array of integers from a buffer
  pcr_buffer_unpack_double - to unpack an array of doubles from a buffer
  pcr_buffer_unpack_char - to unpack an array of characters from a buffer
  pcr_recv_buffer_close - to close a buffer for receiveing data
  pcr_bcast_double - to broadcast an array of doubles
  pcr_bcast_int - to broadcast an array of integers
  pcr_bcast_char - to broadcast an array of characyers
  pcr_allreduce_sum_int - to reduce (sum) and broadcast an array of integers 
  pcr_allreduce_sum_double - to reduce (sum) and broadcast an array of doubles
  pcr_exit_parallel - to exit parallel communication library
  pcr_barrier - to synchronise paralell execution
  pcr_is_parallel_initialized - to determine wheter parallelism is initialized or not
  pcr_my_proc_id - to get my proc id
  pcr_nr_proc - to get number of all processes
------------------------------  			
History:        
10.2012 - Kazimierz Michalik
02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _pch_intf_
#define _pch_intf_

#include "stdint.h"

#ifdef __cplusplus
extern "C"{
#endif

/** @defgroup PCM Parallel Communication
 *
 *  @{
 */

/*** CONSTANTS ***/
extern const int PCC_ANY_PROC; /** wildcard for arbitrary processor id */
extern const int PCC_USE_CURRENT_BUFFER;
extern const int PCC_MASTER_PROC_ID; // We assume that this is the lowest process id.
extern const int PCC_OK; // return code indicating correct execution
extern const int PCC_DEFAULT_BUFFER_SIZE;
extern const int PCC_DEFAULT_COMM;
extern const int PCC_NEIGH_GROUP;
  

/*** GLOBAL VARIABLES ***/
extern int pcv_nr_proc;
extern int pcv_my_proc_id;


/*** FUNCTION DECLARATIONS - headers for external functions ***/

/**--------------------------------------------------------
  pcr_init_parallel - to initialize parallel communication library
---------------------------------------------------------*/
extern int pcr_init_parallel( /** returns: >=0 - success code, <0 - error code */
  int* argc,       /** in: */
  char** argv,     /** in: */
  char* Work_dir,  /** in: */
  char *interactive_output_name, /** in/out: */
  FILE **interactive_output_p, /** out: reset pointer to interactive output */
  int* Nr_pr,      /** out: number of processors/processes */
  int* My_id       /** out: local process ID */
  );

/**--------------------------------------------------------
pcr_print_master - to return the id of read/write master processor
---------------------------------------------------------*/
extern int pcr_print_master(
  /** returns the id of read/write master processor or error code <0 */
);

/**--------------------------------------------------------
  pcr_send_buffer_open - to open a buffer for sending data
---------------------------------------------------------*/
extern int pcr_send_buffer_open( 
  /** returns: >=0 - message buffer ID = success code, <0 - error code */
  int Message_id,  /** in: message ID */
  int Buffer_size  /** in: message buffer size (0 for DEFAULT_BUFFER_SIZE) */
  );

/**--------------------------------------------------------
  pcr_buffer_pack_int - to pack an array of integers to a buffer
---------------------------------------------------------*/
extern int pcr_buffer_pack_int( 
                  /** returns: >=0 - success code, <0 - error code */
  const int Message_id,  /** in: message ID */
  const int Buffer_id,  /** in: buffer ID */
  const int Nr_num,     /** in: number of values to pack to the buffer */
  const int* Numbers    /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
  pcr_buffer_pack_double - to pack an array of doubles to a buffer
---------------------------------------------------------*/
extern int pcr_buffer_pack_double( 
                  /** returns: >=0 - success code, <0 - error code */
  const int Message_id,  /** in: message ID */
  const int Buffer_id,  /** in: message ID */
  const int Nr_num,     /** in: number of values to pack to the buffer */
  const double* Numbers /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
  pcr_buffer_pack_char - to pack an array of chars to a buffer
---------------------------------------------------------*/
int pcr_buffer_pack_char( /** returns: >=0 - success code, <0 - error code */
  const int Message_id,  /** in: message ID */
  const int Buffer_id,  /** in: message ID */
  const int Nr_num,     /** in: number of values to pack to the buffer */
  const char* Numbers /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
  pcr_buffer_send - to send a buffer to destination processor
---------------------------------------------------------*/
extern int pcr_buffer_send( 
                  /** returns: >=0 - success code, <0 - error code */
  int Message_id,     /** in: message ID */
  int Buffer_id,  /** in: buffer ID */
  int Dest_proc   /** in : destination processor ID */
  );


/**-------------------------------------------------------
  pcr_buffer_bcast - to broadcast a buffer between all proces
---------------------------------------------------------*/
int pcr_buffer_bcast( /** returns: >=0 - success code, <0 - error code*/
  int Message_id, /** in: message ID */
  int Buffer_id,  /** in: buffer ID */
  int Sender_proc /** in: sender processor ID */
  );

/**--------------------------------------------------------
  pcr_buffer_receive - to receive a buffer from a particular processor
---------------------------------------------------------*/
extern int pcr_buffer_receive( 
		   /** returns: >= 0 - Buffer_id, <0 - error code */
  int Message_id,   /** in: id of the message containing the data */
  int Sender_proc,  /** in : sender processor ID (or PCC_ANY_PROC) */
  int Buffer_size /** in: message buffer size (0 for DEFAULT_BUFFER_SIZE) */
  );


/**--------------------------------------------------------
  pcr_buffer_unpack_int - to unpack an array of integers from a buffer
---------------------------------------------------------*/
extern int pcr_buffer_unpack_int( 
                  /** returns: >=0 - success code, <0 - error code */
  int Message_id,  /** in: message ID */
  int Buffer_id,  /** in: buffer ID */
  int Nr_num,     /** in: number of values to pack to the buffer */
  int* Numbers    /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
  pcr_buffer_unpack_double - to unpack an array of integers from a buffer
---------------------------------------------------------*/
extern int pcr_buffer_unpack_double( 
                  /** returns: >=0 - success code, <0 - error code */
  int Message_id,  /** in: message ID */
  int Buffer_id,  /** in: buffer ID */
  int Nr_num,     /** in: number of values to unpack from the buffer */
  double* Numbers /** in: array of numbers to unpack */
  );

/**--------------------------------------------------------
  pcr_buffer_unpack_char - to unpack an array of integers from a buffer
---------------------------------------------------------*/
int pcr_buffer_unpack_char( /** returns: >=0 - success code, <0 - error code */
  int Message_id,  /** in: message ID */
  int Buffer_id,  /** in: buffer ID */
  int Nr_num,     /** in: number of values to unpack from the buffer */
  char* Numbers /** in: array of numbers to unpack */
  );

/**--------------------------------------------------------
  pcr_recv_buffer_close - to close a buffer for communication
---------------------------------------------------------*/
extern int pcr_recv_buffer_close( 
                 /** returns: >=0 - success code, <0 - error code */
  int Message_id,  /** in: message ID */
  int Buffer_id  /** in: buffer ID */
  );

/**--------------------------------------------------------
  pcr_send_int - to send an array of integers
---------------------------------------------------------*/
extern int pcr_send_int(/** returns: >=0 - success code, <0 - error code */
  const int Dest_proc_id, /** in : destination processor ID */
  const int Message_id,   /** in: id of the message containing the data */
  const int Nr_num,       /** in: number of values to pack to the buffer */
  const int *Numbers      /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
  pcr_send_double - to send an array of doubles 
---------------------------------------------------------*/
extern int pcr_send_double( 
                    /** returns: >=0 - success code, <0 - error code */
  const int Dest_proc_id, /** in : destination processor ID */
  const int Message_id,   /** in: id of the message containing the data */
  const int Nr_num,       /** in: number of values to send */
  const double* Numbers   /** in: array of numbers to send */
  );

/**--------------------------------------------------------
  pcr_receive_int - to receive an array of integers
---------------------------------------------------------*/
extern int pcr_receive_int(/** returns: >=0 - success code, <0 - error code */
  int Sender_proc_id, /** in : destination processor ID */
  int Message_id,   /** in: id of the message containing the data */
  int Nr_num,       /** in: number of values to pack to the buffer */
  int *Numbers      /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
  pcr_receive_double - to receive an array of doubles 
---------------------------------------------------------*/
extern int pcr_receive_double( 
                    /** returns: >=0 - success code, <0 - error code */
  int Sender_proc_id, /** in : destination processor ID */
  int Message_id,   /** in: id of the message containing the data */
  int Nr_num,       /** in: number of values to pack to the buffer */
  double* Numbers   /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
  pcr_bcast_double - to broadcast an array of doubles
---------------------------------------------------------*/
extern int pcr_bcast_double(int Sender_proc_id, /** in : source  processor ID */
  int Nr_num,       /** in: number of values to pack to the buffer */
  double *Numbers   /** in: array of numbers to pack */
  );


/**--------------------------------------------------------
  pcr_bcast_int - to broadcast an array of integers
---------------------------------------------------------*/
extern int pcr_bcast_int(
  int Sender_proc_id, /** in : source processor ID */
  int Nr_num,       /** in: number of values to pack to the buffer */
  int* Numbers      /** in: array of numbers to pack */
  );


/**-------------------------------------------------------
  pcr_bcast_char - to broadcast an array of charegers
---------------------------------------------------------*/
extern int pcr_bcast_char( /** returns: >=0 - success code, <0 - error code*/
  int Sender_proc_id, /** in : source processor ID */
  int Nr_num,       /** in: number of values to pack to the buffer */
  char* Numbers      /** in: array of numbers to pack */
  );

/**--------------------------------------------------------
pcr_allreduce_sum_int - to reduce by summing and broadcast an array of integers 
---------------------------------------------------------*/
extern int pcr_allreduce_sum_int(
  const int Nr_num,          /** in: number of values to reduce */
  const int* Numbers,        /** in: array of numbers to sum */
  int *Numbers_reduced /** in: array of summed numbers */
  );



/**--------------------------------------------------------
pcr_allreduce_sum_double - to reduce by summing and broadcast an array of doubles
---------------------------------------------------------*/
extern int pcr_allreduce_sum_double(
  const int Nr_num,          /** in: number of values to reduce */
  const double* Numbers,        /** in: array of numbers to sum */
   double* Numbers_reduced /** in: array of summed numbers */
  );

/**--------------------------------------------------------
pcr_allreduce_max_int - to reduce by maxing and broadcast an array of integers 
---------------------------------------------------------*/
extern int pcr_allreduce_max_int(
  const int Nr_num,          /** in: number of values to reduce */
  const int* Numbers,        /** in: array of numbers to max */
  int* Numbers_reduced /** in: array of maxed numbers */
  );



/**--------------------------------------------------------
pcr_allreduce_max_double - to reduce by maxing and broadcast an array of doubles
---------------------------------------------------------*/
extern int pcr_allreduce_max_double(
  int Nr_num,          /** in: number of values to reduce */
  double* Numbers,        /** in: array of numbers to max */
  double* Numbers_reduced /** in: array of maxed numbers */
  );

/**--------------------------------------------------------
  pcr_exit_parallel - to exit parallel communication library
---------------------------------------------------------*/
extern int pcr_exit_parallel( /** returns: >=0 - success code, <0 - error code */
  void
  );

//------------------------------------------------------//
// pcr_barrier - to synchronise parallel execution
//------------------------------------------------------//
extern int pcr_barrier( // returns >=0 - success, <0- error code
  void
);

//------------------------------------------------------//
// pcr_is_parallel_initialized - to determine wheter parallelism is initialized or not
//------------------------------------------------------//
extern int pcr_is_parallel_initialized(  //returns 1 - is initialized, 0 - is NOT initialized
  void
);
  
//------------------------------------------------------//
// pcr_is_this_master - to determine wheter parallelism is initialized or not
//------------------------------------------------------//
extern int pcr_is_this_master( //returns 1 - is master, 0 - is NOT master
  void
);

//------------------------------------------------------//
// pcr_my_proc_id - to get id of current process
//------------------------------------------------------//
extern int pcr_my_proc_id( //returns id of current process
  void
);

//------------------------------------------------------//
// pcr_nr_proc - to get number of all processes  
//------------------------------------------------------//
extern int pcr_nr_proc( //returns number of all process
  void
);

//------------------------------------------------------//
// pcr_allgather_int - gathers together values from a group of processors and distributes to all
//------------------------------------------------------//
extern int pcr_allgather_int( //returns success code
                              const int send_values[],
                              const int n_send_values,
                              int gathered_values[],
                              int n_gathered_values

);

/**--------------------------------------------------------
  pcr_send_bytes - to send an array of integers
---------------------------------------------------------*/
extern int pcr_send_bytes(
                    /** returns: >=0 - success code, <0 - error code */
  const int Dest_proc_id, /** in : destination processor ID */
  const int Message_id,   /** in: id of the message containing the data */
  const int Nr_bytes,       /** in: number of values to send */
  const uint8_t* Bytes   /** in: array of numbers to send */
  );

/**--------------------------------------------------------
  pcr_receive_bytes - to receive an array of integers
---------------------------------------------------------*/
extern int pcr_receive_bytes(
                    /** returns: >=0 - success code, <0 - error code */
  int Sender_proc_id, /** in : destination processor ID */
  int Message_id,   /** in: id of the message containing the data */
  int Nr_bytes,       /** in: number of values to pack to the buffer */
  uint8_t* Bytes      /** in: array of numbers to pack */
  );
 
/**--------------------------------------------------------
  pcr_buffer_source - to get the source of message in buffer
---------------------------------------------------------*/
int pcr_buffer_source_id(int Buffer);

/** @} */ // end of group


#ifdef __cplusplus
}
#endif
  
#endif
