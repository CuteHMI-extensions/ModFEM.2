/*
 * utd_threads.h
 *
 *  Created on: 31-03-2012
 *      Author: Paweł Macioł
 *  Description:
 *  	Simple library for creation and managment of threads.
 */

#ifndef UTD_THREADS_H_
#define UTD_THREADS_H_

/* Defines */
#ifdef _WIN32
  // Windows threads use different data structures
  #include <windows.h>
  typedef HANDLE   utd_thread;
  typedef unsigned (WINAPI *utd_thread_routine)(void *);
  typedef CRITICAL_SECTION utd_critical_section;

  #define ENTERCRITICALSECTION(crit_section) EnterCriticalSection(&(crit_section));
  #define LEAVECRITICALSECTION(crit_section) LeaveCriticalSection(&(crit_section));

#else

  // Includes POSIX thread headers for Linux thread support
  #include <pthread.h>
  typedef pthread_t utd_thread;
  typedef void *(*utd_thread_routine)(void *);
  typedef pthread_mutex_t utd_critical_section;

  #define ENTERCRITICALSECTION(crit_section) pthread_mutex_lock  (&(crit_section));
  #define LEAVECRITICALSECTION(crit_section) pthread_mutex_unlock(&(crit_section));

#endif


#ifdef __cplusplus
    extern "C" {
#endif

    extern int utd_internal_error;

    //Create thread.
    extern utd_thread utd_start_thread(utd_thread_routine, void *data);

    //Wait for thread to finish.
    extern void utd_end_thread(utd_thread thread);

    //Destroy thread.
    extern void utd_destroy_thread(utd_thread thread);

    //Wait for multiple threads.
    extern void utd_wait_for_threads(const utd_thread *threads, int num);

#ifdef __cplusplus
} //extern "C"
#endif


#endif /* UTD_THREADS_H_ */
