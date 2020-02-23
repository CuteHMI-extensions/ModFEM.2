/*
 * utd_threads.c
 *
 *  Created on: 31-03-2012
 *      Author: Paweł Macioł
 * Description:
 * 		Implementation.
 */

#include <modfem/ut_system/utd_threads.h>

int utd_internal_error;

#if _WIN32
	//Create thread
	utd_thread utd_start_thread(utd_thread_routine func, void *data){
		return CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, data, 0, NULL);
	}

	//Wait for thread to finish
	void utd_end_thread(utd_thread thread){
		WaitForSingleObject(thread, INFINITE);
		CloseHandle(thread);
	}

	//Destroy thread
	void utd_destroy_thread(utd_thread thread){
		TerminateThread(thread, 0);
		CloseHandle(thread);
	}

	//Wait for multiple threads
	void utd_wait_for_threads(const utd_thread * threads, int num){
		WaitForMultipleObjects(num, threads, true, INFINITE);

		int i=0; for(; i < num; i++)
			CloseHandle(threads[i]);
	}

#else
	//Create thread
	utd_thread utd_start_thread(utd_thread_routine func, void * data){
		pthread_t thread;
		utd_internal_error = pthread_create(&thread, NULL, func, data);
		return thread;
	}

	//Wait for thread to finish
	void utd_end_thread(utd_thread thread){
		pthread_join(thread, NULL);
	}

	//Destroy thread
	void utd_destroy_thread(utd_thread thread){
		pthread_cancel(thread);
	}

	//Wait for multiple threads
	void utd_wait_for_threads(const utd_thread* threads, int num){
		int i=0; for(; i < num; i++)
			utd_end_thread(threads[i]);
	}

#endif

