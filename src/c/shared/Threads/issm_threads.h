/*!\file: issm_threads.h
 * \brief prototypes for issm_threads.h
 */ 

#ifndef _ISSM_THREADS_H_
#define  _ISSM_THREADS_H_

/*structure that holds the local data for each thread (in the gate), 
 * + the thread specific information (my id + number of threads) : */
typedef struct{
	void* gate;
	int   id;
	int   num;
} pthread_handle;

/*routine that launches "function" in a multi-threaded way if requested, 
 * or just serially if not requested: */
void LaunchThread(void* function(void*), void* gate,int num_threads);
void PartitionRange(int* pi0,int* pi1, int num_el,int num_threads,int my_thread);

#endif //ifndef _ISSM_THREADS_H_
