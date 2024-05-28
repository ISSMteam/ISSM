/*!\file:  LaunchThread.cpp
 * \brief  launch thread in a generic way, covering single and multi-threaded cases
 * This routine attempts to simplify management of multi-threading. When multi-threadeing 
 * is not requested (serial run), LaunchThread will just call the function (provided in argument 
 * list), nothing fancy there.  If multi-threading is requested, LaunchThread will launch the 
 * function on multiple threads (num_threads of them), and provide these functions with the 
 * local data they need (folded in the "gate" structure) + the thread id + the number of threads 
 * All this info is collected in the pthread_handle structure. 
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef _MULTITHREADING_
#include <pthread.h>
#endif

#include "./issm_threads.h"
#include "../MemOps/MemOps.h"
#include "../Exceptions/exceptions.h"

void LaunchThread(void* function(void*), void* gate,int num_threads){

	#ifdef _MULTITHREADING_
	int i;
	int            *status  = NULL;
	pthread_t      *threads = NULL;
	pthread_handle *handles = NULL;

	/*check number of threads*/
	if(num_threads<1)    _error_("number of threads must be at least 1");
	if(num_threads>2000) _error_("number of threads seems to be too high ("<<num_threads<<">2000)");

	/*dynamically allocate: */
	threads=xNew<pthread_t>(num_threads);
	handles=xNew<pthread_handle>(num_threads);

	for(i=0;i<num_threads;i++){
		handles[i].gate=gate;
		handles[i].id=i;
		handles[i].num=num_threads;
	}
	for(i=0;i<num_threads;i++){

		if(pthread_create(threads+i,NULL,function,(void*)(handles+i))){
			_error_("pthread_create error");
		}
	}
	for(i=0;i<num_threads;i++){
		if(pthread_join(threads[i],(void**)&status)){
			_error_("pthread_join error");
		}
	}

	/*Free resources:*/
	xDelete<pthread_t>(threads);
	xDelete<pthread_handle>(handles);

	#else
	pthread_handle handle;
	handle.gate=gate;
	handle.id=0;
	handle.num=1;

	function((void*)&handle);
	#endif
}
