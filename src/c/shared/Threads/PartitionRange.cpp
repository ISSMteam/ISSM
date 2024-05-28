/*!\file:  PartitionRange.cpp
 * \brief: return i0,i1, range of local thread.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <math.h>

void PartitionRange(int* pi0,int* pi1, int num_el,int num_threads,int my_thread){

	/*output: */
	int i0=-1;
	int i1=-1;

	int step;
	int i;

	/*distribute elements across threads :*/
	step=(int)floor((double)num_el/(double)num_threads);
	for(i=0;i<(my_thread+1);i++){
		i0=i*step;
		if(i==(num_threads-1))i1=num_el;
		else i1=i0+step;
	}

	/*Assign output pointers:*/
	*pi0=i0;
	*pi1=i1;
}
