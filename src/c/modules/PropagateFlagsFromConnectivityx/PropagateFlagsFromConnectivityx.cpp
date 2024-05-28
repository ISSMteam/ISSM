/*!\file PropagateFlagsFromConnectivityx
 */

#include "./PropagateFlagsFromConnectivityx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void RecursivePropagation(double* pool, double* connectivity,int index, double* flags);

void PropagateFlagsFromConnectivityx( double* pool, double* connectivity,int index, double* flags){

	/*Call recursive propagation routine: */
	RecursivePropagation(pool, connectivity,index, flags);
}

void RecursivePropagation(double* pool, double* connectivity, int index, double* flags){

	int i;
	int newel;

	/*if this element (index) belongs to the pool already, skip: */
	if(pool[index-1])return;

	/*if this element does not belong to the flags set, skip: */
	if(flags[index-1]==0)return;

	/*put this element (index), which belongs to the flags, into the pool: */
	pool[index-1]=1;

	/*now, propagate recursively using connectivity of this element: */
	for(i=0;i<3;i++){
		newel=(int)*(connectivity+(index-1)*3+i);
		RecursivePropagation(pool, connectivity, newel, flags);
	}
}
