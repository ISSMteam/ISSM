/*!\file:  DistanceToMaskBoundaryxt.cpp
 * \brief  "thread" core code 
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./DistanceToMaskBoundaryx.h"

void* DistanceToMaskBoundaryxt(void* vpthread_handle){

	/*gate variables :*/
	DistanceToMaskBoundaryxThreadStruct *gate        = NULL;
	pthread_handle             *handle      = NULL;
	int  i,j,i0,i1;

	/*recover handle and gate: */
	handle          = (pthread_handle*)vpthread_handle;
	gate            = (DistanceToMaskBoundaryxThreadStruct*)handle->gate;
	int my_thread   = handle->id;
	int num_threads = handle->num;

	/*recover parameters :*/
	int       nods      = gate->nods;
	IssmDouble   *distance    = gate->distance;
	IssmDouble   *x         = gate->x;
	IssmDouble   *y         = gate->y;
	IssmDouble   *mask         = gate->mask;

	/*distribute indices across threads :*/
	PartitionRange(&i0,&i1,nods,num_threads,my_thread);

	/*Loop through vertices: */
	for(i=i0;i<i1;i++){

		IssmDouble d0=1.e+10;
		IssmDouble xi,yi;

		//recover vertex position: 
		xi=x[i];  yi=y[i];

		//figure out if we are inside the mask, or outside: 
		if(mask[i]==1){
			//we are inside, look for nearest vertex that is outside the mask: 
			for(j=0;j<nods;j++){
				if(j==i)continue;
				if (mask[j]==0){
					IssmDouble xj,yj,deltaphi,deltalambda,d;
					xj=x[j]; yj=y[j];
					/*figure  out the distance to xi,yi in lat,long mode, using the greatest circle distance:*/
					deltaphi=fabs(xj-xi); deltalambda=fabs(yj-yi);
					d=2*asin(sqrt(pow(sin(deltaphi/2),2)+cos(xi)*cos(xj)*pow(sin(deltalambda/2),2)));
					if(d<d0)d0=d;
				}
			}
			distance[i]=d0;
		}
		else{
			//we are outside, look for nearest vertex that is inside the mask: 
			for(j=0;j<nods;j++){
				if(j==i)continue;
				if (mask[j]==1){
					IssmDouble xj,yj,deltaphi,deltalambda,d;
					xj=x[j]; yj=y[j];
					/*figure  out the distance to xi,yi in lat,long mode, using the greatest circle distance:*/
					deltaphi=fabs(xj-xi); deltalambda=fabs(yj-yi);
					d=2*asin(sqrt(pow(sin(deltaphi/2),2)+cos(xi)*cos(xj)*pow(sin(deltalambda/2),2)));
					if(d<d0)d0=d;
				}
			}
			distance[i]=d0;
		}
	}

	return NULL;
}
