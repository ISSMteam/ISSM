/*!\file:  PointCloudFindNeighborst.cpp
 * \brief  thread core for PointCloudFindNeighborst code
 */ 

#include "./PointCloudFindNeighborsx.h"
#include "../../shared/shared.h"

void* PointCloudFindNeighborsxt(void* vpthread_handle){

	/*gate variables :*/
	PointCloudFindNeighborsThreadStruct* gate=NULL;
	pthread_handle* handle=NULL;
	int     my_thread;
	int     num_threads;
	double* x;
	double* y;
	int     nods;
	double  mindistance;
	IssmSeqVec<IssmPDouble>*     flags;

	/*recover handle and gate: */
	handle=(pthread_handle*)vpthread_handle;
	gate=(PointCloudFindNeighborsThreadStruct*)handle->gate;
	my_thread=handle->id;
	num_threads=handle->num;

	/*recover parameters :*/
	x=gate->x;
	y=gate->y;
	nods=gate->nods;
	mindistance=gate->mindistance;
	flags=gate->flags;

	/*intermediary: */
	int i,j;
	int i0,i1;
	double distance;
	bool* already=NULL;

	/*allocate: */
	already=xNewZeroInit<bool>(nods);

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,nods,num_threads,my_thread);

	/*Loop over the nodes*/
	for (i=i0;i<i1;i++){

		/*display current iteration*/
		if (my_thread==0 && fmod((double)i,(double)100)==0)
		 _printf_("\r      loop progress: "<<setw(6)<<setprecision(2)<<double(i-i0)/double(i1-i0)*100<<"%   ");

		distance=mindistance+100; //make sure initialization respects min distance criterion.
		for (j=0;j<nods;j++){

			/*skip himself: */
			if (j==i)continue;
			distance=sqrt(pow(x[i]-x[j],2)+ pow(y[i]-y[j],2));

			if(distance<=mindistance){

				/*insert value and go to the next point*/
				if(!already[i]) flags->SetValue(i,1,INS_VAL);
				if(!already[j]) flags->SetValue(j,2,INS_VAL);
				already[i]=true;
				already[j]=true;
				break;
			}
		}
	}
	if (my_thread==0)
	 _printf_("\r      loop progress: "<<fixed<<setw(6)<<setprecision(2)<<100.<<"%  \n");

	/*Free resources:*/
	xDelete<bool>(already);

	return NULL;
}
