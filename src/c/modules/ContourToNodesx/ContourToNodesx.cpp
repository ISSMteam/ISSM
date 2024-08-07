/*! \file  ContourToNodesx.c
 */

#include "./ContourToNodesx.h"

int ContourToNodesx(IssmPDouble** pflags,double* x, double* y, int nods, Contour<IssmPDouble>** contours,int numcontours,int edgevalue){

	/*output: */
	IssmPDouble* flags = xNewZeroInit<IssmPDouble>(nods);

	/*Loop through all contours: */
	for(int i=0;i<numcontours;i++){
		Contour<IssmPDouble>* contouri=*(contours+i);
		int numnodes=contouri->nods;
		IssmPDouble* xc=contouri->x;
		IssmPDouble* yc=contouri->y;
		IsInPoly(flags,xc,yc,numnodes,x,y,0,nods,edgevalue);
	}

	/*Assign output pointers: */
	*pflags=flags;
	return 1;
}

int ContourToNodesx(IssmPDouble** pflags,double* x, double* y, int nods, Contours* contours, int edgevalue){

	/*output: */
	IssmPDouble* flags = xNewZeroInit<IssmPDouble>(nods);

	/*Loop through all contours: */
	if(contours){
		/*initialize thread parameters: */
		ContourToNodesThreadStruct gate;
		gate.contours  = contours;
		gate.nods      = nods;
		gate.edgevalue = edgevalue;
		gate.in_nod    = flags;
		gate.x         = x;
		gate.y         = y;

		/*launch the thread manager with ContourToMeshxt as a core: */
		LaunchThread(ContourToNodesxt,(void*)&gate,_NUMTHREADS_);
	}

	/*Assign output pointers: */
	*pflags=flags;
	return 1;
}

void* ContourToNodesxt(void* vpthread_handle){

	/*gate variables :*/
	ContourToNodesThreadStruct *gate    = NULL;
	pthread_handle             *handle  = NULL;
	int  i1,i0;

	/*recover handle and gate: */
	handle          = (pthread_handle*)vpthread_handle;
	gate            = (ContourToNodesThreadStruct*)handle->gate;
	int my_thread   = handle->id;
	int num_threads = handle->num;

	/*recover parameters :*/
	Contours* contours  = gate->contours;
	int       nods      = gate->nods;
	int       edgevalue = gate->edgevalue;
	double   *in_nod    = gate->in_nod;
	double   *x         = gate->x;
	double   *y         = gate->y;

	/*distribute indices across threads :*/
	PartitionRange(&i0,&i1,nods,num_threads,my_thread);

	/*Loop through all contours: */
	for(Object* & object : contours->objects){
		Contour<double>* contour=(Contour<double>*)object;
		IsInPoly(in_nod,contour->x,contour->y,contour->nods,x,y,i0,i1,edgevalue);
	}

	return NULL;
}
