/*!\file:  ContourToMeshxt.cpp
 * \brief  "thread" core code for interpolating values from a structured grid.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./ContourToMeshx.h"

void* ContourToMeshxt(void* vpthread_handle){

	/*gate variables :*/
	ContourToMeshxThreadStruct *gate        = NULL;
	pthread_handle             *handle      = NULL;
	int  i,i1,i0;

	/*recover handle and gate: */
	handle          = (pthread_handle*)vpthread_handle;
	gate            = (ContourToMeshxThreadStruct*)handle->gate;
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
	for (Object* & object : contours->objects){
		Contour<double>* contour=(Contour<double>*)object;
		IsInPoly(in_nod,contour->x,contour->y,contour->nods,x,y,i0,i1,edgevalue);
	}

	return NULL;
}
