/*!\file:  ExpToLevelSetxt.cpp
 * \brief  "thread" core code for figuring out level set value from a contour and a cloud of points.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Include files*/
#include "./ExpToLevelSetx.h"
double minimum_distance(double x1, double y1, double x2, double y2, double x0, double y0);
void   ContourToLevelSet(double* distance,double* contourx, double* contoury, int contournods, double* x, double* y, int i0, int i10);

void* ExpToLevelSetxt(void* vpthread_handle){

	/*gate variables :*/
	ExpToLevelSetxThreadStruct *gate        = NULL;
	pthread_handle             *handle      = NULL;
	int  i,i1,i0;

	/*recover handle and gate: */
	handle          = (pthread_handle*)vpthread_handle;
	gate            = (ExpToLevelSetxThreadStruct*)handle->gate;
	int my_thread   = handle->id;
	int num_threads = handle->num;

	/*recover parameters :*/
	Contours* contours = gate->contours;
	int       nods     = gate->nods;
	double   *distance = gate->distance;
	double   *x        = gate->x;
	double   *y        = gate->y;

	/*distribute indices across threads :*/
	PartitionRange(&i0,&i1,nods,num_threads,my_thread);

	/*Loop through all contours: */
	for(Object* & object : contours->objects){
		Contour<double>* contour=(Contour<double>*)object;
		ContourToLevelSet(distance,contour->x,contour->y,contour->nods,x,y,i0,i1);
	}

	return NULL;
}

void ContourToLevelSet(double* dist,double* contourx, double* contoury, int contournods, double* x, double* y, int i0, int i1){/*{{{*/

	double x0,y0;
	double x1,y1;
	double x2,y2;
	double mind;

	for(int i=i0;i<i1;i++){

      /*Get current point*/
		x0=x[i]; y0=y[i];

		/*Figure out distance from (x0,y0) to contour: */
		mind=1e+50;
		for(int j=0;j<contournods-1;j++){
         /*Get distance from current segment*/
			x1=contourx[j];   y1=contoury[j];
			x2=contourx[j+1]; y2=contoury[j+1];
			mind=min(mind,minimum_distance(x1,y1,x2,y2,x0,y0));
		}
		dist[i]=min(dist[i],mind);
	}
}
double minimum_distance(double x1, double y1, double x2, double y2, double x0, double y0){
	/* Return minimum distance between line segment [(x1,y1) (x2,y2)] and point (x0,y0)
	 * We use the following notations:
	 * segment: v=(x1,y1), w=(x2,y2)
	 * point:   p=(x0,y0)
	 */

   /*Get segment length square (avoid sqrt) |w-v|^2*/
	double l2 = pow(x2-x1,2)+pow(y2-y1,2);

   /*segment is single point: v == w*/
	if(l2 == 0.) return sqrt(pow(x1-x0,2)+pow(y1-y0,2));

	/*Consider the line extending the segment, parameterized as v + t (w - v).
	We find projection of point p onto the line.  It falls where t = [(p-v) . (w-v)] / |w-v|^2*/
	double t = ((x0-x1)*(x2-x1) + (y0-y1)*(y2-y1)) / l2;
	if(t < 0.){
      /*Beyond the 'v' end of the segment*/
      return sqrt(pow(x1-x0,2)+pow(y1-y0,2));
   }
	else if(t > 1.){
      /*Beyond the 'w' end of the segment*/
      return sqrt(pow(x2-x0,2)+pow(y2-y0,2));
   }

   /*Projection falls on segment*/
	double projx= x1 + t* (x2-x1);
	double projy= y1 + t* (y2-y1);
	return sqrt(pow(projx-x0,2)+pow(projy-y0,2));
}
/*}}}*/
