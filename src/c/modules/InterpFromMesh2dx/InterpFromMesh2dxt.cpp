/*!\file:  InterpFromMesh2dxt.cpp
 * \brief  thread core for InterpFromMesh2dxt code
 */ 

#include "./InterpFromMesh2dx.h"
#include "../../shared/shared.h"

void* InterpFromMesh2dxt(void* vpthread_handle){

	/*intermediary: */
	int     i0,i1,i,j;
	double  area,area_1,area_2,area_3;
	double  data_value;

	/*recover handle and gate: */
	pthread_handle                *handle      = (pthread_handle*)vpthread_handle;
	InterpFromMesh2dxThreadStruct *gate        = (InterpFromMesh2dxThreadStruct*)handle->gate;
	int                            my_thread   = handle->id;
	int                            num_threads = handle->num;

	/*recover parameters :*/
	int     interpolation_type      = gate->interpolation_type;
	bool    debug                   = gate->debug;
	int     nels_data               = gate->nels_data;
	int    *index_data              = gate->index_data;
	double *x_data                  = gate->x_data;
	double *y_data                  = gate->y_data;
	double *data                    = gate->data;
	double  xmin                    = gate->xmin;
	double  xmax                    = gate->xmax;
	double  ymin                    = gate->ymin;
	double  ymax                    = gate->ymax;
	int     nods_prime              = gate->nods_prime;
	IssmSeqVec<IssmPDouble>* data_prime = gate->data_prime;
	double *x_prime                 = gate->x_prime;
	double *y_prime                 = gate->y_prime;
	double *default_values          = gate->default_values;
	int     num_default_values      = gate->num_default_values;
	double *incontour               = gate->incontour;

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,nels_data,num_threads,my_thread);

	/*Loop over the elements*/
	for(i=i0;i<i1;i++){

		/*display current iteration*/
		if (debug && my_thread==0 && fmod((double)i,(double)100)==0)
		 _printf_("\r      interpolation progress: "<<setw(6)<<setprecision(2)<<double(i-i0)/double(i1-i0)*100<<"%   ");

		/*if there is no point inside the domain, go to next iteration*/
		if ( (x_data[index_data[3*i+0]-1]<xmin) && (x_data[index_data[3*i+1]-1]<xmin) && (x_data[index_data[3*i+2]-1]<xmin)) continue;
		if ( (x_data[index_data[3*i+0]-1]>xmax) && (x_data[index_data[3*i+1]-1]>xmax) && (x_data[index_data[3*i+2]-1]>xmax)) continue;
		if ( (y_data[index_data[3*i+0]-1]<ymin) && (y_data[index_data[3*i+1]-1]<ymin) && (y_data[index_data[3*i+2]-1]<ymin)) continue;
		if ( (y_data[index_data[3*i+0]-1]>ymax) && (y_data[index_data[3*i+1]-1]>ymax) && (y_data[index_data[3*i+2]-1]>ymax)) continue;

		/*get area of the current element (Jacobian = 2 * area)*/
		//area =x2 * y3 - y2*x3 + x1 * y2 - y1 * x2 + x3 * y1 - y3 * x1;
		area=x_data[index_data[3*i+1]-1]*y_data[index_data[3*i+2]-1]-y_data[index_data[3*i+1]-1]*x_data[index_data[3*i+2]-1]
		  +  x_data[index_data[3*i+0]-1]*y_data[index_data[3*i+1]-1]-y_data[index_data[3*i+0]-1]*x_data[index_data[3*i+1]-1]
		  +  x_data[index_data[3*i+2]-1]*y_data[index_data[3*i+0]-1]-y_data[index_data[3*i+2]-1]*x_data[index_data[3*i+0]-1];

		/*loop over the prime nodes*/
		for (j=0;j<nods_prime;j++){

			if(incontour[j]){

				/*Get first area coordinate = det(x-x3  x2-x3 ; y-y3   y2-y3)/area*/
				area_1=((x_prime[j]-x_data[index_data[3*i+2]-1])*(y_data[index_data[3*i+1]-1]-y_data[index_data[3*i+2]-1]) 
						-  (y_prime[j]-y_data[index_data[3*i+2]-1])*(x_data[index_data[3*i+1]-1]-x_data[index_data[3*i+2]-1]))/area;
				/*Get second area coordinate =det(x1-x3  x-x3 ; y1-y3   y-y3)/area*/
				area_2=((x_data[index_data[3*i+0]-1]-x_data[index_data[3*i+2]-1])*(y_prime[j]-y_data[index_data[3*i+2]-1]) 
						- (y_data[index_data[3*i+0]-1]-y_data[index_data[3*i+2]-1])*(x_prime[j]-x_data[index_data[3*i+2]-1]))/area;
				/*Get third area coordinate = 1-area1-area2*/
				area_3=1-area_1-area_2;

				/*is the current point in the current element?*/
				if (area_1>=-1.e-8 && area_2>=-1.e-8 && area_3>=-1.e-8){

					/*Yes ! compute the value on the point*/
					if (interpolation_type==1){
						/*nodal interpolation*/
						data_value=area_1*data[index_data[3*i+0]-1]+area_2*data[index_data[3*i+1]-1]+area_3*data[index_data[3*i+2]-1];
					}
					else{
						/*element interpolation*/
						data_value=data[i];
					}
					if (xIsNan<IssmPDouble>(data_value)){
						if(num_default_values==1) data_value=default_values[0];
						else data_value=default_values[j];
					}

					/*insert value and go to the next point*/
					data_prime->SetValue(j,data_value,INS_VAL);
				}
			}
		}
	}
	if(debug && my_thread==0)
	 _printf_("\r      interpolation progress: "<<fixed<<setw(6)<<setprecision(2)<<100.<<"%  \n");
	return NULL;
}
