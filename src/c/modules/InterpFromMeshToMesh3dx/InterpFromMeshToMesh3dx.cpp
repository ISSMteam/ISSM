/*!\file:  InterpFromMeshToMesh3dx.cpp
 * \brief  "c" core code for interpolating values from a structured grid.
 */ 

#include "./InterpFromMeshToMesh3dx.h"
#include "../../shared/shared.h"

int InterpFromMeshToMesh3dx( IssmSeqVec<IssmPDouble>** pdata_prime,double* index_data, double* x_data, double* y_data, double* z_data, int nods_data,int nels_data, double* data, int data_length, double* x_prime, double* y_prime, double* z_prime, int nods_prime,double default_value) {

	/*Output*/
	IssmSeqVec<IssmPDouble>* data_prime=NULL;

	/*Intermediary*/
	int i,j;
	int interpolation_type;
	bool debug;
	double area;
	double area_1,area_2,area_3;
	double zeta,bed,surface;
	double data_value;
	double x_prime_min,x_prime_max;
	double y_prime_min,y_prime_max;
	double x_tria_min,y_tria_min;
	double x_tria_max,y_tria_max;

	/*some checks*/
	if (nels_data<1 || nods_data<6 || nods_prime==0){
		_error_("nothing to be done according to the mesh given in input");
	}

	/*Set debug to 1 if there are lots of elements*/
	debug=(bool)((double)nels_data*(double)nods_prime >= pow((double)10,(double)9));

	/*figure out what kind of interpolation is needed*/
	if (data_length==nods_data){
		interpolation_type=1;
	}
	else if (data_length==nels_data){
		interpolation_type=2;
	}
	else{
		_error_("length of vector data not supported yet. It should be of length (number of nodes) or (number of elements)!");
	}

	/*Get prime mesh extrema coordinates*/
	x_prime_min=x_prime[0]; x_prime_max=x_prime[0];y_prime_min=y_prime[0]; y_prime_max=y_prime[0];
	for (i=1;i<nods_prime;i++){
		if (x_prime[i]<x_prime_min) x_prime_min=x_prime[i];
		if (x_prime[i]>x_prime_max) x_prime_max=x_prime[i];
		if (y_prime[i]<y_prime_min) y_prime_min=y_prime[i];
		if (y_prime[i]>y_prime_max) y_prime_max=y_prime[i];
	}

	/*Initialize output*/
	data_prime=new IssmSeqVec<IssmPDouble>(nods_prime);
	for (i=0;i<nods_prime;i++) data_prime->SetValue(i,default_value,INS_VAL);

	/*Loop over the elements*/
	for (i=0;i<nels_data;i++){

		/*display current iteration*/
		if (debug && fmod((double)i,(double)100)==0)
		 _printf_("\r      interpolation progress: "<<setw(6)<<setprecision(2)<<double(i)/double(nels_data)*100<<"%   ");

		/*Get extrema coordinates of current elements*/
		x_tria_min=x_data[(int)index_data[6*i+0]-1]; x_tria_max=x_tria_min;
		y_tria_min=y_data[(int)index_data[6*i+0]-1]; y_tria_max=y_tria_min;
		for (j=1;j<3;j++){
			if(x_data[(int)index_data[6*i+j]-1]<x_tria_min) x_tria_min=x_data[(int)index_data[6*i+j]-1];
			if(x_data[(int)index_data[6*i+j]-1]>x_tria_max) x_tria_max=x_data[(int)index_data[6*i+j]-1];
			if(y_data[(int)index_data[6*i+j]-1]<y_tria_min) y_tria_min=y_data[(int)index_data[6*i+j]-1];
			if(y_data[(int)index_data[6*i+j]-1]>y_tria_max) y_tria_max=y_data[(int)index_data[6*i+j]-1];
		}

		/*if there is no point inside the domain, go to next iteration*/
		if ( x_prime_max < x_tria_min ) continue; 
		if ( x_prime_min > x_tria_max ) continue; 
		if ( y_prime_max < y_tria_min ) continue; 
		if ( y_prime_min > y_tria_max ) continue; 

		/*get area of the current element (Jacobian = 2 * area)*/
		//area =x2 * y3 - y2*x3 + x1 * y2 - y1 * x2 + x3 * y1 - y3 * x1;
		area=x_data[(int)index_data[6*i+1]-1]*y_data[(int)index_data[6*i+2]-1]-y_data[(int)index_data[6*i+1]-1]*x_data[(int)index_data[6*i+2]-1]
		  +  x_data[(int)index_data[6*i+0]-1]*y_data[(int)index_data[6*i+1]-1]-y_data[(int)index_data[6*i+0]-1]*x_data[(int)index_data[6*i+1]-1]
		  +  x_data[(int)index_data[6*i+2]-1]*y_data[(int)index_data[6*i+0]-1]-y_data[(int)index_data[6*i+2]-1]*x_data[(int)index_data[6*i+0]-1];

		/*loop over the prime nodes*/
		for (j=0;j<nods_prime;j++){

			/*if the current point is not in the triangle, continue*/
			if ( x_prime[j] < x_tria_min ) continue; 
			if ( x_prime[j] > x_tria_max ) continue; 
			if ( y_prime[j] < y_tria_min ) continue; 
			if ( y_prime[j] > y_tria_max ) continue; 

			/*Get first area coordinate = det(x-x3  x2-x3 ; y-y3   y2-y3)/area*/
			area_1=((x_prime[j]-x_data[(int)index_data[6*i+2]-1])*(y_data[(int)index_data[6*i+1]-1]-y_data[(int)index_data[6*i+2]-1]) 
						-  (y_prime[j]-y_data[(int)index_data[6*i+2]-1])*(x_data[(int)index_data[6*i+1]-1]-x_data[(int)index_data[6*i+2]-1]))/area;
			/*Get second area coordinate =det(x1-x3  x-x3 ; y1-y3   y-y3)/area*/
			area_2=((x_data[(int)index_data[6*i+0]-1]-x_data[(int)index_data[6*i+2]-1])*(y_prime[j]-y_data[(int)index_data[6*i+2]-1]) 
						- (y_data[(int)index_data[6*i+0]-1]-y_data[(int)index_data[6*i+2]-1])*(x_prime[j]-x_data[(int)index_data[6*i+2]-1]))/area;
			/*Get third area coordinate = 1-area1-area2*/
			area_3=1-area_1-area_2;

			/*is the current point in the current 2d element?*/
			if (area_1>=0 && area_2>=0 && area_3>=0){

				/*compute bottom and top height of the element at this 2d position*/
				bed    =area_1*z_data[(int)index_data[6*i+0]-1]+area_2*z_data[(int)index_data[6*i+1]-1]+area_3*z_data[(int)index_data[6*i+2]-1];
				surface=area_1*z_data[(int)index_data[6*i+3]-1]+area_2*z_data[(int)index_data[6*i+4]-1]+area_3*z_data[(int)index_data[6*i+5]-1];

				/*Compute zeta*/
				zeta=2*(z_prime[j]-bed)/(surface-bed)-1;

				if (zeta >=-1 && zeta<=1){
					if (interpolation_type==1){
						/*nodal interpolation*/
						data_value=(1-zeta)/2*(area_1*data[(int)index_data[6*i+0]-1]+area_2*data[(int)index_data[6*i+1]-1]+area_3*data[(int)index_data[6*i+2]-1]) 
						  +        (1+zeta)/2*(area_1*data[(int)index_data[6*i+3]-1]+area_2*data[(int)index_data[6*i+4]-1]+area_3*data[(int)index_data[6*i+5]-1]);
					}
					else{
						/*element interpolation*/
						data_value=data[i];
					}
					if (xIsNan<IssmPDouble>(data_value)) data_value=default_value;

					/*insert value and go to the next point*/
					data_prime->SetValue(j,data_value,INS_VAL);
				}
			}
		}
	}
	if (debug)
	 _printf_("\r      interpolation progress: "<<fixed<<setw(6)<<setprecision(2)<<100.<<"%  \n");

	/*Assign output pointers:*/
	*pdata_prime=data_prime;
	return 1;
}
