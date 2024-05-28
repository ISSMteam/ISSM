/*!\file:  InterpFromMeshToGridx.cpp
 * \brief  "c" core code for interpolating values from a structured grid.
 */ 

#include "./InterpFromMeshToGridx.h"
#include "../../shared/shared.h"

void InterpFromMeshToGridx(double** pgriddata,int* index_mesh, double* x_mesh, double* y_mesh, int nods,int nels, double* data_mesh,int data_length,double* x_grid,double* y_grid,int nlines,int ncols,double default_value){

	/*Intermediary*/
	int    i1,i2,j1,j2;
	int    interpolation_type;
	double area;
	double area_1,area_2,area_3;
	double x_tria_min,y_tria_min;
	double x_tria_max,y_tria_max;
	double x_grid_min,y_grid_min;
	double x_grid_max,y_grid_max;
	double data_value;

	/*some checks*/
	if(nels<1 || nods<3 || nlines<2 || ncols<2) _error_("nothing to be done according to the mesh given in input");

	/*figure out what kind of interpolation is needed*/
	if(data_length==nods)      interpolation_type=1;
	else if(data_length==nels) interpolation_type=2;
	else _error_("length of vector data not supported yet. It should be of length (number of nodes) or (number of elements)!");

	/*First, prepare output*/
	double* griddata=xNew<double>(nlines*ncols);
	for(int i=0;i<nlines;i++) for(int j=0;j<ncols; j++) griddata[i*ncols+j]=default_value;

	/*Set debug to "true" if there are lots of elements*/
	bool debug=(bool)((double)ncols*nlines*nels >= 5.e10);

	/*figure out if x or y are flipped*/
	bool xflip = false;
	bool yflip = false;
	if(x_grid[1]-x_grid[0]<0) xflip=true;
	if(y_grid[1]-y_grid[0]<0) yflip=true;

	/*Get min/max coordinates of the grid*/
	double xposting = x_grid[1]-x_grid[0];
	double yposting = y_grid[1]-y_grid[0];
	if(xflip){
		x_grid_min=x_grid[ncols-1];
		x_grid_max=x_grid[0];
	}
	else{
		x_grid_min=x_grid[0];
		x_grid_max=x_grid[ncols-1];
	}
	if(yflip){
		y_grid_min=y_grid[nlines-1];
		y_grid_max=y_grid[0];
	}
	else{
		y_grid_min=y_grid[0];
		y_grid_max=y_grid[nlines-1];
	}

	/*Loop over the elements*/
	for(int n=0;n<nels;n++){

		/*display current iteration*/
		if(debug && n%10000==0)
		 _printf_("\r      interpolation progress: "<<setw(6)<<setprecision(2)<<double(n)/double(nels)*100<<"%   ");

		/*Get extrema coordinates of current elements*/
		x_tria_min=x_mesh[index_mesh[3*n+0]-1]; x_tria_max=x_tria_min;
		y_tria_min=y_mesh[index_mesh[3*n+0]-1]; y_tria_max=y_tria_min;
		for(int i=1;i<3;i++){
			if(x_mesh[index_mesh[3*n+i]-1]<x_tria_min) x_tria_min=x_mesh[index_mesh[3*n+i]-1];
			if(x_mesh[index_mesh[3*n+i]-1]>x_tria_max) x_tria_max=x_mesh[index_mesh[3*n+i]-1];
			if(y_mesh[index_mesh[3*n+i]-1]<y_tria_min) y_tria_min=y_mesh[index_mesh[3*n+i]-1];
			if(y_mesh[index_mesh[3*n+i]-1]>y_tria_max) y_tria_max=y_mesh[index_mesh[3*n+i]-1];
		}

		/*if the current triangle is not in the grid, continue*/
		if( (x_tria_min>x_grid_max) || (x_tria_max<x_grid_min) || (y_tria_min>y_grid_max) || (y_tria_max<y_grid_min) ) continue;

		/*Get indices i and j that form a square around the current triangle*/
		if(yflip){
			i1=max(0,       (int)floor((y_tria_max-y_grid_max)/yposting)-1);
			i2=min(nlines-1,(int)ceil((y_tria_min-y_grid_max)/yposting));
		}
		else{
			i1=max(0,       (int)floor((y_tria_min-y_grid_min)/yposting)-1);
			i2=min(nlines-1,(int)ceil((y_tria_max-y_grid_min)/yposting));
		}
		if(xflip){
			j1=max(0,      (int)floor((x_tria_max-x_grid_max)/xposting)-1);
			j2=min(ncols-1,(int)ceil((x_tria_min-x_grid_max)/xposting));
		}
		else{
			j1=max(0,      (int)floor((x_tria_min-x_grid_min)/xposting)-1);
			j2=min(ncols-1,(int)ceil((x_tria_max-x_grid_min)/xposting));
		}

		/*get area of the current element (Jacobian = 2 * area)*/
		//area =x2 * y3 - y2*x3 + x1 * y2 - y1 * x2 + x3 * y1 - y3 * x1;
		area=x_mesh[index_mesh[3*n+1]-1]*y_mesh[index_mesh[3*n+2]-1]-y_mesh[index_mesh[3*n+1]-1]*x_mesh[index_mesh[3*n+2]-1]
		  +  x_mesh[index_mesh[3*n+0]-1]*y_mesh[index_mesh[3*n+1]-1]-y_mesh[index_mesh[3*n+0]-1]*x_mesh[index_mesh[3*n+1]-1]
		  +  x_mesh[index_mesh[3*n+2]-1]*y_mesh[index_mesh[3*n+0]-1]-y_mesh[index_mesh[3*n+2]-1]*x_mesh[index_mesh[3*n+0]-1];

		/*Go through x_grid and y_grid and interpolate if necessary*/
		for(int i=i1;i<=i2;i++){

			//exit if y not between y_tria_min and y_tria_max
			if((y_grid[i]>y_tria_max) || (y_grid[i]<y_tria_min)) continue;

			for(int j=j1;j<=j2; j++){

				//exit if x not between x_tria_min and x_tria_max
				if((x_grid[j]>x_tria_max) || (x_grid[j]<x_tria_min)) continue;

				/*Get first area coordinate = det(x-x3  x2-x3 ; y-y3   y2-y3)/area*/
				area_1=((x_grid[j]-x_mesh[index_mesh[3*n+2]-1])*(y_mesh[index_mesh[3*n+1]-1]-y_mesh[index_mesh[3*n+2]-1]) 
							-  (y_grid[i]-y_mesh[index_mesh[3*n+2]-1])*(x_mesh[index_mesh[3*n+1]-1]-x_mesh[index_mesh[3*n+2]-1]))/area;
				/*Get second area coordinate =det(x1-x3  x-x3 ; y1-y3   y-y3)/area*/
				area_2=((x_mesh[index_mesh[3*n+0]-1]-x_mesh[index_mesh[3*n+2]-1])*(y_grid[i]-y_mesh[index_mesh[3*n+2]-1]) 
							- (y_mesh[index_mesh[3*n+0]-1]-y_mesh[index_mesh[3*n+2]-1])*(x_grid[j]-x_mesh[index_mesh[3*n+2]-1]))/area;
				/*Get third area coordinate = 1-area1-area2*/
				area_3=1.-area_1-area_2;

				/*is the current point in the current element?*/
				if(area_1>-1e-12 && area_2>-1e-12 && area_3>-1e-12){

					/*Yes ! compute the value on the point*/
					if(interpolation_type==1){
						/*nodal interpolation*/
						data_value=area_1*data_mesh[index_mesh[3*n+0]-1]+area_2*data_mesh[index_mesh[3*n+1]-1]+area_3*data_mesh[index_mesh[3*n+2]-1];
					}
					else{
						/*element interpolation*/
						data_value=data_mesh[n];
					}
					if (xIsNan<IssmDouble>(data_value)) data_value=default_value;

					/*insert value and go to the next point*/
					griddata[i*ncols+j]=data_value;
				}
			}
		}
	}
	if(debug) _printf_("\r      interpolation progress: "<<setw(6)<<setprecision(2)<<fixed<<100.<<"%   \n");

	/*Assign output pointers:*/
	*pgriddata=griddata;
}
