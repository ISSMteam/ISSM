/*!\file:  InterpFromGridToMeshx.cpp
 * \brief  "c" core code for interpolating values from a structured grid.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>
#include "./InterpFromGridToMeshx.h"
#include "../../shared/shared.h"
#include "../../shared/io/io.h"

/*InterpFromGridToMeshx{{{*/
int InterpFromGridToMeshx(IssmSeqVec<IssmPDouble>** pdata_mesh,double* x_in, int x_rows, double* y_in, int y_rows, double* data, int M, int N, double* x_mesh, double* y_mesh, int nods,double default_value,const char* interptype){

	/*output: */
	IssmSeqVec<IssmPDouble>* data_mesh=NULL;

	/*Intermediary*/
	double* x=NULL;
	double* y=NULL;
	int     i;

	/*Some checks on arguments: */
	if ((M<2) || (N<2) || (nods<=0)){
		_error_("nothing to be done according to the dimensions of input matrices and vectors.");
	}
	if (x_in[1]-x_in[0]<0){
		_error_("x coordinate vector should be increasing.\n   use Matlab's command x=flipud(x), also flip the data matrix data=fliplr(data)");
	}
	if (y_in[1]-y_in[0]<0){
		_error_("y coordinate vector should be increasing.\n   use Matlab's command y=flipud(y), also flip the data matrix data=flipud(data)");
	}

	/*Allocate output vector: */
	data_mesh=new IssmSeqVec<IssmPDouble>(nods);

	/*Find out what kind of coordinates (x_in,y_in) have been given is input*/
	if(N==(x_rows-1) && M==(y_rows-1)){

		/*The coordinates given in input describe the contour of each pixel. Take the center of each pixel*/
		x=xNew<double>(N);
		y=xNew<double>(M);
		for(i=0;i<N;i++) x[i]=(x_in[i]+x_in[i+1])/2.;
		for(i=0;i<M;i++) y[i]=(y_in[i]+y_in[i+1])/2.;
		x_rows=x_rows-1;
		y_rows=y_rows-1;
	}
	else if (N==x_rows && M==y_rows){

		/*The coordinates given in input describe the center each pixel. Keep them*/
		x=xNew<double>(N);
		y=xNew<double>(M);
		for(i=0;i<N;i++) x[i]=x_in[i];
		for(i=0;i<M;i++) y[i]=y_in[i];
	}
	else{
		_error_("x and y vectors length should be 1 or 0 more than data number of rows.");
	}

	/*initialize thread parameters: */
	InterpFromGridToMeshxThreadStruct gate;
	gate.x_mesh        = x_mesh;
	gate.y_mesh        = y_mesh;
	gate.x_rows        = x_rows;
	gate.y_rows        = y_rows;
	gate.x             = x;
	gate.y             = y;
	gate.nods          = nods;
	gate.data_mesh     = data_mesh;
	gate.data          = data;
	gate.default_value = default_value;
	gate.interp        = interptype;
	gate.M             = M;
	gate.N             = N;

	/*launch the thread manager with InterpFromGridToMeshxt as a core: */
	LaunchThread(InterpFromGridToMeshxt,(void*)&gate,_NUMTHREADS_);
	//_printf_("\r      interpolation progress: "<<fixed<<setw(6)<<setprecision(2)<<100.<<"%  \n");

	/*Assign output pointers:*/
	*pdata_mesh=data_mesh;
	return 1;
}
/*}}}*/
/*InterpFromGridToMeshxt {{{*/
void* InterpFromGridToMeshxt(void* vpthread_handle){

	/*gate variables :*/
	InterpFromGridToMeshxThreadStruct *gate    = NULL;
	pthread_handle                    *handle  = NULL;
	int my_thread;
	int num_threads;
	int i0,i1;

	/*intermediary: */
	int    i,m,n;
	double x_grid;
	double y_grid;
	double data_value;
	double x1,x2,y1,y2;
	double Q11,Q12,Q21,Q22;

	/*recover handle and gate: */
	handle=(pthread_handle*)vpthread_handle;
	gate=(InterpFromGridToMeshxThreadStruct*)handle->gate;
	my_thread=handle->id;
	num_threads=handle->num;

	/*recover parameters :*/
	double *x_mesh                = gate->x_mesh;
	double *y_mesh                = gate->y_mesh;
	int     x_rows                = gate->x_rows;
	int     y_rows                = gate->y_rows;
	double *x                     = gate->x;
	double *y                     = gate->y;
	int     nods                  = gate->nods;
	IssmSeqVec<IssmPDouble>*data_mesh = gate->data_mesh;
	double *data                  = gate->data;
	double  default_value         = gate->default_value;
	const char* interptype        = gate->interp;
	int     M                     = gate->M;
	int     N                     = gate->N;

	bool debug = M*N>1? true:false;

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,nods,num_threads,my_thread);
	for (i=i0;i<i1;i++) {

		//if(debug && my_thread==0) _printf_("\r      interpolation progress: "<<setw(6)<<setprecision(2)<<double(i-i0)/double(i1-i0)*100<<"%   ");
		x_grid=*(x_mesh+i);
		y_grid=*(y_mesh+i);

		/*Find indices m and n into y and x, for which  y(m)<=y_grids<=y(m+1) and x(n)<=x_grid<=x(n+1)*/
		if(findindices<double>(&n,&m,x,x_rows, y,y_rows, x_grid,y_grid)){

			/*    Q12             Q22
			 * y2 x---------+-----x
			 *    |         |     |
			 *    |         |P    |
			 *    |---------+-----|
			 *    |         |     |
			 *    |         |     |
			 * y1 x---------+-----x Q21
			 *    x1                 x2
			 *
			 */
			x1=x[n]; x2=x[n+1];
			y1=y[m]; y2=y[m+1];
			Q11=data[m*N+n];
			Q12=data[(m+1)*N+n];
			Q21=data[m*N+n+1];
			Q22=data[(m+1)*N+n+1];

			if(strcmp(interptype,"triangle")==0){
				data_value=triangleinterp(x1,x2,y1,y2,Q11,Q12,Q21,Q22,x_grid,y_grid);
			}
			else if(strcmp(interptype,"bilinear")==0){
				data_value=bilinearinterp<double>(x1,x2,y1,y2,Q11,Q12,Q21,Q22,x_grid,y_grid);
			}
			else if(strcmp(interptype,"nearest")==0){
				data_value=nearestinterp(x1,x2,y1,y2, Q11,Q12,Q21,Q22,x_grid,y_grid);
			}
			else{
				_printf_("Interpolation " << interptype << " not supported yet (supported intepolations are: triangle, bilinear and nearest)\n");
				return NULL; /*WARNING: no error because it would blow up the multithreading!*/
			}
			if(xIsNan<IssmPDouble>(data_value)) data_value=default_value;
		}
		else{
			data_value=default_value;
		}

		data_mesh->SetValue(i,data_value,INS_VAL);
	}

	return NULL;
}/*}}}*/
