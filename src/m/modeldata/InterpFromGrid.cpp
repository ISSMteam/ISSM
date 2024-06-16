/*Written by Mathieu Morlighem April 19th 2019*/

/*includes*/
#include <mex.h>
#include <pthread.h>
#include <math.h>   //for isnan
#include <cstring>  // for strcmp
#define f(m,n)\
  data[n*dataM+m] //Warning: matrix is transposed!

/*Inputs{{{*/
#define DATAX   (mxArray*)prhs[0]
#define DATAY   (mxArray*)prhs[1]
#define DATA    (mxArray*)prhs[2]
#define INTERPX (mxArray*)prhs[3]
#define INTERPY (mxArray*)prhs[4]
#define METHOD  (mxArray*)prhs[5]
/*}}}*/
/*Outputs{{{*/
#define INTERP (mxArray**)&plhs[0]
/*}}}*/
/*threading structs{{{*/
typedef struct{
	void* usr;
	int   my_thread;
	int   num_threads;
} pthread_handle;

typedef struct{
	int     dataM;
	int     dataN;
	double* datax;
	double* datay;
	double* data;
	int     interpN;
	double* interpx;
	double* interpy;
	double* interp;
	int     method;
} AppStruct; /*}}}*/
/*Prototypes{{{*/
void  FetchMatrixPointer(double** pmatrix,int *pM,int *pN,const mxArray* dataref);
void  FetchVectorPointer(double** pvector,int *pN,const mxArray* dataref);
void  FetchString(char** pstring,const mxArray* dataref);
void  WriteMatrix(mxArray** pdataref,double* matrix,int M,int N);
void  WriteVector(mxArray** pdataref,double* vector,int N);
void* InterpFromGridt(void* vpthread_handle);
void  LaunchThread(void* function(void*), void* usr,int num_threads);
bool  binary_search_increasing(int* pindex,double target,double* list,int n);
bool  binary_search_decreasing(int* pindex,double target,double* list,int n);
void  dataderivatives(double* A,double* x,double* y,double* data,int M,int N, int m0, int m1,int m2,int m3, int n0, int n1,int n2,int n3);
/*}}}*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){/*{{{*/

	double *datax   = NULL;
	double *datay   = NULL;
	double *data    = NULL;
	int     dataM,dataN;
	double *interpx = NULL;
	double *interpy = NULL;
	double *interp  = NULL;
	int     interpM,interpN;

	int     num_threads = 20;
	int     test1,test2,test3,test4;
	int     method = 1; // 0 = nearest, 1 = bilinear, 2 = bicubic

	/*Check arguments to avoid crash*/
	if(nlhs>1 || (nrhs<5 || nrhs>6)) mexErrMsgTxt("Wrong usage");

	/*Get variables from matlab to C*/
	FetchVectorPointer(&datax,&dataN,DATAX);
	FetchVectorPointer(&datay,&dataM,DATAY);
	FetchMatrixPointer(&data ,&test1,&test2,DATA); 
	FetchMatrixPointer(&interpx,&interpM,&interpN,INTERPX);
	FetchMatrixPointer(&interpy,&test3,&test4,INTERPY);
	if(!dataM*dataN)     mexErrMsgTxt("data is empty");
	if(!interpM*interpN) mexErrMsgTxt("no interpolation requested");
	if(test1!=dataM)     mexErrMsgTxt("x should have as many elements as there are columns in the data");
	if(test2!=dataN)     mexErrMsgTxt("y should have as many elements as there are lines in the data");
	if(test3!=interpM)   mexErrMsgTxt("interpolation locations (x,y) do not have the same size");
	if(test4!=interpN)   mexErrMsgTxt("interpolation locations (x,y) do not have the same size");
	if(nrhs==6){
		char* method_string = NULL;
		FetchString(&method_string,METHOD);
		if(strcmp(method_string,"nearest")==0)      method = 0;
		else if(strcmp(method_string,"linear")==0)  method = 1;
		else if(strcmp(method_string,"cubic")==0)   method = 2;
		else{
			mexErrMsgTxt("Method not supported yet");
		}

		mxFree(method_string);
	}

	/*Check inputs*/
	if(true){
		for(int i=0;i<interpM*interpN;i++){
			if(isnan(interpx[i])) mexErrMsgTxt("NaN found in interpx");
			if(isnan(interpy[i])) mexErrMsgTxt("NaN found in interpy");
		}
	}
	if(method==3){
		if(datax[1]-datax[0]<0) mexErrMsgTxt("x needs to be increasing for cubic interpolation");
		if(datay[1]-datay[0]<0) mexErrMsgTxt("y needs to be increasing for cubic interpolation");
	}

	/*Allocate output*/
	interp=(double*)mxMalloc(interpM*interpN*sizeof(double));

	/*Multithreaded core*/
	AppStruct usr;
	usr.dataM   = dataM;
	usr.dataN   = dataN;
	usr.datax   = datax;
	usr.datay   = datay;
	usr.data    = data;
	usr.interpN = interpM*interpN;
	usr.interpx = interpx;
	usr.interpy = interpy;
	usr.interp  = interp;
	usr.method  = method;
	LaunchThread(InterpFromGridt,(void*)&usr,num_threads);

	/*Write output vector*/
	WriteMatrix(INTERP,interp,interpM,interpN);

	/*Clean-up and return*/
	/*Do not erase pointers!*/
	return;
}/*}}}*/

/*InterpFromGridt{{{*/
void* InterpFromGridt(void* vpthread_handle){

	/*recover this thread info*/
	pthread_handle *handle = (pthread_handle*)vpthread_handle;
	int my_thread   = handle->my_thread;
	int num_threads = handle->num_threads;

	/*Recover struct*/
	AppStruct *usr = (AppStruct*)handle->usr;
	int     dataM   = usr->dataM;
	int     dataN   = usr->dataN;
	double *datax   = usr->datax;
	double *datay   = usr->datay;
	double *data    = usr->data;
	int     interpN = usr->interpN;
	double *interpx = usr->interpx;
	double *interpy = usr->interpy;
	double *interp  = usr->interp;
	int     method = usr->method;

	/*Intermediary*/
	double xprime,yprime;
	double x,y,x0,x1,x2,x3,y0,y1,y2,y3;
	double Q11,Q12;
	double Q21,Q22;
	double A[16];
	int    m,n,m0,m1,m2,m3,n0,n1,n2,n3;
	int    oldm=-1,oldn=-1;

	/*Is our matrix inverted?*/
	bool invertx = (datax[1]-datax[0])<0 ? true:false;
	bool inverty = (datay[1]-datay[0])<0 ? true:false;

	for(int idx=my_thread;idx<interpN;idx+=num_threads){

		x=interpx[idx];
		y=interpy[idx];

		/*Find indices m and n into y and x, for which  y(m)<=y_grids<=y(m+1) and x(n)<=x_grid<=x(n+1)*/
		if(invertx) binary_search_decreasing(&n,x,datax,dataN);
		else        binary_search_increasing(&n,x,datax,dataN);
		if(inverty) binary_search_decreasing(&m,y,datay,dataM);
		else        binary_search_increasing(&m,y,datay,dataM);

		if(n>=0 && n<dataN && m>=0 && m<dataM){

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
			if(invertx){
				n1=n+1; n2=n;
			}
			else{
				n1=n; n2=n+1;
			}
			if(inverty){
				m1=m+1; m2=m;
			}
			else{
				m1=m; m2=m+1;
			}

			x1 = datax[n1]; x2 = datax[n2];
			y1 = datay[m1]; y2 = datay[m2];

			if(method==0){
				/*Nearest neighbor interpolation*/
				if(x > (x1+x2)/2.){
					if(y > (y1+y2)/2.)
						interp[idx] = f(m2,n2);
					else
						interp[idx] = f(m1,n2);
					}
				else{
					if(y > (y1+y2)/2.)
						interp[idx] = f(m2,n1);
					else
						interp[idx] = f(m1,n1);
				}
				continue;
			}
			else if(method==1){
				/*Bilinear interpolation*/
				if(Q11==-9999 || Q12==-9999 || Q21==-9999 || Q22==-9999){
					interp[idx] = -9999;
					continue;
				}

				interp[idx] =
				  +f(m1,n1)*(x2-x)*(y2-y)/((x2-x1)*(y2-y1))
				  +f(m1,n2)*(x-x1)*(y2-y)/((x2-x1)*(y2-y1))
				  +f(m2,n1)*(x2-x)*(y-y1)/((x2-x1)*(y2-y1))
				  +f(m2,n2)*(x-x1)*(y-y1)/((x2-x1)*(y2-y1));
			}
			else{
				/*Bicubic interpolation*/
				if(invertx){n0=n+2; n3=n-1;}
				else{ n0=n-1; n3=n+2; }
				if(inverty){ m0=m+2; m3=m-1; }
				else{ m0=m-1; m3=m+2; }

				if(n0<0 || n3>=dataN || m0<0 || m3>=dataM){
					interp[idx] = -9999.;
					continue;
				}

				/*Local coordinates (between 0 and 1)*/
				xprime = (x - datax[n1])/(datax[n2]-datax[n1]);
				yprime = (y - datay[m1])/(datay[m2]-datay[m1]);

				/*Get derivatives at current pixel*/
				if(oldm!=m || oldn!=n){
					dataderivatives(&A[0],datax,datay,data,dataM,dataN,m0,m1,m2,m3,n0,n1,n2,n3);
					oldm = m;
					oldn = n;
				}

				double a00 = A[0];
				double a10 = A[4];
				double a20 = -3*A[0]+3*A[1]-2*A[4]-A[5];
				double a30 = 2*A[0]-2*A[1]+A[4]+A[5];
				double a01 = A[8];
				double a11 = A[12];
				double a21 = -3*A[8]+3*A[9]-2*A[12]-A[13];
				double a31 = 2*A[8]-2*A[9]+A[12]+A[13];
				double a02 = -3*A[0]+3*A[2]-2*A[8]-A[10];
				double a12 = -3*A[4]+3*A[6]-2*A[12]-A[14];
				double a22 = 9*A[0]-9*A[1]-9*A[2]+9*A[3]+6*A[4]+3*A[5]-6*A[6]-3*A[7]+6*A[8]-6*A[9]+3*A[10]-3*A[11]+4*A[12]+2*A[13]+2*A[14]+A[15];
				double a32 =-6*A[0]+6*A[1]+6*A[2]-6*A[3]-3*A[4]-3*A[5]+3*A[6]+3*A[7]-4*A[8]+4*A[9]-2*A[10]+2*A[11]-2*A[12]-2*A[13]-A[14]-A[15];
				double a03 = 2*A[0]-2*A[2]+A[8]+A[10];
				double a13 = 2*A[4]-2*A[6]+A[12]+A[14];
				double a23 =-6*A[0]+6*A[1]+6*A[2]-6*A[3]-4*A[4]-2*A[5]+4*A[6]+2*A[7]-3*A[8]+3*A[9]-3*A[10]+3*A[11]-2*A[12]-A[13]-2*A[14]-A[15] ;
				double a33 = 4*A[0]-4*A[1]-4*A[2]+4*A[3]+2*A[4]+2*A[5]-2*A[6]-2*A[7]+2*A[8]-2*A[9]+2*A[10]-2*A[11]+A[12]+A[13]+A[14]+A[15];

				x1= xprime;
				x2= x1*x1;
				x3= x2*x1;
				y1= yprime;
				y2= y1*y1;
				y3= y2*y1;
				interp[idx] = (a00+a01*y1+a02*y2+a03*y3)+(a10+a11*y1+a12*y2+a13*y3)*x1+(a20+a21*y1+a22*y2+a23*y3)*x2+(a30+a31*y1+a32*y2+a33*y3)*x3;
			}
		}
		else{
			interp[idx] = -9999.;
		}
	}
	//if(my_thread==0) printf("\r   interpolation progress = %5.1f%%\n",100.);

	return NULL;
}/*}}}*/
/*binary_search_increasing {{{*/
bool binary_search_increasing(int* pindex,double target,double* list,int n){

	/*output*/
	int  index;       //index, if found
	bool found=false; //found=0 if target is not found, 1 otherwise.

	/*intermediary*/
	int n0 = 0;
	int n1 = int(n/2);
	//int n1 = int((target-list[0])/(list[1]-list[0]));
	int n2 = n-1;

	if(target<list[n0]){
		found  = true;
		index  = -1;
	}
	else if(target>list[n2]){
		found  = true;
		index  = n;
	}
	else{
		while(!found){
			/*did we find the target?*/
			if(list[n1]<=target && list[n1+1]>=target){
				found = true;
				index = n1;
				break;
			}
			if(target < list[n1]){
				n2 = n1;
				n1 = n0 + int((n2-n0)/2);
			}
			else{
				n0 = n1;
				n1 = n0 + int((n2-n0)/2);
			}
		}
	}

	/*Assign output pointers:*/
	*pindex=index;

	/*Return result: */
	return found;
}/*}}}*/
/*binary_search_decreasing{{{*/
bool binary_search_decreasing(int* pindex,double target,double* list,int n){

	/*output*/
	int  index;       //index, if found
	bool found=false; //found=0 if target is not found, 1 otherwise.

	/*intermediary*/
	int n0 = 0;
	int n1 = int(n/2);
	//int n1 = int((target-list[0])/(list[0]-list[1]));
	int n2 = n-1;

	if (target>list[n0]){
		found  = true;
		index  = -1;
	}
	else if(target<list[n2]){
		found  = true;
		index  = n;
	}
	else{
		while(!found){
			/*did we find the target?*/
			if(list[n1]>=target && list[n1+1]<=target){
				found = true;
				index = n1;
				break;
			}
			if(target > list[n1]){
				n2 = n1;
				n1 = n0 + int((n2-n0)/2);
			}
			else{
				n0 = n1;
				n1 = n0 + int((n2-n0)/2);
			}
		}
	}

	/*Assign output pointers:*/
	*pindex=index;

	/*Return result: */
	return found;
}/*}}}*/
/*dataderivatives{{{*/
void  dataderivatives(double* A,double* x,double* y,double* data,int dataM,int dataN,
			int m0, int m1,int m2,int m3, int n0, int n1,int n2,int n3){

   /* i+1 +  +-------+ f(1,1)
    *     |  |       |
    *     |  |f(0,0) |
    *   i +  +-------+ f(1,0)
    *     +--+-------+-----> x
    *        j       j+1
	 */

   /*Function at corners*/
   A[0] = f(m1,n1); // f(0,0)
   A[1] = f(m1,n2); // f(1,0)
   A[2] = f(m2,n1); // f(0,1)
   A[3] = f(m2,n2); // f(1,1)

   /*x component of the gradient*/
   A[4] = .5*(f(m1,n2) - f(m1,n0));///(x[n2]-x[n0]); // dfdx(0,0)
   A[5] = .5*(f(m1,n3) - f(m1,n1));///(x[n3]-x[n1]); // dfdx(1,0)
   A[6] = .5*(f(m2,n2) - f(m2,n0));///(x[n2]-x[n0]); // dfdx(0,1)
   A[7] = .5*(f(m2,n3) - f(m2,n1));///(x[n3]-x[n1]); // dfdx(1,1)

   /*y component of the gradient*/
   A[ 8] = .5*(f(m2,n1) - f(m0,n1));///(y[m2]-y[m0]); // dfdy(0,0)
   A[ 9] = .5*(f(m2,n2) - f(m0,n2));///(y[m2]-y[m0]); // dfdy(1,0)
   A[10] = .5*(f(m3,n1) - f(m1,n1));///(y[m3]-y[m1]); // dfdy(0,1)
   A[11] = .5*(f(m3,n2) - f(m1,n2));///(y[m3]-y[m1]); // dfdy(1,1)

   /*cross-component of the gradient*/
   A[12] = .25*( (f(m2,n2) - f(m2,n0)) - (f(m0,n2) - f(m0,n0)) );///( (x[n2]-x[n0])*(y[m2]-y[m0]) ); // d2f/dxdy (0,0)
   A[13] = .25*( (f(m2,n3) - f(m2,n1)) - (f(m0,n3) - f(m0,n1)) );///( (x[n3]-x[n1])*(y[m2]-y[m0]) ); // d2f/dxdy (1,0)
   A[14] = .25*( (f(m3,n2) - f(m3,n0)) - (f(m1,n2) - f(m1,n0)) );///( (x[n2]-x[n0])*(y[m3]-y[m1]) ); // d2f/dxdy (0,1)
   A[15] = .25*( (f(m3,n3) - f(m3,n1)) - (f(m1,n3) - f(m1,n1)) );///( (x[n3]-x[n1])*(y[m3]-y[m1]) ); // d2f/dxdy (1,1)
}/*}}}*/
/*LaunchThread{{{*/
void LaunchThread(void* function(void*), void* usr,int num_threads){

	int i;
	int            *status  = NULL;
	pthread_t      *threads = NULL;
	pthread_handle *handles = NULL;

	/*dynamically allocate: */
	threads=(pthread_t*)mxMalloc(num_threads*sizeof(pthread_t));
	handles=(pthread_handle*)mxMalloc(num_threads*sizeof(pthread_handle));

	for(i=0;i<num_threads;i++){
		handles[i].usr=usr;
		handles[i].my_thread  =i;
		handles[i].num_threads=num_threads;
	}

	if(num_threads==1){
		function(handles);
	}
	else{
		for(i=0;i<num_threads;i++){
			if(pthread_create(threads+i,NULL,function,(void*)(handles+i))){
				mexErrMsgTxt("pthread_create error");
			}
		}
		for(i=0;i<num_threads;i++){
			if(pthread_join(threads[i],(void**)&status)){
				mexErrMsgTxt("pthread_join error");
			}
		}
	}

	/*Free resources:*/
	mxFree(threads);
	mxFree(handles);
}/*}}}*/
/*FetchMatrixPointer {{{*/
void FetchMatrixPointer(double** pmatrix,int *pM,int *pN,const mxArray* dataref){

	double *matrix=NULL;
	double *values=NULL;
	int     N,M;

	if(mxIsEmpty(dataref) ){
		M=N=0;
		matrix=NULL;
	}
	else if (mxIsDouble(dataref) ){
		M=mxGetM(dataref);
		N=mxGetN(dataref);
		matrix=(double*)mxGetPr(dataref);
	}
	else{
		mexErrMsgTxt("matrix type not supported");
	}

	*pmatrix=matrix;
	if (pN)*pN=N;
	if (pM)*pM=M;
}/*}}}*/
/*FetchVectorPointer {{{*/
void FetchVectorPointer(double** pvector,int *pN,const mxArray* dataref){

	double *vector=NULL;
	double *values=NULL;
	int     N;

	if(mxIsEmpty(dataref) ){
		N=0;
		vector=NULL;
	}
	else if (mxIsDouble(dataref) ){
		if(mxGetM(dataref)!=1 && mxGetN(dataref)!=1){
			mexErrMsgTxt("input is a matrix and not a vector");
		}
		N=mxGetN(dataref)*mxGetM(dataref);
		vector=(double*)mxGetPr(dataref);
	}
	else{
		mexErrMsgTxt("vector type not supported");
	}

	*pvector=vector;
	if (pN)*pN=N;
}/*}}}*/
/*FetchString{{{*/
void FetchString(char** pstring,const mxArray* dataref){

	char* outstring=NULL;

	/*Ok, the string should be coming directly from the matlab workspace: */
	if (!mxIsClass(dataref,"char")){
		mexErrMsgTxt("input data_type is not a string!");
	}
	else{
		/*Recover the string:*/
		int stringlen;

		stringlen = mxGetM(dataref)*mxGetN(dataref)+1;
		outstring = (char*)mxMalloc(stringlen*sizeof(char));
		mxGetString(dataref,outstring,stringlen);
	}

	/*Assign output pointers:*/
	*pstring=outstring;
	return;
}/*}}}*/
/*WriteMatrix {{{*/
void WriteMatrix(mxArray** pdataref,double* matrix,int M,int N){

	mxArray* dataref=NULL;

	if(matrix){
		/*data is a double* pointer. set pointer and invert sizes*/
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)M); 
		mxSetN(dataref,(mwSize)N);
		mxSetPr(dataref,(double*)matrix);
	}
	else{
		dataref = mxCreateDoubleScalar(0.0);
	}
	*pdataref=dataref;
}
/*}}}*/
/*WriteVector {{{*/
void WriteVector(mxArray** pdataref,double* vector,int N){

	mxArray* dataref=NULL;

	if(vector){
		/*data is a double* pointer. Copy into a vector: */
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)N);
		mxSetN(dataref,(mwSize)1);
		mxSetPr(dataref,(double*)vector);
	}
	else{
		dataref = mxCreateDoubleScalar(0.0);
	}
	*pdataref=dataref;
}
/*}}}*/
