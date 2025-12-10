/*!\file:  Kriging.cpp
 * \brief  "c" core code for Kriging
 */ 

#include "./Krigingx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../modules.h"
#ifdef _HAVE_GSL_
#include <gsl/gsl_linalg.h>
#endif
int   Krigingx(double** ppredictions,double **perror,double* obs_x, double* obs_y, double* obs_list, int obs_length,double* x_interp,double* y_interp,int n_interp,Options* options){/*{{{*/

	/*output*/
	double *predictions = NULL;
	double *error       = NULL;

	/*Intermediaries*/
	int           mindata,maxdata;
	double        dmindata,dmaxdata,dnumthreads; //FIXME (Options come as double but we want to retrive integers)
	double        radius;
	char         *output       = NULL;
	Variogram    *variogram    = NULL;
	Observations *observations = NULL;

	/*threading: */
	KrigingxThreadStruct gate;
	int num = _NUMTHREADS_;

	/*Get Variogram from Options*/
	ProcessVariogram(&variogram,options);
	options->Get(&radius,"searchradius",0.);
	options->Get(&dmindata,"mindata",1.);  mindata=int(dmindata);//FIXME (Options come as double but we want to retrive integers)
	options->Get(&dmaxdata,"maxdata",50.); maxdata=int(dmaxdata);//FIXME (Options come as double but we want to retrive integers)
	options->Get(&dnumthreads,"numthreads",double(num)); num=int(dnumthreads);//FIXME (Options come as double but we want to retrive integers)

	/*Process observation dataset*/
	observations=new Observations(obs_list,obs_x,obs_y,obs_length,options);

	/*Allocate output*/
	predictions =xNewZeroInit<double>(n_interp);
	error       =xNewZeroInit<double>(n_interp);

	/*Get output*/
	options->Get(&output,"output",(char*)"prediction");

	if(strcmp(output,"quadtree")==0){
		observations->QuadtreeColoring(predictions,x_interp,y_interp,n_interp);
	}
	else if(strcmp(output,"variomap")==0){
		observations->Variomap(predictions,x_interp,n_interp);
	}
	else if(strcmp(output,"distance")==0){
		/*initialize thread parameters: */
		gate.n_interp     = n_interp;
		gate.x_interp     = x_interp;
		gate.y_interp     = y_interp;
		gate.radius       = radius;
		gate.mindata      = mindata;
		gate.maxdata      = maxdata;
		gate.variogram    = variogram;
		gate.observations = observations;
		gate.predictions  = predictions;
		gate.error        = error;
		gate.numdone      = xNewZeroInit<int>(num);

		/*launch the thread manager with Krigingxt as a core: */
		LaunchThread(Distancest,(void*)&gate,num);
		xDelete<int>(gate.numdone);
	}
	else if(strcmp(output,"delaunay")==0){

		#ifdef _HAVE_BAMG_
		int nobs,nel;
		double *x     = NULL;
		double *y     = NULL;
		double *data  = NULL;
		int    *index = NULL;

		observations->ObservationList(&x,&y,&data,&nobs);

		_printf_("Generation Delaunay Triangulation\n");
		BamgTriangulatex(&index,&nel,x,y,nobs);

		_printf_("Interpolating\n");
		xDelete<double>(predictions);
		InterpFromMeshToMesh2dx(&predictions,index,x,y,nobs,nel,data,nobs,1,x_interp,y_interp,n_interp,options);
		xDelete<double>(x);
		xDelete<double>(y);
		xDelete<double>(data);
		xDelete<int>(index);
		#else
		_error_("you did not compile ISSM with bamg");
		#endif
	}
	else if(strcmp(output,"nearestneighbor")==0){
		/*initialize thread parameters: */
		gate.n_interp     = n_interp;
		gate.x_interp     = x_interp;
		gate.y_interp     = y_interp;
		gate.radius       = radius;
		gate.mindata      = mindata;
		gate.maxdata      = maxdata;
		gate.variogram    = variogram;
		gate.observations = observations;
		gate.predictions  = predictions;
		gate.error        = error;
		gate.numdone      = xNewZeroInit<int>(num);

		/*launch the thread manager with Krigingxt as a core: */
		LaunchThread(NearestNeighbort,(void*)&gate,num);
		printf("\r      interpolation progress: 100%%   \n");
		xDelete<int>(gate.numdone);
	}
	else if(strcmp(output,"idw")==0){ //Inverse distance weighting
		double power;
		options->Get(&power,"power",2.);
		/*initialize thread parameters: */
		gate.n_interp     = n_interp;
		gate.x_interp     = x_interp;
		gate.y_interp     = y_interp;
		gate.radius       = radius;
		gate.mindata      = mindata;
		gate.maxdata      = maxdata;
		gate.variogram    = variogram;
		gate.observations = observations;
		gate.predictions  = predictions;
		gate.error        = error;
		gate.numdone      = xNewZeroInit<int>(num);
		gate.power        = power;

		/*launch the thread manager with Krigingxt as a core: */
		LaunchThread(idwt,(void*)&gate,num);
		printf("\r      interpolation progress: 100%%   \n");
		xDelete<int>(gate.numdone);
	}
	else if(strcmp(output,"v4")==0){ //Inverse distance weighting
#if !defined(_HAVE_GSL_)
		_error_("GSL is required for v4 interpolation");
#endif
		/*initialize thread parameters: */
		gate.n_interp     = n_interp;
		gate.x_interp     = x_interp;
		gate.y_interp     = y_interp;
		gate.radius       = radius;
		gate.mindata      = mindata;
		gate.maxdata      = maxdata;
		gate.variogram    = variogram;
		gate.observations = observations;
		gate.predictions  = predictions;
		gate.error        = error;
		gate.numdone      = xNewZeroInit<int>(num);

		/*launch the thread manager with Krigingxt as a core: */
		LaunchThread(v4t,(void*)&gate,num);
		printf("\r      interpolation progress: 100%%   \n");
		xDelete<int>(gate.numdone);
	}
	else if(strcmp(output,"prediction")==0){
#if !defined(_HAVE_GSL_)
		_error_("GSL is required for v4 interpolation");
#endif

		/*initialize thread parameters: */
		gate.n_interp     = n_interp;
		gate.x_interp     = x_interp;
		gate.y_interp     = y_interp;
		gate.radius       = radius;
		gate.mindata      = mindata;
		gate.maxdata      = maxdata;
		gate.variogram    = variogram;
		gate.observations = observations;
		gate.predictions  = predictions;
		gate.error        = error;
		gate.numdone      = xNewZeroInit<int>(num);

		/*launch the thread manager with Krigingxt as a core: */
		LaunchThread(Krigingxt,(void*)&gate,num);
		printf("\r      interpolation progress: 100%%   \n");
		xDelete<int>(gate.numdone);
	}
	else{
		_error_("output '" << output << "' not supported yet");
	}

	/*clean-up and Assign output pointer*/
	delete variogram;
	delete observations;
	xDelete<char>(output);
	*ppredictions = predictions;
	*perror       = error;
	return 1;
}/*}}}*/
void* Krigingxt(void* vpthread_handle){/*{{{*/

	/*gate variables :*/
	KrigingxThreadStruct *gate        = NULL;
	pthread_handle       *handle      = NULL;
	int my_thread;
	int num_threads;
	int i0,i1;

	/*recover handle and gate: */
	handle      = (pthread_handle*)vpthread_handle;
	gate        = (KrigingxThreadStruct*)handle->gate;
	my_thread   = handle->id;
	num_threads = handle->num;

	/*recover parameters :*/
	int           n_interp     = gate->n_interp;
	double       *x_interp     = gate->x_interp;
	double       *y_interp     = gate->y_interp;
	double        radius       = gate->radius;
	int           mindata      = gate->mindata;
	int           maxdata      = gate->maxdata;
	Variogram    *variogram    = gate->variogram;
	Observations *observations = gate->observations;
	double       *predictions  = gate->predictions;
	double       *error        = gate->error;
	int          *numdone      = gate->numdone;

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,n_interp,num_threads,my_thread);
	for(int idx=i0;idx<i1;idx++){

		/*Print info*/
		numdone[my_thread]=idx-i0;
		if(my_thread==0){
			int alldone=numdone[0];
			for(int i=1;i<num_threads;i++) alldone+=numdone[i];
			printf("\r      interpolation progress: %6.2g%%   ",double(alldone)/double(n_interp)*100.);
		}

		/*Kriging interpolation*/
		observations->InterpolationKriging(&predictions[idx],&error[idx],x_interp[idx],y_interp[idx],radius,mindata,maxdata,variogram);
	}

	return NULL;
}/*}}}*/
void* NearestNeighbort(void* vpthread_handle){/*{{{*/

	/*gate variables :*/
	KrigingxThreadStruct *gate        = NULL;
	pthread_handle       *handle      = NULL;
	int my_thread;
	int num_threads;
	int i0,i1;

	/*recover handle and gate: */
	handle      = (pthread_handle*)vpthread_handle;
	gate        = (KrigingxThreadStruct*)handle->gate;
	my_thread   = handle->id;
	num_threads = handle->num;

	/*recover parameters :*/
	int           n_interp     = gate->n_interp;
	double       *x_interp     = gate->x_interp;
	double       *y_interp     = gate->y_interp;
	double        radius       = gate->radius;
	int           mindata      = gate->mindata;
	int           maxdata      = gate->maxdata;
	Variogram    *variogram    = gate->variogram;
	Observations *observations = gate->observations;
	double       *predictions  = gate->predictions;
	double       *error        = gate->error;
	int          *numdone      = gate->numdone;

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,n_interp,num_threads,my_thread);
	for(int idx=i0;idx<i1;idx++){

		/*Print info*/
		numdone[my_thread]=idx-i0;
		if(my_thread==0){
			int alldone=numdone[0];
			for(int i=1;i<num_threads;i++) alldone+=numdone[i];
			printf("\r      interpolation progress: %6.2g%%   ",double(alldone)/double(n_interp)*100.);
		}

		observations->InterpolationNearestNeighbor(&predictions[idx],x_interp[idx],y_interp[idx],radius);
	}

	return NULL;
}/*}}}*/
void* idwt(void* vpthread_handle){/*{{{*/

	/*gate variables :*/
	KrigingxThreadStruct *gate        = NULL;
	pthread_handle       *handle      = NULL;
	int my_thread;
	int num_threads;
	int i0,i1;

	/*recover handle and gate: */
	handle      = (pthread_handle*)vpthread_handle;
	gate        = (KrigingxThreadStruct*)handle->gate;
	my_thread   = handle->id;
	num_threads = handle->num;

	/*recover parameters :*/
	int           n_interp     = gate->n_interp;
	double       *x_interp     = gate->x_interp;
	double       *y_interp     = gate->y_interp;
	double        radius       = gate->radius;
	int           mindata      = gate->mindata;
	int           maxdata      = gate->maxdata;
	Variogram    *variogram    = gate->variogram;
	Observations *observations = gate->observations;
	double       *predictions  = gate->predictions;
	double       *error        = gate->error;
	int          *numdone      = gate->numdone;
	double        power        = gate->power;

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,n_interp,num_threads,my_thread);
	for(int idx=i0;idx<i1;idx++){

		/*Print info*/
		numdone[my_thread]=idx-i0;
		if(my_thread==0){
			int alldone=numdone[0];
			for(int i=1;i<num_threads;i++) alldone+=numdone[i];
			printf("\r      interpolation progress: %6.2g%%   ",double(alldone)/double(n_interp)*100.);
		}

		observations->InterpolationIDW(&predictions[idx],x_interp[idx],y_interp[idx],radius,mindata,maxdata,power);
	}
	return NULL;
}/*}}}*/
void* v4t(void* vpthread_handle){/*{{{*/

	/*gate variables :*/
	KrigingxThreadStruct *gate        = NULL;
	pthread_handle       *handle      = NULL;
	int my_thread;
	int num_threads;
	int i0,i1;

	/*recover handle and gate: */
	handle      = (pthread_handle*)vpthread_handle;
	gate        = (KrigingxThreadStruct*)handle->gate;
	my_thread   = handle->id;
	num_threads = handle->num;

	/*recover parameters :*/
	int           n_interp     = gate->n_interp;
	double       *x_interp     = gate->x_interp;
	double       *y_interp     = gate->y_interp;
	double        radius       = gate->radius;
	int           mindata      = gate->mindata;
	int           maxdata      = gate->maxdata;
	Variogram    *variogram    = gate->variogram;
	Observations *observations = gate->observations;
	double       *predictions  = gate->predictions;
	double       *error        = gate->error;
	int          *numdone      = gate->numdone;

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,n_interp,num_threads,my_thread);
	for(int idx=i0;idx<i1;idx++){

		/*Print info*/
		numdone[my_thread]=idx-i0;
		if(my_thread==0){
			int alldone=numdone[0];
			for(int i=1;i<num_threads;i++) alldone+=numdone[i];
			printf("\r      interpolation progress: %6.2g%%   ",double(alldone)/double(n_interp)*100.);
		}

		observations->InterpolationV4(&predictions[idx],x_interp[idx],y_interp[idx],radius,mindata,maxdata);
	}
	return NULL;
}/*}}}*/
void* Distancest(void* vpthread_handle){/*{{{*/

	/*gate variables :*/
	KrigingxThreadStruct *gate        = NULL;
	pthread_handle       *handle      = NULL;
	int my_thread;
	int num_threads;
	int i0,i1;

	/*recover handle and gate: */
	handle      = (pthread_handle*)vpthread_handle;
	gate        = (KrigingxThreadStruct*)handle->gate;
	my_thread   = handle->id;
	num_threads = handle->num;

	/*recover parameters :*/
	int           n_interp     = gate->n_interp;
	double       *x_interp     = gate->x_interp;
	double       *y_interp     = gate->y_interp;
	double        radius       = gate->radius;
	int           mindata      = gate->mindata;
	int           maxdata      = gate->maxdata;
	Variogram    *variogram    = gate->variogram;
	Observations *observations = gate->observations;
	double       *predictions  = gate->predictions;
	double       *error        = gate->error;
	int          *numdone      = gate->numdone;

	/*partition loop across threads: */
	PartitionRange(&i0,&i1,n_interp,num_threads,my_thread);
	observations->Distances(&predictions[i0],&x_interp[i0],&y_interp[i0],i1-i0,radius);
	return NULL;
}/*}}}*/

void ProcessVariogram(Variogram **pvariogram,Options* options){/*{{{*/

	/*Intermediaries*/
	Variogram* variogram = NULL;
	char      *model     = NULL;

	if(options->GetOption("model")){
		options->Get(&model,"model");
		if     (strcmp(model,"gaussian")==0)    variogram = new GaussianVariogram(options);
		else if(strcmp(model,"exponential")==0) variogram = new ExponentialVariogram(options);
		else if(strcmp(model,"spherical")==0)   variogram = new SphericalVariogram(options);
		else if(strcmp(model,"power")==0)       variogram = new PowerVariogram(options);
		else _error_("variogram " << model << " not supported yet (list of supported variogram: gaussian, exponential, spherical and power)");
	}
	else variogram = new GaussianVariogram(options);

	/*Assign output pointer*/
	xDelete<char>(model);
	*pvariogram = variogram;
}/*}}}*/
