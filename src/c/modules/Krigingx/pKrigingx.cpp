/*!\file:  Kriging.cpp
 * \brief  "c" core code for Kriging
 */ 

#include "./Krigingx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/io/io.h"

int pKrigingx(double** ppredictions,double **perror,double* obs_x, double* obs_y, double* obs_list, int obs_length,double* x_interp,double* y_interp,int n_interp,Options* options){/*{{{*/

#ifdef _HAVE_MPI_
	int num_procs;
	int my_rank;

	/*output*/
	double *predictions = NULL;
	double *error       = NULL;

	/*Intermediaries*/
	int           mindata,maxdata;
	double        dmindata,dmaxdata;
	double        radius;
	char         *output       = NULL;
	Variogram    *variogram    = NULL;
	Observations *observations = NULL;

	/*timing*/
	double   start, finish;
	double   start_core, finish_core;
	double   start_init, finish_init;

	/*Get my_rank: */
	my_rank=IssmComm::GetRank();
	num_procs=IssmComm::GetSize();

	/*Get some Options*/
	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD); start=ISSM_MPI_Wtime();
	options->Get(&radius,"searchradius",0.);

	options->Get(&dmindata,"mindata",1.);  mindata=int(dmindata);//FIXME (Options come as double but we want to retrive integers)
	options->Get(&dmaxdata,"maxdata",50.); maxdata=int(dmaxdata);//FIXME (Options come as double but we want to retrive integers)

	/*Process observation dataset*/
	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD); start_init=ISSM_MPI_Wtime();
	observations=new Observations(obs_list,obs_x,obs_y,obs_length,options);
	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD); finish_init=ISSM_MPI_Wtime();

	/*Allocate output*/
	predictions =xNewZeroInit<double>(n_interp);
	error       =xNewZeroInit<double>(n_interp);

	/*Get output*/
	options->Get(&output,"output",(char*)"prediction");

	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD); start_core=ISSM_MPI_Wtime( );
	if(strcmp(output,"quadtree")==0){
		observations->QuadtreeColoring(predictions,x_interp,y_interp,n_interp);
	}
	else if(strcmp(output,"variomap")==0){
		observations->Variomap(predictions,x_interp,n_interp);
	}
	else if(strcmp(output,"prediction")==0){

		/*Process Variogram*/
		ProcessVariogram2(&variogram,options);

		/*partition loop across threads: */
		for(int idx=my_rank;idx<n_interp;idx+=num_procs){
			_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<double(idx)/double(n_interp)*100.<<"%  \n");
			observations->InterpolationKriging(&predictions[idx],&error[idx],x_interp[idx],y_interp[idx],radius,mindata,maxdata,variogram);
		}
		_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<100.<<"%  \n");

		double *sumpredictions =xNew<double>(n_interp);
		double *sumerror       =xNew<double>(n_interp);
		ISSM_MPI_Allreduce(predictions,sumpredictions,n_interp,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		ISSM_MPI_Allreduce(error,sumerror,n_interp,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		xDelete<double>(error); error=sumerror;
		xDelete<double>(predictions); predictions=sumpredictions;
	}
	else if(strcmp(output,"v4")==0){

		/*partition loop across threads: */
		for(int idx=my_rank;idx<n_interp;idx+=num_procs){
			_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<double(idx)/double(n_interp)*100.<<"%  \n");
			observations->InterpolationV4(&predictions[idx],x_interp[idx],y_interp[idx],radius,mindata,maxdata);
		}
		_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<100.<<"%  \n");

		double *sumpredictions =xNew<double>(n_interp);
		ISSM_MPI_Allreduce(predictions,sumpredictions,n_interp,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		xDelete<double>(predictions); predictions=sumpredictions;
	}
	else if(strcmp(output,"nearestneighbor")==0){

		/*partition loop across threads: */
		for(int idx=my_rank;idx<n_interp;idx+=num_procs){
			_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<double(idx)/double(n_interp)*100.<<"%  \n");
			observations->InterpolationNearestNeighbor(&predictions[idx],x_interp[idx],y_interp[idx],radius);
		}
		_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<100.<<"%  \n");

		double *sumpredictions =xNew<double>(n_interp);
		ISSM_MPI_Allreduce(predictions,sumpredictions,n_interp,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		xDelete<double>(predictions); predictions=sumpredictions;
	}
	else if(strcmp(output,"distance")==0){

		/*partition loop across threads: */
		for(int idx=my_rank;idx<n_interp;idx+=num_procs){
			_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<double(idx)/double(n_interp)*100.<<"%  \n");
			observations->Distances(&predictions[idx],&x_interp[idx],&y_interp[idx],1,radius);
		}
		_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<100.<<"%  \n");

		double *sumpredictions =xNew<double>(n_interp);
		ISSM_MPI_Allreduce(predictions,sumpredictions,n_interp,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		xDelete<double>(predictions); predictions=sumpredictions;
	}
	else if(strcmp(output,"idw")==0){
		double power;
		options->Get(&power,"power",2.);

		/*partition loop across threads: */
		for(int idx=my_rank;idx<n_interp;idx+=num_procs){
			_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<double(idx)/double(n_interp)*100.<<"%  \n");
			observations->InterpolationIDW(&predictions[idx],x_interp[idx],y_interp[idx],radius,mindata,maxdata,power);
		}
		_printf0_("      interpolation progress: "<<fixed<<setw(6)<<setprecision(4)<<100.<<"%  \n");

		double *sumpredictions =xNew<double>(n_interp);
		ISSM_MPI_Allreduce(predictions,sumpredictions,n_interp,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		xDelete<double>(predictions); predictions=sumpredictions;
	}
	else{
		_error_("output '" << output << "' not supported yet");
	}
	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD); finish_core=ISSM_MPI_Wtime( );

	/*clean-up and Assign output pointer*/
	delete variogram;
	delete observations;
	xDelete<char>(output);
	*ppredictions = predictions;
	*perror       = error;

	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD); finish=ISSM_MPI_Wtime( );
	_printf0_("\n   " << setw(34) << left << "Observation fitering elapsed time: " << finish_init-start_init << " seconds  \n\n");
	_printf0_("   " << setw(34) << left << "Kriging prediction elapsed time: " << finish_core-start_core << " seconds  \n\n");
	_printf0_("\n   " << "Total elapsed time " << int((finish-start)/3600) << " hrs " << int(int(finish-start)%3600/60) << " min " << int(finish-start)%60 << " sec\n\n\n");
	return 1;
#else
	_error_("MPI not available");
#endif
}/*}}}*/
void ProcessVariogram2(Variogram **pvariogram,Options* options){/*{{{*/

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
