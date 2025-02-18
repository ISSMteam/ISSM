/*!\file:  kriging.cpp
 * \brief: kriging main parallel program
 */

#include "./issm.h"

/*Local prototypes*/
void ProcessArguments2(char** pbinfilename,char** poutbinfilename,char** plockfilename,char** prootpath,int argc,char **argv);
void ProcessInputfile(double **px,double **py,double **pdata,int *pnobs,double **px_interp,double **py_interp,int *pninterp,Options **poptions,FILE* fid);

int main(int argc,char **argv){

	/*I/O: */
	FILE *output_fid = NULL;
	FILE *input_fid  = NULL;

	/*File names*/
	char *lockfilename   = NULL;
	char *binfilename    = NULL;
	char *outbinfilename = NULL;
	char *rootpath       = NULL;

	/*Input*/
	int         ninterp,nobs;
	double *x        = NULL;
	double *y        = NULL;
	double *data     = NULL;
	double *x_interp = NULL;
	double *y_interp = NULL;
	Options    *options  = NULL;

	/*Output*/
	double *predictions = NULL;
	double *error       = NULL;

	/*Initialize exception trapping: */
	ExceptionTrapBegin();

	/*Initialize environment (MPI, PETSC, MUMPS, etc ...)*/
	ISSM_MPI_Comm comm=EnvironmentInit(argc,argv);
	IssmComm::SetComm(comm);

	ProcessArguments2(&binfilename,&outbinfilename,&lockfilename,&rootpath,argc,argv);

	/*Process input files*/
	input_fid=pfopen(binfilename,"rb");
	ProcessInputfile(&x,&y,&data,&nobs,&x_interp,&y_interp,&ninterp,&options,input_fid);
	pfclose(input_fid,binfilename);

	_printf0_("call computational core:\n");
	pKrigingx(&predictions,&error,x,y,data,nobs,x_interp,y_interp,ninterp,options);

	_printf0_("write results to disk:\n");
	Results *results = new Results();
	if(IssmComm::GetRank()==0){
		output_fid=pfopen0(outbinfilename,"wb");
		results->AddObject(new GenericExternalResult<double*>(results->Size()+1,"predictions",predictions,ninterp,1,1,0));
		results->AddObject(new GenericExternalResult<double*>(results->Size()+1,"error",error,ninterp,1,1,0));
		for(Object* & object :results->objects){
			ExternalResult* result=xDynamicCast<ExternalResult*>(object);
			result->WriteData(output_fid,1);
		}
		pfclose(output_fid,outbinfilename);
	}

	/*Close output and toolkits options file and write lock file if requested*/
	_printf0_("write lock file:\n");
	WriteLockFile(lockfilename);

	/*Free resources: */
	xDelete<char>(lockfilename);
	xDelete<char>(binfilename);
	xDelete<char>(outbinfilename);
	xDelete<char>(rootpath);
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<double>(data);
	xDelete<double>(x_interp);
	xDelete<double>(y_interp);
	xDelete<double>(predictions);
	xDelete<double>(error);
	delete options;
	delete results;

	/*Finalize environment:*/
	EnvironmentFinalize();

	/*Finalize exception trapping: */
	ExceptionTrapEnd();

	return 0; //unix success return;
}

void ProcessArguments2(char** pbinfilename,char** poutbinfilename,char** plockfilename,char** prootpath,int argc,char **argv){

	if(argc<1)_error_("Usage error: no execution path provided");
	if(argc<2)_error_("Usage error: missing model name");

	/*Create rootpath from argument*/
	char* rootpatharg=argv[1];
	int   rootpath_len = strlen(rootpatharg)+2;
	char* rootpath = xNew<char>(rootpath_len);
	snprintf(rootpath, rootpath_len,"%s/",rootpatharg);

	/*Create all file paths*/
	char* modelname = argv[2];
	int base_length = strlen(rootpath)+strlen(modelname);
	int binfilename_len    = base_length+strlen(".bin")     +1;
	int outbinfilename_len = base_length+strlen(".outbin")  +1;
	int lockfilename_len   = base_length+strlen(".lock")    +1;
	char* binfilename    = xNew<char>(binfilename_len);      snprintf(binfilename,      binfilename_len,     "%s%s%s",rootpath,modelname,".bin");
	char* outbinfilename = xNew<char>(outbinfilename_len);   snprintf(outbinfilename,   outbinfilename_len,  "%s%s%s",rootpath,modelname,".outbin");
	char* lockfilename   = xNew<char>(lockfilename_len);     snprintf(lockfilename,     lockfilename_len,    "%s%s%s",rootpath,modelname,".lock");

	/*Clean up and assign output pointer*/
	*pbinfilename    = binfilename;
	*poutbinfilename = outbinfilename;
	*plockfilename   = lockfilename;
	*prootpath       = rootpath;
}

void ProcessInputfile(double **px,double **py,double **pdata,int *pnobs,double **px_interp,double **py_interp,int *pninterp,Options **poptions,FILE* fid){

	int     ninterp,nobs;
	double *x        = NULL;
	double *y        = NULL;
	double *data     = NULL;
	double *x_interp = NULL;
	double *y_interp = NULL;
	Options    *options  = NULL;

	int      M,N;
	IoModel* iomodel = new IoModel();
	iomodel->fid=fid;
	iomodel->CheckFile();
	iomodel->FetchData(&x,&M,&N,"md.x");               nobs=M*N;
	iomodel->FetchData(&y,&M,&N,"md.y");               _assert_(M*N==nobs);
	iomodel->FetchData(&data,&M,&N,"md.data");         _assert_(M*N==nobs);
	iomodel->FetchData(&x_interp,&M,&N,"md.x_interp"); ninterp=M*N;
	iomodel->FetchData(&y_interp,&M,&N,"md.y_interp"); _assert_(M*N==ninterp);

	/*Read options*/
	options = new Options();
	iomodel->FetchData(options,"md.y_interp");

	/*Assign output pointer*/
	*px        = x;
	*py        = y;
	*pdata     = data;
	*pnobs     = nobs;
	*px_interp = x_interp;
	*py_interp = y_interp;
	*pninterp  = ninterp;
	*poptions  = options;
	delete iomodel;
}
