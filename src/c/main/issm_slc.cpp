/*!\file:  issm_slc.cpp
 * \brief: ISSM SLR main program. 
 */ 

#include "./issm.h"
#include <stdlib.h>

int main(int argc,char **argv){

	/*diverse:*/
	int    nummodels;
	int*   commsizes=NULL;
	int*   rankzeros=NULL;
	char** dirnames=NULL;
	char** modelnames=NULL;
	int    modelid; 
	int    earthid; 
	int    my_rank;
	int    count=0;
	ISSM_MPI_Comm worldcomm;
	ISSM_MPI_Comm modelcomm;
	ISSM_MPI_Comm toearthcomm;
	ISSM_MPI_Comm* fromicecomms=NULL;

	/*Initialize exception trapping: */
	ExceptionTrapBegin();

	/*Initialize environment (MPI, PETSC, MUMPS, etc ...)*/
	worldcomm=EnvironmentInit(argc,argv);

	/*What is my rank?:*/
	ISSM_MPI_Comm_rank(worldcomm,&my_rank);

	/*How many models are we going to run (along with description and number of dedicated cores):{{{*/
	nummodels=(int) strtol(argv[4], (char **)NULL, 10);
	commsizes=xNew<int>(nummodels);
	dirnames=xNew<char*>(nummodels);
	modelnames=xNew<char*>(nummodels);
	rankzeros=xNew<int>(nummodels);
	for(int i=0;i<nummodels;i++){
		char* string=NULL;

		string=xNew<char>(strlen(argv[5+3*i])+1);
		xMemCpy<char>(string,argv[5+3*i],strlen(argv[5+3*i])+1);
		dirnames[i]=string;

		string=xNew<char>(strlen(argv[5+3*i+1])+1);
		xMemCpy<char>(string,argv[5+3*i+1],strlen(argv[5+3*i+1])+1);
		modelnames[i]=string;

		commsizes[i]=(int) strtol(argv[5+3*i+2], (char **)NULL, 10);
	}

	/*Figure out which model each cpu will belong to: */
	count=0;
	for(int i=0;i<nummodels;i++){
		if(my_rank>=count && my_rank<(count+commsizes[i])){
			modelid=i;
			break;
		}
		count+=commsizes[i];
	} 
	/*Buil array of who is rank 0 of their own group:*/
	count=0;
	for(int i=0;i<nummodels;i++){
		rankzeros[i]=count;
		count+=commsizes[i];
	}
	/*}}}*/

	/*Split world into sub-communicators for each and every model:*/
	ISSM_MPI_Comm_split(worldcomm,modelid, my_rank, &modelcomm);

	/*Build inter communicators:*/
	earthid=nummodels-1; //last model to be provided in the argument list if the earth model.
	if(modelid==earthid){
		fromicecomms=xNew<ISSM_MPI_Comm>(nummodels-1);
		for(int i=0;i<earthid;i++){
			ISSM_MPI_Intercomm_create( modelcomm, 0, worldcomm, rankzeros[i], i, fromicecomms+i); //communicate from local erth comm 9rank 0) to ice comm (rank 0) using modelid tag.
		}
	}
	else{
		ISSM_MPI_Intercomm_create( modelcomm, 0, worldcomm, rankzeros[earthid], modelid, &toearthcomm); //communicate from local ice comm (rank 0) to earth comm (rank 0) using modelid tag.
	}

	/*Supply specific argc and argv for each sub-communicator (corresponding to each  model specificatiions):{{{*/
	char** arguments=xNew<char*>(4);
	arguments[0]=xNew<char>(strlen(argv[0])+1); xMemCpy<char>(arguments[0],argv[0],strlen(argv[0])+1); //executable name
	arguments[1]=xNew<char>(strlen(argv[1])+1); xMemCpy<char>(arguments[1],argv[1],strlen(argv[1])+1); //solution name
	arguments[2]=xNew<char>(strlen(argv[5+3*modelid])+1); xMemCpy<char>(arguments[2],argv[5+3*modelid],strlen(argv[5+3*modelid])+1); //directory name
	arguments[3]=xNew<char>(strlen(argv[5+3*modelid+1])+1); xMemCpy<char>(arguments[3],argv[5+3*modelid+1],strlen(argv[5+3*modelid+1])+1); //model name
	/*}}}*/

	/*Initialize femmodel from arguments provided command line: */
	FemModel *femmodel = new FemModel(4,arguments,modelcomm);
	xDelete<char>(arguments[0]);
	xDelete<char>(arguments[1]);
	xDelete<char>(arguments[2]);
	xDelete<char>(arguments[3]);
	xDelete<char*>(arguments);

	/*Now that the models are initialized, keep communicator information in the parameters datasets of each model: */
	femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(worldcomm,WorldCommEnum));
	femmodel->parameters->AddObject(new IntParam(NumModelsEnum,nummodels));
	femmodel->parameters->AddObject(new IntParam(ModelIdEnum,modelid));
	femmodel->parameters->AddObject(new IntParam(EarthIdEnum,earthid));
	femmodel->parameters->AddObject(new IntParam(IsSlcCouplingEnum,1));
	if(modelid==earthid) femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm*>(fromicecomms,IcecapToEarthCommEnum));
	else femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(toearthcomm,IcecapToEarthCommEnum));

	/*Solve: */
	femmodel->Solve();

	/*Output results: */
	OutputResultsx(femmodel);

	/*Wrap up: */
	femmodel->CleanUp();
	delete femmodel;

	/*Delete communicators */
	ISSM_MPI_Comm_free(&modelcomm);
	if(modelid==earthid){
		for(int i=0;i<earthid;i++) ISSM_MPI_Comm_free(&fromicecomms[i]);
		xDelete<ISSM_MPI_Comm>(fromicecomms);
	}
	else{
		ISSM_MPI_Comm_free(&toearthcomm);
	}
	xDelete<int>(rankzeros);
	xDelete<int>(commsizes);

	/*Finalize environment:*/
	EnvironmentFinalize();

	/*Finalize exception trapping: */
	ExceptionTrapEnd();

	/*Free resources:*/
	for(int i=0;i<nummodels;i++){
		char* string=NULL;
		string=dirnames[i]; xDelete<char>(string);
		string=modelnames[i]; xDelete<char>(string);
	}
	xDelete<char*>(dirnames);
	xDelete<char*>(modelnames);

	/*Return unix success: */
	return 0; 
}
