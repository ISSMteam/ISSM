/*!\file:  issm_slc.cpp
 * \brief: ISSM SLR main program. 
 */ 

#include "./issm.h"
#include <stdlib.h>

int main(int argc,char **argv){

	/*diverse:*/
	int*   commsizes=NULL;
	int*   rankzeros=NULL;
	int    my_rank;
	ISSM_MPI_Comm  worldcomm,modelcomm,toearthcomm;
	ISSM_MPI_Comm* fromicecomms=NULL;

	/*Initialize exception trapping: */
	ExceptionTrapBegin();

	/*Initialize environment (MPI, PETSC, MUMPS, etc ...)*/
	worldcomm=EnvironmentInit(argc,argv);

	/*What is my rank?:*/
	ISSM_MPI_Comm_rank(worldcomm,&my_rank);

	/*How many models are we going to run (along with description and number of dedicated cores)*/
	int nummodels=(int) strtol(argv[4], (char **)NULL, 10);
	commsizes  = xNew<int>(nummodels);
	rankzeros  = xNew<int>(nummodels);
	for(int i=0;i<nummodels;i++){
		commsizes[i]=(int) strtol(argv[5+3*i+2], (char **)NULL, 10);
	}

	/*Figure out which model each cpu will belong to*/
	int count   = 0;
	int modelid = -1;
	for(int i=0;i<nummodels;i++){
		if(my_rank>=count && my_rank<(count+commsizes[i])){
			modelid=i;
			break;
		}
		count+=commsizes[i];
	} 
	_assert_(modelid>=0);

	/*Buil array of who is rank 0 of their own group:*/
	count = 0;
	for(int i=0;i<nummodels;i++){
		rankzeros[i]=count;
		count+=commsizes[i];
	}

	/*Split world into sub-communicators for each and every model:*/
	ISSM_MPI_Comm_split(worldcomm, modelid, my_rank, &modelcomm);

	/*Build inter communicators:*/
	int earthid=nummodels-1; //last model to be provided in the argument list if the earth model.
	if(modelid==earthid){
		fromicecomms=xNew<ISSM_MPI_Comm>(nummodels-1);
		for(int i=0;i<earthid;i++){
         /*communicate from local erth comm 9rank 0) to ice comm (rank 0) using modelid tag.*/
			ISSM_MPI_Intercomm_create( modelcomm, 0, worldcomm, rankzeros[i], i, fromicecomms+i); 
		}
	}
	else{
      /*communicate from local ice comm (rank 0) to earth comm (rank 0) using modelid tag.*/
		ISSM_MPI_Intercomm_create( modelcomm, 0, worldcomm, rankzeros[earthid], modelid, &toearthcomm);
	}

	/*Supply specific argc and argv for each sub-communicator (corresponding to each  model specifications)*/
	char** arguments=xNew<char*>(4);
	arguments[0]=xNew<char>(strlen(argv[0])+1); xMemCpy<char>(arguments[0],argv[0],strlen(argv[0])+1); //executable name
	arguments[1]=xNew<char>(strlen(argv[1])+1); xMemCpy<char>(arguments[1],argv[1],strlen(argv[1])+1); //solution name
	arguments[2]=xNew<char>(strlen(argv[5+3*modelid])+1); xMemCpy<char>(arguments[2],argv[5+3*modelid],strlen(argv[5+3*modelid])+1); //directory name
	arguments[3]=xNew<char>(strlen(argv[5+3*modelid+1])+1); xMemCpy<char>(arguments[3],argv[5+3*modelid+1],strlen(argv[5+3*modelid+1])+1); //model name

	/*Initialize femmodel from arguments provided command line: */
	FemModel *femmodel = new FemModel(4,arguments,modelcomm);
	for(int i=0;i<4;i++) xDelete<char>(arguments[i]);
	xDelete<char*>(arguments);

	/*Now that the models are initialized, keep communicator information in the parameters datasets of each model: */
	femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(worldcomm,WorldCommEnum));
	femmodel->parameters->AddObject(new IntParam(NumModelsEnum,nummodels));
	femmodel->parameters->AddObject(new IntParam(ModelIdEnum,modelid));
	femmodel->parameters->AddObject(new IntParam(EarthIdEnum,earthid));
	femmodel->parameters->AddObject(new IntParam(IsSlcCouplingEnum,1));
	if(modelid==earthid){
      femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm*>(fromicecomms,IcecapToEarthCommEnum));
   }
	else{
      femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(toearthcomm,IcecapToEarthCommEnum));
   }

	/*Solve: */
	femmodel->Solve();

	/*Output results: */
	OutputResultsx(femmodel);

   /*Free resources:*/
	femmodel->CleanUp();
	delete femmodel;
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

	/*Return unix success: */
	return 0; 
}
