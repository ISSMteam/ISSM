/*!\file:  issm_ocean.cpp
 * \brief: ISSM OCEAN main program. 
 */ 

#include "./issm.h"
#include <stdlib.h>

int main(int argc,char **argv){

	/*diverse:*/
	int    icecommsize;
	int    my_rank,my_local_rank,my_size,my_local_size;
	ISSM_MPI_Comm worldcomm;
	ISSM_MPI_Comm modelcomm;
	ISSM_MPI_Comm frommitgcm;
	ISSM_MPI_Comm tomitgcmcomm;
	ISSM_MPI_Status status;

	/*Initialize exception trapping: */
	ExceptionTrapBegin();

	/*Initialize environment (MPI, PETSC, MUMPS, etc ...)*/
	worldcomm=EnvironmentInit(argc,argv);

	/*What is my rank?:*/
	ISSM_MPI_Comm_rank(worldcomm,&my_rank);
	ISSM_MPI_Comm_size(worldcomm,&my_size);

	/*First model is ice, second is ocean*/
	/*ice comm size: */
	icecommsize=(int) strtol(argv[2], (char **)NULL, 10);

	/*Split world into sub-communicators for each and every model:*/
	ISSM_MPI_Comm_split(worldcomm,0, my_rank, &modelcomm);
	ISSM_MPI_Comm_rank(modelcomm,&my_local_rank);
	ISSM_MPI_Comm_size(modelcomm,&my_local_size);

	ISSM_MPI_Intercomm_create( modelcomm, 0, worldcomm, my_local_size, 0, &tomitgcmcomm); 

	FemModel *femmodel = new FemModel(argc,argv,modelcomm);

	/*Now that the models are initialized, keep communicator information in the parameters datasets of each model: */
	femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(worldcomm,WorldCommEnum));
	femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(tomitgcmcomm,ToMITgcmCommEnum));

	/*Solve: */
	femmodel->Solve();

	/*Output results: */
	OutputResultsx(femmodel);

	/*Wrap up: */
	femmodel->CleanUp();

	/*Delete Model: */
	delete femmodel;

	/*Finalize environment:*/
	EnvironmentFinalize();

	/*Finalize exception trapping: */
	ExceptionTrapEnd();

	/*Free resources:*/

	/*Return unix success: */
	return 0; 
}
