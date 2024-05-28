/*!\file:  issm.cpp
 * \brief: ISSM main program
 */ 

#include "../../c/main/issm.h"

extern "C" { 
	int IssmModule(char** poutput,int* poutputsize, double* buffer, int buffersize, char* toolkits,char* solution,char* modelname){

		/*output variables:*/
		char* output=NULL;
		size_t size;

		/*Initialize exception trapping: */
		ExceptionTrapBegin();

		/*Initialize environment: */
		ISSM_MPI_Comm comm_init=EnvironmentInit(0,NULL);

		/*Initialize femmodel from arguments provided command line: */
		FemModel *femmodel = new FemModel(buffer,buffersize,toolkits,solution,modelname,comm_init);

		/*Solve: */
		femmodel->Solve();

		/*Output results: */
		OutputResultsx(femmodel);

		/*Wrap up: */
		femmodel->CleanUpJs(&output,&size);

		/*Delete Model: */
		delete femmodel;
		
		/*Finalize environment:*/
		EnvironmentFinalize();

		/*Finalize exception trapping: */
		ExceptionTrapEnd();

		/*Assign output pointers:*/
		*poutputsize=(int)size;
		*poutput=output;

		/*Return output stream: */
		return 0 ;

	} 
} //extern "C" 
