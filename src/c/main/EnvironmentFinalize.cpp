/*!\file:  EnvironmentFinalize.cpp
 * \brief: finalize Petsc, MPI, you name it
 */ 
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "../toolkits/toolkits.h"
#include "../shared/shared.h"

void EnvironmentFinalize(void){

	int my_rank;

	/*Make sure we are all here*/
	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD);

	/*Print closing statement*/
	ISSM_MPI_Comm_rank(ISSM_MPI_COMM_WORLD,&my_rank);

	/*Finalize: */
	//if(!my_rank) printf("closing MPI\n");
	ISSM_MPI_Finalize();
}
