/*!\file:  EnvironmentInit.cpp
 * \brief: initialize Petsc, MPI, you name it
 */ 
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <stdio.h>
#include "../toolkits/toolkits.h"

ISSM_MPI_Comm EnvironmentInit(int argc,char** argv){

	/*Output*/
	ISSM_MPI_Comm comm = 0;

	/*Initialize MPI environment: */
	#if defined(_HAVE_MPI_)
	ISSM_MPI_Init(&argc,&argv);
	comm = ISSM_MPI_COMM_WORLD;
	#else
	comm = 1; //bogus number for comm, which does not exist anyway.
	#endif

	/*Print Banner*/
	int my_rank = 0;
	ISSM_MPI_Comm_rank(comm,&my_rank);
	if(!my_rank) printf("\n");
	if(!my_rank) printf("──────────────────────────────────────────────────────────────────────\n");
	if(!my_rank) printf("%s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
	if(!my_rank) printf("          GitHub: %s\n", PACKAGE_BUGREPORT);
	if(!my_rank) printf("   Documentation: %s\n", PACKAGE_URL);
	if(!my_rank) printf("──────────────────────────────────────────────────────────────────────\n");

	/*Return communicator*/
	return comm;
}
