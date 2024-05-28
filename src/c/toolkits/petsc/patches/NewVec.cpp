/*!\file NewVec.cpp
 * \brief: create distributed Petsc vector.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscksp.h>

#include "./petscpatches.h"
#include "../../mpi/issmmpi.h"

Vec NewVec(int size,ISSM_MPI_Comm comm,bool fromlocalsize){

	int local_size;

	/*output: */
	Vec vector=NULL;

	/*determine local size of vector: */
	if(fromlocalsize){
		local_size=size;
	}
	else{
		local_size=DetermineLocalSize(size,comm);
	}

	VecCreate(comm,&vector); 
	VecSetSizes(vector,local_size,PETSC_DECIDE); 
	VecSetFromOptions(vector); 
	VecSetOption(vector,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

	return vector;
}
