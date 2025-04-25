/*!\file:  MatFree.cpp
 * \brief wrapper to MatDestroy
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscksp.h>

#include "./petscpatches.h"

void MatFree(PMat* pmat){

	#if PETSC_VERSION_LT(3,2,0)
	if(*pmat)MatDestroy(*pmat);
	*pmat=NULL;
	#else
	if(*pmat)MatDestroy(pmat);
	*pmat=NULL;
	#endif

}
