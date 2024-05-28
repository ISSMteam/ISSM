/*!\file:  VecFree.cpp
 * \brief wrapper to VecDestroy
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscksp.h>

void VecFree(Vec* pvec){

	#if PETSC_VERSION_LT(3,2,0)
	if(*pvec)VecDestroy(*pvec);
	#else
	if(*pvec)VecDestroy(pvec);
	#endif
	*pvec=NULL;

}
