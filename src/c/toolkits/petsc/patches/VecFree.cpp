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

#include "./petscpatches.h"

template<typename vectype>
void VecFree(vectype* pvec){

	#if PETSC_VERSION_LT(3,2,0)
	if(*pvec)VecDestroy(*pvec);
	#else
	if(*pvec)VecDestroy(pvec);
	#endif
	*pvec=NULL;

}

template void VecFree<PVec>(PVec*);
#if _HAVE_CODIPACK_
template void VecFree<Vec>(Vec*);
#endif
