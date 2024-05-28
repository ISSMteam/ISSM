/*!\file:  KSPFree.cpp
 * \brief wrapper to KSPDestroy
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscksp.h>

void KSPFree(KSP* pksp){

	#if PETSC_VERSION_LT(3,2,0)
	if(*pksp)KSPDestroy(*pksp);
	*pksp=NULL;
	#else
	if(*pksp)KSPDestroy(pksp);
	*pksp=NULL;
	#endif

}
