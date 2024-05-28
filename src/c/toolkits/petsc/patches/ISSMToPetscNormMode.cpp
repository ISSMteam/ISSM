/* \file ISSMToPetscNormMode.cpp
 * \brief: convert NormMode from ISSM to Petsc
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscksp.h>

/*ISSM includes: */
#include "../../toolkitsenums.h"
#include "../../../shared/shared.h"

NormType ISSMToPetscNormMode(NormMode mode){

	switch(mode){
		case NORM_INF:  
			return NORM_INFINITY;
			break;
		case NORM_TWO:  
			return NORM_2;
			break;
		default: 
			_error_("unknown norm !");
			break;
	}
}
