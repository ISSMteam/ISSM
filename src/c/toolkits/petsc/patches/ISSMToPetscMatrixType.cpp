/* \file ISSMToPetscMatrixType.cpp
 * \brief: convert MatrixType from ISSM to Petsc
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

MatType ISSMToPetscMatrixType(MatrixType type){

	switch(type){
		case DENSE_SEQUENTIAL:  
			return MATSEQDENSE;
			break;
		case SPARSE_SEQUENTIAL:  
			return MATSEQAIJ;
			break;
		default: 
			_error_("unknown matrix type !");
			break;
	}
}
