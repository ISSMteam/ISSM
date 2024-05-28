/* \file ISSMToPetscInsertMode.cpp
 * \brief: convert InsertMode from ISSM to Petsc
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

InsertMode ISSMToPetscInsertMode(InsMode mode){

	switch(mode){
		case ADD_VAL:  
			return ADD_VALUES;
			break;
		case INS_VAL:
			return INSERT_VALUES;
			break;
		default: 
			_error_("unknown insert mode!");
			break;
	}
}
