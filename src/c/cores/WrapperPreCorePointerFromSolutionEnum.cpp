/*!\file:  WrapperPreCorePointerFromSolutionEnum.cpp
 * \brief: return solution core that is carried out once only for Dakota runs.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void WrapperPreCorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype){

	/*output: */
	void (*solutioncore)(FemModel*)=NULL;

	switch(solutiontype){

		case TransientSolutionEnum:
			solutioncore=&transient_precore;
			break;
		default:
			break;
	}

	/*Assign output pointer:*/
	*psolutioncore=solutioncore;

}
