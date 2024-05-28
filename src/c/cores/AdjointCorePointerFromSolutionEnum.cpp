/*!\file:  AdjointCorePointerFromSolutionEnum.cpp
 * \brief: return type of analyses, number of analyses and core solution function.
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

void AdjointCorePointerFromSolutionEnum(void (**padjointcore)(FemModel*),int solutiontype){

	/*output: */
	void (*adjointcore)(FemModel*)=NULL;

	switch(solutiontype){

		case StressbalanceSolutionEnum:
			adjointcore=&adjointstressbalance_core;
			break;
		case SteadystateSolutionEnum:
			adjointcore=&adjointstressbalance_core;
			break;
		case BalancethicknessSolutionEnum:
			adjointcore=&adjointbalancethickness_core;
			break;
		case Balancethickness2SolutionEnum:
			adjointcore=&adjointbalancethickness2_core;
			break;
		case BalancethicknessSoftSolutionEnum:
			adjointcore=&dummy_core;
			break;
		default:
			_error_("No adjoint has been implemented for solution " << EnumToStringx(solutiontype) << " yet");
			break;
	}

	/*Assign output pointer:*/
	_assert_(padjointcore);
	*padjointcore=adjointcore;

}
