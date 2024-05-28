/*!\file:  CorePointerFromSolutionEnum.cpp
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

void CorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype){

	/*output: */
	void (*solutioncore)(FemModel*)=NULL;

	switch(solutiontype){

		case StressbalanceSolutionEnum:
			solutioncore=&stressbalance_core;
			break;
		case SteadystateSolutionEnum:
			solutioncore=&steadystate_core;
			break;
		case ThermalSolutionEnum:
			solutioncore=&thermal_core;
			break;
		case BalancethicknessSolutionEnum:
			solutioncore=&balancethickness_core;
			break;
		case Balancethickness2SolutionEnum:
			solutioncore=&balancethickness2_core;
			break;
		case BalancethicknessSoftSolutionEnum:
			solutioncore=&dummy_core;
			break;
		case BalancevelocitySolutionEnum:
			solutioncore=&balancevelocity_core;
			break;
		case HydrologySolutionEnum:
			solutioncore=&hydrology_core;
			break;
		case SurfaceSlopeSolutionEnum:
			solutioncore=&surfaceslope_core;
			break;
		case BedSlopeSolutionEnum:
			solutioncore=&bedslope_core;
			break;
		case TransientSolutionEnum:
			solutioncore=&transient_core;
			break;
		case MasstransportSolutionEnum:
			solutioncore=&masstransport_core;
			break;
		case OceantransportSolutionEnum:
			solutioncore=&oceantransport_core;
			break;
		case EsaSolutionEnum:
			solutioncore=&esa_core;
			break;
		case DamageEvolutionSolutionEnum:
			solutioncore=&damage_core;
			break;
		case LoveSolutionEnum:
			#if _HAVE_LOVE_
			solutioncore=&love_core;
			#else
			_error_("ISSM not compiled with Love capability");
			#endif
			break;
		case SamplingSolutionEnum:
			solutioncore=&sampling_core;
			break;

		default:
			_error_("solution type: " << EnumToStringx(solutiontype) << " not supported yet!");
			break;
	}

	/*Assign output pointer:*/
	_assert_(psolutioncore);
	*psolutioncore=solutioncore;
}
