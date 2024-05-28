/*!\file:  WrapperCorePointerFromSolutionEnum.cpp
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

void WrapperCorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype,bool nodakotacore){

	/*output: */
	void (*solutioncore)(FemModel*)=NULL;

	/*parameters: */
	bool control_analysis;
	bool dakota_analysis;
	int  inversiontype;

	/* retrieve some parameters that tell us whether wrappers are allowed, or whether we return 
	 * a pure core. Wrappers can be dakota_core (which samples many solution_cores) or control_core (which 
	 * carries out adjoint based inversion on a certain core: */
	parameters->FindParam(&dakota_analysis,QmuIsdakotaEnum);
	parameters->FindParam(&control_analysis,InversionIscontrolEnum);
	parameters->FindParam(&inversiontype,InversionTypeEnum);

	if(nodakotacore)dakota_analysis=false;

	if(dakota_analysis){
		#ifdef _HAVE_DAKOTA_
		solutioncore=dakota_core;
		#else
		_error_("ISSM was not compiled with dakota support, cannot carry out dakota analysis!");
		#endif
	}
	else if(control_analysis){
		switch(inversiontype){
			case 0: solutioncore=control_core; break;
			case 1: solutioncore=controltao_core; break;
			case 2: solutioncore=controlm1qn3_core; break;
			case 3: solutioncore=controlvalidation_core; break;
			case 4: solutioncore=controladm1qn3_core; break;
			default: _error_("control type not supported");
		}
	}
	else CorePointerFromSolutionEnum(&solutioncore,parameters,solutiontype);  /*This means we retrieve a core solution that is not a wrapper*/

	/*Assign output pointer:*/
	_assert_(psolutioncore);
	*psolutioncore=solutioncore;

}
