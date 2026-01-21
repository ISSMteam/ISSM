/*!\file: balancethickness_core.cpp
 * \brief: core of the balancethickness solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void balancethickness_core(FemModel* femmodel){

	/*recover parameters: */
	bool save_results;
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	/*Depth average velocities if necessary*/
	int domaintype;
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum){
		femmodel->parameters->SetParam(VxEnum,InputToDepthaverageInEnum);
		femmodel->parameters->SetParam(VxAverageEnum,InputToDepthaverageOutEnum);
		depthaverage_core(femmodel);
		if(domaintype==Domain3DEnum){
			femmodel->parameters->SetParam(VyEnum,InputToDepthaverageInEnum);
			femmodel->parameters->SetParam(VyAverageEnum,InputToDepthaverageOutEnum);
			depthaverage_core(femmodel);
		}
	}

	if(VerboseSolution()) _printf0_("computing balance thickness\n");
	femmodel->SetCurrentConfiguration(BalancethicknessAnalysisEnum);
	solutionsequence_linear(femmodel);

	if(save_results){
		int outputs = ThicknessEnum;
		femmodel->RequestedOutputsx(&femmodel->results,&outputs,1);
	}

}
