/*!\file: balancevelocity_core.cpp
 * \brief: core of the balancevelocity solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void balancevelocity_core(FemModel* femmodel){

	/*parameters: */
	bool        save_results;
	IssmDouble  l = 8.;

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	if(VerboseSolution()) _printf0_("computing smooth driving stress:\n");
	femmodel->parameters->SetParam(l,SmoothThicknessMultiplierEnum);
	femmodel->SetCurrentConfiguration(SmoothAnalysisEnum);
	femmodel->parameters->SetParam(DrivingStressXEnum,InputToSmoothEnum);
	solutionsequence_linear(femmodel);
	femmodel->parameters->SetParam(DrivingStressYEnum,InputToSmoothEnum);
	solutionsequence_linear(femmodel);

	if(VerboseSolution()) _printf0_("computing balance velocities\n");
	femmodel->SetCurrentConfiguration(BalancevelocityAnalysisEnum);
	solutionsequence_linear(femmodel);

	if(save_results){
		int outputs[3] = {DrivingStressXEnum,DrivingStressYEnum,VelEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],3);
	}

}
