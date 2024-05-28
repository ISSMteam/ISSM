/*!\file: balancethickness_core.cpp
 * \brief: core of the balancethickness solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void balancethickness2_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;
	//IssmDouble  l = 3.;

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	//if(VerboseSolution()) _printf0_("computing smooth surface slopes:\n");
	//femmodel->parameters->SetParam(l,SmoothThicknessMultiplierEnum);
	//femmodel->SetCurrentConfiguration(SmoothAnalysisEnum);
	//femmodel->parameters->SetParam(SurfaceSlopeXEnum,InputToSmoothEnum);
	//solutionsequence_linear(femmodel);
	//femmodel->parameters->SetParam(SurfaceSlopeYEnum,InputToSmoothEnum);
	//solutionsequence_linear(femmodel);
	//surfaceslope_core(femmodel);

	femmodel->SetCurrentConfiguration(Balancethickness2AnalysisEnum);
	solutionsequence_linear(femmodel);
	//solutionsequence_nonlinear(femmodel,false);

	if(save_results){
		const int numoutputs = 1;
		int outputs[numoutputs] = {ThicknessEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],numoutputs);
	}

}
