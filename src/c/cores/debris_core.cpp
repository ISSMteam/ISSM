/*!\file: debris_core.cpp
 * \brief: core of the debris solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
#include "../classes/Inputs/TransientInput.h"

void debris_core(FemModel* femmodel){ /*{{{*/

	/*Start profiler*/
	femmodel->profiler->Start(DEBRISCORE);

	/*parameters: */
	int    numoutputs,domaintype;
	bool   save_results;
	int    solution_type,stabilization;
	char** requested_outputs = NULL;
	DebrisAnalysis * debris_analysis = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(DebrisAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,DebrisNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,DebrisRequestedOutputsEnum);

	if(VerboseSolution()) _printf0_("   computing debris transport\n");

	// We need surface slopes for removal model
	surfaceslope_core(femmodel);

	/*Transport Debris*/
	femmodel->inputs->DuplicateInput(VxEnum,VxDebrisEnum);
	if(domaintype!=Domain2DverticalEnum){
		femmodel->inputs->DuplicateInput(VyEnum,VyDebrisEnum);	
	}
	femmodel->parameters->SetParam(VxEnum,InputToDepthaverageInEnum);
	femmodel->parameters->SetParam(VxAverageEnum,InputToDepthaverageOutEnum);
	depthaverage_core(femmodel);
	if(domaintype!=Domain2DverticalEnum){
		femmodel->parameters->SetParam(VyEnum,InputToDepthaverageInEnum);
		femmodel->parameters->SetParam(VyAverageEnum,InputToDepthaverageOutEnum);
		depthaverage_core(femmodel);
	}

	debris_analysis = new DebrisAnalysis();
	debris_analysis->Core(femmodel);
	delete debris_analysis;	

	femmodel->parameters->SetParam(DebrisThicknessEnum,InputToExtrudeEnum);
	extrudefromtop_core(femmodel);	

	if(save_results) femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	if(solution_type==DebrisSolutionEnum)femmodel->RequestedDependentsx();

	/*Free resources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*profiler*/
	femmodel->profiler->Stop(DEBRISCORE);
} /*}}}*/
