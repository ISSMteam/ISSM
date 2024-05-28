/*!\file: thermal_core.cpp
 * \brief: core of the thermal solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../analyses/analyses.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void thermal_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(THERMALCORE);

	/*intermediary*/
	bool   save_results,isenthalpy;
	bool   dakota_analysis;
	int    solution_type,numoutputs;
	char** requested_outputs = NULL;
	EnthalpyAnalysis * enthalpy_analysis = NULL;

	/*first recover parameters common to all solutions*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
	femmodel->parameters->FindParam(&numoutputs,ThermalNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,ThermalRequestedOutputsEnum);

	/*Calculate geothermalflux*/
	GeothermalFluxx(femmodel);

	if(isenthalpy){
		femmodel->InputMakeDiscontinuous(BasalforcingsGroundediceMeltingRateEnum);
		enthalpy_analysis = new EnthalpyAnalysis();
		enthalpy_analysis->Core(femmodel);
		delete enthalpy_analysis;
	}
	else{
		if(VerboseSolution()) _printf0_("   computing temperatures\n");
		femmodel->SetCurrentConfiguration(ThermalAnalysisEnum);
		solutionsequence_thermal_nonlinear(femmodel);

		if(VerboseSolution()) _printf0_("   computing melting\n");
		femmodel->SetCurrentConfiguration(MeltingAnalysisEnum);
		solutionsequence_linear(femmodel);
	}

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	/*Free resources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*End profiler*/
        femmodel->profiler->Stop(THERMALCORE);
}
