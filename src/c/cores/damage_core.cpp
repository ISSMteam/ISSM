/*
 * \brief: damage_core.cpp: core for the damage solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void damage_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(DAMAGECORE);

	/*intermediary*/
	bool   save_results;
	bool   dakota_analysis     = false;
	int    solution_type,stabilization;
	int    numoutputs          = 0;
	char   **requested_outputs = NULL;

	//first recover parameters common to all solutions
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,DamageEvolutionNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,DamageEvolutionRequestedOutputsEnum);
	femmodel->parameters->FindParam(&stabilization,DamageStabilizationEnum);

	if(VerboseSolution()) _printf0_("   computing damage\n");
	Damagex(femmodel); /* optionally calculate damage analytically first */
	femmodel->SetCurrentConfiguration(DamageEvolutionAnalysisEnum);
	if(stabilization==4){
		solutionsequence_fct(femmodel);
	}
	else{
		solutionsequence_linear(femmodel);
	}

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	/*Free resources:*/
	if(numoutputs){
		for(int i=0;i<numoutputs;i++){
			xDelete<char>(requested_outputs[i]);
		}
		xDelete<char*>(requested_outputs);
	}

	/*End profiler*/
	femmodel->profiler->Stop(DAMAGECORE);
}
