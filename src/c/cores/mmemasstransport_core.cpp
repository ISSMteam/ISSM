/*!\file: mmemasstransport_core.cpp
 * \brief: core of the mmemasstransport solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void SolidEarthMmemasstransportUpdates(FemModel* femmodel);
void mmemasstransport_core(FemModel* femmodel){ /*{{{*/

	/*Start profiler*/
	femmodel->profiler->Start(MMEMASSTRANSPORTCORE);

	/*parameters: */
	int    numoutputs;
	bool   save_results;
	bool   dakota_analysis;
	int    solution_type;
	Vector<IssmDouble>*  ug  = NULL;
	char** requested_outputs = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(MmemasstransportAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&numoutputs,MmemasstransportNumRequestedOutputsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,MmemasstransportRequestedOutputsEnum);

	if(VerboseSolution()) _printf0_("   computing MME mass transport\n");

	/*save current thickness before updating, and make sure it's P0:*/
	femmodel->InputToP0(ThicknessEnum,ThicknessOldEnum);

	/*grab thickness and ice/ocean levelsets from MmemasstransportThicknessEnum inputs in each element, assemble into a vector and feed to 
	 * InputUpdateFromSolutionx which will deal with accumulating such inputs:*/
	GetSolutionFromInputsx(&ug,femmodel); 
	InputUpdateFromSolutionx(femmodel,ug); 

	SolidEarthMmemasstransportUpdates(femmodel);

	if(save_results){
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}
	if(solution_type==MmemasstransportSolutionEnum)femmodel->RequestedDependentsx();

	/*Free ressources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*profiler*/
	femmodel->profiler->Stop(MMEMASSTRANSPORTCORE);

	/*free ressources:*/
	delete ug;

} /*}}}*/
void SolidEarthMmemasstransportUpdates(FemModel* femmodel){ /*{{{*/

	IssmDouble time;
	int frequency,count;
	
	if(VerboseModule()) _printf0_("   SolidEarth update for Mme mass transport solution");

	/*retrieve parameters:*/
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);
	femmodel->parameters->FindParam(&count,SealevelchangeRunCountEnum);

	/* From old and new thickness, create delta ice thicknes, and accumulate:*/
	femmodel->inputs->ZAXPY(-1, ThicknessOldEnum,ThicknessEnum,DeltaIceThicknessEnum);
	femmodel->inputs->AXPY(+1, DeltaIceThicknessEnum,AccumulatedDeltaIceThicknessEnum);

	/* Compute total ice thickness change between two sea-level solver time steps, ie. every frequency*dt. */
	if(count==frequency){
		femmodel->inputs->ZAXPY(-1, OldAccumulatedDeltaIceThicknessEnum,AccumulatedDeltaIceThicknessEnum,DeltaIceThicknessEnum);
		femmodel->inputs->DuplicateInput(AccumulatedDeltaIceThicknessEnum,OldAccumulatedDeltaIceThicknessEnum);
	}
	return;
}/*}}}*/
