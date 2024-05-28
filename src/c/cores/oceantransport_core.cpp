/*!\file: oceantransport_core.cpp
 * \brief: core of the oceantransport solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void SolidEarthOceanUpdates(FemModel* femmodel);
void oceantransport_core(FemModel* femmodel){ /*{{{*/

	/*Start profiler*/
	femmodel->profiler->Start(OCEANTRANSPORTCORE);

	/*parameters: */
	int    numoutputs;
	bool   save_results;
	bool   dakota_analysis;
	int    solution_type;
	Vector<IssmDouble>*  ug  = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(OceantransportAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);

	if(VerboseSolution()) _printf0_("   computing ocean mass transport\n");

	femmodel->SetCurrentConfiguration(OceantransportAnalysisEnum);

	/*save current bottom pressures before updating:*/
	InputDuplicatex(femmodel,BottomPressureEnum,BottomPressureOldEnum);
	InputDuplicatex(femmodel,DslEnum,DslOldEnum);
	InputDuplicatex(femmodel,StrEnum,StrOldEnum);

	/*grab bottom pressures, dsl and str from OceantransportSpcbottompressure, OceantransportSpcdslEnum 
	 * and OceantransportSpcstrEnum  inputs in each element, assemble into a vector and feed to 
	 * InputUpdateFromSolutionx which will deal with accumulating such inputs:*/
	GetSolutionFromInputsx(&ug,femmodel); 
	InputUpdateFromSolutionx(femmodel,ug); 

	SolidEarthOceanUpdates(femmodel);

	if(solution_type==OceantransportSolutionEnum)femmodel->RequestedDependentsx();

	/*profiler*/
	femmodel->profiler->Stop(OCEANTRANSPORTCORE);

	/*Free resources:*/
	delete ug;

} /*}}}*/
void SolidEarthOceanUpdates(FemModel* femmodel){ /*{{{*/

	int isgrd;
	int grdmodel;
	IssmDouble time;
	int frequency,count;

	/*retrieve parameters:*/
	femmodel->parameters->FindParam(&isgrd,SolidearthSettingsGRDEnum);
	femmodel->parameters->FindParam(&grdmodel,GrdModelEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);
	femmodel->parameters->FindParam(&count,SealevelchangeRunCountEnum);

	/*early return?:*/
	if(!isgrd)return;

	/* From old and new bottom pressures, create delta bottom pressure, delta dsl and delta str. 
	 * Accumulate delta bottom pressure: */
	femmodel->inputs->ZAXPY(-1, BottomPressureOldEnum,BottomPressureEnum,DeltaBottomPressureEnum);
	femmodel->inputs->AXPY(+1, DeltaBottomPressureEnum,AccumulatedDeltaBottomPressureEnum);

	femmodel->inputs->ZAXPY(-1, DslOldEnum,DslEnum,DeltaDslEnum);
	femmodel->inputs->ZAXPY(-1, StrOldEnum,StrEnum,DeltaStrEnum);

	/* Compute total bottom pressure change between two sea-level solver time steps, ie. every frequency*dt. */
	if(count==frequency){
		femmodel->inputs->ZAXPY(-1, OldAccumulatedDeltaBottomPressureEnum,AccumulatedDeltaBottomPressureEnum,DeltaBottomPressureEnum);
		femmodel->inputs->DuplicateInput(AccumulatedDeltaBottomPressureEnum,OldAccumulatedDeltaBottomPressureEnum);
	}
	return;
}/*}}}*/
