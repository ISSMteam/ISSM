/*!\file: levelsetfunctionslope_core.cpp
 * \brief: core of the slope solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

void levelsetfunctionslope_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;
	int  domaintype;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);

	if(VerboseSolution()) _printf0_("   computing slope of levelset function...\n");

	/*Call on core computations: */
	femmodel->SetCurrentConfiguration(L2ProjectionBaseAnalysisEnum);

	femmodel->parameters->SetParam(LevelsetfunctionSlopeXEnum,InputToL2ProjectEnum);
	solutionsequence_linear(femmodel);

	if(domaintype!=Domain2DverticalEnum){
		femmodel->parameters->SetParam(LevelsetfunctionSlopeYEnum,InputToL2ProjectEnum);
		solutionsequence_linear(femmodel);
	}
	if(domaintype==Domain2DverticalEnum){
		femmodel->parameters->SetParam(LevelsetfunctionSlopeXEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
	}
}
