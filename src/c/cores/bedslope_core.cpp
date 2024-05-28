/*!\file: bedslope_core.cpp
 * \brief: core of the slope solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void bedslope_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;
	int  domaintype;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);

	if(VerboseSolution()) _printf0_("   computing slope\n");

	/*Call on core computations: */
	femmodel->SetCurrentConfiguration(L2ProjectionBaseAnalysisEnum);

	femmodel->parameters->SetParam(BedSlopeXEnum,InputToL2ProjectEnum);
	solutionsequence_linear(femmodel);

	if(domaintype!=Domain2DverticalEnum){
		femmodel->parameters->SetParam(BedSlopeYEnum,InputToL2ProjectEnum);
		solutionsequence_linear(femmodel);
	}

	if(save_results){
		if(domaintype!=Domain2DverticalEnum){
			int outputs[2] = {BedSlopeXEnum,BedSlopeYEnum};
			femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],2);
		}
		else{
			int outputs[1] = {BedSlopeXEnum};
			femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],1);
		}
	}

}
