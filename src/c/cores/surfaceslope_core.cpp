/*!\file: surfaceslope_core.cpp
 * \brief: core of the slope solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

void surfaceslope_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;
	int  domaintype;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);

	if(VerboseSolution()) _printf0_("computing slope...\n");

	/*Call on core computations: */
	femmodel->SetCurrentConfiguration(L2ProjectionBaseAnalysisEnum);

	femmodel->parameters->SetParam(SurfaceSlopeXEnum,InputToL2ProjectEnum);
	solutionsequence_linear(femmodel);

	if(domaintype!=Domain2DverticalEnum){
		femmodel->parameters->SetParam(SurfaceSlopeYEnum,InputToL2ProjectEnum);
		solutionsequence_linear(femmodel);
	}
	if(domaintype==Domain2DverticalEnum){
		femmodel->parameters->SetParam(SurfaceSlopeXEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
	}

	if(save_results){
		if(domaintype!=Domain2DverticalEnum){
			int outputs[2] = {SurfaceSlopeXEnum,SurfaceSlopeYEnum};
			femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],2);

		}
		else{
			int outputs = SurfaceSlopeXEnum;
			femmodel->RequestedOutputsx(&femmodel->results,&outputs,1);
		}
	}

}
