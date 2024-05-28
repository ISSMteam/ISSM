/*!\file: bmb_core.cpp
 * \brief: core of the bmb (Basal mass balance) solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void bmb_core(FemModel* femmodel){

	/*First, get BMB model from parameters*/
	int  basalforcing_model;
	bool isplume = false;
	femmodel->parameters->FindParam(&basalforcing_model,BasalforcingsEnum);

	if(VerboseSolution()) _printf0_("   computing basal mass balance\n");

	/*In some cases we need to run additional analyses to get the required input data*/
	if(basalforcing_model==BasalforcingsPicoEnum){
		femmodel->parameters->FindParam(&isplume,BasalforcingsPicoIsplumeEnum);
		if(isplume){
			femmodel->SetCurrentConfiguration(L2ProjectionBaseAnalysisEnum);
			femmodel->parameters->SetParam(BaseSlopeXEnum,InputToL2ProjectEnum);
			solutionsequence_linear(femmodel);
			femmodel->parameters->SetParam(BaseSlopeYEnum,InputToL2ProjectEnum);
			solutionsequence_linear(femmodel);
			femmodel->SetCurrentConfiguration(GLheightadvectionAnalysisEnum);
			solutionsequence_linear(femmodel);
		}
	}

	/*Call module now*/
	FloatingiceMeltingRatex(femmodel);

	/*Extrude basal melt if not default melting rate (which may be a transient input that can't be extruded)*/
	if(basalforcing_model!=FloatingMeltRateEnum){
		femmodel->parameters->SetParam(BasalforcingsFloatingiceMeltingRateEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
	}
}
