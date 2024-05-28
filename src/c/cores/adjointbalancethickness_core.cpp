/*!\file:  adjointbalancethickness_core.cpp
 * \brief compute inverse method adjoint state
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void adjointbalancethickness_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;

	/*retrieve parameters:*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	/*compute thickness */
	if(VerboseSolution()) _printf0_("   computing thickness\n");
	femmodel->SetCurrentConfiguration(BalancethicknessAnalysisEnum);
	solutionsequence_linear(femmodel);

	/*Call SurfaceAreax, because some it might be needed by PVector*/
	SurfaceAreax(NULL,femmodel);

	/*compute adjoint*/
	if(VerboseSolution()) _printf0_("   computing adjoint\n");
	femmodel->SetCurrentConfiguration(BalancethicknessAnalysisEnum,AdjointBalancethicknessAnalysisEnum);
	solutionsequence_adjoint_linear(femmodel);

	/*Save results*/
	if(save_results){
		int outputs[1] = {AdjointEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],1);
	}
}
