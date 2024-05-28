/*!\file:  adjointbalancethickness2_core.cpp
 * \brief compute inverse method adjoint state
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void adjointbalancethickness2_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;

	/*retrieve parameters:*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	/*compute thickness2 */
	if(VerboseSolution()) _printf0_("   computing thickness2\n");
	femmodel->SetCurrentConfiguration(Balancethickness2AnalysisEnum);
	solutionsequence_linear(femmodel);

	/*Call SurfaceAreax, because some it might be needed by PVector*/
	//SurfaceAreax(NULL,femmodel);

	/*compute adjoint*/
	if(VerboseSolution()) _printf0_("   computing adjoint\n");
	femmodel->SetCurrentConfiguration(Balancethickness2AnalysisEnum,AdjointBalancethickness2AnalysisEnum);
	solutionsequence_adjoint_linear(femmodel);

	/*Save results*/
	if(save_results || true){
		int outputs[1] = {AdjointEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],1);
	}
}
