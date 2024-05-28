/*!\file: groundingline_core.cpp
 * \brief: core of the groundingline solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void groundingline_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(GROUNDINGLINECORE);

	/* intermediaries */
	int numoutputs;
	bool save_results;
	char** requested_outputs = NULL;

	/* recover parameters */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&numoutputs,GroundinglineNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,GroundinglineRequestedOutputsEnum);

	/*Move grounding line*/
	if(VerboseSolution()) _printf0_("   computing new grounding line position\n");
	GroundinglineMigrationx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);

	/*Update geometry and mask accordingly*/
	femmodel->parameters->SetParam(MaskOceanLevelsetEnum,InputToExtrudeEnum);
	extrudefrombase_core(femmodel);
	femmodel->parameters->SetParam(BaseEnum,InputToExtrudeEnum);
	extrudefrombase_core(femmodel);
	femmodel->parameters->SetParam(SurfaceEnum,InputToExtrudeEnum);
	extrudefrombase_core(femmodel);

	/*Save results*/
	if(save_results) femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);

	/*Free resources:*/
   if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}

	/*Stop profiler*/
	femmodel->profiler->Stop(GROUNDINGLINECORE);
}
