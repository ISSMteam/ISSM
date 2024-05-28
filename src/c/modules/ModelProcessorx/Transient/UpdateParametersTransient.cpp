/*
 * UpdateParametersTransient.cpp:
 */

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx.h"

void UpdateParametersTransient(Parameters* parameters,IoModel* iomodel){/*{{{*/
	
	bool isgroundingline;
	int numoutputs;
	char** requestedoutputs;
	parameters->FindParam(&isgroundingline,TransientIsgroundinglineEnum);
	if(isgroundingline){
		iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.groundingline.requested_outputs");
		parameters->AddObject(new IntParam(GroundinglineNumRequestedOutputsEnum,numoutputs));
      if(numoutputs)parameters->AddObject(new StringArrayParam(GroundinglineRequestedOutputsEnum,requestedoutputs,numoutputs));
      iomodel->DeleteData(&requestedoutputs,numoutputs,"md.groundingline.requested_outputs");
	}
}/*}}}*/


