#include "./HydrologyTwsAnalysis.h"
#include <math.h>
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../classes/Inputs/TransientInput.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyTwsAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void HydrologyTwsAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void HydrologyTwsAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	::CreateNodes(nodes,iomodel,HydrologyTwsAnalysisEnum,P1Enum);
}/*}}}*/
int  HydrologyTwsAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyTwsAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int nature=0;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Plug inputs into element:*/
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.spcwatercolumn", HydrologyTwsSpcEnum);

	/*Initialize sea level cumulated sea level loads :*/
	iomodel->ConstantToInput(inputs,elements,0.,AccumulatedDeltaTwsEnum,P1Enum);
	iomodel->ConstantToInput(inputs,elements,0.,OldAccumulatedDeltaTwsEnum,P1Enum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",WatercolumnEnum);

}/*}}}*/
void HydrologyTwsAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Tws?*/
	if(hydrology_model!=HydrologyTwsEnum) return;
	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
	parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyTwsAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyTwsAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyTwsAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyTwsAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyTwsAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* HydrologyTwsAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           HydrologyTwsAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           HydrologyTwsAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyTwsAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Only update if on base*/
	if(!element->IsOnBase()) return;

	/*Fetch dof list and allocate solution vector*/
	int *doflist = NULL;
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);

	int numnodes = element->GetNumberOfNodes();
	IssmDouble* watercolumn = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		watercolumn[i]=solution[doflist[i]];
		/*Check solution*/
		if(xIsNan<IssmDouble>(watercolumn[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(watercolumn[i])) _error_("Inf found in solution vector");
	}
	element->AddBasalInput(WatercolumnEnum,watercolumn,element->GetElementType());

	xDelete<int>(doflist);
	xDelete<IssmDouble>(watercolumn);

}/*}}}*/
void           HydrologyTwsAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
