#include "./MmemasstransportAnalysis.h"
#include <math.h>
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../classes/Inputs/TransientInput.h"
#include "../classes/Inputs/TriaInput.h"
#include "../classes/gauss/Gauss.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void MmemasstransportAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void MmemasstransportAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void MmemasstransportAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	::CreateNodes(nodes,iomodel,MmemasstransportAnalysisEnum,P1Enum);
}/*}}}*/
int  MmemasstransportAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 3;
}/*}}}*/
void MmemasstransportAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int  nature=0;
	bool isdakota=0;

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
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MmemasstransportMaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MmemasstransportMaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mmemasstransport.thickness", MmemasstransportThicknessEnum);

	/*Initialize sea level cumulated sea level loads :*/
	iomodel->ConstantToInput(inputs,elements,0.,AccumulatedDeltaIceThicknessEnum,P0Enum);
	iomodel->ConstantToInput(inputs,elements,0.,OldAccumulatedDeltaIceThicknessEnum,P0Enum);
	iomodel->ConstantToInput(inputs,elements,0.,DeltaIceThicknessEnum,P0Enum); 

}/*}}}*/
void MmemasstransportAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	int nids,npart,nel;
	IssmDouble* ids=NULL; 
	IssmDouble* partition = NULL;

	iomodel->FetchData(&nel,"md.mesh.numberofelements");
	iomodel->FetchData(&ids,&nids,NULL,"md.mmemasstransport.ids");
	//_printf_("nids: " << nids << "\n"); for(int i=0;i<nids;i++)_printf_(ids[i] << "|");  _printf_("\n");
	parameters->AddObject(new DoubleMatParam(MmemasstransportModelidsEnum,ids,nids,1));
	iomodel->FetchData(&partition,&npart,NULL,"md.mmemasstransport.partition");
	if (npart!=nel)_error_("MmemasstransportAnalysis:UpdateParameters: partition vector should be distributed over elements, not vertices!");
	parameters->AddObject(new DoubleMatParam(MmemasstransportPartitionEnum,partition,nel,1));
	
	xDelete<IssmDouble>(ids);
	xDelete<IssmDouble>(partition);
	
	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.mmemasstransport.requested_outputs");
	parameters->AddObject(new IntParam(MmemasstransportNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(MmemasstransportRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.mmemasstransport.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           MmemasstransportAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           MmemasstransportAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* MmemasstransportAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* MmemasstransportAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* MmemasstransportAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* MmemasstransportAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           MmemasstransportAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	/*do nothing:*/
	return;
}/*}}}*/
void           MmemasstransportAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           MmemasstransportAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	IssmDouble time;
	IssmDouble ice[3];
	IssmDouble ocean[3];
	IssmDouble height;
	int numnodes = element->GetNumberOfNodes();
	
	element->parameters->FindParam(&time,TimeEnum);

	TriaInput* h_input=xDynamicCast<TriaInput*>(element->GetInput(MmemasstransportThicknessEnum,time)); _assert_(h_input);
	TriaInput* i_input=xDynamicCast<TriaInput*>(element->GetInput(MmemasstransportMaskIceLevelsetEnum,time)); _assert_(i_input);
	TriaInput* o_input=xDynamicCast<TriaInput*>(element->GetInput(MmemasstransportMaskOceanLevelsetEnum,time)); _assert_(o_input);

	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<3;iv++){
		gauss->GaussVertex(iv);
		i_input->GetInputValue(&ice[iv],gauss);
		o_input->GetInputValue(&ocean[iv],gauss);
	}
	h_input->GetInputAverage(&height);
	delete gauss;

	/*Add thickness ice and ocean levelsets as inputs to the tria element: */
	element->AddInput(ThicknessEnum,&height,P0Enum); //very important, do not change the type of ThicknessEnum to P1 when it should be P0! 
	element->AddInput(MaskIceLevelsetEnum,ice,P1Enum);
	element->AddInput(MaskOceanLevelsetEnum,ocean,P1Enum);

}/*}}}*/
void           MmemasstransportAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
