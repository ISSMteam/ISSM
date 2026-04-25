#include "./HydrologyPrescribeAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyPrescribeAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*Nothing to be done*/
}/*}}}*/
void HydrologyPrescribeAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*Nothing to be done*/
}/*}}}*/
void HydrologyPrescribeAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Prescribe?*/
	if(hydrology_model!=HydrologyprescribeEnum) return;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,HydrologyPrescribeAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");

}/*}}}*/
int  HydrologyPrescribeAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyPrescribeAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model,frictionlaw;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Prescribe?*/
	if(hydrology_model!=HydrologyprescribeEnum) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Add input to elements*/
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.head",HydrologyHeadEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);

	/*Initialize requested outputs in case they are not defined later for this partition*/
	iomodel->ConstantToInput(inputs,elements,0.,EffectivePressureEnum,P1Enum);
}/*}}}*/
void HydrologyPrescribeAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Prescribe?*/
	if(hydrology_model!=HydrologyprescribeEnum) return;
	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
	parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");

	/*Nothing else to add for now*/
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyPrescribeAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPrescribeAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyPrescribeAnalysis::CreateDVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementMatrix* HydrologyPrescribeAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyPrescribeAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyPrescribeAnalysis::CreatePVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPrescribeAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPrescribeAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyPrescribeAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPrescribeAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/

/*Additional methods*/
void HydrologyPrescribeAnalysis::UpdateEffectivePressure(FemModel* femmodel){/*{{{*/

	/*Loop over each element to compute Subglacial Water Pressure at vertices*/
   for(Object* &object:femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		UpdateEffectivePressure(element);
	}

}/*}}}*/
void HydrologyPrescribeAnalysis::UpdateEffectivePressure(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(!element->IsOnBase()) return;

	/*Intermediaries*/
	IssmDouble bed,thickness,head;

	/* Fetch number of nodes and allocate output*/
   int numnodes = element->GetNumberOfNodes();
   IssmDouble* N = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	IssmDouble  g          = element->FindParam(ConstantsGEnum);
	IssmDouble  rho_ice    = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water  = element->FindParam(MaterialsRhoFreshwaterEnum);
	Input* head_input      = element->GetInput(HydrologyHeadEnum); _assert_(head_input);
	Input* thickness_input = element->GetInput(ThicknessEnum);     _assert_(thickness_input);
	Input* base_input      = element->GetInput(BaseEnum);          _assert_(base_input);

   Gauss* gauss=element->NewGauss();
   for (int i=0;i<numnodes;i++){
      gauss->GaussVertex(i);

		base_input->GetInputValue(&bed,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		head_input->GetInputValue(&head,gauss);

		N[i] = rho_ice*g*thickness - rho_water*g*(head-bed);
	}

	/*Set to 0 if inactive element*/
   if(element->IsAllFloating() || !element->IsIceInElement()){
      for(int iv=0;iv<numnodes;iv++) N[iv] = 0.;
      element->AddInput(EffectivePressureEnum,N,P1Enum);
      xDelete<IssmDouble>(N);
      return;
   }

	/*Add new gap as an input*/
	element->AddInput(EffectivePressureEnum,N,P1Enum);

	/*Clean up and return*/
   xDelete<IssmDouble>(N);
	delete gauss;
}/*}}}*/

