#include "./HydrologyPismAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyPismAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	return;

}/*}}}*/
void HydrologyPismAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	return;

}/*}}}*/
void HydrologyPismAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	return;
}/*}}}*/
int  HydrologyPismAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 0;
}/*}}}*/
void HydrologyPismAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model,frictionlaw;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Pism?*/
	if(hydrology_model!=HydrologypismEnum) return;

	/*Add input to elements*/
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.drainage_rate",HydrologyDrainageRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.watercolumn_max",HydrologyWatercolumnMaxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",WatercolumnEnum,0.);
}/*}}}*/
void HydrologyPismAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Pism?*/
	if(hydrology_model!=HydrologypismEnum) return;
	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
	parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");

	/*Nothing else to add for now*/
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyPismAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPismAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyPismAnalysis::CreateDVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementMatrix* HydrologyPismAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyPismAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyPismAnalysis::CreatePVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPismAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPismAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyPismAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyPismAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/

/*Additional methods*/
void HydrologyPismAnalysis::UpdateWaterColumn(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		this->UpdateWaterColumn(element);
	}

}/*}}}*/
void HydrologyPismAnalysis::UpdateWaterColumn(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating()) return;

	/*Intermediaries */
	IssmDouble  dt,drainage_rate,water_column;

	/*Retrieve all inputs and parameters*/
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);

	/*Get water column and drainage rate*/
	const int  numvertices= element->GetNumberOfVertices();
	IssmDouble* watercolumn  = xNew<IssmDouble>(numvertices);
	IssmDouble* drainagerate = xNew<IssmDouble>(numvertices);
	IssmDouble* meltingrate  = xNew<IssmDouble>(numvertices);
 	IssmDouble* watercolumn_max  = xNew<IssmDouble>(numvertices);
	element->GetInputListOnVertices(&watercolumn[0],WaterColumnOldEnum);
	element->GetInputListOnVertices(&drainagerate[0],HydrologyDrainageRateEnum);
	element->GetInputListOnVertices(&meltingrate[0],BasalforcingsGroundediceMeltingRateEnum);
	element->GetInputListOnVertices(&watercolumn_max[0],HydrologyWatercolumnMaxEnum);

	/*Add water*/
	for(int i=0;i<numvertices;i++){
		watercolumn[i] += (meltingrate[i]/rho_ice*rho_water-drainagerate[i])*dt;
		/*Check that water column height is within 0 and upper bound, correct if needed*/
 		if(watercolumn[i]>watercolumn_max[i]) watercolumn[i]=watercolumn_max[i];
 		if(watercolumn[i]<0) watercolumn[i]=0.;
	}

	/* Divide by connectivity, add degree of channelization as an input */
	/*FIXME: should be changed to P1, this is due to the NR, IsAllFloating will return 0 on this element, but it should not be DG*/
	element->AddInput(WatercolumnEnum,&watercolumn[0],P1DGEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(watercolumn);
	xDelete<IssmDouble>(meltingrate);
	xDelete<IssmDouble>(watercolumn_max);
	xDelete<IssmDouble>(drainagerate);
}/*}}}*/
