#include "./HydrologyDCEfficientAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
int  HydrologyDCEfficientAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyDCEfficientAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	bool        isefficientlayer;
	int         hydrology_model;
	int         eplflip_lock;
	int         eplthickcomp;
	IssmDouble  eplinitthick;
	IssmDouble  eplcolapsethick;
	IssmDouble  eplmaxthick;
	IssmDouble  eplcond;

	/*retrieve some parameters: */
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");

	/*If not return*/
	if(!isefficientlayer) return;

	/*If yes, get the parameters*/
	iomodel->FetchData(&eplflip_lock,"md.hydrology.eplflip_lock");
	iomodel->FetchData(&eplthickcomp,"md.hydrology.epl_thick_comp");

	parameters->AddObject(new IntParam(HydrologydcEplflipLockEnum,eplflip_lock));
	parameters->AddObject(new IntParam(HydrologydcEplThickCompEnum,eplthickcomp));

	iomodel->FetchData(&eplinitthick,"md.hydrology.epl_initial_thickness");
	iomodel->FetchData(&eplcolapsethick,"md.hydrology.epl_colapse_thickness");
	iomodel->FetchData(&eplmaxthick,"md.hydrology.epl_max_thickness");
	iomodel->FetchData(&eplcond,"md.hydrology.epl_conductivity");
	parameters->AddObject(new DoubleParam(HydrologydcEplInitialThicknessEnum,eplinitthick));
	parameters->AddObject(new DoubleParam(HydrologydcEplColapseThicknessEnum,eplcolapsethick));
	parameters->AddObject(new DoubleParam(HydrologydcEplMaxThicknessEnum,eplmaxthick));
	parameters->AddObject(new DoubleParam(HydrologydcEplConductivityEnum,eplcond));

}/*}}}*/
void HydrologyDCEfficientAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	bool   isefficientlayer;
	int    hydrology_model;

	/*Now, do we really want DC?*/
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.epl_head",EplHeadSubstepEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sediment_head",SedimentHeadSubstepEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.epl_thickness",HydrologydcEplThicknessSubstepEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.basal_moulin_input",HydrologydcBasalMoulinInputEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}/*}}}*/
void HydrologyDCEfficientAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Now, do we really want DC?*/
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	bool isefficientlayer;
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	}
	::CreateNodes(nodes,iomodel,HydrologyDCEfficientAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
void HydrologyDCEfficientAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Do we really want DC?*/
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	bool isefficientlayer;
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spcepl_head",HydrologyDCEfficientAnalysisEnum,P1Enum);
}/*}}}*/
void HydrologyDCEfficientAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Do we really want DC?*/
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	bool isefficientlayer;
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	/*Fetch parameters: */
	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(1,"md.mesh.vertexonbase");
	}

	//Add moulin inputs as loads
	CreateSingleNodeToElementConnectivity(iomodel);
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");
}/*}}}*/
void HydrologyDCEfficientAnalysis::InitZigZagCounter(FemModel* femmodel){/*{{{*/

	int*   eplzigzag_counter =NULL;
	eplzigzag_counter=xNewZeroInit<int>(femmodel->nodes->Size());
	femmodel->parameters->AddObject(new IntVecParam(EplZigZagCounterEnum,eplzigzag_counter,femmodel->nodes->Size()));
	xDelete<int>(eplzigzag_counter);
}/*}}}*/
void HydrologyDCEfficientAnalysis::ResetCounter(FemModel* femmodel){/*{{{*/

	int*     eplzigzag_counter=NULL;
	femmodel->parameters->FindParam(&eplzigzag_counter,NULL,EplZigZagCounterEnum);
	for(int i=0;i<femmodel->nodes->Size();i++){
		eplzigzag_counter[i]=0;
	}
	femmodel->parameters->SetParam(eplzigzag_counter,femmodel->nodes->Size(),EplZigZagCounterEnum);
	xDelete<int>(eplzigzag_counter);
}/*}}}*/

/*Finite Element Analysis*/
void HydrologyDCEfficientAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void HydrologyDCEfficientAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyDCEfficientAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyDCEfficientAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyDCEfficientAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	bool     active_element;
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);
	/*Check that all nodes are active, else return empty matrix*/
	if(!active_element) {
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		return NULL;
	}

	/* Intermediaries */
	IssmDouble  D_scalar,Jdet,dt;
	IssmDouble  transfer;
	IssmDouble  epl_transmitivity;
	IssmDouble  epl_storing;
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* epl_thick_input = basalelement->GetInput(HydrologydcEplThicknessSubstepEnum); _assert_(epl_thick_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = basalelement->NewGauss(2);
	while(gauss->next()){
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		epl_transmitivity = EplTransmitivity(basalelement,gauss,epl_thick_input);
		epl_storing			= EplStoring(basalelement,gauss,epl_thick_input);

		/*Diffusivity*/
		D_scalar=epl_transmitivity*gauss->weight*Jdet;
		if(dt!=0.) D_scalar=D_scalar*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D_scalar*(dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]);
			}
		}
		/*Transient*/
		if(dt!=0.){
			D_scalar=epl_storing*gauss->weight*Jdet;
			for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[j]*basis[i];
			/*Transfer EPL part*/
			transfer=GetHydrologyKMatrixTransfer(basalelement);
			D_scalar=dt*transfer*gauss->weight*Jdet;
			for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[j]*basis[i];
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* HydrologyDCEfficientAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	bool     active_element;
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);

	/*Check that all nodes are active, else return empty matrix*/
	if(!active_element) {
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		return NULL;
	}
	/*Intermediaries */
	int        smb_model,smb_averaging;
	int        smbsubstepping,hydrologysubstepping;
	IssmDouble dt,scalar,water_head;
	IssmDouble water_load,transfer,runoff_value;
	IssmDouble epl_storing;  //,epl_transmitivity;
	IssmDouble Jdet,time;
	IssmDouble residual,connectivity;
	IssmDouble active_node;

	IssmDouble *xyz_list             = NULL;
	Input     *old_wh_input         = NULL;
	Input     *dummy_input          = NULL;
	Input     *surface_runoff_input = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = basalelement->GetNumberOfNodes();
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize Element vector*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);


	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement ->FindParam(&smb_model,SmbEnum);

	Input*	epl_thick_input  = basalelement->GetInput(HydrologydcEplThicknessSubstepEnum); _assert_(epl_thick_input);
	Input*	sed_head_input   = basalelement->GetInput(SedimentHeadSubstepEnum); _assert_(sed_head_input);
	Input*	basal_melt_input = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(basal_melt_input);
	Input*	residual_input   = basalelement->GetInput(SedimentHeadResidualEnum); _assert_(residual_input);

	if(dt!= 0.){
		old_wh_input = basalelement->GetInput(EplHeadOldEnum);            _assert_(old_wh_input);
	}
	if(smb_model==SMBgradientscomponentsEnum){
		basalelement->FindParam(&time,TimeEnum);
		basalelement->FindParam(&smbsubstepping,SmbStepsPerStepEnum);
		basalelement->FindParam(&hydrologysubstepping,HydrologyStepsPerStepEnum);

		if(smbsubstepping==1){
			//no substeping for the smb we take the result from there
			dummy_input = basalelement->GetInput(SmbRunoffEnum); _assert_(dummy_input);
		}
		else if(smbsubstepping>1 && smbsubstepping<=hydrologysubstepping){
			//finer hydro stepping, we take the value at the needed time
			dummy_input = basalelement->GetInput(SmbRunoffTransientEnum, time); _assert_(dummy_input);
		}
		else{
			//finer stepping in smb, we average the runoff from transient input
			basalelement->FindParam(&smb_averaging,SmbAveragingEnum);
			dummy_input = basalelement->GetInput(SmbRunoffTransientEnum,time-dt,time,smb_averaging); _assert_(dummy_input);
		}
		surface_runoff_input=xDynamicCast<Input*>(dummy_input); _assert_(surface_runoff_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = basalelement->NewGauss(2);
	while(gauss->next()){
		basalelement ->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement ->NodalFunctions(basis,gauss);
		//epl_transmitivity = EplTransmitivity(basalelement,gauss,epl_thick_input);
		epl_storing	= EplStoring(basalelement,gauss,epl_thick_input);

		/*Loading term*/
		basal_melt_input->GetInputValue(&water_load,gauss);
		if(surface_runoff_input) surface_runoff_input->GetInputValue(&runoff_value,gauss);
		else                     runoff_value = 0.;
		scalar = Jdet*gauss->weight*(water_load+runoff_value);
		if(dt!=0.) scalar = scalar*dt;
		for(int i=0;i<numnodes;i++){
			//This is the original
			pe->values[i]+=scalar*basis[i];
			//This is the noded version
			/* basalelement->GetInputValue(&active_node,basalelement->nodes[i],HydrologydcMaskEplactiveNodeEnum); */
			/* if(!reCast<bool>(active_node)){ */
			/* 	pe->values[i]+=scalar*basis[i]; */
			//}
			//if(basalelement->nodes[i]->Sid()==42)_printf_("EPL uni Input "<<scalar*basis[i]<<"\n");
		}
		/*Transient and transfer terms*/
		if(dt!=0.){
			old_wh_input->GetInputValue(&water_head,gauss);
			/*Dealing with the epl part of the transfer term*/
			transfer=GetHydrologyPVectorTransfer(basalelement,gauss,sed_head_input);
			scalar = Jdet*gauss->weight*((water_head*epl_storing)+(dt*transfer));
			for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
		}
	}
	delete gauss;

	/*	Add residual if necessary*/
	gauss = basalelement->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);
		//epl_transmitivity = EplTransmitivity(basalelement,gauss,epl_thick_input);
		connectivity = IssmDouble(basalelement->VertexConnectivity(iv));
		residual_input->GetInputValue(&residual,gauss);
		pe->values[iv]+=residual/connectivity;
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void HydrologyDCEfficientAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,EplHeadSubstepEnum);
}/*}}}*/
void HydrologyDCEfficientAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void HydrologyDCEfficientAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	/*Intermediaries*/
	int      domaintype;
	Element* basalelement=NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	/*Intermediary*/
	int* doflist = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	basalelement->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* eplHeads = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	/*If the EPL is not active we revert to the bedrock elevation when deactivating*/
	for(int i=0;i<numnodes;i++){
		eplHeads[i]=solution[doflist[i]];
		if(xIsNan<IssmDouble>(eplHeads[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(eplHeads[i])) _error_("Inf found in solution vector");
	}
	/*Add input to the element: */
	element->AddBasalInput(EplHeadSubstepEnum,eplHeads,P1Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(eplHeads);
	xDelete<int>(doflist);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
} /*}}}*/
void HydrologyDCEfficientAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/

/*Intermediaries*/
IssmDouble HydrologyDCEfficientAnalysis::EplStoring(Element* element,Gauss* gauss, Input* epl_thick_input){/*{{{*/
	IssmDouble epl_storing;
	IssmDouble epl_thickness;
	IssmDouble porewater_mass           = element->FindParam(HydrologydcEplPoreWaterMassEnum);
	IssmDouble layer_compressibility    = element->FindParam(HydrologydcEplLayerCompressibilityEnum);

	epl_thick_input->GetInputValue(&epl_thickness,gauss);
	epl_storing=porewater_mass*epl_thickness*layer_compressibility;

	return epl_storing;
}/*}}}*/
IssmDouble HydrologyDCEfficientAnalysis::EplTransmitivity(Element* element,Gauss* gauss, Input* epl_thick_input){/*{{{*/
	IssmDouble epl_transmitivity;
	IssmDouble epl_thickness;
	IssmDouble epl_conductivity      = element->FindParam(HydrologydcEplConductivityEnum);
	epl_thick_input->GetInputValue(&epl_thickness,gauss);

	epl_transmitivity=epl_conductivity*epl_thickness;
	return epl_transmitivity;
}/*}}}*/
void HydrologyDCEfficientAnalysis::GetHydrologyDCInefficientHmax(IssmDouble* ph_max,Element* element, Node* innode){/*{{{*/

	int        hmax_flag;
	IssmDouble h_max;
	IssmDouble rho_ice,rho_water;
	IssmDouble thickness,bed;
	/*Get the flag to the limitation method*/
	element->FindParam(&hmax_flag,HydrologydcSedimentlimitFlagEnum);

	/*Switch between the different cases*/
	switch(hmax_flag){
	case 0:
		h_max=1.0e+10;
		break;
	case 1:
		element->FindParam(&h_max,HydrologydcSedimentlimitEnum);
		break;
	case 2:
		/*Compute max*/
		rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
		rho_ice   = element->FindParam(MaterialsRhoIceEnum);
		element-> GetInputValue(&thickness,innode,ThicknessEnum);
		element-> GetInputValue(&bed,innode,BaseEnum);
		h_max=((rho_ice*thickness)/rho_water)+bed;
		break;
	case 3:
		_error_("Using normal stress  not supported yet");
		break;
	default:
		_error_("no case higher than 3 for SedimentlimitFlag");
	}
	/*Assign output pointer*/
	*ph_max=h_max;
}
/*}}}*/
IssmDouble HydrologyDCEfficientAnalysis::GetHydrologyKMatrixTransfer(Element* element){/*{{{*/

	int transfermethod;
	IssmDouble leakage,transfer;

	element->FindParam(&transfermethod,HydrologydcTransferFlagEnum);
	/*Switch between the different transfer methods cases*/
	switch(transfermethod){
	case 0:
		/*Just keepping the transfer to zero*/
		transfer=0.0;
		break;
	case 1:
		element->FindParam(&leakage,HydrologydcLeakageFactorEnum);
		transfer=+leakage;
		break;
	default:
		_error_("no case higher than 1 for the Transfer method");
	}
	return transfer;
}/*}}}*/
IssmDouble HydrologyDCEfficientAnalysis::GetHydrologyPVectorTransfer(Element* element, Gauss* gauss, Input* sed_head_input){/*{{{*/

	int transfermethod;
	IssmDouble sediment_head;
	IssmDouble leakage,transfer;

	element->FindParam(&transfermethod,HydrologydcTransferFlagEnum);
	/*Switch between the different transfer methods cases*/
	switch(transfermethod){
	case 0:
		/*Just keepping the transfer to zero*/
		transfer=0.0;
		break;
	case 1:
		_assert_(sed_head_input);
		sed_head_input->GetInputValue(&sediment_head,gauss);
		element->FindParam(&leakage,HydrologydcLeakageFactorEnum);
		transfer=+sediment_head*leakage;
		break;
	default:
		_error_("no case higher than 1 for the Transfer method");
	}

	return transfer;
}/*}}}*/
void HydrologyDCEfficientAnalysis::ComputeEPLThickness(FemModel* femmodel){/*{{{*/

	int iseplthickcomp;

	/*Skip if we don't want to compute thicknesses*/
	femmodel->parameters->FindParam(&iseplthickcomp,HydrologydcEplThickCompEnum);
	if(iseplthickcomp==0) return;

	/*Get Parameters*/
	int        domaintype;
	IssmDouble max_thick;
	IssmDouble epl_conductivity;
	IssmDouble latentheat;
	IssmDouble rho_ice;
	IssmDouble rho_water;
	IssmDouble gravity;
	IssmDouble dt;

	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&max_thick,HydrologydcEplMaxThicknessEnum);
	femmodel->parameters->FindParam(&epl_conductivity,HydrologydcEplConductivityEnum);
	femmodel->parameters->FindParam(&latentheat,MaterialsLatentheatEnum);
	femmodel->parameters->FindParam(&rho_ice,MaterialsRhoIceEnum);
	femmodel->parameters->FindParam(&rho_water,MaterialsRhoFreshwaterEnum);
	femmodel->parameters->FindParam(&gravity,ConstantsGEnum);
	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);

		/*skip element if 3d and not on base*/
		if(domaintype==Domain3DEnum && !element->IsOnBase()) continue;

		int         numnodes      = element->GetNumberOfNodes();
		IssmDouble* thickness     = xNew<IssmDouble>(numnodes);
		IssmDouble* B             = xNew<IssmDouble>(numnodes);
		IssmDouble* n             = xNew<IssmDouble>(numnodes);
		IssmDouble* eplhead       = xNew<IssmDouble>(numnodes);
		IssmDouble* epl_slopeX    = xNew<IssmDouble>(numnodes);
		IssmDouble* epl_slopeY    = xNew<IssmDouble>(numnodes);
		IssmDouble* old_thickness = xNew<IssmDouble>(numnodes);
		IssmDouble* ice_thickness = xNew<IssmDouble>(numnodes);
		IssmDouble* bed           = xNew<IssmDouble>(numnodes);

		bool       active_element;
		IssmDouble init_thick;
		element->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);

		/* Intermiedaries */
		IssmDouble  A;
		IssmDouble  EPLgrad2;
		IssmDouble  EPL_N;
		IssmDouble  opening,closing;

		if(!active_element){
			init_thick = element->FindParam(HydrologydcEplInitialThicknessEnum);
			/*Keeping thickness to initial value if EPL is not active*/
			for(int i=0;i<numnodes;i++){
				thickness[i]=init_thick;
			}
		}
		else{
			switch(domaintype){
				case Domain2DhorizontalEnum: element->GetInputListOnVertices(&B[0],MaterialsRheologyBbarEnum); break;
				case Domain3DEnum:           element->GetInputListOnVertices(&B[0],MaterialsRheologyBEnum); break;
				default: _error_("not Implemented Yet");
			}
			element->GetInputListOnVertices(&eplhead[0],EplHeadSubstepEnum);
			element->GetInputListOnVertices(&epl_slopeX[0],EplHeadSlopeXEnum);
			element->GetInputListOnVertices(&epl_slopeY[0],EplHeadSlopeYEnum);
			element->GetInputListOnVertices(&old_thickness[0],HydrologydcEplThicknessOldEnum);
			element->GetInputListOnVertices(&ice_thickness[0],ThicknessEnum);
			element->GetInputListOnVertices(&bed[0],BaseEnum);
			element->GetInputListOnVertices(&n[0],MaterialsRheologyNEnum);

			for(int i=0;i<numnodes;i++){
				A=pow(B[i],-n[i]);
				/*Compute first the effective pressure in the EPL*/
				EPL_N=gravity*((rho_ice*ice_thickness[i])-(rho_water*(eplhead[i]-bed[i])));
				if(EPL_N<0.0)EPL_N=0.0;
				/*Get then the square of the gradient of EPL heads*/
				EPLgrad2 = (epl_slopeX[i]*epl_slopeX[i])+(epl_slopeY[i]*epl_slopeY[i]);
				/*And proceed to the real thing*/
				opening=(rho_water*gravity*epl_conductivity*EPLgrad2*dt)/(rho_ice*latentheat);
				closing=(2.0*A*dt*pow(EPL_N,n[i]))/(pow(n[i],n[i]));
				/*implicit*/
				thickness[i] = old_thickness[i]/(1.0-opening+closing);
				/*explicit*/
				//thickness[i] = old_thickness[i]*(1.0+opening-closing);
				/*centered*/
				//thickness[i] = old_thickness[i]*(1.0+opening-closing)/(1.0-opening+closing);
				/*Take care of otherthikening*/
				if(thickness[i]>max_thick){
					thickness[i] = max_thick;
				}
			}
		}
		element->AddInput(HydrologydcEplThicknessSubstepEnum,thickness,element->GetElementType());
		xDelete<IssmDouble>(thickness);
		xDelete<IssmDouble>(eplhead);
		xDelete<IssmDouble>(epl_slopeX);
		xDelete<IssmDouble>(epl_slopeY);
		xDelete<IssmDouble>(old_thickness);
		xDelete<IssmDouble>(ice_thickness);
		xDelete<IssmDouble>(bed);
		xDelete<IssmDouble>(B);
		xDelete<IssmDouble>(n);
	}
}/*}}}*/
void  HydrologyDCEfficientAnalysis::HydrologyEPLGetMask(Vector<IssmDouble>* vec_mask, Vector<IssmDouble>* recurence, Element* element){/*{{{*/

	bool        active_element;
	int         domaintype;
	int         current_mask;
	IssmDouble  h_max;
	IssmDouble  sedheadmin;
	Element*    basalelement=NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	/*Intermediaries*/
	int         numnodes      =basalelement->GetNumberOfNodes();
	IssmDouble* epl_thickness =xNew<IssmDouble>(numnodes);
	IssmDouble* old_active    =xNew<IssmDouble>(numnodes);
	IssmDouble* sedhead       =xNew<IssmDouble>(numnodes);
	IssmDouble* eplhead       =xNew<IssmDouble>(numnodes);
	IssmDouble* residual      =xNew<IssmDouble>(numnodes);

	IssmDouble init_thick    =basalelement->FindParam(HydrologydcEplInitialThicknessEnum);
	IssmDouble colapse_thick =basalelement->FindParam(HydrologydcEplColapseThicknessEnum);

	basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);

	basalelement-> GetInputListOnVertices(&sedhead[0],SedimentHeadSubstepEnum);
	basalelement-> GetInputListOnVertices(&old_active[0],HydrologydcMaskEplactiveNodeEnum);
	basalelement-> GetInputListOnVertices(&residual[0],SedimentHeadResidualEnum);
	basalelement-> GetInputListOnVertices(&epl_thickness[0],HydrologydcEplThicknessSubstepEnum);
	basalelement-> GetInputListOnVertices(&eplhead[0],EplHeadSubstepEnum);

	/*Get minimum sediment head of the element*/
	sedheadmin=sedhead[0];
	for(int i=1;i<numnodes;i++) if(sedhead[i]<=sedheadmin)sedheadmin=sedhead[i];
	for(int i=0;i<numnodes;i++){
		current_mask=reCast<int>(old_active[i]);
		/*If mask was already one, keep one or colapse*/
		if(old_active[i]>0.){
			vec_mask->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
			current_mask=1;
			/* If epl thickness gets under colapse thickness, close the layer if there is no residual*/
			if(epl_thickness[i]<colapse_thick && residual[i]<=0.){
				vec_mask->SetValue(basalelement->nodes[i]->Sid(),0.,INS_VAL);
				current_mask=0;
				recurence->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
			}
		}
		else if (old_active[i]==0.){
			/*Activate if we have a residual from sediment*/
			if(residual[i]>0.){
				vec_mask->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
				current_mask=1;
				recurence->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
			}
			else{
				/*If node is now closed bring its thickness back to initial*/
				epl_thickness[i]=init_thick;
				vec_mask->SetValue(basalelement->nodes[i]->Sid(),0.,INS_VAL);
				current_mask=0;
			}
		}
		/*Increase of the efficient system is needed if the epl head reach the maximum value (sediment max value for now) we check that only if the epl is active here*/
		if(current_mask>0){
			GetHydrologyDCInefficientHmax(&h_max,basalelement,basalelement->nodes[i]);
			if(eplhead[i]>=h_max && active_element){
				for(int j=0;j<numnodes;j++){
					/*Increase of the domain is on the downstream node in term of sediment head*/
					if((sedhead[j] == sedheadmin) && (i!=j)){
						vec_mask->SetValue(basalelement->nodes[j]->Sid(),1.,INS_VAL);
						if(old_active[j]==0.){
							recurence->SetValue(basalelement->nodes[j]->Sid(),1.,INS_VAL);
						}
					}
				}
			}
		}
	}
	element->AddBasalInput(HydrologydcEplThicknessSubstepEnum,epl_thickness,basalelement->GetElementType());

	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(epl_thickness);
	xDelete<IssmDouble>(old_active);
	xDelete<IssmDouble>(sedhead);
	xDelete<IssmDouble>(eplhead);
	xDelete<IssmDouble>(residual);
}
/*}}}*/
void HydrologyDCEfficientAnalysis::HydrologyEPLGetActive(Vector<IssmDouble>* active_vec, Element* element){/*{{{*/
	/*Constants*/
	int      domaintype;
	Element*   basalelement=NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	const int   numnodes = basalelement->GetNumberOfNodes();
	IssmDouble  flag     = 0.;
	IssmDouble* active   = xNew<IssmDouble>(numnodes);
	bool active_element;

	/*Pass the activity mask from elements to nodes*/
	basalelement->GetInputListOnVertices(&active[0],HydrologydcMaskEplactiveNodeEnum);
	basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);

	for(int i=0;i<numnodes;i++) flag+=active[i];

	/*If any node is active all the node in the element are active*/
	if(flag>0.){
		for(int i=0;i<numnodes;i++){
			active_vec->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
		}
	}
	/*If the element is active all its nodes are active*/
	else if(active_element){
		for(int i=0;i<numnodes;i++){
			active_vec->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
		}
	}
	else{
		/*Do not do anything: at least one node is active for this element but this element is not solved for*/
	}
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(active);
}
/*}}}*/
