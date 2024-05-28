#include "./HydrologyDCInefficientAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../classes/Node.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/TransientInput.h"

/*Model processing*/
int  HydrologyDCInefficientAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyDCInefficientAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int         hydrology_model;
	int         sedimentlimit_flag;
	int         transfer_flag;
	int         unconfined_flag;
	int         penalty_lock;
	int         hydro_maxiter;
	int         hydroslices;
	int         averaging_method;
	int         numoutputs;
	bool        isefficientlayer;
	bool        sliceadapt;
	IssmDouble  penalty_factor;
	IssmDouble  rel_tol;
	IssmDouble  leakagefactor;
	IssmDouble  sedimentlimit;
	IssmDouble  sed_poro;
	IssmDouble  sed_thick;
	char**      requestedoutputs = NULL;

	/*retrieve some parameters: */
	bool   issmb;
	iomodel->FindConstant(&issmb,"md.transient.issmb");
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	iomodel->FetchData(&sedimentlimit_flag, "md.hydrology.sedimentlimit_flag" );
	iomodel->FetchData(&transfer_flag,      "md.hydrology.transfer_flag" );
	iomodel->FetchData(&unconfined_flag,    "md.hydrology.unconfined_flag" );
	iomodel->FetchData(&penalty_lock,       "md.hydrology.penalty_lock" );
	iomodel->FetchData(&hydro_maxiter,      "md.hydrology.max_iter" );
	iomodel->FetchData(&hydroslices,        "md.hydrology.steps_per_step");
	iomodel->FetchData(&sliceadapt,         "md.hydrology.step_adapt");
	iomodel->FetchData(&averaging_method,   "md.hydrology.averaging");
	iomodel->FetchData(&isefficientlayer,   "md.hydrology.isefficientlayer");
	iomodel->FetchData(&penalty_factor,     "md.hydrology.penalty_factor" );
	iomodel->FetchData(&rel_tol,            "md.hydrology.rel_tol" );

	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
	parameters->AddObject(new IntParam(HydrologydcSedimentlimitFlagEnum,sedimentlimit_flag));
	parameters->AddObject(new IntParam(HydrologydcTransferFlagEnum,transfer_flag));
	parameters->AddObject(new IntParam(HydrologydcUnconfinedFlagEnum,unconfined_flag));
	parameters->AddObject(new IntParam(HydrologydcPenaltyLockEnum,penalty_lock));
	parameters->AddObject(new IntParam(HydrologydcMaxIterEnum,hydro_maxiter));
	parameters->AddObject(new IntParam(HydrologyStepsPerStepEnum,hydroslices));
	parameters->AddObject(new BoolParam(HydrologyStepAdaptEnum,sliceadapt));
	parameters->AddObject(new IntParam(HydrologyAveragingEnum,averaging_method));
	parameters->AddObject(new BoolParam(HydrologydcIsefficientlayerEnum,isefficientlayer));
	parameters->AddObject(new DoubleParam(HydrologydcPenaltyFactorEnum,penalty_factor));
	parameters->AddObject(new DoubleParam(HydrologydcRelTolEnum,rel_tol));

	iomodel->FetchData(&sed_poro,  "md.hydrology.sediment_porosity" );
	iomodel->FetchData(&sed_thick, "md.hydrology.sediment_thickness" );

	parameters->AddObject(new DoubleParam(HydrologydcSedimentPorosityEnum,sed_poro));
	parameters->AddObject(new DoubleParam(HydrologydcSedimentThicknessEnum,sed_thick));

	if(transfer_flag==1){
		iomodel->FetchData(&leakagefactor,"md.hydrology.leakage_factor");
		parameters->AddObject(new DoubleParam(HydrologydcLeakageFactorEnum,leakagefactor));
	}
	if(sedimentlimit_flag==1){
		iomodel->FetchData(&sedimentlimit,"md.hydrology.sedimentlimit");
		parameters->AddObject(new DoubleParam(HydrologydcSedimentlimitEnum,sedimentlimit));
	}
	if(!issmb){
		parameters->AddObject(iomodel->CopyConstantObject("md.smb.model",SmbEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.smb.averaging",SmbAveragingEnum));
	}

  /*Requested outputs*/
  iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
  parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
  if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
  iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");
}/*}}}*/
void HydrologyDCInefficientAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	bool   isefficientlayer;
	int    hydrology_model;

	/*Fetch data needed: */
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	/*Fetch data needed: */
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");

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
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.basal_moulin_input",HydrologydcBasalMoulinInputEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sediment_head",SedimentHeadSubstepEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sediment_transmitivity",HydrologydcSedimentTransmitivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.mask_thawed_node",HydrologydcMaskThawedNodeEnum);


	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	if(isefficientlayer){
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.mask_eplactive_node",HydrologydcMaskEplactiveNodeEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.epl_head",EplHeadSubstepEnum);
	}

}/*}}}*/
void HydrologyDCInefficientAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	}
	::CreateNodes(nodes,iomodel,HydrologyDCInefficientAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
void HydrologyDCInefficientAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spcsediment_head",HydrologyDCInefficientAnalysisEnum,P1Enum);
}/*}}}*/
void HydrologyDCInefficientAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(1,"md.mesh.vertexonbase");
	}
	//create penalties for nodes: no node can have water above the max
	CreateSingleNodeToElementConnectivity(iomodel);
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Pengrid(i+1,i,iomodel));
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Pengrid(i+1,i,iomodel));
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyDCInefficientAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyDCInefficientAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyDCInefficientAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyDCInefficientAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyDCInefficientAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	bool     thawed_element;
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

	basalelement->GetInputValue(&thawed_element,HydrologydcMaskThawedEltEnum);

	/*Check that all nodes are active, else return empty matrix*/
	if(!thawed_element) {
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		return NULL;
	}

	/*Intermediaries */
	bool        active_element,isefficientlayer;
	IssmDouble  D_scalar,Jdet,dt;
	IssmDouble  sediment_transmitivity;
	IssmDouble  transfer,sediment_storing;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement ->GetVerticesCoordinates(&xyz_list);
	basalelement ->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement ->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);
	Input* SedTrans_input = basalelement->GetInput(HydrologydcSedimentTransmitivityEnum); _assert_(SedTrans_input);
	Input* sed_head_input = basalelement->GetInput(SedimentHeadSubstepEnum);
	Input* base_input     = basalelement->GetInput(BaseEnum);

	/*Transfer related Inputs*/
	if(isefficientlayer){
		basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);
	}
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		sediment_transmitivity = SedimentTransmitivity(basalelement,gauss,sed_head_input,base_input,SedTrans_input);
		sediment_storing       = SedimentStoring(basalelement,gauss,sed_head_input,base_input);

		/*Diffusivity*/
		D_scalar=sediment_transmitivity*gauss->weight*Jdet;
		if(dt!=0.) D_scalar=D_scalar*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D_scalar*(dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]);
			}
		}
		/*Transient*/
		if(dt!=0.){
			D_scalar=sediment_storing*gauss->weight*Jdet;
			for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[j]*basis[i];
			/*Transfer EPL part*/
			if(isefficientlayer){
				if(active_element){
					transfer=GetHydrologyKMatrixTransfer(basalelement);
					D_scalar=dt*transfer*gauss->weight*Jdet;
					for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[j]*basis[i];
				}
			}
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
ElementVector* HydrologyDCInefficientAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	bool		thawed_element;
	int		domaintype;
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

	basalelement->GetInputValue(&thawed_element,HydrologydcMaskThawedEltEnum);

	/*Check that all nodes are active, else return empty matrix*/
	if(!thawed_element) {
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		return NULL;
	}

	/*Intermediaries */
	bool       active_element,isefficientlayer;
	int        smb_model,smb_averaging;
	int        smbsubstepping, hydrologysubstepping;
	IssmDouble dt,scalar,water_head;
	IssmDouble sediment_storing,sediment_transmitivity;
	IssmDouble water_load,runoff_value,transfer;
	IssmDouble Jdet,time;
	IssmDouble active_node;

	IssmDouble *xyz_list             = NULL;
	Input     *active_element_input = NULL;
	Input     *old_wh_input         = NULL;
	Input     *dummy_input          = NULL;
	Input     *surface_runoff_input = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);
	basalelement->FindParam(&smb_model,SmbEnum);

	Input*	sed_head_input   = basalelement->GetInput(SedimentHeadSubstepEnum);
	Input*	epl_head_input   = basalelement->GetInput(EplHeadSubstepEnum);
	Input*	base_input       = basalelement->GetInput(BaseEnum);
	Input*	basal_melt_input = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(basal_melt_input);
	Input*	SedTrans_input   = basalelement->GetInput(HydrologydcSedimentTransmitivityEnum); _assert_(SedTrans_input);

	if(dt!= 0.){
		old_wh_input = basalelement->GetInput(SedimentHeadOldEnum); _assert_(old_wh_input);
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

	/*Transfer related Inputs*/
	if(isefficientlayer){
		basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		sediment_transmitivity = SedimentTransmitivity(basalelement,gauss,sed_head_input,base_input,SedTrans_input);

		/*Loading term*/
		if(!isefficientlayer){
			basal_melt_input->GetInputValue(&water_load,gauss);
			if(surface_runoff_input) surface_runoff_input->GetInputValue(&runoff_value,gauss);
			else                     runoff_value = 0.;
			scalar = Jdet*gauss->weight*(water_load+runoff_value);
			if(dt!=0.) scalar = scalar*dt;
			for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
		}
		else{
			/*if EPL is present and active input is there not here*/
			if(!active_element){
				basal_melt_input->GetInputValue(&water_load,gauss);
				if(surface_runoff_input)surface_runoff_input->GetInputValue(&runoff_value,gauss);
				else runoff_value = 0.;
				scalar = Jdet*gauss->weight*(water_load+runoff_value);
				if(dt!=0.) scalar = scalar*dt;
				for(int i=0;i<numnodes;i++){
					//This is the original
					pe->values[i]+=scalar*basis[i];
					//This is the noded version
					/* basalelement->GetInputValue(&active_node,basalelement->nodes[i],HydrologydcMaskEplactiveNodeEnum); */
					/* if(!reCast<bool>(active_node)){ */
					/* 	pe->values[i]+=scalar*basis[i]; */
					/* 	//if(basalelement->nodes[i]->Sid()==42)_printf_("IDS uni Input "<<scalar*basis[i]<<"\n"); */
					//}
				}
			}
		}

		/*Transient and transfer terms*/
		if(dt!=0.){
			old_wh_input->GetInputValue(&water_head,gauss);
			sediment_storing = SedimentStoring(basalelement,gauss,sed_head_input,base_input);
			if(isefficientlayer){
				/*Dealing with the sediment part of the transfer term*/
				if(active_element){
					transfer=GetHydrologyPVectorTransfer(basalelement,gauss,epl_head_input);
				}
				else{
					transfer=0.0;
				}
				scalar = Jdet*gauss->weight*((water_head*sediment_storing)+(dt*transfer));
				for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
			}
			else{
				scalar = Jdet*gauss->weight*(water_head*sediment_storing);
				for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
			}
		}
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void HydrologyDCInefficientAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,SedimentHeadSubstepEnum);
}/*}}}*/
void HydrologyDCInefficientAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void HydrologyDCInefficientAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

 	/*Intermediaries*/
	bool	 converged;
	int*     doflist = NULL;
	int	 domaintype;
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
	/*Fetch number of nodes for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	basalelement->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values   = xNew<IssmDouble>(numnodes);
	IssmDouble* pressure = xNew<IssmDouble>(numnodes);
	IssmDouble* residual = xNew<IssmDouble>(numnodes);

	for(int i=0;i<numnodes;i++){
	  	values[i] =solution[doflist[i]];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*If converged keep the residual in mind, also compute effective pressure*/
	basalelement->GetInputValue(&converged,ConvergedEnum);

	/*Get inputs*/
	if(converged){
		IssmDouble  penalty_factor,kmax,kappa,h_max;
		IssmDouble* thickness = xNew<IssmDouble>(numnodes);
		IssmDouble* base = xNew<IssmDouble>(numnodes);
		IssmDouble* transmitivity = xNew<IssmDouble>(numnodes);

		basalelement->FindParam(&kmax,HydrologySedimentKmaxEnum);
		basalelement->FindParam(&penalty_factor,HydrologydcPenaltyFactorEnum);
		IssmDouble rho_freshwater = basalelement->FindParam(MaterialsRhoFreshwaterEnum);
		IssmDouble rho_ice        = basalelement->FindParam(MaterialsRhoIceEnum);
		IssmDouble g              = basalelement->FindParam(ConstantsGEnum);

		basalelement->GetInputListOnVertices(&thickness[0],ThicknessEnum);
		basalelement->GetInputListOnVertices(&base[0],BaseEnum);
		basalelement->GetInputListOnVertices(&transmitivity[0],HydrologydcSedimentTransmitivityEnum);

		kappa=kmax*pow(10.,penalty_factor);

		for(int i=0;i<numnodes;i++){
			GetHydrologyDCInefficientHmax(&h_max,basalelement,basalelement->GetNode(i));
			if(values[i]>h_max) {
				residual[i] = kappa*(values[i]-h_max);
				//residual[i] = kappa*(values[i]-h_max)*transmitivity[i];
			}
			else{
				residual[i] = 0.;
			}
			pressure[i]=(rho_ice*g*thickness[i])-(rho_freshwater*g*(values[i]-base[i]));
		}
		xDelete<IssmDouble>(thickness);
		xDelete<IssmDouble>(base);
		xDelete<IssmDouble>(transmitivity);
	}

	/*Add input to the element: */
	element->AddBasalInput(SedimentHeadSubstepEnum,values,P1Enum);
	element->AddBasalInput(EffectivePressureSubstepEnum,pressure,P1Enum);
	element->AddBasalInput(SedimentHeadResidualEnum,residual,P1Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(residual);
	xDelete<IssmDouble>(pressure);
	xDelete<int>(doflist);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void HydrologyDCInefficientAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
/*Intermediaries*/
IssmDouble HydrologyDCInefficientAnalysis::SedimentStoring(Element* element,Gauss* gauss,Input* sed_head_input, Input* base_input){/*{{{*/
	int unconf_scheme;
	IssmDouble expfac;
	IssmDouble sediment_storing;
	IssmDouble storing,yield;
	IssmDouble base_elev,prestep_head,water_sheet;
	IssmDouble porewater_mass           = element->FindParam(HydrologydcSedimentPoreWaterMassEnum);
	IssmDouble layer_compressibility    = element->FindParam(HydrologydcSedimentLayerCompressibilityEnum);
	IssmDouble sediment_thickness       = element->FindParam(HydrologydcSedimentThicknessEnum);
	element->FindParam(&unconf_scheme,HydrologydcUnconfinedFlagEnum);
	switch(unconf_scheme){
	case 0:
		sediment_storing=porewater_mass*sediment_thickness*layer_compressibility;
		break;
	case 1:
		yield = element->FindParam(HydrologydcSedimentPorosityEnum);
		base_input->GetInputValue(&base_elev,gauss);
		sed_head_input->GetInputValue(&prestep_head,gauss);

		water_sheet=max(0.0,(prestep_head-(base_elev-sediment_thickness)));
		storing=porewater_mass*sediment_thickness*layer_compressibility;
		//using logistic function for heavyside approximation
		expfac=10.;
		sediment_storing=yield+(storing-yield)/(1+exp(-2*expfac*(water_sheet-0.99*sediment_thickness)));
		break;
	default:
		_error_("UnconfinedFlag is 0 or 1");
	}
	return sediment_storing;
}/*}}}*/
IssmDouble HydrologyDCInefficientAnalysis::SedimentTransmitivity(Element* element,Gauss* gauss,Input* sed_head_input, Input* base_input,Input* SedTrans_input){/*{{{*/
	int unconf_scheme;
	IssmDouble sediment_transmitivity;
	IssmDouble FullLayer_transmitivity;
	IssmDouble base_elev,prestep_head,water_sheet;
	IssmDouble sediment_thickness       = element->FindParam(HydrologydcSedimentThicknessEnum);

	element->FindParam(&unconf_scheme,HydrologydcUnconfinedFlagEnum);
	SedTrans_input->GetInputValue(&FullLayer_transmitivity,gauss);

	switch(unconf_scheme){
	case 0:
		sediment_transmitivity=FullLayer_transmitivity;
		break;
	case 1:
		base_input->GetInputValue(&base_elev,gauss);
		sed_head_input->GetInputValue(&prestep_head,gauss);
		water_sheet=max(0.0,(prestep_head-(base_elev-sediment_thickness)));

		//min definition of the if test
		sediment_transmitivity=FullLayer_transmitivity*min(water_sheet/sediment_thickness,1.);
		if (sediment_transmitivity==0){
			sediment_transmitivity=1.0e-20;
		}

		break;
	default:
		_error_("UnconfinedFlag is 0 or 1");
	}
	return sediment_transmitivity;
}/*}}}*/
void  HydrologyDCInefficientAnalysis::GetHydrologyDCInefficientHmax(IssmDouble* ph_max,Element* element, Node* innode){/*{{{*/

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
		element->GetInputValue(&thickness,innode,ThicknessEnum);
		element->GetInputValue(&bed,innode,BaseEnum);
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
IssmDouble HydrologyDCInefficientAnalysis::GetHydrologyKMatrixTransfer(Element* element){/*{{{*/

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
IssmDouble HydrologyDCInefficientAnalysis::GetHydrologyPVectorTransfer(Element* element, Gauss* gauss, Input* epl_head_input){/*{{{*/

	int transfermethod;
	IssmDouble epl_head;
	IssmDouble leakage,transfer;
	element->FindParam(&transfermethod,HydrologydcTransferFlagEnum);

	/*Switch between the different transfer methods cases*/
	switch(transfermethod){
	case 0:
		/*Just keepping the transfer to zero*/
		transfer=0.0;
		break;
	case 1:
		_assert_(epl_head_input);
		epl_head_input->GetInputValue(&epl_head,gauss);
		element->FindParam(&leakage,HydrologydcLeakageFactorEnum);
		transfer=+epl_head*leakage;
		break;
	default:
		_error_("no case higher than 1 for the Transfer method");
	}
	return transfer;
}/*}}}*/
void HydrologyDCInefficientAnalysis::ElementizeEplMask(FemModel* femmodel){/*{{{*/

	bool     element_active;
	Element* element=NULL;
	int      elementssize=femmodel->elements->Size();
	for(Object* & object : femmodel->elements->objects){
		element = xDynamicCast<Element*>(object);
		Input* input=element->GetInput(HydrologydcMaskEplactiveNodeEnum); _assert_(input);
		if(input->GetInputMax()>0.){
			element_active = true;
		}
		else{
			element_active = false;
		}
		element->SetBoolInput(element->inputs,HydrologydcMaskEplactiveEltEnum,element_active);
	}
}/*}}}*/
void  HydrologyDCInefficientAnalysis::HydrologyIDSGetMask(Vector<IssmDouble>* vec_mask, Element* element){/*{{{*/
	bool        active_element;
	int         domaintype;
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
	int				numnodes    =	basalelement->GetNumberOfNodes();
	IssmDouble*		meltingrate =	xNew<IssmDouble>(numnodes);
	IssmDouble*		groundedice =	xNew<IssmDouble>(numnodes);

	basalelement->GetInputListOnVertices(&meltingrate[0],BasalforcingsGroundediceMeltingRateEnum);
	basalelement->GetInputListOnVertices(&groundedice[0],MaskOceanLevelsetEnum);

	/*if melting rate is not positive and node is not floating, deactivate*/
	for(int i=0;i<numnodes;i++){
		if ((meltingrate[i]<=0.0) && (groundedice[i]>0)){
			vec_mask->SetValue(basalelement->nodes[i]->Sid(),0.,INS_VAL);
		}
		else{
			vec_mask->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
		}
	}

	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(meltingrate);
	xDelete<IssmDouble>(groundedice);
}/*}}}*/
void HydrologyDCInefficientAnalysis::ElementizeIdsMask(FemModel* femmodel){/*{{{*/

	bool     element_active;
	Element* element=NULL;
	int elementssize = femmodel->elements->Size();
	for(Object* & object : femmodel->elements->objects){
		element = xDynamicCast<Element*>(object);

		Input* input=element->GetInput(HydrologydcMaskThawedNodeEnum); _assert_(input);
		if(input->GetInputMax()>0.){
			element_active = true;
		}
		else{
			element_active = false;
		}
		element->SetBoolInput(element->inputs,HydrologydcMaskThawedEltEnum,element_active);
	}
}/*}}}*/
void HydrologyDCInefficientAnalysis::HydrologyIdsGetActive(Vector<IssmDouble>* active_vec, Element* element){/*{{{*/
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
	basalelement->GetInputListOnVertices(&active[0],HydrologydcMaskThawedNodeEnum);
	basalelement->GetInputValue(&active_element,HydrologydcMaskThawedEltEnum);

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
