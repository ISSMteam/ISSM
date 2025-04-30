#include "./HydrologyShaktiAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Define 2 hardcoded parameters*/
#define OMEGA 0.001    // parameter controlling transition to nonlinear resistance in basal system (dimensionless)
#define NU    1.787e-6 //kinematic water viscosity m^2/s
#define CT    7.5e-8  // Clapeyron slope (K/Pa) 
#define CW    4.22e3   // specific heat capacity of water (J/kg/K)

/*Model processing*/
void HydrologyShaktiAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologyshaktiEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spchead",HydrologyShaktiAnalysisEnum,P1Enum);

}/*}}}*/
void HydrologyShaktiAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Shakti?*/
	if(hydrology_model!=HydrologyshaktiEnum) return;

	/*Create discrete loads for Moulins*/
	CreateSingleNodeToElementConnectivity(iomodel);
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(1,"md.mesh.vertexonbase");
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

	/*Deal with Neumann BC*/
	int M,N;
	int *segments = NULL;
	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(&segments,&M,&N,"md.mesh.segments2d");
	}
	else if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchData(&segments,&M,&N,"md.mesh.segments");
	}
	else{
		_error_("mesh type not supported yet");
	}

	/*Check that the size seem right*/
	_assert_(N==3); _assert_(M>=3);

	for(int i=0;i<M;i++){
		if(iomodel->my_elements[segments[i*3+2]-1]){
			loads->AddObject(new Neumannflux(i+1,i,iomodel,segments));
		}
	}
	xDelete<int>(segments);

}/*}}}*/
void HydrologyShaktiAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Shakti?*/
	if(hydrology_model!=HydrologyshaktiEnum) return;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,HydrologyShaktiAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  HydrologyShaktiAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyShaktiAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model,frictionlaw;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Shakti?*/
	if(hydrology_model!=HydrologyshaktiEnum) return;

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
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.geothermalflux",BasalforcingsGeothermalfluxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.head",HydrologyHeadEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.gap_height",HydrologyGapHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.englacial_input",HydrologyEnglacialInputEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.moulin_input",HydrologyMoulinInputEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_spacing",HydrologyBumpSpacingEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.reynolds",HydrologyReynoldsEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.neumannflux",HydrologyNeumannfluxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.storage",HydrologyStorageEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxBaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyBaseEnum);
	}

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);
}/*}}}*/
void HydrologyShaktiAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Shakti?*/
	if(hydrology_model!=HydrologyshaktiEnum) return;

	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
   parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.relaxation",HydrologyRelaxationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.gap_height_min",HydrologyGapHeightMinEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.gap_height_max",HydrologyGapHeightMaxEnum));

  /*Requested outputs*/
  iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
  parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
  if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
  iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");

	/*Friction*/
	FrictionUpdateParameters(parameters, iomodel);
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyShaktiAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyShaktiAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyShaktiAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyShaktiAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyShaktiAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  gap,bed,thickness,head,g,rho_ice,rho_water,A,B,n,storage;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);

	/*Get conductivity from inputs*/
	IssmDouble conductivity = GetConductivity(basalelement);

	/*Get Params*/
	IssmDouble dt;
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);

	/*Get all inputs and parameters*/
	basalelement->FindParam(&rho_water,MaterialsRhoFreshwaterEnum);
	basalelement->FindParam(&rho_ice,MaterialsRhoIceEnum);
	basalelement->FindParam(&g,ConstantsGEnum);
	Input* B_input = basalelement->GetInput(MaterialsRheologyBEnum);      _assert_(B_input);
	Input* n_input = basalelement->GetInput(MaterialsRheologyNEnum);      _assert_(n_input);
	Input* gap_input = basalelement->GetInput(HydrologyGapHeightEnum);    _assert_(gap_input);
	Input* thickness_input = basalelement->GetInput(ThicknessEnum);       _assert_(thickness_input);
	Input* head_input = basalelement->GetInput(HydrologyHeadEnum);        _assert_(head_input);
	Input* base_input = basalelement->GetInput(BaseEnum);                 _assert_(base_input);
	Input* storage_input = basalelement->GetInput(HydrologyStorageEnum);  _assert_(storage_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(1);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		base_input->GetInputValue(&bed,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		gap_input->GetInputValue(&gap,gauss);
		head_input->GetInputValue(&head,gauss);
		storage_input->GetInputValue(&storage, gauss);

		/*Get ice A parameter*/
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		A=pow(B,-n);

		/*Get water and ice pressures*/
		IssmDouble pressure_ice   = rho_ice*g*thickness;    _assert_(pressure_ice>0.);
		IssmDouble pressure_water = rho_water*g*(head-bed);
		if(pressure_water>pressure_ice) pressure_water = pressure_ice;

		IssmDouble factor = conductivity*gauss->weight*Jdet;
		IssmDouble factor2 = gauss->weight*Jdet*storage/dt + gauss->weight*Jdet*A*(n)*(pow(fabs(pressure_ice-pressure_water),(n-1))*rho_water*g)*gap;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += factor*(dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j])
				  + factor2*basis[i]*basis[j];
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
ElementVector* HydrologyShaktiAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble  Jdet,meltrate,G,dh[2],B,A,n;
	IssmDouble  gap,bed,thickness,head,ieb,head_old,storage;
	IssmDouble  lr,br,vx,vy,beta,lc;
	IssmDouble  alpha2,frictionheat;
   IssmDouble  PMPheat,dissipation,dpressure_water[2],dbed[2];	
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	IssmDouble  latentheat      = basalelement->FindParam(MaterialsLatentheatEnum);
	IssmDouble  g               = basalelement->FindParam(ConstantsGEnum);
	IssmDouble  rho_ice         = basalelement->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water       = basalelement->FindParam(MaterialsRhoFreshwaterEnum);
	Input* geothermalflux_input = basalelement->GetInput(BasalforcingsGeothermalfluxEnum);_assert_(geothermalflux_input);
	Input* head_input           = basalelement->GetInput(HydrologyHeadEnum);              _assert_(head_input);
	Input* gap_input            = basalelement->GetInput(HydrologyGapHeightEnum);         _assert_(gap_input);
	Input* thickness_input      = basalelement->GetInput(ThicknessEnum);                  _assert_(thickness_input);
	Input* base_input           = basalelement->GetInput(BaseEnum);                       _assert_(base_input);
	Input* B_input              = basalelement->GetInput(MaterialsRheologyBEnum);         _assert_(B_input);
	Input* n_input              = basalelement->GetInput(MaterialsRheologyNEnum);         _assert_(n_input);
	Input* englacial_input      = basalelement->GetInput(HydrologyEnglacialInputEnum);    _assert_(englacial_input);
	Input* lr_input             = basalelement->GetInput(HydrologyBumpSpacingEnum);       _assert_(lr_input);
	Input* br_input             = basalelement->GetInput(HydrologyBumpHeightEnum);        _assert_(br_input);
   Input* headold_input        = basalelement->GetInput(HydrologyHeadOldEnum);           _assert_(headold_input);
	Input* storage_input        = basalelement->GetInput(HydrologyStorageEnum);           _assert_(storage_input);

	/*Get conductivity from inputs*/
	IssmDouble conductivity = GetConductivity(basalelement);

	/*Get Params*/
	IssmDouble dt;
   basalelement->FindParam(&dt,TimesteppingTimeStepEnum);

	/*Build friction basalelement, needed later: */
	Friction* friction=new Friction(basalelement,2);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		geothermalflux_input->GetInputValue(&G,gauss);
		base_input->GetInputValue(&bed,gauss);
		base_input->GetInputDerivativeValue(&dbed[0],xyz_list,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		gap_input->GetInputValue(&gap,gauss);
		head_input->GetInputValue(&head,gauss);
		head_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
		englacial_input->GetInputValue(&ieb,gauss);
		lr_input->GetInputValue(&lr,gauss);
		br_input->GetInputValue(&br,gauss);
		headold_input->GetInputValue(&head_old,gauss);
		storage_input->GetInputValue(&storage,gauss);

		/*Get ice A parameter*/
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		A=pow(B,-n);

		/*Compute beta term*/
		if(gap<br)
		 beta = (br-gap)/lr;
		else
		 beta = 0.;

		/*Compute frictional heat flux*/
		friction->GetAlpha2(&alpha2,gauss);
		friction->GetBasalSlidingSpeeds(&vx, &vy, gauss);
		frictionheat=alpha2*(vx*vx+vy*vy);

		/*Get water and ice pressures*/
		IssmDouble pressure_ice   = rho_ice*g*thickness;    _assert_(pressure_ice>0.); 
		IssmDouble pressure_water = rho_water*g*(head-bed);
		if(pressure_water>pressure_ice) pressure_water = pressure_ice;

		/*Get water pressure from previous time step to use in lagged creep term*/
		IssmDouble pressure_water_old = rho_water*g*(head_old-bed);
		if(pressure_water_old>pressure_ice) pressure_water_old = pressure_ice;

		/*Compute change in sensible heat due to changes in pressure melting point*/
		dpressure_water[0] = rho_water*g*(dh[0] - dbed[0]);
		dpressure_water[1] = rho_water*g*(dh[1] - dbed[1]);

		meltrate = 1/latentheat*(G+frictionheat+rho_water*g*conductivity*(dh[0]*dh[0]+dh[1]*dh[1]));

		IssmDouble factor = Jdet*gauss->weight*
		 (meltrate*(1/rho_water-1/rho_ice)
		  +A*pow(fabs(pressure_ice - pressure_water),n-1)*(pressure_ice + rho_water*g*bed)*gap
		  +(n-1)*A*pow(fabs(pressure_ice - pressure_water),n-1)*(rho_water*g*head)*gap
		  -beta*sqrt(vx*vx+vy*vy)
		  +ieb
		  +storage*head_old/dt);
		for(int i=0;i<numnodes;i++) pe->values[i]+=factor*basis[i];

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete friction;
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void           HydrologyShaktiAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,HydrologyHeadEnum);
}/*}}}*/
void           HydrologyShaktiAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyShaktiAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Only update if on base*/
	if(!element->IsOnBase()) return;

	/*Intermediary*/
	IssmDouble dh[3];
	int* doflist = NULL;
	IssmDouble* xyz_list = NULL;

	/*Get gravity from parameters*/
	IssmDouble  g = element->FindParam(ConstantsGEnum);

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numnodes);

	/*Get thickness and base on nodes to apply cap on water head*/
	IssmDouble* thickness = xNew<IssmDouble>(numnodes);
	IssmDouble* bed       = xNew<IssmDouble>(numnodes);
	IssmDouble  rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	element->GetInputListOnNodes(&thickness[0],ThicknessEnum);
	element->GetInputListOnNodes(&bed[0],BaseEnum);

	/*Get head from previous time-step and under-relaxation coefficient to use in under-relaxation for nonlinear convergence*/
   IssmDouble* head_old  = xNew<IssmDouble>(numnodes); 
	element->GetInputListOnNodes(&head_old[0],HydrologyHeadEnum);
   IssmDouble relaxation; 
	element->FindParam(&relaxation,HydrologyRelaxationEnum);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];

		/*Under-relaxation*/
	   values[i] = head_old[i] - relaxation*(head_old[i]-values[i]);

		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Add input to the element: */
	element->AddBasalInput(HydrologyHeadEnum,values,element->GetElementType());

	/*Update reynolds number according to new solution*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* head_input = element->GetInput(HydrologyHeadEnum);_assert_(head_input);
	IssmDouble conductivity = GetConductivity(element);

	/*Get gap height derivatives at the center of the element*/
	Gauss* gauss=element->NewGauss(1);
	head_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
	delete gauss;

	IssmDouble reynolds = conductivity*sqrt(dh[0]*dh[0]+dh[1]*dh[1])/NU;
	element->AddBasalInput(HydrologyReynoldsEnum,&reynolds,P0Enum);

   /*Compute new effective pressure*/
   this->UpdateEffectivePressure(element);

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(bed);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
   xDelete<IssmDouble>(head_old);
}/*}}}*/
void           HydrologyShaktiAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Update active elements based on ice levelset and ocean levelset*/
	GetMaskOfIceVerticesLSMx(femmodel,true);
	SetActiveNodesLSMx(femmodel,true);

	IssmDouble rho_ice   = femmodel->parameters->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = femmodel->parameters->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = femmodel->parameters->FindParam(ConstantsGEnum);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){

		/*Get current element and return if not on base*/
		Element *element  = xDynamicCast<Element*>(object);
		if(!element->IsOnBase()) continue;

		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *bed       = xNew<IssmDouble>(numnodes);
		IssmDouble *thickness = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&bed[0],BaseEnum);
		element->GetInputListOnNodes(&thickness[0],ThicknessEnum);
		element->GetInputListOnNodes(&ls_active[0],HydrologyMaskNodeActivationEnum);

		//for(int in=0;in<numnodes;in++){ //
		for(int in=0;in<3;in++){ //
			Node* node=element->GetNode(in);
			if(mask[in]>0. && ls_active[in]==1.){
				node->Activate(); //Not sure if we need this!
			}
			else{
				IssmDouble phi =  rho_ice*g*thickness[in] + rho_water*g*bed[in]; //FIXME this is correct!
				node->Deactivate();// Not sure if we need this
				node->ApplyConstraint(0,phi);
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(bed);
		xDelete<IssmDouble>(thickness);
		xDelete<IssmDouble>(ls_active);
	}

	return;
}/*}}}*/

/*Additional methods*/
IssmDouble HydrologyShaktiAnalysis::GetConductivity(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble gap,reynolds;

	/*Get gravity from parameters*/
	IssmDouble  g = element->FindParam(ConstantsGEnum);

	/*Get Reynolds and gap average values*/
	Input* reynolds_input = element->GetInput(HydrologyReynoldsEnum);  _assert_(reynolds_input);
	Input* gap_input      = element->GetInput(HydrologyGapHeightEnum); _assert_(gap_input);
	reynolds_input->GetInputAverage(&reynolds);
	gap_input->GetInputAverage(&gap);

	/*Compute conductivity*/
	IssmDouble conductivity = pow(gap,3)*g/(12.*NU*(1+OMEGA*reynolds));
	_assert_(conductivity>0);

	/*Clean up and return*/
	return conductivity;
}/*}}}*/
void HydrologyShaktiAnalysis::UpdateGapHeight(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		UpdateGapHeight(element);
	}

}/*}}}*/
void HydrologyShaktiAnalysis::UpdateGapHeight(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
	if(!element->IsOnBase()) return;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble  newgap = 0.;
	IssmDouble  Jdet,meltrate,G,dh[3],B,A,n,dt;
	IssmDouble  gap,bed,thickness,head;
	IssmDouble  lr,br,vx,vy,beta,lc;
	IssmDouble  alpha2,frictionheat;
	IssmDouble* xyz_list = NULL;
	IssmDouble  dpressure_water[3],dbed[3],PMPheat,dissipation;
	IssmDouble q = 0.;
	IssmDouble channelization = 0.;

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  latentheat      = basalelement->FindParam(MaterialsLatentheatEnum);
	IssmDouble  g               = basalelement->FindParam(ConstantsGEnum);
	IssmDouble  rho_ice         = basalelement->FindParam(MaterialsRhoIceEnum);
	IssmDouble  rho_water       = basalelement->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble  gap_height_min  = basalelement->FindParam(HydrologyGapHeightMinEnum);
	IssmDouble  gap_height_max  = basalelement->FindParam(HydrologyGapHeightMaxEnum);
	Input* geothermalflux_input = basalelement->GetInput(BasalforcingsGeothermalfluxEnum);_assert_(geothermalflux_input);
	Input* head_input           = basalelement->GetInput(HydrologyHeadEnum);              _assert_(head_input);
	Input* gap_input            = basalelement->GetInput(HydrologyGapHeightEnum);         _assert_(gap_input);
	Input* thickness_input      = basalelement->GetInput(ThicknessEnum);                  _assert_(thickness_input);
	Input* base_input           = basalelement->GetInput(BaseEnum);                       _assert_(base_input);
	Input* B_input              = basalelement->GetInput(MaterialsRheologyBEnum);         _assert_(B_input);
	Input* n_input              = basalelement->GetInput(MaterialsRheologyNEnum);         _assert_(n_input);
	Input* lr_input             = basalelement->GetInput(HydrologyBumpSpacingEnum);       _assert_(lr_input);
	Input* br_input             = basalelement->GetInput(HydrologyBumpHeightEnum);        _assert_(br_input);

	/*Get conductivity from inputs*/
	IssmDouble conductivity = GetConductivity(basalelement);

	/*Build friction basalelement, needed later: */
	Friction* friction=new Friction(basalelement,2);

	/*Keep track of weights*/
	IssmDouble totalweights=0.;

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		geothermalflux_input->GetInputValue(&G,gauss);
		base_input->GetInputValue(&bed,gauss);
		base_input->GetInputDerivativeValue(&dbed[0],xyz_list,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		gap_input->GetInputValue(&gap,gauss);
		head_input->GetInputValue(&head,gauss);
		head_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
		lr_input->GetInputValue(&lr,gauss);
		br_input->GetInputValue(&br,gauss);

		/*Get ice A parameter*/
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		A=pow(B,-n);

		/*Compute beta term*/
		if(gap<br)
		 beta = (br-gap)/lr;
		else
		 beta = 0.;

		/*Compute frictional heat flux*/
		friction->GetAlpha2(&alpha2,gauss);
		friction->GetBasalSlidingSpeeds(&vx, &vy, gauss);
		frictionheat=alpha2*(vx*vx+vy*vy);

		/*Get water and ice pressures*/
		IssmDouble pressure_ice   = rho_ice*g*thickness;    _assert_(pressure_ice>0.); 
		IssmDouble pressure_water = rho_water*g*(head-bed);
		if(pressure_water>pressure_ice) pressure_water = pressure_ice;

		/* Compute change in sensible heat due to changes in pressure melting point*/
		dpressure_water[0] = rho_water*g*(dh[0] - dbed[0]);
		dpressure_water[1] = rho_water*g*(dh[1] - dbed[1]);
		dissipation=rho_water*g*conductivity*(dh[0]*dh[0]+dh[1]*dh[1]);

		meltrate = 1/latentheat*(G+frictionheat+rho_water*g*conductivity*(dh[0]*dh[0]+dh[1]*dh[1]));

		element->AddBasalInput(HydrologyMeltRateEnum,&meltrate,P0Enum);
		element->AddBasalInput(HydrologyFrictionHeatEnum,&frictionheat,P0Enum);
		element->AddBasalInput(HydrologyDissipationEnum,&dissipation,P0Enum);
		element->AddBasalInput(HydrologyPmpHeatEnum,&PMPheat,P0Enum);

		newgap += gauss->weight*Jdet*(gap+dt*(
						meltrate/rho_ice
						-A*pow(fabs(pressure_ice-pressure_water),n-1)*(pressure_ice-pressure_water)*gap
						+beta*sqrt(vx*vx+vy*vy)
						));

		totalweights +=gauss->weight*Jdet;

		/* Compute basal water flux */
		q += gauss->weight*Jdet*(conductivity*sqrt(dh[0]*dh[0]+dh[1]*dh[1]));

		/* Compute "degree of channelization" (ratio of melt opening to opening by sliding) */
		channelization += gauss->weight*Jdet*(meltrate/rho_ice/(meltrate/rho_ice+beta*sqrt(vx*vx+vy*vy)));
	}

	/*Divide by connectivity and constrain gap height*/
	newgap = newgap/totalweights;
	if(newgap<gap_height_min) newgap=gap_height_min;
	if(newgap>gap_height_max) newgap=gap_height_max;

	/*Add new gap as an input*/
	element->AddBasalInput(HydrologyGapHeightEnum,&newgap,P0Enum);

	/*Divide by connectivity, add basal flux as an input*/
	q = q/totalweights;
	element->AddBasalInput(HydrologyBasalFluxEnum,&q,P0Enum);

	/* Divide by connectivity, add degree of channelization as an input */
	channelization = channelization/totalweights;
	element->AddBasalInput(DegreeOfChannelizationEnum,&channelization,P0Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	delete friction;
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void HydrologyShaktiAnalysis::UpdateEffectivePressure(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		UpdateEffectivePressure(element);
	}

}/*}}}*/
void HydrologyShaktiAnalysis::UpdateEffectivePressure(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return;
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
      gauss->GaussNode(element->GetElementType(),i);

		base_input->GetInputValue(&bed,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		head_input->GetInputValue(&head,gauss);

		N[i] = rho_ice*g*thickness - rho_water*g*(head-bed);
	}

	/*Add new gap as an input*/
	element->AddBasalInput(EffectivePressureEnum,N,element->GetElementType());

	/*Clean up and return*/
   xDelete<IssmDouble>(N);
	delete gauss;
}/*}}}*/
