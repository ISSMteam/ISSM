#include "./EnthalpyAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
#include "../cores/cores.h"

/*Model processing*/
void EnthalpyAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Intermediary*/
	int        count;
	int        M,N;
	bool       spcpresent = false;
	int        finiteelement;
	IssmDouble heatcapacity;
	IssmDouble referencetemperature;

	/*Output*/
	IssmDouble *spcvector  = NULL;
	IssmDouble *spcvectorstatic  = NULL;
	IssmDouble* times=NULL;
	IssmDouble* values=NULL;
	IssmDouble* issurface = NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&heatcapacity,"md.materials.heatcapacity");
	iomodel->FindConstant(&referencetemperature,"md.constants.referencetemperature");
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");

	/*return if 2d mesh*/
	if(iomodel->domaintype==Domain2DhorizontalEnum) return;

	/*Fetch data: */
	iomodel->FetchData(&issurface,&M,&N,"md.mesh.vertexonsurface"); _assert_(N>0); _assert_(M==iomodel->numberofvertices);
	iomodel->FetchData(&spcvector,&M,&N,"md.thermal.spctemperature");
	iomodel->FetchData(&spcvectorstatic,&M,&N,"md.thermal.spctemperature");

	/*Specific case for PDD, we want the constaints to be updated by the PDD scheme itself*/
	bool isdynamic = false;
	if (iomodel->solution_enum==TransientSolutionEnum){
		int smb_model;
		iomodel->FindConstant(&smb_model,"md.smb.model");
		if(smb_model==SMBpddEnum)				isdynamic=true;
		if(smb_model==SMBd18opddEnum)			isdynamic=true;
		if(smb_model==SMBpddSicopolisEnum)	isdynamic=true;
	}

	/*Convert spcs from temperatures to enthalpy*/
	_assert_(N>0); _assert_(M>=iomodel->numberofvertices);
	for(int i=0;i<iomodel->numberofvertices;i++){
		for(int j=0;j<N;j++){
			if (isdynamic){
				if (issurface[i]==1){
					spcvector[i*N+j] = heatcapacity*(spcvector[i*N+j]-referencetemperature);
					spcvectorstatic[i*N+j] = NAN;
				}
				else{
					spcvector[i*N+j] = NAN;
					spcvectorstatic[i*N+j] = heatcapacity*(spcvectorstatic[i*N+j]-referencetemperature);
				}
			}
			else{
				spcvector[i*N+j] = heatcapacity*(spcvector[i*N+j]-referencetemperature);
			}
		}
	}

	if(isdynamic){
		IoModelToDynamicConstraintsx(constraints,iomodel,spcvector,iomodel->numberofvertices,1,EnthalpyAnalysisEnum,finiteelement);
		IoModelToConstraintsx(constraints,iomodel,spcvectorstatic,M,N,EnthalpyAnalysisEnum,finiteelement);
	}
	else{
		IoModelToConstraintsx(constraints,iomodel,spcvector,M,N,EnthalpyAnalysisEnum,finiteelement);
	}

	/*Free resources:*/
	iomodel->DeleteData(spcvector,"md.thermal.spctemperature");
	iomodel->DeleteData(spcvectorstatic,"md.thermal.spctemperature");
	iomodel->DeleteData(issurface,"md.mesh.vertexonsurface");
	xDelete<IssmDouble>(times);
	xDelete<IssmDouble>(values);
}/*}}}*/
void EnthalpyAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads */
}/*}}}*/
void EnthalpyAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,EnthalpyAnalysisEnum,finiteelement);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  EnthalpyAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void EnthalpyAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	bool dakota_analysis,ismovingfront,isenthalpy;
	int  basalforcing_model,materialstype;

	/*Now, is the model 3d? otherwise, do nothing: */
	if(iomodel->domaintype==Domain2DhorizontalEnum)return;

	/*Is enthalpy requested?*/
	iomodel->FindConstant(&isenthalpy,"md.thermal.isenthalpy");
	if(!isenthalpy) return;

	/*Fetch data needed: */
	iomodel->FetchData(3,"md.initialization.temperature","md.initialization.waterfraction","md.initialization.pressure");

	/*Finite element type*/
	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");
	iomodel->FindConstant(&materialstype,"md.materials.type");

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.temperature",TemperatureEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.waterfraction",WaterfractionEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.enthalpy",EnthalpyEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",WatercolumnEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VxMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VyMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VzMeshEnum);
	if(ismovingfront){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum); // required for updating active nodes
	}

	/*Basal forcings variables*/
	iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
	switch(basalforcing_model){
		case MantlePlumeGeothermalFluxEnum:
			break;
		default:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.geothermalflux",BasalforcingsGeothermalfluxEnum);
			break;
	}

	/*Rheology type*/
	iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
	switch(materialstype){
		case MatenhancediceEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_E",MaterialsRheologyEEnum);
			break;
		case MatdamageiceEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			break;
		case MatestarEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_Ec",MaterialsRheologyEcEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_Es",MaterialsRheologyEsEnum);
			break;
		case MaticeEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			break;
		default:
			_error_("not supported");
	}

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);

	/*Free data: */
	iomodel->DeleteData(3,"md.initialization.temperature","md.initialization.waterfraction","md.initialization.pressure");

}/*}}}*/
void EnthalpyAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.stabilization",ThermalStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.maxiter",ThermalMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.reltol",ThermalReltolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.isenthalpy",ThermalIsenthalpyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.isdynamicbasalspc",ThermalIsdynamicbasalspcEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.isdrainicecolumn",ThermalIsdrainicecolumnEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.watercolumn_upperlimit",ThermalWatercolumnUpperlimitEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.thermal.requested_outputs");
	parameters->AddObject(new IntParam(ThermalNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(ThermalRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.thermal.requested_outputs");

	/*Friction*/
	FrictionUpdateParameters(parameters, iomodel);
}/*}}}*/

/*Finite Element Analysis*/
void           EnthalpyAnalysis::ApplyBasalConstraints(IssmDouble* local_spc,Element* element){/*{{{*/

	/* Do not check if ice in element, this may lead to inconsistencies between cpu partitions */
	/* Only update constraints at the base. */
	if(!(element->IsOnBase())) return;

	/*Intermediary*/
	bool        isdynamicbasalspc;
	int         numindices;
	int        *indices = NULL;
	IssmDouble	pressure;

	/*Check wether dynamic basal boundary conditions are activated */
	element->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);
	if(!isdynamicbasalspc) return;

	/*Get parameters and inputs: */
	Input* pressure_input = element->GetInput(PressureEnum); _assert_(pressure_input);

	/*Fetch indices of basal & surface nodes for this finite element*/
	Penta *penta =  (Penta *) element; // TODO: add Basal-/SurfaceNodeIndices to element.h, and change this to Element*
	penta->BasalNodeIndices(&numindices,&indices,element->GetElementType());

	GaussPenta* gauss=new GaussPenta();
	for(int i=0;i<numindices;i++){
		gauss->GaussNode(element->GetElementType(),indices[i]);

		pressure_input->GetInputValue(&pressure,gauss);

		/*apply or release spc*/
		Node* node=element->GetNode(indices[i]);
		if(!node->IsActive()) continue;
		if(local_spc[node->Lid()]==1.){
			pressure_input->GetInputValue(&pressure, gauss);
			node->ApplyConstraint(0,PureIceEnthalpy(element,pressure));
		}
		else {
			node->DofInFSet(0);
		}
	}

	/*Free resources:*/
	xDelete<int>(indices);
	delete gauss;
}/*}}}*/
void           EnthalpyAnalysis::ComputeBasalMeltingrate(FemModel* femmodel){/*{{{*/
	/*Compute basal melting rates: */
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		ComputeBasalMeltingrate(element);
	}

	/*extrude inputs*/
	femmodel->parameters->SetParam(BasalforcingsGroundediceMeltingRateEnum,InputToExtrudeEnum);
	extrudefrombase_core(femmodel);
}/*}}}*/
void           EnthalpyAnalysis::ComputeBasalMeltingrate(Element* element){/*{{{*/
	/*Calculate the basal melt rates of the enthalpy model after Aschwanden 2012*/
	/* melting rate is positive when melting, negative when refreezing*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return;

	/* Only compute melt rates at the base of grounded ice*/
	if(!element->IsOnBase() || element->IsAllFloating()) return;

	/* Intermediaries */
	bool			converged;
	const int   dim=3;
	int         i,is,state;
	int			nodedown,nodeup,numnodes,numsegments;
	int			enthalpy_enum;
	IssmDouble  vec_heatflux[dim],normal_base[dim],d1enthalpy[dim],d1pressure[dim];
	IssmDouble  basalfriction,alpha2,geothermalflux,heatflux;
	IssmDouble  dt,yts;
	IssmDouble  melting_overshoot,lambda;
	IssmDouble  vx,vy,vz;
	IssmDouble *xyz_list      = NULL;
	IssmDouble *xyz_list_base = NULL;
	int        *pairindices   = NULL;

	/*Fetch parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GetInputValue(&converged,ConvergedEnum);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&yts, ConstantsYtsEnum);

	if(dt==0. && !converged) enthalpy_enum=EnthalpyPicardEnum;
	else enthalpy_enum=EnthalpyEnum;

	IssmDouble latentheat = element->FindParam(MaterialsLatentheatEnum);
	IssmDouble rho_ice    = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water  = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble beta		 = element->FindParam(MaterialsBetaEnum);
	IssmDouble kappa		 = EnthalpyDiffusionParameterVolume(element,enthalpy_enum);     _assert_(kappa>=0.);
	IssmDouble kappa_mix;

	/*retrieve inputs*/
	Input* enthalpy_input       = element->GetInput(enthalpy_enum);                   _assert_(enthalpy_input);
	Input* pressure_input       = element->GetInput(PressureEnum);                    _assert_(pressure_input);
	Input* geothermalflux_input = element->GetInput(BasalforcingsGeothermalfluxEnum); _assert_(geothermalflux_input);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,3);

	/******** MELTING RATES  ************************************//*{{{*/
	element->NormalBase(&normal_base[0],xyz_list_base);
	element->VerticalSegmentIndicesBase(&pairindices,&numsegments);
	IssmDouble* meltingrate_enthalpy = xNew<IssmDouble>(numsegments);
	IssmDouble* heating = xNew<IssmDouble>(numsegments);

	numnodes=element->GetNumberOfNodes();
	IssmDouble* enthalpies = xNew<IssmDouble>(numnodes);
	IssmDouble* pressures = xNew<IssmDouble>(numnodes);
	IssmDouble* watercolumns = xNew<IssmDouble>(numnodes);
	IssmDouble* basalmeltingrates = xNew<IssmDouble>(numnodes);
	element->GetInputListOnNodes(enthalpies,enthalpy_enum);
	element->GetInputListOnNodes(pressures,PressureEnum);
	element->GetInputListOnNodes(watercolumns,WatercolumnEnum);
	element->GetInputListOnNodes(basalmeltingrates,BasalforcingsGroundediceMeltingRateEnum);

	IssmDouble watercolumnupperlimit = element->FindParam(ThermalWatercolumnUpperlimitEnum);

	Gauss* gauss=element->NewGauss();
	for(is=0;is<numsegments;is++){
		nodedown = pairindices[is*2+0];
		nodeup   = pairindices[is*2+1];
		gauss->GaussNode(element->GetElementType(),nodedown);

		state=GetThermalBasalCondition(element, enthalpies[nodedown], enthalpies[nodeup], pressures[nodedown], pressures[nodeup], watercolumns[nodedown], basalmeltingrates[nodedown]);
		switch (state) {
			case 0:
				// cold, dry base: apply basal surface forcing
				for(i=0;i<3;i++) vec_heatflux[i]=0.;
				break;
			case 1: case 2: case 3:
				// case 1 : cold, wet base: keep at pressure melting point
				// case 2: temperate, thin refreezing base: release spc
				// case 3: temperate, thin melting base: set spc
				enthalpy_input->GetInputDerivativeValue(&d1enthalpy[0],xyz_list,gauss);
				for(i=0;i<3;i++) vec_heatflux[i]=-kappa*d1enthalpy[i];
				break;
			case 4:
				// temperate, thick melting base: set grad H*n=0
				kappa_mix=GetWetIceConductivity(element, enthalpies[nodedown], pressures[nodedown]);
				pressure_input->GetInputDerivativeValue(&d1pressure[0],xyz_list,gauss);
				for(i=0;i<3;i++) vec_heatflux[i]=kappa_mix*beta*d1pressure[i];
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}
		if(state==0) meltingrate_enthalpy[is]=0.;
		else{
			/*heat flux along normal*/
			heatflux=0.;
			for(i=0;i<3;i++) heatflux+=(vec_heatflux[i])*normal_base[i];

			/*basal friction*/
			friction->GetAlpha2(&alpha2,gauss);
			friction->GetBasalSlidingSpeeds(&vx, &vy, &vz, gauss);
			basalfriction=alpha2*(vx*vx + vy*vy + vz*vz);
			geothermalflux_input->GetInputValue(&geothermalflux,gauss);
			/* -Mb= Fb-(q-q_geo)/((1-w)*L*rho), and (1-w)*rho=rho_ice, cf Aschwanden 2012, eqs.1, 2, 66*/
			heating[is]=(heatflux+basalfriction+geothermalflux);
			meltingrate_enthalpy[is]=heating[is]/(latentheat*rho_ice); // m/s water equivalent
		}
	}/*}}}*/

	/******** UPDATE MELTINGRATES AND WATERCOLUMN **************//*{{{*/
	for(is=0;is<numsegments;is++){
		nodedown = pairindices[is*2+0];
		nodeup   = pairindices[is*2+1];
		if(dt!=0.){
			if(watercolumns[nodedown]+meltingrate_enthalpy[is]*dt<0.){	// prevent too much freeze on
				lambda = -watercolumns[nodedown]/(dt*meltingrate_enthalpy[is]); _assert_(lambda>=0.); _assert_(lambda<1.);
				watercolumns[nodedown]=0.;
				basalmeltingrates[nodedown]=lambda*meltingrate_enthalpy[is]; // restrict freeze on only to size of watercolumn
				enthalpies[nodedown]+=(1.-lambda)*dt/yts*meltingrate_enthalpy[is]*latentheat*rho_ice; // use rest of energy to cool down base: dE=L*m, m=(1-lambda)*meltingrate*rho_ice
			}
			else{
				basalmeltingrates[nodedown]=meltingrate_enthalpy[is];
				watercolumns[nodedown]+=dt*meltingrate_enthalpy[is];
			}
			if(watercolumns[nodedown]>watercolumnupperlimit) watercolumns[nodedown]=watercolumnupperlimit;
		}
		else{
			basalmeltingrates[nodedown]=meltingrate_enthalpy[is];
			if(watercolumns[nodedown]+meltingrate_enthalpy[is]<0.)
				watercolumns[nodedown]=0.;
			else
				watercolumns[nodedown]+=meltingrate_enthalpy[is];
		}
		basalmeltingrates[nodedown]*=rho_water/rho_ice; // convert meltingrate from water to ice equivalent
		_assert_(watercolumns[nodedown]>=0.);
	}/*}}}*/

	/*feed updated variables back into model*/
	int finite_element = element->GetElementType(); if(finite_element==P1Enum) finite_element = P1DGEnum;
	if(dt!=0.){
		element->AddInput(enthalpy_enum,enthalpies,finite_element);
		element->AddInput(WatercolumnEnum,watercolumns,finite_element);
	}
	element->AddInput(BasalforcingsGroundediceMeltingRateEnum,basalmeltingrates,P1DGEnum);

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<int>(pairindices);
	xDelete<IssmDouble>(enthalpies);
	xDelete<IssmDouble>(pressures);
	xDelete<IssmDouble>(watercolumns);
	xDelete<IssmDouble>(basalmeltingrates);
	xDelete<IssmDouble>(meltingrate_enthalpy);
	xDelete<IssmDouble>(heating);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_base);
}/*}}}*/
void           EnthalpyAnalysis::Core(FemModel* femmodel){/*{{{*/

	IssmDouble dt;
	bool isdynamicbasalspc;

	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	femmodel->parameters->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);

	if(VerboseSolution()) _printf0_("   computing enthalpy\n");
	femmodel->SetCurrentConfiguration(EnthalpyAnalysisEnum);
	if((dt>0.) && isdynamicbasalspc)	UpdateBasalConstraints(femmodel);
	solutionsequence_thermal_nonlinear(femmodel);

	/*transfer enthalpy to enthalpy picard for the next step: */
	InputDuplicatex(femmodel,EnthalpyEnum,EnthalpyPicardEnum);

	PostProcessing(femmodel);

}/*}}}*/
void           EnthalpyAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* EnthalpyAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixVolume(element);
	ElementMatrix* Ke2=CreateKMatrixShelf(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateKMatrixVolume(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         stabilization;
	IssmDouble  Jdet,dt,u,v,w,um,vm,wm,vel;
	IssmDouble  h,hx,hy,hz,vx,vy,vz;
	IssmDouble  tau_parameter,diameter,factor;
	IssmDouble  tau_parameter_anisotropic[2],tau_parameter_hor,tau_parameter_ver;
	IssmDouble  D_scalar;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke       = element->NewElementMatrix();
	IssmDouble*    basis    = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis   = xNew<IssmDouble>(3*numnodes);
	IssmDouble     K[3][3];

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,ThermalStabilizationEnum);
	IssmDouble  rho_water           = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble  rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  gravity             = element->FindParam(ConstantsGEnum);
	IssmDouble  heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble  thermalconductivity = element->FindParam(MaterialsThermalconductivityEnum);
	Input* vx_input  = element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input  = element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input  = element->GetInput(VzEnum);     _assert_(vz_input);
	Input* vxm_input = element->GetInput(VxMeshEnum); _assert_(vxm_input);
	Input* vym_input = element->GetInput(VyMeshEnum); _assert_(vym_input);
	Input* vzm_input = element->GetInput(VzMeshEnum); _assert_(vzm_input);

	/*Enthalpy diffusion parameter*/
	IssmDouble kappa=this->EnthalpyDiffusionParameterVolume(element,EnthalpyPicardEnum); _assert_(kappa>=0.);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		D_scalar=gauss->weight*Jdet;
		if(dt!=0.) D_scalar=D_scalar*dt;

		/*Conduction: */
		factor = D_scalar*kappa/rho_ice;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += factor*(
							dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[2*numnodes+j]*dbasis[2*numnodes+i]);
			}
		}

		/*Advection: */
		vx_input->GetInputValue(&u,gauss); vxm_input->GetInputValue(&um,gauss); vx=u-um;
		vy_input->GetInputValue(&v,gauss); vym_input->GetInputValue(&vm,gauss); vy=v-vm;
		vz_input->GetInputValue(&w,gauss); vzm_input->GetInputValue(&wm,gauss); vz=w-wm;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D_scalar*(
							vx*dbasis[0*numnodes+j]*basis[i] + vy*dbasis[1*numnodes+j]*basis[i] +vz*dbasis[2*numnodes+j]*basis[i]
							);
			}
		}

		/*Transient: */
		if(dt!=0.){
			D_scalar=gauss->weight*Jdet;
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += D_scalar*basis[j]*basis[i];
				}
			}
			D_scalar=D_scalar*dt;
		}

		/*Artificial diffusivity*/
		if(stabilization==1){
			element->ElementSizes(&hx,&hy,&hz);
			vel=sqrt(vx*vx + vy*vy + vz*vz)+1.e-14;
			h=sqrt( pow(hx*vx/vel,2) + pow(hy*vy/vel,2) + pow(hz*vz/vel,2));
			K[0][0]=h/(2.*vel)*fabs(vx*vx);  K[0][1]=h/(2.*vel)*fabs(vx*vy); K[0][2]=h/(2.*vel)*fabs(vx*vz);
			K[1][0]=h/(2.*vel)*fabs(vy*vx);  K[1][1]=h/(2.*vel)*fabs(vy*vy); K[1][2]=h/(2.*vel)*fabs(vy*vz);
			K[2][0]=h/(2.*vel)*fabs(vz*vx);  K[2][1]=h/(2.*vel)*fabs(vz*vy); K[2][2]=h/(2.*vel)*fabs(vz*vz);
			for(int i=0;i<3;i++) for(int j=0;j<3;j++) K[i][j] = D_scalar*K[i][j];

			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += (
								dbasis[0*numnodes+i] *(K[0][0]*dbasis[0*numnodes+j] + K[0][1]*dbasis[1*numnodes+j]+ K[0][2]*dbasis[2*numnodes+j]) +
								dbasis[1*numnodes+i] *(K[1][0]*dbasis[0*numnodes+j] + K[1][1]*dbasis[1*numnodes+j]+ K[1][2]*dbasis[2*numnodes+j]) +
								dbasis[2*numnodes+i] *(K[2][0]*dbasis[0*numnodes+j] + K[2][1]*dbasis[1*numnodes+j]+ K[2][2]*dbasis[2*numnodes+j])
								);
				}
			}
		}
		/*SUPG*/
		else if(stabilization==2){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			diameter=element->MinEdgeLength(xyz_list);
			tau_parameter=element->StabilizationParameter(u-um,v-vm,w-wm,diameter,kappa/rho_ice);
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=tau_parameter*D_scalar*
					  ((u-um)*dbasis[0*numnodes+i]+(v-vm)*dbasis[1*numnodes+i]+(w-wm)*dbasis[2*numnodes+i])*
					  ((u-um)*dbasis[0*numnodes+j]+(v-vm)*dbasis[1*numnodes+j]+(w-wm)*dbasis[2*numnodes+j]);
				}
			}
			if(dt!=0.){
				D_scalar=gauss->weight*Jdet;
				factor = tau_parameter*D_scalar;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=factor*basis[j]*((u-um)*dbasis[0*numnodes+i]+(v-vm)*dbasis[1*numnodes+i]+(w-wm)*dbasis[2*numnodes+i]);
					}
				}
			}
		}
		/*anisotropic SUPG*/
		else if(stabilization==3){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			element->ElementSizes(&hx,&hy,&hz);
			element->StabilizationParameterAnisotropic(&tau_parameter_anisotropic[0],u-um,v-vm,w-wm,hx,hy,hz,kappa/rho_ice);
			tau_parameter_hor=tau_parameter_anisotropic[0];
			tau_parameter_ver=tau_parameter_anisotropic[1];
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=D_scalar*
						(sqrt(tau_parameter_hor)*(u-um)*dbasis[0*numnodes+i]+sqrt(tau_parameter_hor)*(v-vm)*dbasis[1*numnodes+i]+sqrt(tau_parameter_ver)*(w-wm)*dbasis[2*numnodes+i])*
						(sqrt(tau_parameter_hor)*(u-um)*dbasis[0*numnodes+j]+sqrt(tau_parameter_hor)*(v-vm)*dbasis[1*numnodes+j]+sqrt(tau_parameter_ver)*(w-wm)*dbasis[2*numnodes+j]);
				}
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateKMatrixShelf(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Initialize Element matrix and return if necessary*/
	if(!element->IsOnBase() || !element->IsAllFloating()) return NULL;

	/*Intermediaries*/
	IssmDouble  dt,Jdet,D;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementMatrix* Ke    = element->NewElementMatrix();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  gravity             = element->FindParam(ConstantsGEnum);
	IssmDouble  rho_water           = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble  rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble  mixed_layer_capacity= element->FindParam(MaterialsMixedLayerCapacityEnum);
	IssmDouble  thermal_exchange_vel= element->FindParam(MaterialsThermalExchangeVelocityEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		D=gauss->weight*Jdet*rho_water*mixed_layer_capacity*thermal_exchange_vel/(heatcapacity*rho_ice);
		if(reCast<bool,IssmDouble>(dt)) D=dt*D;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D*basis[i]*basis[j];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return Ke;
}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorVolume(element);
	ElementVector* pe2=CreatePVectorSheet(element);
	ElementVector* pe3=CreatePVectorShelf(element);
	ElementVector* pe =new ElementVector(pe1,pe2,pe3);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	delete pe3;
	return pe;
}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVectorVolume(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         i, stabilization;
	IssmDouble  Jdet,phi,dt;
	IssmDouble  enthalpy, Hpmp;
	IssmDouble  enthalpypicard, d1enthalpypicard[3];
	IssmDouble  pressure, d1pressure[3], d2pressure;
	IssmDouble  waterfractionpicard;
	IssmDouble  kappa,tau_parameter,diameter,hx,hy,hz,kappa_w;
	IssmDouble  tau_parameter_anisotropic[2],tau_parameter_hor,tau_parameter_ver;
	IssmDouble  u,v,w;
	IssmDouble  scalar_def, scalar_sens ,scalar_transient;
	IssmDouble* xyz_list = NULL;
	IssmDouble  d1H_d1P, d1P2;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = element->GetNumberOfNodes();
	int numvertices = element->GetNumberOfVertices();

	/*Initialize Element vector*/
	ElementVector* pe             = element->NewElementVector();
	IssmDouble*    basis          = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis         = xNew<IssmDouble>(3*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble  rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble  thermalconductivity = element->FindParam(MaterialsThermalconductivityEnum);
	IssmDouble  temperateiceconductivity = element->FindParam(MaterialsTemperateiceconductivityEnum);
	IssmDouble  beta                = element->FindParam(MaterialsBetaEnum);
	IssmDouble  latentheat          = element->FindParam(MaterialsLatentheatEnum);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,ThermalStabilizationEnum);
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=element->GetInput(VzEnum); _assert_(vz_input);
	Input* enthalpypicard_input=element->GetInput(EnthalpyPicardEnum); _assert_(enthalpypicard_input);
	Input* pressure_input=element->GetInput(PressureEnum); _assert_(pressure_input);
	Input* enthalpy_input=NULL;
	if(dt>0.){
		enthalpy_input = element->GetInput(EnthalpyEnum); _assert_(enthalpy_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/*viscous dissipation*/
		element->ViscousHeating(&phi,xyz_list,gauss,vx_input,vy_input,vz_input);

		scalar_def=phi/rho_ice*Jdet*gauss->weight;
		if(dt!=0.) scalar_def=scalar_def*dt;

		for(i=0;i<numnodes;i++) pe->values[i]+=scalar_def*basis[i];

		/*sensible heat flux in temperate ice*/
		enthalpypicard_input->GetInputValue(&enthalpypicard,gauss);
		pressure_input->GetInputValue(&pressure,gauss);
		Hpmp=this->PureIceEnthalpy(element, pressure);

		if(enthalpypicard>=Hpmp){
			enthalpypicard_input->GetInputDerivativeValue(&d1enthalpypicard[0],xyz_list,gauss);
			pressure_input->GetInputDerivativeValue(&d1pressure[0],xyz_list,gauss);
			d2pressure=0.; // for linear elements, 2nd derivative is zero

			d1H_d1P=0.;
			for(i=0;i<3;i++) d1H_d1P+=d1enthalpypicard[i]*d1pressure[i];
			d1P2=0.;
			for(i=0;i<3;i++) d1P2+=pow(d1pressure[i],2.);

			scalar_sens=-beta*((temperateiceconductivity - thermalconductivity)/latentheat*(d1H_d1P + beta*heatcapacity*d1P2))/rho_ice;
			if(dt!=0.) scalar_sens=scalar_sens*dt;
			for(i=0;i<numnodes;i++) pe->values[i]+=scalar_sens*basis[i];
		}

		/* Build transient now */
		if(dt>0.){
			enthalpy_input->GetInputValue(&enthalpy, gauss);
			scalar_transient=enthalpy*Jdet*gauss->weight;
			for(i=0;i<numnodes;i++) pe->values[i]+=scalar_transient*basis[i];
		}

		/* SUPG */
		if(stabilization==2){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			diameter=element->MinEdgeLength(xyz_list);
			kappa=this->EnthalpyDiffusionParameterVolume(element,EnthalpyPicardEnum); _assert_(kappa>=0.);
			vx_input->GetInputValue(&u,gauss);
			vy_input->GetInputValue(&v,gauss);
			vz_input->GetInputValue(&w,gauss);
			tau_parameter=element->StabilizationParameter(u,v,w,diameter,kappa/rho_ice);

			for(i=0;i<numnodes;i++) pe->values[i]+=tau_parameter*scalar_def*(u*dbasis[0*numnodes+i]+v*dbasis[1*numnodes+i]+w*dbasis[2*numnodes+i]);

			if(dt!=0.){
				for(i=0;i<numnodes;i++) pe->values[i]+=tau_parameter*scalar_transient*(u*dbasis[0*numnodes+i]+v*dbasis[1*numnodes+i]+w*dbasis[2*numnodes+i]);
			}
		}
		/* anisotropic SUPG */
		else if(stabilization==3){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			element->ElementSizes(&hx,&hy,&hz);
			kappa=this->EnthalpyDiffusionParameterVolume(element,EnthalpyPicardEnum); _assert_(kappa>=0.);
			vx_input->GetInputValue(&u,gauss);
			vy_input->GetInputValue(&v,gauss);
			vz_input->GetInputValue(&w,gauss);
			element->StabilizationParameterAnisotropic(&tau_parameter_anisotropic[0],u,v,w,hx,hy,hz,kappa/rho_ice);
			tau_parameter_hor=tau_parameter_anisotropic[0];
			tau_parameter_ver=tau_parameter_anisotropic[1];

			for(i=0;i<numnodes;i++) pe->values[i]+=scalar_def*(tau_parameter_hor*u*dbasis[0*numnodes+i]+tau_parameter_hor*v*dbasis[1*numnodes+i]+tau_parameter_ver*w*dbasis[2*numnodes+i]);
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;

}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVectorSheet(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/* implementation of the basal condition decision chart of Aschwanden 2012, Fig.5 */
	if(!element->IsOnBase() || element->IsAllFloating()) return NULL;

	bool converged, isdynamicbasalspc;
	int i, state;
	int enthalpy_enum;
	IssmDouble  dt,Jdet,scalar;
	IssmDouble	enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate;
	IssmDouble	vx,vy,vz;
	IssmDouble  alpha2,basalfriction,geothermalflux,heatflux;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);
	element->GetInputValue(&converged,ConvergedEnum);
	if(dt==0. && !converged) enthalpy_enum=EnthalpyPicardEnum; // use enthalpy from last iteration
	else enthalpy_enum=EnthalpyEnum; // use enthalpy from last time step
	Input* enthalpy_input		 = element->GetInput(enthalpy_enum);					 _assert_(enthalpy_input);
	Input* pressure_input		 = element->GetInput(PressureEnum);							 _assert_(pressure_input);
	Input* watercolumn_input	 = element->GetInput(WatercolumnEnum);							 _assert_(watercolumn_input);
	Input* meltingrate_input	 = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);							 _assert_(meltingrate_input);
	Input* geothermalflux_input = element->GetInput(BasalforcingsGeothermalfluxEnum); _assert_(geothermalflux_input);
	IssmDouble  rho_ice			 = element->FindParam(MaterialsRhoIceEnum);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,3);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	Gauss* gaussup=element->NewGaussTop(4);
	while(gauss->next() && gaussup->next()){

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		if(isdynamicbasalspc){
			enthalpy_input->GetInputValue(&enthalpy,gauss);
			enthalpy_input->GetInputValue(&enthalpyup,gaussup);
			pressure_input->GetInputValue(&pressure,gauss);
			pressure_input->GetInputValue(&pressureup,gaussup);
			watercolumn_input->GetInputValue(&watercolumn,gauss);
			meltingrate_input->GetInputValue(&meltingrate,gauss);
			state=GetThermalBasalCondition(element, enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate);
		}
		else
			state=0;

		switch (state) {
			case 0: case 1: case 2: case 3:
				// cold, dry base; cold, wet base; refreezing temperate base; thin temperate base:
				// Apply basal surface forcing.
				// Interpolated values of enthalpy on gauss nodes may indicate cold base,
				// although one node might have become temperate. So keep heat flux switched on.
				geothermalflux_input->GetInputValue(&geothermalflux,gauss);
				friction->GetAlpha2(&alpha2,gauss);
				friction->GetBasalSlidingSpeeds(&vx, &vy, &vz, gauss);
				basalfriction=alpha2*(vx*vx+vy*vy+vz*vz);
				heatflux=(basalfriction+geothermalflux)/(rho_ice);
				scalar=gauss->weight*Jdet*heatflux;
				if(dt!=0.) scalar=dt*scalar;
				for(i=0;i<numnodes;i++)
					pe->values[i]+=scalar*basis[i];
				break;
			case 4:
				// temperate, thick melting base: set grad H*n=0
				for(i=0;i<numnodes;i++)
					pe->values[i]+=0.;
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}
	}

	/*Clean up and return*/
	delete gauss;
	delete gaussup;
	delete friction;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;

}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVectorShelf(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Get basal element*/
	if(!element->IsOnBase() || !element->IsAllFloating()) return NULL;

	IssmDouble  Hpmp,dt,Jdet,scalar_ocean,pressure;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	Input*      pressure_input=element->GetInput(PressureEnum); _assert_(pressure_input);
	IssmDouble  gravity             = element->FindParam(ConstantsGEnum);
	IssmDouble  rho_water           = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble  rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble  mixed_layer_capacity= element->FindParam(MaterialsMixedLayerCapacityEnum);
	IssmDouble  thermal_exchange_vel= element->FindParam(MaterialsThermalExchangeVelocityEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		pressure_input->GetInputValue(&pressure,gauss);
		Hpmp=element->PureIceEnthalpy(pressure);

		scalar_ocean=gauss->weight*Jdet*rho_water*mixed_layer_capacity*thermal_exchange_vel*Hpmp/(heatcapacity*rho_ice);
		if(reCast<bool,IssmDouble>(dt)) scalar_ocean=dt*scalar_ocean;

		for(int i=0;i<numnodes;i++) pe->values[i]+=scalar_ocean*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;
}/*}}}*/
void           EnthalpyAnalysis::DrainWaterfraction(FemModel* femmodel){/*{{{*/
	/*Drain excess water fraction in ice column: */
	ComputeWaterfractionDrainage(femmodel);
	DrainageUpdateWatercolumn(femmodel);
	DrainageUpdateEnthalpy(femmodel);
}/*}}}*/
void				EnthalpyAnalysis::ComputeWaterfractionDrainage(FemModel* femmodel){/*{{{*/

	int k,numnodes;
	IssmDouble dt;
	Element* element= NULL;

	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	for(Object* & object : femmodel->elements->objects){
		element=xDynamicCast<Element*>(object);
		numnodes=element->GetNumberOfNodes();
		IssmDouble* waterfractions= xNew<IssmDouble>(numnodes);
		IssmDouble* drainage= xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(waterfractions,WaterfractionEnum);
		for(k=0; k<numnodes;k++){
			drainage[k]=DrainageFunctionWaterfraction(waterfractions[k], dt);
		}
		int finite_element = element->GetElementType(); if(finite_element==P1Enum) finite_element = P1DGEnum;
		element->AddInput(WaterfractionDrainageEnum,drainage,finite_element);

		xDelete<IssmDouble>(waterfractions);
		xDelete<IssmDouble>(drainage);
	}
}/*}}}*/
void				EnthalpyAnalysis::DrainageUpdateWatercolumn(FemModel* femmodel){/*{{{*/

	int k,numnodes, numbasalnodes;
	IssmDouble dt;
	int* basalnodeindices=NULL;
	Element* element= NULL;

	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	/*depth-integrate the drained water fraction */
	femmodel->parameters->SetParam(WaterfractionDrainageEnum,InputToDepthaverageInEnum);
	femmodel->parameters->SetParam(WaterfractionDrainageIntegratedEnum,InputToDepthaverageOutEnum);
	depthaverage_core(femmodel);
	femmodel->parameters->SetParam(WaterfractionDrainageIntegratedEnum,InputToExtrudeEnum);
	extrudefrombase_core(femmodel);
	/*multiply depth-average by ice thickness*/
	for(Object* & object : femmodel->elements->objects){
		element=xDynamicCast<Element*>(object);
		numnodes=element->GetNumberOfNodes();
		IssmDouble* drainage_int= xNew<IssmDouble>(numnodes);
		IssmDouble* thicknesses= xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(drainage_int,WaterfractionDrainageIntegratedEnum);
		element->GetInputListOnNodes(thicknesses,ThicknessEnum);
		for(k=0;k<numnodes;k++){
			drainage_int[k]*=thicknesses[k];
		}
		int finite_element = element->GetElementType(); if(finite_element==P1Enum) finite_element = P1DGEnum;
		element->AddInput(WaterfractionDrainageIntegratedEnum, drainage_int,finite_element);

		xDelete<IssmDouble>(drainage_int);
		xDelete<IssmDouble>(thicknesses);
	}

	/*update water column*/
	for(Object* & object : femmodel->elements->objects){
		element=xDynamicCast<Element*>(object);
		/* Check if ice in element */
		if(!element->IsIceInElement()) continue;
		if(!element->IsOnBase()) continue;

		numnodes=element->GetNumberOfNodes();
		IssmDouble* watercolumn= xNew<IssmDouble>(numnodes);
		IssmDouble* drainage_int= xNew<IssmDouble>(numnodes);
		element->GetInputListOnNodes(watercolumn,WatercolumnEnum);
		element->GetInputListOnNodes(drainage_int,WaterfractionDrainageIntegratedEnum);

		element->BasalNodeIndices(&numbasalnodes,&basalnodeindices,element->GetElementType());
		for(k=0;k<numbasalnodes;k++){
			watercolumn[basalnodeindices[k]]+=dt*drainage_int[basalnodeindices[k]];
		}
		int finite_element = element->GetElementType(); if(finite_element==P1Enum) finite_element = P1DGEnum;
		element->AddInput(WatercolumnEnum, watercolumn,finite_element);

		xDelete<IssmDouble>(watercolumn);
		xDelete<IssmDouble>(drainage_int);
		xDelete<int>(basalnodeindices);
	}
}/*}}}*/
void				EnthalpyAnalysis::DrainageUpdateEnthalpy(FemModel* femmodel){/*{{{*/

	int k,numnodes;
	IssmDouble dt;
	Element* element= NULL;
	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	for(Object* & object : femmodel->elements->objects){
		element=xDynamicCast<Element*>(object);
		numnodes=element->GetNumberOfNodes();
		IssmDouble* enthalpies= xNew<IssmDouble>(numnodes);
		IssmDouble* pressures= xNew<IssmDouble>(numnodes);
		IssmDouble* temperatures= xNew<IssmDouble>(numnodes);
		IssmDouble* waterfractions= xNew<IssmDouble>(numnodes);
		IssmDouble* drainage= xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(pressures,PressureEnum);
		element->GetInputListOnNodes(temperatures,TemperatureEnum);
		element->GetInputListOnNodes(waterfractions,WaterfractionEnum);
		element->GetInputListOnNodes(drainage,WaterfractionDrainageEnum);

		for(k=0;k<numnodes;k++){
			if(dt==0.)
				waterfractions[k]-=drainage[k];
			else
				waterfractions[k]-=dt*drainage[k];

			element->ThermalToEnthalpy(&enthalpies[k], temperatures[k], waterfractions[k], pressures[k]);
		}
		int finite_element = element->GetElementType(); if(finite_element==P1Enum) finite_element = P1DGEnum;
		element->AddInput(WaterfractionEnum,waterfractions,finite_element);
		element->AddInput(EnthalpyEnum,enthalpies,finite_element);

		xDelete<IssmDouble>(enthalpies);
		xDelete<IssmDouble>(pressures);
		xDelete<IssmDouble>(temperatures);
		xDelete<IssmDouble>(waterfractions);
		xDelete<IssmDouble>(drainage);
	}
}/*}}}*/
IssmDouble     EnthalpyAnalysis::EnthalpyDiffusionParameter(Element* element,IssmDouble enthalpy,IssmDouble pressure){/*{{{*/

	IssmDouble heatcapacity             = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble temperateiceconductivity = element->FindParam(MaterialsTemperateiceconductivityEnum);
	IssmDouble thermalconductivity      = element->FindParam(MaterialsThermalconductivityEnum);

	if(enthalpy < PureIceEnthalpy(element,pressure))
		return thermalconductivity/heatcapacity;
	else
		return temperateiceconductivity/heatcapacity;
}/*}}}*/
IssmDouble     EnthalpyAnalysis::EnthalpyDiffusionParameterVolume(Element* element,int enthalpy_enum){/*{{{*/

	int         iv;
	IssmDouble  lambda;                   /* fraction of cold ice    */
	IssmDouble  kappa,kappa_c,kappa_t; /* enthalpy conductivities */
	IssmDouble  Hc,Ht;

	/*Get pressures and enthalpies on vertices*/
	int         numvertices = element->GetNumberOfVertices();
	int         effectiveconductivity_averaging;
	IssmDouble* pressures   = xNew<IssmDouble>(numvertices);
	IssmDouble* enthalpies  = xNew<IssmDouble>(numvertices);
	IssmDouble* PIE         = xNew<IssmDouble>(numvertices);
	IssmDouble* dHpmp       = xNew<IssmDouble>(numvertices);
	element->GetInputListOnVertices(pressures,PressureEnum);
	element->GetInputListOnVertices(enthalpies,enthalpy_enum);
	element->FindParam(&effectiveconductivity_averaging,MaterialsEffectiveconductivityAveragingEnum);

	for(iv=0;iv<numvertices;iv++){
		PIE[iv]   = PureIceEnthalpy(element,pressures[iv]);
		dHpmp[iv] = enthalpies[iv]-PIE[iv];
	}

	bool allequalsign = true;
	if(dHpmp[0]<0.){
		for(iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]<0.));
	}
	else{
		for(iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]>=0.));
	}

	if(allequalsign){
		kappa = EnthalpyDiffusionParameter(element,enthalpies[0],pressures[0]);
	}
	else{
		kappa_c = EnthalpyDiffusionParameter(element,PureIceEnthalpy(element,0.)-1.,0.);
		kappa_t = EnthalpyDiffusionParameter(element,PureIceEnthalpy(element,0.)+1.,0.);

		Hc=0.; Ht=0.;
		for(iv=0; iv<numvertices;iv++){
			if(enthalpies[iv]<PIE[iv])
			 Hc+=(PIE[iv]-enthalpies[iv]);
			else
			 Ht+=(enthalpies[iv]-PIE[iv]);
		}
		_assert_((Hc+Ht)>0.);
		lambda = Hc/(Hc+Ht);
		_assert_(lambda>=0.);
		_assert_(lambda<=1.);

		if(effectiveconductivity_averaging==0){
			/* return arithmetic mean (volume average) of thermal conductivities, weighted by fraction of cold/temperate ice */
			kappa=kappa_c*lambda+(1.-lambda)*kappa_t;
		}
		else if(effectiveconductivity_averaging==1){
			/* return harmonic mean (reciprocal avarage) of thermal conductivities, weighted by fraction of cold/temperate ice, cf Patankar 1980, pp44 */
			kappa=kappa_c*kappa_t/(lambda*kappa_t+(1.-lambda)*kappa_c);
		}
		else if(effectiveconductivity_averaging==2){
			/* return geometric mean (power law) of thermal conductivities, weighted by fraction of cold/temperate ice */
			kappa=pow(kappa_c,lambda)*pow(kappa_t,1.-lambda);
		}
		else{
			_error_("effectiveconductivity_averaging not supported yet");
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(PIE);
	xDelete<IssmDouble>(dHpmp);
	xDelete<IssmDouble>(pressures);
	xDelete<IssmDouble>(enthalpies);
	return kappa;
}/*}}}*/
void           EnthalpyAnalysis::GetBasalConstraints(Vector<IssmDouble>* vec_spc,Element* element){/*{{{*/

	/*Intermediary*/
	bool        isdynamicbasalspc;
	IssmDouble	dt;

	/*Check wether dynamic basal boundary conditions are activated */
	element->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);
	if(!isdynamicbasalspc) return;

	element->FindParam(&dt,TimesteppingTimeStepEnum);
	if(dt==0.){
		GetBasalConstraintsSteadystate(vec_spc,element);
	}
	else{
		GetBasalConstraintsTransient(vec_spc,element);
	}
}/*}}}*/
void           EnthalpyAnalysis::GetBasalConstraintsSteadystate(Vector<IssmDouble>* vec_spc,Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return;

	/* Only update constraints at the base.
	 * Floating ice is not affected by basal BC decision chart. */
	if(!(element->IsOnBase()) || element->IsAllFloating()) return;

	/*Intermediary*/
	int         numindices, numindicesup, state;
	int        *indices = NULL, *indicesup = NULL;
	IssmDouble	enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate;

	/*Get parameters and inputs: */
	Input* enthalpy_input		 = element->GetInput(EnthalpyPicardEnum);					 _assert_(enthalpy_input);
	Input* pressure_input		 = element->GetInput(PressureEnum);							 _assert_(pressure_input);
	Input* watercolumn_input	 = element->GetInput(WatercolumnEnum);							 _assert_(watercolumn_input);
	Input* meltingrate_input	 = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);							 _assert_(meltingrate_input);

	/*Fetch indices of basal & surface nodes for this finite element*/
	Penta *penta =  (Penta *) element; // TODO: add Basal-/SurfaceNodeIndices to element.h, and change this to Element*
	penta->BasalNodeIndices(&numindices,&indices,element->GetElementType());
	penta->SurfaceNodeIndices(&numindicesup,&indicesup,element->GetElementType());	_assert_(numindices==numindicesup);

	GaussPenta* gauss=new GaussPenta();
	GaussPenta* gaussup=new GaussPenta();
	for(int i=0;i<numindices;i++){
		gauss->GaussNode(element->GetElementType(),indices[i]);
		gaussup->GaussNode(element->GetElementType(),indicesup[i]);

		enthalpy_input->GetInputValue(&enthalpy,gauss);
		enthalpy_input->GetInputValue(&enthalpyup,gaussup);
		pressure_input->GetInputValue(&pressure,gauss);
		pressure_input->GetInputValue(&pressureup,gaussup);
		watercolumn_input->GetInputValue(&watercolumn,gauss);
		meltingrate_input->GetInputValue(&meltingrate,gauss);

		state=GetThermalBasalCondition(element, enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate);
		switch (state) {
			case 0:
				// cold, dry base: apply basal surface forcing
				vec_spc->SetValue(element->nodes[i]->Pid(),0.,INS_VAL);
				break;
			case 1:
				// cold, wet base: keep at pressure melting point
				vec_spc->SetValue(element->nodes[i]->Pid(),1.,INS_VAL);
				break;
			case 2:
				// temperate, thin refreezing base:
				vec_spc->SetValue(element->nodes[i]->Pid(),1.,INS_VAL);
				break;
			case 3:
				// temperate, thin melting base: set spc
				vec_spc->SetValue(element->nodes[i]->Pid(),1.,INS_VAL);
				break;
			case 4:
				// temperate, thick melting base:
				vec_spc->SetValue(element->nodes[i]->Pid(),1.,INS_VAL);
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}
	}

	/*Free resources:*/
	xDelete<int>(indices);
	xDelete<int>(indicesup);
	delete gauss;
	delete gaussup;
}/*}}}*/
void           EnthalpyAnalysis::GetBasalConstraintsTransient(Vector<IssmDouble>* vec_spc,Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return;

	/* Only update constraints at the base.
	 * Floating ice is not affected by basal BC decision chart.*/
	if(!(element->IsOnBase()) || element->IsAllFloating()) return;

	/*Intermediary*/
	int         numindices, numindicesup, state;
	int        *indices = NULL, *indicesup = NULL;
	IssmDouble	enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate;

	/*Get parameters and inputs: */
	Input* enthalpy_input    = element->GetInput(EnthalpyEnum);                            _assert_(enthalpy_input); //TODO: check EnthalpyPicard?
	Input* pressure_input    = element->GetInput(PressureEnum);                            _assert_(pressure_input);
	Input* watercolumn_input = element->GetInput(WatercolumnEnum);                         _assert_(watercolumn_input);
	Input* meltingrate_input = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(meltingrate_input);

	/*Fetch indices of basal & surface nodes for this finite element*/
	Penta *penta =  (Penta *) element; // TODO: add Basal-/SurfaceNodeIndices to element.h, and change this to Element*
	penta->BasalNodeIndices(&numindices,&indices,element->GetElementType());
	penta->SurfaceNodeIndices(&numindicesup,&indicesup,element->GetElementType());	_assert_(numindices==numindicesup);

	GaussPenta* gauss=new GaussPenta();
	GaussPenta* gaussup=new GaussPenta();

	for(int i=0;i<numindices;i++){
		gauss->GaussNode(element->GetElementType(),indices[i]);
		gaussup->GaussNode(element->GetElementType(),indicesup[i]);

		enthalpy_input->GetInputValue(&enthalpy,gauss);
		enthalpy_input->GetInputValue(&enthalpyup,gaussup);
		pressure_input->GetInputValue(&pressure,gauss);
		pressure_input->GetInputValue(&pressureup,gaussup);
		watercolumn_input->GetInputValue(&watercolumn,gauss);
		meltingrate_input->GetInputValue(&meltingrate,gauss);

		state=GetThermalBasalCondition(element, enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate);

		switch (state) {
			case 0:
				// cold, dry base: apply basal surface forcing
				vec_spc->SetValue(element->nodes[i]->Pid(),0.,INS_VAL);
				break;
			case 1:
				// cold, wet base: keep at pressure melting point
				vec_spc->SetValue(element->nodes[i]->Pid(),1.,INS_VAL);
				break;
			case 2:
				// temperate, thin refreezing base: release spc
				vec_spc->SetValue(element->nodes[i]->Pid(),0.,INS_VAL);
				break;
			case 3:
				// temperate, thin melting base: set spc
				vec_spc->SetValue(element->nodes[i]->Pid(),1.,INS_VAL);
				break;
			case 4:
				// temperate, thick melting base: set grad H*n=0
				vec_spc->SetValue(element->nodes[i]->Pid(),0.,INS_VAL);
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}

	}

	/*Free resources:*/
	xDelete<int>(indices);
	xDelete<int>(indicesup);
	delete gauss;
	delete gaussup;
}/*}}}*/
void           EnthalpyAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,EnthalpyEnum);
}/*}}}*/
int            EnthalpyAnalysis::GetThermalBasalCondition(Element* element, IssmDouble enthalpy, IssmDouble enthalpyup, IssmDouble pressure, IssmDouble pressureup, IssmDouble watercolumn, IssmDouble meltingrate){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return -1;

	/* Only update Constraints at the base of grounded ice*/
	if(!(element->IsOnBase())) return -1;


	/*Get parameters and inputs: */
	IssmDouble  dt = element->FindParam(TimesteppingTimeStepEnum);

	/*Determine basal state*/
	int state=-1;
	if(enthalpy<PureIceEnthalpy(element,pressure)){
		if(watercolumn<=0.) state=0; // cold, dry base
		else state=1; // cold, wet base (refreezing)
	}
	else{
		if(enthalpyup<PureIceEnthalpy(element,pressureup)){
			if((dt==0.) && (meltingrate<0.)) state=2;	// refreezing temperate base (non-physical, only for steadystate solver)
			else	state=3; // temperate base, but no temperate layer
		}
		else state=4; // temperate layer with positive thickness
	}

	_assert_(state>=0);
	return state;
}/*}}}*/
IssmDouble     EnthalpyAnalysis::GetWetIceConductivity(Element* element, IssmDouble enthalpy, IssmDouble pressure){/*{{{*/

	IssmDouble temperature, waterfraction;
	IssmDouble kappa_w = 0.6; // thermal conductivity of water (in W/m/K)
	IssmDouble kappa_i = element->FindParam(MaterialsThermalconductivityEnum);
	element->EnthalpyToThermal(&temperature, &waterfraction, enthalpy, pressure);

	return (1.-waterfraction)*kappa_i + waterfraction*kappa_w;
}/*}}}*/
void           EnthalpyAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           EnthalpyAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	bool        converged;
	int         i,rheology_law;
	IssmDouble  B_average,s_average,T_average=0.,P_average=0.;
	int        *doflist   = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values        = xNew<IssmDouble>(numnodes);
	IssmDouble* pressure      = xNew<IssmDouble>(numnodes);
	IssmDouble* surface       = xNew<IssmDouble>(numnodes);
	IssmDouble* B             = xNew<IssmDouble>(numnodes);
	IssmDouble* temperature   = xNew<IssmDouble>(numnodes);
	IssmDouble* waterfraction = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];

		/*Check solution*/
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Get all inputs and parameters*/
	element->GetInputValue(&converged,ConvergedEnum);
	element->GetInputListOnNodes(&pressure[0],PressureEnum);
	int finite_element = element->GetElementType(); if(finite_element==P1Enum) finite_element = P1DGEnum;
	if(converged){
		for(i=0;i<numnodes;i++){
			element->EnthalpyToThermal(&temperature[i],&waterfraction[i],values[i],pressure[i]);
			if(waterfraction[i]<0.) _error_("Negative water fraction found in solution vector");
			//if(waterfraction[i]>1.) _error_("Water fraction >1 found in solution vector");
		}
		element->AddInput(EnthalpyEnum,values,finite_element);
		element->AddInput(WaterfractionEnum,waterfraction,finite_element);
		element->AddInput(TemperatureEnum,temperature,finite_element);

		IssmDouble* n = xNew<IssmDouble>(numnodes);
		if(element->material->ObjectEnum()==MatestarEnum){
			for(i=0;i<numnodes;i++) n[i]=3.;
		}
		else{
			element->GetInputListOnNodes(&n[0],MaterialsRheologyNEnum);
		}

		/*Update Rheology only if converged (we must make sure that the temperature is below melting point
		 * otherwise the rheology could be negative*/
		element->FindParam(&rheology_law,MaterialsRheologyLawEnum);
		element->GetInputListOnNodes(&surface[0],SurfaceEnum);
		switch(rheology_law){
			case NoneEnum:
				/*Do nothing: B is not temperature dependent*/
				break;
			case BuddJackaEnum:
				for(i=0;i<numnodes;i++) B[i]=BuddJacka(temperature[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
			case CuffeyEnum:
				for(i=0;i<numnodes;i++) B[i]=Cuffey(temperature[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
			case CuffeyTemperateEnum:
				for(i=0;i<numnodes;i++) B[i]=CuffeyTemperate(temperature[i], waterfraction[i],n[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
			case PatersonEnum:
				for(i=0;i<numnodes;i++) B[i]=Paterson(temperature[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
			case NyeH2OEnum:
				for(i=0;i<numnodes;i++) B[i]=NyeH2O(values[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
			case NyeCO2Enum:
				for(i=0;i<numnodes;i++) B[i]=NyeCO2(values[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
			case ArrheniusEnum:{
				element->GetVerticesCoordinates(&xyz_list);
				for(i=0;i<numnodes;i++) B[i]=Arrhenius(temperature[i],surface[i]-xyz_list[i*3+2],n[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
				}
			case LliboutryDuvalEnum:{
				for(i=0;i<numnodes;i++) B[i]=LliboutryDuval(values[i],pressure[i],n[i],element->FindParam(MaterialsBetaEnum),element->FindParam(ConstantsReferencetemperatureEnum),element->FindParam(MaterialsHeatcapacityEnum),element->FindParam(MaterialsLatentheatEnum));
				element->AddInput(MaterialsRheologyBEnum,&B[0],finite_element);
				break;
				}
			default: _error_("Rheology law " << EnumToStringx(rheology_law) << " not supported yet");
		}
		xDelete<IssmDouble>(n);
	}
	else{
		element->AddInput(EnthalpyPicardEnum,values,finite_element);
	}

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(temperature);
	xDelete<IssmDouble>(waterfraction);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
}/*}}}*/
void           EnthalpyAnalysis::PostProcessing(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	bool computebasalmeltingrates=true;
	bool isdrainicecolumn;
	IssmDouble dt;

	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	femmodel->parameters->FindParam(&isdrainicecolumn,ThermalIsdrainicecolumnEnum);

	if(isdrainicecolumn){
		DrainWaterfraction(femmodel);
	}
	if(computebasalmeltingrates){
		ComputeBasalMeltingrate(femmodel);
	}

}/*}}}*/
IssmDouble     EnthalpyAnalysis::PureIceEnthalpy(Element* element,IssmDouble pressure){/*{{{*/

	IssmDouble heatcapacity         = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble referencetemperature = element->FindParam(ConstantsReferencetemperatureEnum);

	return heatcapacity*(TMeltingPoint(element,pressure)-referencetemperature);
}/*}}}*/
IssmDouble     EnthalpyAnalysis::TMeltingPoint(Element* element,IssmDouble pressure){/*{{{*/

	IssmDouble meltingpoint = element->FindParam(MaterialsMeltingpointEnum);
	IssmDouble beta         = element->FindParam(MaterialsBetaEnum);

	return meltingpoint-beta*pressure;
}/*}}}*/
void           EnthalpyAnalysis::UpdateBasalConstraints(FemModel* femmodel){/*{{{*/

	/*Update basal dirichlet BCs for enthalpy: */
	int numnodes            = femmodel->nodes->NumberOfNodes();
   int localmasters        = femmodel->nodes->NumberOfNodesLocal();
	Vector<IssmDouble>* spc = new Vector<IssmDouble>(localmasters,numnodes);

	/*First create a vector to figure out what elements should be constrained*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		GetBasalConstraints(spc,element);
	}

	/*Assemble*/
	spc->Assemble();

	/*Get local vector with both masters and slaves:*/
	IssmDouble *local_spc = NULL;
	femmodel->GetLocalVectorWithClonesNodes(&local_spc,spc);
	delete spc;

	/*Then update basal constraints nodes accordingly*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		ApplyBasalConstraints(local_spc,element);
	}

	femmodel->UpdateConstraintsx();

	/*Delete*/
	xDelete<IssmDouble>(local_spc);
}/*}}}*/
void           EnthalpyAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/
