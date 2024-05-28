#include "./ThermalAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void ThermalAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Intermediary*/
	int        M,N;
	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");
	_assert_(finiteelement==P1Enum);

	/*Output*/
	IssmDouble *spcvector  = NULL;
	IssmDouble *spcvectorstatic  = NULL;
	IssmDouble *issurface = NULL;

	/*Only 3d mesh supported*/
	if(iomodel->domaintype!=Domain3DEnum) _error_("not supported yet");

	/*Specific case for PDD, we want the constaints to be updated by the PDD scheme itself*/
	bool isdynamic = false;
	if(iomodel->solution_enum==ThermalSolutionEnum){
		/*No PDD scheme, keep default*/
	}
	else if(iomodel->solution_enum==SteadystateSolutionEnum){
		/*No PDD scheme, keep default*/
	}
	else if (iomodel->solution_enum==TransientSolutionEnum){
		int smb_model;
		iomodel->FindConstant(&smb_model,"md.smb.model");
		if(smb_model==SMBpddEnum)				isdynamic=true;
		if(smb_model==SMBd18opddEnum)			isdynamic=true;
		if(smb_model==SMBpddSicopolisEnum)  isdynamic=true;
	}
	else{
		_error_("Solution "<<EnumToStringx(iomodel->solution_enum)<<" not supported yet");
	}

	/*Fetch data: */
	iomodel->FetchData(&issurface,&M,&N,"md.mesh.vertexonsurface"); _assert_(N>0); _assert_(M==iomodel->numberofvertices);
	iomodel->FetchData(&spcvector,&M,&N,"md.thermal.spctemperature");
	iomodel->FetchData(&spcvectorstatic,&M,&N,"md.thermal.spctemperature");

	/*Convert spcs from temperatures to enthalpy*/
	_assert_(N>0); _assert_(M>=iomodel->numberofvertices);
	for(int i=0;i<iomodel->numberofvertices;i++){
		for(int j=0;j<N;j++){
			if (isdynamic){
				if (issurface[i]==1)spcvectorstatic[i*N+j] = NAN;
				else spcvector[i*N+j] = NAN;
			}
		}
	}

	if(isdynamic){
		IoModelToDynamicConstraintsx(constraints,iomodel,spcvector,iomodel->numberofvertices,1,ThermalAnalysisEnum,finiteelement);
		IoModelToConstraintsx(constraints,iomodel,spcvectorstatic,M,N,ThermalAnalysisEnum,finiteelement);
	}
	else{
		IoModelToConstraintsx(constraints,iomodel,spcvector,M,N,ThermalAnalysisEnum,finiteelement);
	}

	/*Free resources:*/
	iomodel->DeleteData(spcvector,"md.thermal.spctemperature");
	iomodel->DeleteData(spcvectorstatic,"md.thermal.spctemperature");
	iomodel->DeleteData(issurface,"md.mesh.vertexonsurface");

}/*}}}*/
void ThermalAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	if(iomodel->domaintype==Domain2DhorizontalEnum) _error_("2d meshes not supported yet");

	/*create penalties for nodes: no node can have a temperature over the melting point*/
	iomodel->FetchData(1,"md.thermal.spctemperature");
	CreateSingleNodeToElementConnectivity(iomodel);

	for(int i=0;i<iomodel->numberofvertices;i++){

		/*keep only this partition's nodes:*/
		if(iomodel->my_vertices[i]){
			if (xIsNan<IssmDouble>(iomodel->Data("md.thermal.spctemperature")[i])){ //No penalty applied on spc nodes!
				loads->AddObject(new Pengrid(i+1,i,iomodel));
			}
		}
	}
	iomodel->DeleteData(1,"md.thermal.spctemperature");

}/*}}}*/
void ThermalAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");
	_assert_(finiteelement==P1Enum);

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,ThermalAnalysisEnum,finiteelement);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  ThermalAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void ThermalAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int basalforcing_model,materialstype;

	/*Now, is the model 3d? otherwise, do nothing: */
	if(iomodel->domaintype==Domain2DhorizontalEnum)return;

	/*Update elements: */
	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");
	_assert_(finiteelement==P1Enum);
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	bool dakota_analysis, ismovingfront;
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");
	iomodel->FindConstant(&materialstype,"md.materials.type");

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.temperature",TemperatureEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VxMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VyMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VzMeshEnum);

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

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);
}/*}}}*/
void ThermalAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.maxiter",ThermalMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.stabilization",ThermalStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.penalty_factor",ThermalPenaltyFactorEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.penalty_threshold",ThermalPenaltyThresholdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.penalty_lock",ThermalPenaltyLockEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.isenthalpy",ThermalIsenthalpyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.isdynamicbasalspc",ThermalIsdynamicbasalspcEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.thermal.requested_outputs");
	parameters->AddObject(new IntParam(ThermalNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(ThermalRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.thermal.requested_outputs");

	/*Friction*/
	FrictionUpdateParameters(parameters, iomodel);
}/*}}}*/

/*Finite Element Analysis*/
void           ThermalAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           ThermalAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* ThermalAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* ThermalAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* ThermalAnalysis::CreateKMatrix(Element* element){/*{{{*/

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
ElementMatrix* ThermalAnalysis::CreateKMatrixShelf(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Initialize Element matrix and return if necessary*/
	if(!element->IsOnBase() || !element->IsAllFloating()) return NULL;

	IssmDouble  dt,Jdet,D;
	IssmDouble *xyz_list_base = NULL;

	/*Get basal element*/
	if(!element->IsOnBase() || !element->IsAllFloating()) return NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementMatrix* Ke    = element->NewElementMatrix();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble gravity             = element->FindParam(ConstantsGEnum);
	IssmDouble rho_water           = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble mixed_layer_capacity= element->FindParam(MaterialsMixedLayerCapacityEnum);
	IssmDouble thermal_exchange_vel= element->FindParam(MaterialsThermalExchangeVelocityEnum);

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
ElementMatrix* ThermalAnalysis::CreateKMatrixVolume(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         stabilization;
	IssmDouble  Jdet,dt,u,v,w,um,vm,wm,vel;
	IssmDouble  h,hx,hy,hz,vx,vy,vz,D_scalar;
	IssmDouble  tau_parameter,diameter;
	IssmDouble  tau_parameter_anisotropic[2],tau_parameter_hor,tau_parameter_ver;	
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(3*numnodes);
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
	IssmDouble  kappa = thermalconductivity/(rho_ice*heatcapacity);
	Input* vx_input  = element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input  = element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input  = element->GetInput(VzEnum);     _assert_(vz_input);
	Input* vxm_input = element->GetInput(VxMeshEnum); _assert_(vxm_input);
	Input* vym_input = element->GetInput(VyMeshEnum); _assert_(vym_input);
	Input* vzm_input = element->GetInput(VzMeshEnum); _assert_(vzm_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		D_scalar=gauss->weight*Jdet;
		if(dt!=0.) D_scalar=D_scalar*dt;

		/*Conduction: */
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D_scalar*kappa*(
							dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[2*numnodes+j]*dbasis[2*numnodes+i]
							);
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
			for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[j]*basis[i];
			D_scalar=D_scalar*dt;
		}

		/*Artifficial diffusivity*/
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
		else if(stabilization==2){
			diameter=element->MinEdgeLength(xyz_list);
			tau_parameter=element->StabilizationParameter(u-um,v-vm,w-wm,diameter,kappa);
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=tau_parameter*D_scalar*
					  ((u-um)*dbasis[0*numnodes+i]+(v-vm)*dbasis[1*numnodes+i]+(w-wm)*dbasis[2*numnodes+i])*
					  ((u-um)*dbasis[0*numnodes+j]+(v-vm)*dbasis[1*numnodes+j]+(w-wm)*dbasis[2*numnodes+j]);
				}
			}
			if(dt!=0.){
				D_scalar=gauss->weight*Jdet;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=tau_parameter*D_scalar*basis[j]*((u-um)*dbasis[0*numnodes+i]+(v-vm)*dbasis[1*numnodes+i]+(w-wm)*dbasis[2*numnodes+i]);
					}
				}
			}
		}
		/*anisotropic SUPG*/
		else if(stabilization==3){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			element->ElementSizes(&hx,&hy,&hz);
			element->StabilizationParameterAnisotropic(&tau_parameter_anisotropic[0],u-um,v-vm,w-wm,hx,hy,hz,kappa);
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
ElementVector* ThermalAnalysis::CreatePVector(Element* element){/*{{{*/

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
ElementVector* ThermalAnalysis::CreatePVectorSheet(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/* Geothermal flux on ice sheet base and basal friction */
	if(!element->IsOnBase() || element->IsAllFloating()) return NULL;

	IssmDouble  dt,Jdet,geothermalflux,vx,vy,vz;
	IssmDouble  alpha2,scalar,basalfriction,heatflux;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* vx_input             = element->GetInput(VxEnum);                          _assert_(vx_input);
	Input* vy_input             = element->GetInput(VyEnum);                          _assert_(vy_input);
	Input* vz_input             = element->GetInput(VzEnum);                          _assert_(vz_input);
	Input* geothermalflux_input = element->GetInput(BasalforcingsGeothermalfluxEnum); _assert_(geothermalflux_input);
	IssmDouble  rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,3);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss   = element->NewGaussBase(4);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		geothermalflux_input->GetInputValue(&geothermalflux,gauss);
		friction->GetAlpha2(&alpha2,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vz_input->GetInputValue(&vz,gauss);
		vz = 0.;//FIXME
		basalfriction = alpha2*(vx*vx + vy*vy + vz*vz);
		heatflux      = (basalfriction+geothermalflux)/(rho_ice*heatcapacity);

		scalar = gauss->weight*Jdet*heatflux;
		if(dt!=0.) scalar=dt*scalar;

		for(int i=0;i<numnodes;i++) pe->values[i]+=scalar*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;
}/*}}}*/
ElementVector* ThermalAnalysis::CreatePVectorShelf(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	IssmDouble  t_pmp,dt,Jdet,scalar_ocean,pressure;
	IssmDouble *xyz_list_base = NULL;

	/*Get basal element*/
	if(!element->IsOnBase() || !element->IsAllFloating()) return NULL;

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
		t_pmp=element->TMeltingPoint(pressure);

		scalar_ocean=gauss->weight*Jdet*rho_water*mixed_layer_capacity*thermal_exchange_vel*(t_pmp)/(heatcapacity*rho_ice);
		if(reCast<bool,IssmDouble>(dt)) scalar_ocean=dt*scalar_ocean;

		for(int i=0;i<numnodes;i++) pe->values[i]+=scalar_ocean*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;
}/*}}}*/
ElementVector* ThermalAnalysis::CreatePVectorVolume(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         stabilization;
	IssmDouble  Jdet,phi,dt;
	IssmDouble  temperature;
	IssmDouble  tau_parameter,diameter,hx,hy,hz;
	IssmDouble  tau_parameter_anisotropic[2],tau_parameter_hor,tau_parameter_ver;
	IssmDouble  u,v,w;
	IssmDouble  scalar_def,scalar_transient;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe     = element->NewElementVector();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(3*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble  rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble  thermalconductivity = element->FindParam(MaterialsThermalconductivityEnum);
	IssmDouble  kappa = thermalconductivity/(rho_ice*heatcapacity);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,ThermalStabilizationEnum);
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=element->GetInput(VzEnum); _assert_(vz_input);
	Input* temperature_input = NULL;
	if(reCast<bool,IssmDouble>(dt)){temperature_input = element->GetInput(TemperatureEnum); _assert_(temperature_input);}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->ViscousHeating(&phi,xyz_list,gauss,vx_input,vy_input,vz_input);

		scalar_def=phi/(rho_ice*heatcapacity)*Jdet*gauss->weight;
		if(reCast<bool,IssmDouble>(dt)) scalar_def=scalar_def*dt;

		for(int i=0;i<numnodes;i++) pe->values[i]+=scalar_def*basis[i];

		/* Build transient now */
		if(reCast<bool,IssmDouble>(dt)){
			temperature_input->GetInputValue(&temperature, gauss);
			scalar_transient=temperature*Jdet*gauss->weight;
			for(int i=0;i<numnodes;i++) pe->values[i]+=scalar_transient*basis[i];
		}

		if(stabilization==2){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			diameter=element->MinEdgeLength(xyz_list);
			vx_input->GetInputValue(&u,gauss);
			vy_input->GetInputValue(&v,gauss);
			vz_input->GetInputValue(&w,gauss);

			tau_parameter=element->StabilizationParameter(u,v,w,diameter,kappa);

			for(int i=0;i<numnodes;i++) pe->values[i]+=tau_parameter*scalar_def*(u*dbasis[0*numnodes+i]+v*dbasis[1*numnodes+i]+w*dbasis[2*numnodes+i]);
			if(reCast<bool,IssmDouble>(dt)){
				for(int i=0;i<numnodes;i++) pe->values[i]+=tau_parameter*scalar_transient*(u*dbasis[0*numnodes+i]+v*dbasis[1*numnodes+i]+w*dbasis[2*numnodes+i]);
			}
		}
		/* anisotropic SUPG */
		else if(stabilization==3){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			element->ElementSizes(&hx,&hy,&hz);
			vx_input->GetInputValue(&u,gauss);
			vy_input->GetInputValue(&v,gauss);
			vz_input->GetInputValue(&w,gauss);
			element->StabilizationParameterAnisotropic(&tau_parameter_anisotropic[0],u,v,w,hx,hy,hz,kappa);
			tau_parameter_hor=tau_parameter_anisotropic[0];
			tau_parameter_ver=tau_parameter_anisotropic[1];

			for(int i=0;i<numnodes;i++) pe->values[i]+=scalar_def*(tau_parameter_hor*u*dbasis[0*numnodes+i]+tau_parameter_hor*v*dbasis[1*numnodes+i]+tau_parameter_ver*w*dbasis[2*numnodes+i]);
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;

}/*}}}*/
void           ThermalAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,TemperatureEnum);
}/*}}}*/
void           ThermalAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           ThermalAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	bool        converged;
	int         i,rheology_law;
	int        *doflist   = NULL;
	IssmDouble *xyz_list  = NULL;
	bool        hack      = false;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values    = xNew<IssmDouble>(numnodes);
	IssmDouble* surface   = xNew<IssmDouble>(numnodes);
	IssmDouble* B         = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];

		/*Check solution*/
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
		//if(values[i]<0)      _printf_("temperature < 0°K found in solution vector\n");
		//if(values[i]>275)    _printf_("temperature > 275°K found in solution vector (Paterson's rheology associated is negative)\n");
	}

	/*Force temperature between [Tpmp-50 Tpmp] to disable penalties*/
	if(hack){
		IssmDouble* pressure = xNew<IssmDouble>(numnodes);
		element->GetInputListOnNodes(&pressure[0],PressureEnum);
		for(i=0;i<numnodes;i++){
			if(values[i]>element->TMeltingPoint(pressure[i]))     values[i]=element->TMeltingPoint(pressure[i]);
			if(values[i]<element->TMeltingPoint(pressure[i])-50.) values[i]=element->TMeltingPoint(pressure[i])-50.;
		}
		xDelete<IssmDouble>(pressure);
	}

	/*Get all inputs and parameters*/
	element->GetInputValue(&converged,ConvergedEnum);
	if(converged){
		element->AddInput(TemperatureEnum,values,element->GetElementType());

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
				for(i=0;i<numnodes;i++) B[i]=BuddJacka(values[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case CuffeyEnum:
				for(i=0;i<numnodes;i++) B[i]=Cuffey(values[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case PatersonEnum:
				for(i=0;i<numnodes;i++) B[i]=Paterson(values[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case NyeH2OEnum:
				for(i=0;i<numnodes;i++) B[i]=NyeH2O(values[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break; 
			case NyeCO2Enum:
				for(i=0;i<numnodes;i++) B[i]=NyeCO2(values[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case ArrheniusEnum:{
				element->GetVerticesCoordinates(&xyz_list);
				for(i=0;i<numnodes;i++) B[i]=Arrhenius(values[i],surface[i]-xyz_list[i*3+2],n[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
				}
			default:
				_error_("Rheology law " << EnumToStringx(rheology_law) << " not supported yet");
		}
		xDelete<IssmDouble>(n);
	}
	else{
		element->AddInput(TemperaturePicardEnum,values,element->GetElementType());
	}

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
}/*}}}*/
void           ThermalAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/
