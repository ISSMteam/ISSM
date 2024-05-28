#include "./AgeAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void AgeAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	IoModelToConstraintsx(constraints,iomodel,"md.age.spcage",AgeAnalysisEnum,P1Enum);

}/*}}}*/
void AgeAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	if(iomodel->domaintype==Domain2DhorizontalEnum) _error_("2d meshes not supported yet");

}/*}}}*/
void AgeAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int finiteelement = P1Enum;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,AgeAnalysisEnum,finiteelement);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  AgeAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void AgeAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Now, is the model 3d? otherwise, do nothing: */
	if(iomodel->domaintype==Domain2DhorizontalEnum)return;

	/*Update elements: */
	int finiteelement = P1Enum;
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	bool dakota_analysis, ismovingfront;
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.age",AgeEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VxMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VyMeshEnum);
	InputUpdateFromConstantx(inputs,elements,0.,VzMeshEnum);
}/*}}}*/
void AgeAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.age.stabilization",AgeStabilizationEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.age.requested_outputs");
	parameters->AddObject(new IntParam(AgeNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(AgeRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.age.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           AgeAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           AgeAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* AgeAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* AgeAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* AgeAnalysis::CreateKMatrix(Element* element){/*{{{*/

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
	element->FindParam(&stabilization,AgeStabilizationEnum);
	IssmDouble  rho_water           = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble  rho_ice             = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  gravity             = element->FindParam(ConstantsGEnum);
	IssmDouble  heatcapacity        = element->FindParam(MaterialsHeatcapacityEnum);
	IssmDouble  thermalconductivity = 1.;
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
ElementVector* AgeAnalysis::CreatePVector(Element* element){/*{{{*/

	_error_("STOP");

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         stabilization;
	IssmDouble  Jdet,dt;
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
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,AgeStabilizationEnum);
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

		scalar_def=1.*Jdet*gauss->weight;
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

			tau_parameter=element->StabilizationParameter(u,v,w,diameter,1.e-15); //assume very small conductivity to get tau

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
			element->StabilizationParameterAnisotropic(&tau_parameter_anisotropic[0],u,v,w,hx,hy,hz,1.e-15); //assume very small conductivity to get tau
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
void           AgeAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,AgeEnum);
}/*}}}*/
void           AgeAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           AgeAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	element->InputUpdateFromSolutionOneDof(solution,AgeEnum);

}/*}}}*/
void           AgeAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
	_error_("Should also automatically constrain surface/basal nodes where we have inflow");

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);

		if(element->IsOnSurface()){
			element = element->SpawnTopElement();
			int         numnodes  = element->GetNumberOfNodes();
			IssmDouble *mask      = xNew<IssmDouble>(numnodes);
			IssmDouble *bed       = xNew<IssmDouble>(numnodes);
			IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

			element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
			element->GetInputListOnNodes(&bed[0],BaseEnum);
			element->GetInputListOnNodes(&ls_active[0],IceMaskNodeActivationEnum);

			for(int in=0;in<numnodes;in++){
				Node* node=element->GetNode(in);
				if(mask[in]<0. && ls_active[in]==1.){
					node->Activate();
				}
				else{
					node->Deactivate();
					node->ApplyConstraint(0,bed[in]);
				}
			}
			xDelete<IssmDouble>(mask);
			xDelete<IssmDouble>(bed);
			xDelete<IssmDouble>(ls_active);
		}
		else if(element->IsOnBase()){
			element = element->SpawnBasalElement();
			_error_("not supported");
		}
	}
}/*}}}*/
