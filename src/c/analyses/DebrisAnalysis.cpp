#include "./DebrisAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/TransientInput.h"
#include "../solutionsequences/solutionsequences.h"
#include "../cores/cores.h"
#include <math.h>

#define FINITEELEMENT P1Enum
#define EPS 1e-14

/*Model processing*/
void DebrisAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*add constraints*/
	IoModelToConstraintsx(constraints,iomodel,"md.debris.spcthickness",DebrisAnalysisEnum,FINITEELEMENT);

}/*}}}*/
void DebrisAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int penpair_ids[2];
	int count=0;
	int numvertex_pairing;

	/*Create Penpair for vertex_pairing: */
	IssmDouble *vertex_pairing=NULL;
	IssmDouble *nodeonsurface=NULL;
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.debris.vertex_pairing");
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(&nodeonsurface,NULL,NULL,"md.mesh.vertexonsurface");
	for(int i=0;i<numvertex_pairing;i++){

		if(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+0])-1]){

			/*In debugging mode, check that the second node is in the same cpu*/
			_assert_(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+1])-1]);

			/*Skip if one of the two is not on the bed*/
			if(iomodel->domaintype!=Domain2DhorizontalEnum){
				if(!(reCast<bool>(nodeonsurface[reCast<int>(vertex_pairing[2*i+0])-1])) || !(reCast<bool>(nodeonsurface[reCast<int>(vertex_pairing[2*i+1])-1]))) continue;
			}

			/*Get node ids*/
			penpair_ids[0]=reCast<int>(vertex_pairing[2*i+0]);
			penpair_ids[1]=reCast<int>(vertex_pairing[2*i+1]);

			/*Create Load*/
			loads->AddObject(new Penpair(count+1,&penpair_ids[0]));
			count++;
		}
	}

	/*free resources: */
	iomodel->DeleteData(vertex_pairing,"md.debris.vertex_pairing");
	iomodel->DeleteData(nodeonsurface,"md.mesh.vertexonsurface");
}/*}}}*/
void DebrisAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,DebrisAnalysisEnum,FINITEELEMENT,isamr);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");

}/*}}}*/
int  DebrisAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void DebrisAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Now, is the model 3d? otherwise, do nothing: */
	if (iomodel->domaintype==Domain2DhorizontalEnum)return;

	int smb_model;

	/*Fetch data needed: */
	iomodel->FindConstant(&smb_model,"md.smb.model");

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,FINITEELEMENT);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.initialization.debris",DebrisThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxDebrisEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyDebrisEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
	}

	switch(smb_model){
		case SMBforcingEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.mass_balance",SmbMassBalanceEnum,0.);
			break;
		default:
			/*Nothing for now*/
			;
	}

}/*}}}*/
void DebrisAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.debris.stabilization",DebrisStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debris.removalmodel",DebrisRemovalmodelEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debris.displacementmodel",DebrisDisplacementmodelEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debris.min_thickness",DebrisMinThicknessEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debris.packingfraction",DebrisPackingFractionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debris.removal_slope_threshold",DebrisRemovalSlopeThresholdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debris.removal_stress_threshold",DebrisRemovalStressThresholdEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.debris.requested_outputs");
	parameters->AddObject(new IntParam(DebrisNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(DebrisRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.debris.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           DebrisAnalysis::Core(FemModel* femmodel){/*{{{*/

	//PreProcessing(femmodel);
	//femmodel->parameters->SetParam(VxDebrisEnum,InputToExtrudeEnum);
	//extrudefromtop_core(femmodel);
	femmodel->SetCurrentConfiguration(DebrisAnalysisEnum);        
	solutionsequence_linear(femmodel);
	PostProcessing(femmodel);

}/*}}}*/
void           DebrisAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* DebrisAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* DebrisAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* DebrisAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int        stabilization,domaintype,dim;
	IssmDouble Jdet,D_scalar,dt,h;
	IssmDouble vel,vx,vy,dvxdx,dvydy;
	IssmDouble yts=31536000.;
	IssmDouble tau;
	IssmDouble dvx[2],dvy[2];
	Element*    topelement = NULL;
	IssmDouble* xyz_list = NULL;

	/*Get top element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			topelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = topelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = topelement->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*   	basis  = xNew<IssmDouble>(numnodes);
	IssmDouble* 	dbasis = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*	D      = xNew<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	topelement->GetVerticesCoordinates(&xyz_list);
	topelement->FindParam(&dt,TimesteppingTimeStepEnum);
	topelement->FindParam(&stabilization,DebrisStabilizationEnum);
	Input* vx_input=topelement->GetInput(VxDebrisEnum); _assert_(vx_input);
	Input* vy_input=NULL;
	if(dim>1){vy_input = topelement->GetInput(VyDebrisEnum); _assert_(vy_input);}
	h=topelement->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	while(gauss->next()){

		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		topelement->NodalFunctions(basis,gauss);
		topelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		vx_input->GetInputValue(&vx,gauss);
		if(dim==2) vy_input->GetInputValue(&vy,gauss);

		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j];

		/*Advection terms*/
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		D_scalar=dt*gauss->weight*Jdet;
		if(dim==2){
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			dvxdx=dvx[0];
			dvydy=dvy[1];
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*\phi_i \phi_j \nabla\cdot v*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy);
					/*\phi_i v\cdot\nabla\phi_j*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
				}
			}
		}else{
			dvxdx=dvx[0];
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += D_scalar*dvxdx*basis[i]*basis[j];
					Ke->values[i*numnodes+j] += D_scalar*vx*dbasis[0*numnodes+j]*basis[i];
				}
			}
		}

		IssmDouble rho;
		if(FINITEELEMENT==P1Enum){
			rho=2.;
		}else if(FINITEELEMENT==P2Enum){
			rho=4.;
		}

		for(int i=0;i<(dim*dim);i++) D[i]=0.;
		if(stabilization==1){
			/*SSA*/
			if(dim==1){
				vx_input->GetInputValue(&vx,gauss);
				D[0]=h/rho*fabs(vx);
			}else{
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				vel=sqrt(vx*vx+vy*vy);
				D[0*dim+0]=h/rho*fabs(vx);
				D[1*dim+1]=h/rho*fabs(vy);
			}
		}else if(stabilization==2){  
			/*Streamline upwinding*/
			if(dim==1){
				vx_input->GetInputValue(&vx,gauss);
				vel=fabs(vx)+EPS;
				D[0] = h/(rho*vel)*vx*vx;
			}else{
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				vel=sqrt(vx*vx+vy*vy)+EPS;
				D[0*dim+0]=h/(rho*vel)*vx*vx;
				D[1*dim+0]=h/(rho*vel)*vy*vx;
				D[0*dim+1]=h/(rho*vel)*vx*vy;
				D[1*dim+1]=h/(rho*vel)*vy*vy;
			}		
		}else if(stabilization==3){  
			/*SUPG*/
			if(dim==1){
				vx_input->GetInputValue(&vx,gauss);
				tau=h/(rho*max(fabs(vx),EPS));
			}else{
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				tau=h/(rho*sqrt(vx*vx+vy*vy)+EPS);
			}
		}

		if(stabilization==1 || stabilization==2){
			for(int i=0;i<dim*dim;i++) D[i]=D_scalar*D[i];
			if(dim==2){
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += (
								dbasis[0*numnodes+i] *(D[0*dim+0]*dbasis[0*numnodes+j] + D[0*dim+1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1*dim+0]*dbasis[0*numnodes+j] + D[1*dim+1]*dbasis[1*numnodes+j]));
					}   
				}
			}else{
				for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += dbasis[0*numnodes+i]*D[0]*dbasis[0*numnodes+j];
			}
		}else if(stabilization==3){ 
			if(dim==1){
				/*Mass matrix - part 2*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=gauss->weight*Jdet*tau*basis[j]*vx*dbasis[0*numnodes+i];
					}
				}
				/*Mass matrix - part 3*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=gauss->weight*Jdet*tau*basis[j]*basis[i]*dvxdx;
					}
				}

				/*Advection matrix - part 2, A*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(vx*dbasis[0*numnodes+j])*(vx*dbasis[0*numnodes+i]);
					}
				}

				/*Advection matrix - part 3, A*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(vx*dbasis[0*numnodes+j])*(basis[i]*dvxdx);
					}
				}

				/*Advection matrix - part 2, B*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(basis[j]*dvxdx)*(vx*dbasis[0*numnodes+i]);
					}
				}

				/*Advection matrix - part 3, B*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(basis[j]*dvxdx)*(basis[i]*dvxdx);;
					}
				}
			}else if(dim==2){
				/*Mass matrix - part 2*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=gauss->weight*Jdet*tau*basis[j]*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
					}
				}
				/*Mass matrix - part 3*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=gauss->weight*Jdet*tau*basis[j]*(basis[i]*dvxdx+basis[i]*dvydy);
					}
				}

				/*Advection matrix - part 2, A*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j])*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
					}
				}
				/*Advection matrix - part 3, A*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j])*(basis[i]*dvxdx+basis[i]*dvydy);
					}
				}

				/*Advection matrix - part 2, B*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(basis[j]*dvxdx+basis[j]*dvydy)*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
					}
				}
				/*Advection matrix - part 3, B*/
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=dt*gauss->weight*Jdet*tau*(basis[j]*dvxdx+basis[j]*dvydy)*(basis[i]*dvxdx+basis[i]*dvydy);
					}
				}
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	return Ke;
}/*}}}*/
ElementVector* DebrisAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int	stabilization,dim,domaintype;
	IssmDouble  Jdet,dt;
	IssmDouble  smb,thickness,psi;
	IssmDouble  vx,vy,vel,dvxdx,dvydy,h,tau,pf;
	IssmDouble yts=31536000.;
	IssmDouble  dvx[2],dvy[2];
	IssmDouble* xyz_list = NULL;
	Element*    topelement = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			topelement = element;
			dim = 2; 
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 1; 
			break;
		case Domain3DEnum:           
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 2; break;
		default: 
			_error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = topelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = topelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis= xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	topelement->GetVerticesCoordinates(&xyz_list);
	topelement->FindParam(&dt,TimesteppingTimeStepEnum);
	topelement->FindParam(&pf,DebrisPackingFractionEnum);
	topelement->FindParam(&stabilization,DebrisStabilizationEnum);
	Input* smb_input      	= topelement->GetInput(SmbMassBalanceEnum);  _assert_(smb_input);
	Input* thickness_input  = topelement->GetInput(DebrisThicknessEnum); _assert_(thickness_input);
	Input* vx_input  = topelement->GetInput(VxDebrisEnum);  _assert_(vx_input);
	Input* vy_input=NULL;
	if(dim==2){
		vy_input=topelement->GetInput(VyDebrisEnum); _assert_(vy_input);
	}
	h=topelement->CharacteristicLength();

	IssmDouble rho;
	if(FINITEELEMENT==P1Enum){
		rho=2.;
	}else if(FINITEELEMENT==P2Enum){
		rho=4.;
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	while(gauss->next()){

		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		topelement->NodalFunctions(basis,gauss);

		smb_input->GetInputValue(&smb,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		if(smb>0.){
			psi=thickness-0.*dt*smb*pf;
		}else{
			psi=thickness-dt*smb*pf;
		}

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*psi*basis[i]; 

		if(stabilization==3){
			/*SUPG*/
			topelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			dvxdx=dvx[0];
			if(dim==1){
				vx_input->GetInputValue(&vx,gauss);
				tau=h/(rho*max(fabs(vx),EPS));

				/*Force vector - part 2*/
				for(int i=0;i<numnodes;i++){
					pe->values[i]+=Jdet*gauss->weight*psi*tau*vx*dbasis[0*numnodes+i];
				}
				/*Force vector - part 3*/
				for(int i=0;i<numnodes;i++){
					pe->values[i]+=Jdet*gauss->weight*psi*tau*basis[i]*dvxdx;
				}

			}else if(dim==2){
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
				vel=sqrt(vx*vx+vy*vy);
				dvydy=dvy[1];
				tau=h/(rho*vel+EPS);

				/*Force vector - part 2*/
				for(int i=0;i<numnodes;i++){
					pe->values[i]+=Jdet*gauss->weight*psi*tau*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
				}
				/*Force vector - part 3*/
				for(int i=0;i<numnodes;i++){
					pe->values[i]+=Jdet*gauss->weight*psi*tau*(basis[i]*dvxdx+basis[i]*dvydy);
				}
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	return pe;
}/*}}}*/
void           DebrisAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           DebrisAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           DebrisAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int *ddoflist=NULL;

	int numnodes = element->GetNumberOfNodes();
	IssmDouble* newthickness  = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	IssmDouble minthickness = element->FindParam(DebrisMinThicknessEnum);
	element->GetDofListLocal(&ddoflist,NoneApproximationEnum,GsetEnum);

	for(int i=0;i<numnodes;i++){
		newthickness[i] = solution[ddoflist[i]];
		if(xIsNan<IssmDouble>(newthickness[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(newthickness[i])) _error_("Inf found in solution vector");

		// check for thickness<minthickness
		if(newthickness[i]<minthickness) newthickness[i]=minthickness;
	}

	// update inputs
	element->AddInput(DebrisThicknessEnum,newthickness,P1Enum);

	// Free resources
	xDelete<IssmDouble>(newthickness);
	xDelete<int>(ddoflist);
	//*/
}/*}}}*/
void           DebrisAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	//SetActiveNodesLSMx(femmodel);

	// Update active elements based on ice levelset and ocean levelset*/
	GetMaskOfIceVerticesLSMx(femmodel,false,true); //FIXME?
	SetActiveNodesLSMx(femmodel,false,true); //FIXME?

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);
		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&ls_active[0],DebrisMaskNodeActivationEnum);

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]>0. && ls_active[in]==1.){
				// Do nothing
				node->Activate(); //Not sure if we need this!
			}
			else{
				IssmDouble phi=0;
				node->Deactivate();// Not sure if we need this
				node->ApplyConstraint(0,phi);
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(ls_active);
	}
	//*/
	return;	
}/*}}}*/
void           DebrisAnalysis::PostProcessing(FemModel* femmodel){/*{{{*/

	if(VerboseSolution()) _printf0_("   Debris postprocessing\n");

	/*Intermediaries*/
	int removalmodel;
	int k,numnodes;
	int domaintype,dim;
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&removalmodel,DebrisRemovalmodelEnum);
	Element* element= NULL;
	Element*    topelement = NULL;

	if(removalmodel==0){
		// no removal, do nothing
	}else{
		// slope or driving stress removal

		for(Object* & object : femmodel->elements->objects){
			element=xDynamicCast<Element*>(object);

			numnodes=element->GetNumberOfNodes();
			element->FindParam(&domaintype,DomainTypeEnum);
			IssmDouble* icethickness   = xNew<IssmDouble>(numnodes);
			IssmDouble* debristhickness= xNew<IssmDouble>(numnodes);
			IssmDouble* slopex	   = xNew<IssmDouble>(numnodes);
			IssmDouble* slopey	   = xNew<IssmDouble>(numnodes); 
			IssmDouble* onsurface	   = xNew<IssmDouble>(numnodes); 
			IssmDouble* ls_active      = xNew<IssmDouble>(numnodes); 
			element->GetInputListOnNodes(debristhickness,DebrisThicknessEnum);
			element->GetInputListOnNodes(icethickness,ThicknessEnum);
			element->GetInputListOnNodes(onsurface,MeshVertexonsurfaceEnum);
			element->GetInputListOnNodes(ls_active,DebrisMaskNodeActivationEnum);

			dim=1;
			element->GetInputListOnNodes(slopex,SurfaceSlopeXEnum);
			if(domaintype!=Domain2DverticalEnum){
				element->GetInputListOnNodes(slopey,SurfaceSlopeYEnum);
				dim=2;
			}
			IssmDouble slope,rad2deg=180./M_PI; //=57.2958
			bool isminthicknessinelement=false;
			bool remove_debris=false;
			bool isactive=false;

			IssmDouble iceminthickness=element->FindParam(MasstransportMinThicknessEnum);                        

			switch(removalmodel){
				case 1:{
					       IssmDouble slope_threshold=element->FindParam(DebrisRemovalSlopeThresholdEnum);
					       int kk=0;
					       for(k=0; k<numnodes;k++){
						       if(icethickness[k]<=(iceminthickness+0.00001)) isminthicknessinelement=true;
						       if(icethickness[k]<=(iceminthickness+0.00001)) kk++;
					       }
					       isminthicknessinelement=true;
					       if(kk<numnodes && isminthicknessinelement){
						       for(k=0; k<numnodes;k++){
							       slope=fabs(slopex[k]);
							       if(dim==2) slope=pow(pow(slopex[k],2)+pow(slopey[k],2),0.5);
							       //slope_mean=slope_mean+slope;
							       if((atan(slope)*rad2deg)>slope_threshold) remove_debris=true;
							       //if((atan(slope)*rad2deg)>slope_threshold) debristhickness[k]=0.;
						       }
						       //if((atan(slope_mean)*rad2deg)>slope_threshold) remove_debris=true;
						       if(remove_debris){
							       for(k=0; k<numnodes;k++){
								       debristhickness[k]=0.;
							       }
						       }
					       }
					       element->AddInput(DebrisThicknessEnum,debristhickness,P1Enum);

					       xDelete<IssmDouble>(debristhickness);
					       xDelete<IssmDouble>(icethickness);
					       xDelete<IssmDouble>(slopex);
					       xDelete<IssmDouble>(slopey);
					       break;
				       }
				case 2:{
					       IssmDouble stress_threshold = element->FindParam(DebrisRemovalStressThresholdEnum);
					       IssmDouble gravity = element->FindParam(ConstantsGEnum);
					       IssmDouble stress,rhod=1900.;
					       int kk=0;
					       for(k=0; k<numnodes;k++){
						       if(icethickness[k]<=(iceminthickness+0.00001)) isminthicknessinelement=true;
						       if(icethickness[k]<=(iceminthickness+0.00001)) kk++;
					       }
					       isminthicknessinelement=true;
					       if(kk<numnodes && isminthicknessinelement){
						       //stress=0;
						       IssmDouble stress_sum=0.;
						       for(k=0; k<numnodes;k++){
							       slope=fabs(slopex[k]);
							       if(dim==2) slope=pow(pow(slopex[k],2)+pow(slopey[k],2),0.5);
							       stress=rhod*gravity*debristhickness[k]*slope;//pow(slope*slope/(slope*slope+1),0.5);//sin(slope/rad2deg);
							       //stress_sum=stress_sum+stress;
							       if(stress>stress_threshold) remove_debris=true;
							       //if(stress>stress_threshold) debristhickness[k]=0.;
						       }
						       //if((stress_sum/double(kk))>stress_threshold) remove_debris=true;
						       if(remove_debris){
							       for(k=0; k<numnodes;k++){
								       debristhickness[k]=0.;
							       }
						       }
					       }
					       element->AddInput(DebrisThicknessEnum,debristhickness,P1Enum);

					       xDelete<IssmDouble>(debristhickness);
					       xDelete<IssmDouble>(icethickness);
					       xDelete<IssmDouble>(slopex);
					       xDelete<IssmDouble>(slopey);
					       xDelete<IssmDouble>(ls_active);
					       break;
				       }
				default: _error_("removalmodel "<<EnumToStringx(removalmodel)<<" not implemented yet");
			}
		}
	}

}/*}}}*/
void DebrisAnalysis::PreProcessing(FemModel* femmodel){/*{{{*/

	if(VerboseSolution()) _printf0_("   Debris preprocessing\n");

	Element* element= NULL;
	for(Object* & object : femmodel->elements->objects){
		element=xDynamicCast<Element*>(object);

		int numvertices = element->GetNumberOfVertices();

		IssmDouble* vx = xNew<IssmDouble>(numvertices);
		IssmDouble* debristhickness = xNew<IssmDouble>(numvertices);
		IssmDouble* slopex         = xNew<IssmDouble>(numvertices);
		IssmDouble* onsurface      = xNew<IssmDouble>(numvertices);
		IssmDouble* icethickness      = xNew<IssmDouble>(numvertices);

		element->GetInputListOnVertices(&debristhickness[0],DebrisThicknessEnum);
		element->GetInputListOnVertices(&vx[0],VxDebrisEnum);
		element->GetInputListOnVertices(&slopex[0],SurfaceSlopeXEnum);
		element->GetInputListOnVertices(&onsurface[0],MeshVertexonsurfaceEnum);
		element->GetInputListOnVertices(&icethickness[0],ThicknessEnum);

		IssmDouble slope,rad2deg=180./M_PI; //=57.2958
		IssmDouble vslipx,rhod=1900.;
		IssmDouble gravity=element->FindParam(ConstantsGEnum);
		IssmDouble slope_threshold=element->FindParam(DebrisRemovalSlopeThresholdEnum);
		IssmDouble iceminthickness=element->FindParam(MasstransportMinThicknessEnum);

		int step;
		IssmDouble dt, maxv;
		IssmDouble yts=31536000.;
		femmodel->parameters->FindParam(&step,StepEnum);
		femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

		bool isminthicknessinelement=false;
		for(int i=0;i<numvertices;i++){
			if(icethickness[i]<=(iceminthickness+0.01)) isminthicknessinelement=true;
		}
		if(isminthicknessinelement){
			//do nothing
			for(int i=0;i<numvertices;i++){
				if(icethickness[i]<=(iceminthickness+0.01)) vx[i]=0.;
			}
		}else{
			for(int i=0;i<numvertices;i++){
				//if(onsurface[i]>.5){
				slope=fabs(slopex[i]);
				//if((atan(slope)*rad2deg)>25.){
				//if(debristhickness[i]>0.01){
				vslipx=1.0/yts;
				//maxv=10.0/2./dt;
				//vslipx=-slope_threshold*rhod*gravity*debristhickness[i]*slopex[i]/yts;
				vx[i]=vx[i]+vslipx;
				//debristhickness[i]=debristhickness[i];
				//if(vx[i]>maxv) vx[i]=maxv;
				//}
				//} 
				//}
			}
		}
		//if(step%100==0)   
		element->AddInput(VxDebrisEnum,vx,P1Enum);

		/* Free resources */
		xDelete<IssmDouble>(debristhickness);
		xDelete<IssmDouble>(icethickness);
		xDelete<IssmDouble>(vx);
		xDelete<IssmDouble>(slopex);
		xDelete<IssmDouble>(onsurface);
	}

}/*}}}*/
