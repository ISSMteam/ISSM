#include "./BalancethicknessAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/DatasetInput.h"

/*Model processing*/
void BalancethicknessAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int    stabilization;	
	iomodel->FindConstant(&stabilization,"md.balancethickness.stabilization");

	/*Do not add constraints in DG*/
	if(stabilization!=3){
		IoModelToConstraintsx(constraints,iomodel,"md.balancethickness.spcthickness",BalancethicknessAnalysisEnum,P1Enum);
	}

}/*}}}*/
void BalancethicknessAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediary*/
	int element;
	int stabilization;

	/*Fetch parameters: */
	iomodel->FindConstant(&stabilization,"md.balancethickness.stabilization");

	/*Loads only in DG*/
	if (stabilization==3){

		/*Get faces and elements*/
		CreateFaces(iomodel);
		iomodel->FetchData(1,"md.geometry.thickness");

		/*First load data:*/
		for(int i=0;i<iomodel->numberoffaces;i++){

			/*Get left and right elements*/
			element=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]

			/*Now, if this element is not in the partition, pass: */
			if(!iomodel->my_elements[element]) continue;

			/* Add load */
			loads->AddObject(new Numericalflux(i+1,i,i,iomodel));
		}

		/*Free data: */
		iomodel->DeleteData(1,"md.geometry.thickness");
	}
}/*}}}*/
void BalancethicknessAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int  stabilization;
	iomodel->FindConstant(&stabilization,"md.balancethickness.stabilization");

	/*Check in 3d*/
	if(stabilization==3 && iomodel->domaintype==Domain3DEnum) _error_("DG 3d not implemented yet");

	/*First fetch data: */
	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	if(stabilization!=3){
		::CreateNodes(nodes,iomodel,BalancethicknessAnalysisEnum,P1Enum);
	}
	else{
		::CreateNodes(nodes,iomodel,BalancethicknessAnalysisEnum,P1DGEnum);
	}
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  BalancethicknessAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void BalancethicknessAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    stabilization,finiteelement;

	/*Fetch data needed: */
	iomodel->FindConstant(&stabilization,"md.balancethickness.stabilization");

	/*Finite element type*/
	finiteelement = P1Enum;
	if(stabilization==3){
		finiteelement = P1DGEnum;
	}

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.smb.mass_balance",SmbMassBalanceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.balancethickness.thickening_rate",BalancethicknessThickeningRateEnum);

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}/*}}}*/
void BalancethicknessAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	parameters->AddObject(iomodel->CopyConstantObject("md.balancethickness.stabilization",BalancethicknessStabilizationEnum));
}/*}}}*/

/*Finite Element Analysis*/
void           BalancethicknessAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BalancethicknessAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BalancethicknessAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BalancethicknessAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BalancethicknessAnalysis::CreateKMatrix(Element* element){/*{{{*/

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementMatrix* Ke = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			Ke = CreateKMatrixCG(basalelement);
			break;
		case P1DGEnum:
			Ke = CreateKMatrixDG(basalelement);
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementMatrix* BalancethicknessAnalysis::CreateKMatrixCG(Element* element){/*{{{*/

	/*Intermediaries */
	int        stabilization;
	int        domaintype;
	IssmDouble Jdet,D_scalar,h;
	IssmDouble vel,vx,vy,dvxdx,dvydy;
	IssmDouble dvx[2],dvy[2];
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*		dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble     D[2][2];

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&stabilization,BalancethicknessStabilizationEnum);
	Input* vxaverage_input=NULL;
	Input* vyaverage_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vxaverage_input=element->GetInput(VxEnum); _assert_(vxaverage_input);
		vyaverage_input=element->GetInput(VyEnum); _assert_(vyaverage_input);
	}
	else{
		vxaverage_input=element->GetInput(VxAverageEnum); _assert_(vxaverage_input);
		vyaverage_input=element->GetInput(VyAverageEnum); _assert_(vyaverage_input);
	}
	h = element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		vxaverage_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vyaverage_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		dvxdx=dvx[0];
		dvydy=dvy[1];
		D_scalar=gauss->weight*Jdet;

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/*\phi_i \phi_j \nabla\cdot v*/
				Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy);
				/*\phi_i v\cdot\nabla\phi_j*/
				Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
			}
		}

		if(stabilization==1){
			/*Streamline upwinding*/
			vel=sqrt(vx*vx+vy*vy);
			D[0][0]=h/(2*vel)*vx*vx;
			D[1][0]=h/(2*vel)*vy*vx;
			D[0][1]=h/(2*vel)*vx*vy;
			D[1][1]=h/(2*vel)*vy*vy;
		}
		else if(stabilization==2){
			/*SSA*/
			vxaverage_input->GetInputAverage(&vx);
			vyaverage_input->GetInputAverage(&vy);
			D[0][0]=h/2.0*fabs(vx);
			D[0][1]=0.;
			D[1][0]=0.;
			D[1][1]=h/2.0*fabs(vy);
		}
		if(stabilization==1 || stabilization==2){
			D[0][0]=D_scalar*D[0][0];
			D[1][0]=D_scalar*D[1][0];
			D[0][1]=D_scalar*D[0][1];
			D[1][1]=D_scalar*D[1][1];

			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += (
								dbasis[0*numnodes+i] *(D[0][0]*dbasis[0*numnodes+j] + D[0][1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1][0]*dbasis[0*numnodes+j] + D[1][1]*dbasis[1*numnodes+j]) 
								);
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
ElementMatrix* BalancethicknessAnalysis::CreateKMatrixDG(Element* element){/*{{{*/

	/*Intermediaries */
	int        domaintype;
	IssmDouble Jdet,D_scalar,vx,vy,dvxdx,dvydy,vel;
	IssmDouble dvx[2],dvy[2];
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*		dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble     D[2][2];

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&domaintype,DomainTypeEnum);
	Input* vxaverage_input=NULL;
	Input* vyaverage_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vxaverage_input=element->GetInput(VxEnum); _assert_(vxaverage_input);
		vyaverage_input=element->GetInput(VyEnum); _assert_(vyaverage_input);
	}
	else{
		vxaverage_input=element->GetInput(VxAverageEnum); _assert_(vxaverage_input);
		vyaverage_input=element->GetInput(VyAverageEnum); _assert_(vyaverage_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		vxaverage_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vyaverage_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		D_scalar=gauss->weight*Jdet;

		/*WARNING: inverted compared to CG*/
		D_scalar = - gauss->weight*Jdet;
		D[0][0]  = D_scalar*vx;
		D[0][1]  = 0.;
		D[1][0]  = 0.;
		D[1][1]  = D_scalar*vy;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j]*(dvxdx+dvydy);
				Ke->values[i*numnodes+j] += D_scalar*basis[j]*(vx*dbasis[0*numnodes+i] + vy*dbasis[1*numnodes+i]);
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
ElementVector* BalancethicknessAnalysis::CreatePVector(Element* element){/*{{{*/

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementVector* pe = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			pe = CreatePVectorCG(basalelement);
			break;
		case P1DGEnum:
			pe = CreatePVectorDG(basalelement);
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
ElementVector* BalancethicknessAnalysis::CreatePVectorCG(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  dhdt,mb,ms,Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* mb_input   = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);       _assert_(mb_input);
	Input* ms_input   = element->GetInput(SmbMassBalanceEnum);     _assert_(ms_input);
	Input* dhdt_input = element->GetInput(BalancethicknessThickeningRateEnum); _assert_(dhdt_input);

	/*Initialize mb_correction to 0, do not forget!:*/
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		ms_input->GetInputValue(&ms,gauss);
		mb_input->GetInputValue(&mb,gauss);
		dhdt_input->GetInputValue(&dhdt,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(ms-mb-dhdt)*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* BalancethicknessAnalysis::CreatePVectorDG(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  dhdt,mb,ms,Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* mb_input   = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);       _assert_(mb_input);
	Input* ms_input   = element->GetInput(SmbMassBalanceEnum);     _assert_(ms_input);
	Input* dhdt_input = element->GetInput(BalancethicknessThickeningRateEnum); _assert_(dhdt_input);

	/*Initialize mb_correction to 0, do not forget!:*/
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		ms_input->GetInputValue(&ms,gauss);
		mb_input->GetInputValue(&mb,gauss);
		dhdt_input->GetInputValue(&dhdt,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(ms-mb-dhdt)*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           BalancethicknessAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           BalancethicknessAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	/* WARNING: this gradient is valid for Soft balance thickness only */

	/*If on water, grad = 0: */
	if(!element->IsIceInElement()) return;

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble thickness,thicknessobs,dH[3],dp[3];
	IssmDouble  vx,vy,vel,dvx[2],dvy[2],dhdt,basal_melting,surface_mass_balance;
	IssmDouble *xyz_list= NULL;

	/*Get list of cost functions*/
	int *responses = NULL;
	int  num_responses,resp,solution;
	element->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	element->FindParam(&responses,NULL,InversionCostFunctionsEnum);
	element->FindParam(&solution,SolutionTypeEnum);
	if(solution!=BalancethicknessSoftSolutionEnum) _error_("not implemented yet");
	if(control_type!=ThicknessEnum)                _error_("Control "<<EnumToStringx(control_type)<<" not supported");

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* dbasis        = xNew<IssmDouble>(2*numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	DatasetInput* weights_input       = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum);  _assert_(weights_input);
	Input* thickness_input            = element->GetInput(ThicknessEnum);                           _assert_(thickness_input);
	Input* thicknessobs_input         = element->GetInput(InversionThicknessObsEnum);               _assert_(thicknessobs_input);
	Input* vx_input                   = element->GetInput(VxEnum);                                  _assert_(vx_input);
	Input* vy_input                   = element->GetInput(VyEnum);                                  _assert_(vy_input);
	Input* surface_mass_balance_input = element->GetInput(SmbMassBalanceEnum);          _assert_(surface_mass_balance_input);
	Input* basal_melting_input        = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(basal_melting_input);
	Input* dhdt_input                 = element->GetInput(BalancethicknessThickeningRateEnum);      _assert_(dhdt_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		thickness_input->GetInputValue(&thickness, gauss);
		thickness_input->GetInputDerivativeValue(&dH[0],xyz_list,gauss);
		thicknessobs_input->GetInputValue(&thicknessobs, gauss);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);
		element->NodalFunctionsP1Derivatives(dbasis,xyz_list,gauss);

		/*Deal with first part (partial derivative a J with respect to k)*/
		for(resp=0;resp<num_responses;resp++){

			weights_input->GetInputValue(&weight,gauss,responses[resp]);

			switch(responses[resp]){
				case ThicknessAbsMisfitEnum:
					for(int i=0;i<numvertices;i++) ge[i]+= (thicknessobs-thickness)*weight*Jdet*gauss->weight*basis[i];
					break;
				case ThicknessAbsGradientEnum:
					for(int i=0;i<numvertices;i++) ge[i]+= - weight*dH[0]*dbasis[0*numvertices+i]*Jdet*gauss->weight;
					for(int i=0;i<numvertices;i++) ge[i]+= - weight*dH[1]*dbasis[1*numvertices+i]*Jdet*gauss->weight;
					break;
				case ThicknessAlongGradientEnum:
					vx_input->GetInputValue(&vx,gauss);
					vy_input->GetInputValue(&vy,gauss);
					vel = sqrt(vx*vx+vy*vy);
					vx  = vx/(vel+1.e-9);
					vy  = vy/(vel+1.e-9);
					for(int i=0;i<numvertices;i++) ge[i]+= - weight*(dH[0]*vx+dH[1]*vy)*(dbasis[0*numvertices+i]*vx+dbasis[1*numvertices+i]*vy)*Jdet*gauss->weight;
					break;
				case ThicknessAcrossGradientEnum:
					vx_input->GetInputValue(&vx,gauss);
					vy_input->GetInputValue(&vy,gauss);
					vel = sqrt(vx*vx+vy*vy);
					vx  = vx/(vel+1.e-9);
					vy  = vy/(vel+1.e-9);
					for(int i=0;i<numvertices;i++) ge[i]+= - weight*(dH[0]*(-vy)+dH[1]*vx)*(dbasis[0*numvertices+i]*(-vy)+dbasis[1*numvertices+i]*vx)*Jdet*gauss->weight;
					break;
				case BalancethicknessMisfitEnum:
					surface_mass_balance_input->GetInputValue(&surface_mass_balance,gauss);
					basal_melting_input->GetInputValue(&basal_melting,gauss);
					dhdt_input->GetInputValue(&dhdt,gauss);
					vx_input->GetInputValue(&vx,gauss);
					vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
					vy_input->GetInputValue(&vy,gauss);
					vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
					for(int i=0;i<numvertices;i++){
						ge[i]+= - weight*Jdet*gauss->weight*(
							(vx*dH[0]+vy*dH[1] + thickness*(dvx[0]+dvy[1]))*(vx*dbasis[0*numvertices+i]+ vy*dbasis[1*numvertices+i] + basis[i]*(dvx[0]+dvy[1]))
							-(surface_mass_balance-basal_melting-dhdt)*(vx*dbasis[0*numvertices+i]+ vy*dbasis[1*numvertices+i] + basis[i]*(dvx[0]+dvy[1]))
							);
					}
					break;
				default:
					_error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
			}
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	xDelete<int>(responses);
	delete gauss;

}/*}}}*/
void           BalancethicknessAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,ThicknessEnum);
			break;
		case Domain3DEnum:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,ThicknessEnum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           BalancethicknessAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/
