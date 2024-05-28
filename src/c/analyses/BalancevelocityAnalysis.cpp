#include "./BalancevelocityAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void BalancevelocityAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*No constraints for now*/
	//IoModelToConstraintsx(constraints,iomodel,"md.balancethickness.spcthickness",BalancevelocityAnalysisEnum,P1Enum);
}/*}}}*/
void BalancevelocityAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads*/
}/*}}}*/
void BalancevelocityAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Check in 3d*/
	if(iomodel->domaintype==Domain3DEnum) _error_("DG 3d not implemented yet");

	/*First fetch data: */
	::CreateNodes(nodes,iomodel,BalancevelocityAnalysisEnum,P1Enum);
}/*}}}*/
int  BalancevelocityAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void BalancevelocityAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

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
void BalancevelocityAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           BalancevelocityAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BalancevelocityAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BalancevelocityAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BalancevelocityAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BalancevelocityAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  dhdt,mb,ms,Jdet;
	IssmDouble  h,gamma,thickness;
	IssmDouble  hnx,hny,dhnx[2],dhny[2];
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    HNx    = xNew<IssmDouble>(numnodes);
	IssmDouble*    HNy    = xNew<IssmDouble>(numnodes);
	IssmDouble*    H      = xNew<IssmDouble>(numnodes);
	IssmDouble*    Nx     = xNew<IssmDouble>(numnodes);
	IssmDouble*    Ny     = xNew<IssmDouble>(numnodes);

	/*Retrieve all Inputs and parameters: */
	element->GetVerticesCoordinates(&xyz_list);
	Input* H_input = element->GetInput(ThicknessEnum); _assert_(H_input);
	h = element->CharacteristicLength();

	/*Get vector N for all nodes and build HNx and HNy*/
	element->GetInputListOnNodes(Nx,DrivingStressXEnum);
	element->GetInputListOnNodes(Ny,DrivingStressYEnum);
	element->GetInputListOnNodes(H,ThicknessEnum);
	for(int i=0;i<numnodes;i++){
		IssmDouble norm=sqrt(Nx[i]*Nx[i]+Ny[i]*Ny[i]+1.e-10);
		HNx[i] = -H[i]*Nx[i]/norm;
		HNy[i] = -H[i]*Ny[i]/norm;
	}

	/*Start looping on the number of gaussian points:*/
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		H_input->GetInputValue(&thickness,gauss);
		if(thickness<50.) thickness=50.;
		element->ValueP1DerivativesOnGauss(&dhnx[0],HNx,xyz_list,gauss);
		element->ValueP1DerivativesOnGauss(&dhny[0],HNy,xyz_list,gauss);
		element->ValueP1OnGauss(&hnx,HNx,gauss);
		element->ValueP1OnGauss(&hny,HNy,gauss);

		gamma=h/(2.*thickness+1.e-10);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += gauss->weight*Jdet*(
							(basis[i]+gamma*(basis[i]*(dhnx[0]+dhny[1]) + dbasis[0*numnodes+i]*hnx + dbasis[1*numnodes+i]*hny))*
							(basis[j]*(dhnx[0]+dhny[1])  + dbasis[0*numnodes+j]*hnx + dbasis[1*numnodes+j]*hny)
							);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(H);
	xDelete<IssmDouble>(Nx);
	xDelete<IssmDouble>(Ny);
	xDelete<IssmDouble>(HNx);
	xDelete<IssmDouble>(HNy);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* BalancevelocityAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
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

	/*Intermediaries */
	IssmDouble dhdt,mb,ms,Jdet;
	IssmDouble gamma,thickness;
	IssmDouble hnx,hny,dhnx[2],dhny[2];
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(numnodes*2);
	IssmDouble*    H      = xNew<IssmDouble>(numnodes);
	IssmDouble*    Nx     = xNew<IssmDouble>(numnodes);
	IssmDouble*    Ny     = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	Input* ms_input   = basalelement->GetInput(SmbMassBalanceEnum);          _assert_(ms_input);
	Input* mb_input   = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(mb_input);
	Input* dhdt_input = basalelement->GetInput(BalancethicknessThickeningRateEnum);      _assert_(dhdt_input);
	Input* H_input    = basalelement->GetInput(ThicknessEnum);                           _assert_(H_input);
	IssmDouble h = basalelement->CharacteristicLength();

	/*Get vector N for all nodes*/
	basalelement->GetInputListOnNodes(Nx,DrivingStressXEnum);
	basalelement->GetInputListOnNodes(Ny,DrivingStressYEnum);
	basalelement->GetInputListOnNodes(H,ThicknessEnum);
	for(int i=0;i<numnodes;i++){
		IssmDouble norm=sqrt(Nx[i]*Nx[i]+Ny[i]*Ny[i]+1.e-10);
		Nx[i] = -H[i]*Nx[i]/norm;
		Ny[i] = -H[i]*Ny[i]/norm;
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		element->ValueP1DerivativesOnGauss(&dhnx[0],Nx,xyz_list,gauss);
		element->ValueP1DerivativesOnGauss(&dhny[0],Ny,xyz_list,gauss);
		element->ValueP1OnGauss(&hnx,Nx,gauss);
		element->ValueP1OnGauss(&hny,Ny,gauss);

		ms_input->GetInputValue(&ms,gauss);
		mb_input->GetInputValue(&mb,gauss);
		dhdt_input->GetInputValue(&dhdt,gauss);
		H_input->GetInputValue(&thickness,gauss);
		if(thickness<50.) thickness=50.;

		gamma=h/(2.*thickness+1.e-10);

		for(int i=0;i<numnodes;i++){
			pe->values[i]+=Jdet*gauss->weight*(ms-mb-dhdt)*( basis[i] + gamma*(basis[i]*(dhnx[0]+dhny[1])+hnx*dbasis[0*numnodes+i] + hny*dbasis[1*numnodes+i]));
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(H);
	xDelete<IssmDouble>(Nx);
	xDelete<IssmDouble>(Ny);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void           BalancevelocityAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           BalancevelocityAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           BalancevelocityAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,VelEnum);
			break;
		case Domain3DEnum:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,VelEnum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           BalancevelocityAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
