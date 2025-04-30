#include "./UzawaPressureAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void UzawaPressureAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	return;
}/*}}}*/
void UzawaPressureAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	return;
}/*}}}*/
void UzawaPressureAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int finiteelement;
	int fe_FS;

	iomodel->FindConstant(&fe_FS,"md.flowequation.fe_FS");
	if(fe_FS==LATaylorHoodEnum) finiteelement = P1Enum;
	else if(fe_FS==LACrouzeixRaviartEnum) finiteelement = P1DGEnum;
	else _error_("solution not supported yet");

	::CreateNodes(nodes,iomodel,UzawaPressureAnalysisEnum,finiteelement);
}/*}}}*/
int  UzawaPressureAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void UzawaPressureAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Update elements: */
	int finiteelement;
	int counter=0;
	int fe_FS;

	iomodel->FindConstant(&fe_FS,"md.flowequation.fe_FS");
	if(fe_FS==LATaylorHoodEnum) finiteelement = P1Enum;
	else if(fe_FS==LACrouzeixRaviartEnum) finiteelement = P1DGEnum;
	else _error_("solution not supported yet");

	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum,0.);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum,0.);
	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum,0.);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum,0.);
	InputUpdateFromConstantx(inputs,elements,0.,SigmaNNEnum);
}/*}}}*/
void UzawaPressureAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.augmented_lagrangian_rhop",AugmentedLagrangianRhopEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.augmented_lagrangian_rholambda",AugmentedLagrangianRholambdaEnum));
}/*}}}*/

/*Finite Element Analysis*/
void           UzawaPressureAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           UzawaPressureAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* UzawaPressureAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* UzawaPressureAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* UzawaPressureAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	IssmDouble  D_scalar,Jdet;
	IssmDouble *xyz_list = NULL;
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke   = element->NewElementMatrix();
	IssmDouble*    M    = xNew<IssmDouble>(numnodes);

	IssmDouble connectivity;
	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	int numvertices = element->GetNumberOfVertices();

	for(int iv=0;iv<numvertices;iv++){
		connectivity=(IssmDouble)element->VertexConnectivity(iv);
		Ke->values[iv*numvertices+iv]=1./connectivity;
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(M);
	return Ke;
}/*}}}*/
ElementVector* UzawaPressureAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	int          dim;
	IssmDouble   Jdet,rhop,divu;
	IssmDouble   *xyz_list = NULL;
	int numnodes = element->GetNumberOfNodes();

	/*Retrieve all inputs and parameters*/
	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&rhop,AugmentedLagrangianRhopEnum);
	element->GetVerticesCoordinates(&xyz_list);

	/*Initialize Element matrix and vectors*/
	ElementVector *pe    = element->NewElementVector();
	IssmDouble    *basis = xNew<IssmDouble>(numnodes);
	IssmDouble     dvx[3];
	IssmDouble     dvy[3];
	IssmDouble     dvz[3];

	Input* vx_input=element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis, gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		if(dim==3){
			vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
		}

		divu=dvx[0]+dvy[1];
		if (dim==3) divu=divu+dvz[2];

		for(int i=0;i<numnodes;i++){
			pe->values[i] += - rhop * divu * Jdet * gauss->weight * basis[i];
		}
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	return pe;
}/*}}}*/
void           UzawaPressureAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           UzawaPressureAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           UzawaPressureAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int        dim;
	int        *doflist       = NULL;
	IssmDouble rholambda,un,vx,vy,vz,sigmann;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes      = element->GetNumberOfNodes();
	int numnodessigma;
	if(element->element_type==P1Enum) numnodessigma=element->GetNumberOfNodes(P2Enum);
	else if(element->element_type==P1DGEnum) numnodessigma=element->GetNumberOfNodes(P2Enum);
	else _error_("finite element not supported yet");

	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch dof list and allocate solution vector*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values        = xNew<IssmDouble>(numnodes);
	IssmDouble* valueslambda  = xNewZeroInit<IssmDouble>(numnodessigma);
	IssmDouble* pressure      = xNew<IssmDouble>(numnodes);
	Input* vx_input           = element->GetInput(VxEnum);      _assert_(vx_input);
	Input* vy_input           = element->GetInput(VyEnum);      _assert_(vy_input);
	Input* vz_input           = NULL;
	if(dim==3){vz_input       = element->GetInput(VzEnum);      _assert_(vz_input);}
	element->GetInputListOnNodes(&pressure[0],PressureEnum);

	/*Update pressure enum first*/
	for(int i=0;i<numnodes;i++){
		values[i]  = pressure[i] + solution[doflist[i]];
	}
	element->AddInput(PressureEnum,values,element->GetElementType());

	/*Now compute sigmann if on base*/
	if(element->IsOnBase() && 0){ 
		Input* sigmann_input      = element->GetInput(SigmaNNEnum); _assert_(sigmann_input);
		if(dim==3) _error_("not implemented yet");

		int baselist[3];
		int onbase=0;
		IssmDouble  Jdet;
		IssmDouble bed_normal[3];
		IssmDouble  Jlambda[3][3]  = {0.0};
		IssmDouble  Cuk[3]         = {0.0};
		IssmDouble  deltalambda[3] = {0.0};
		IssmDouble* vertexonbase  = xNew<IssmDouble>(numnodessigma);
		Input* vertexonbase_input = element->GetInput(MeshVertexonbaseEnum); _assert_(vertexonbase_input);
		Gauss* gauss = element->NewGauss();

		IssmDouble* basis = xNewZeroInit<IssmDouble>(numnodessigma);
		element->GetVerticesCoordinatesBase(&xyz_list_base);
		element->NormalBase(&bed_normal[0],xyz_list_base);
		element->FindParam(&rholambda,AugmentedLagrangianRholambdaEnum);

		for(int i=0;i<numnodessigma;i++){
			gauss->GaussNode(P2Enum,i);
			vertexonbase_input->GetInputValue(&vertexonbase[i], gauss);
			if(vertexonbase[i]==1){ 
				baselist[onbase]=i;
				onbase += 1;
			}
		}
		if(onbase!=3) _error_("basal nodes of element not found");

		delete gauss;
		gauss = element->NewGaussBase(3);
		while(gauss->next()){

			/*Compute Jlambda*/
			element->NodalFunctionsP2(basis,gauss);
			element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					Jlambda[i][j] += Jdet*gauss->weight*basis[baselist[i]]*basis[baselist[j]];
				}
			}

			/*Compute rho_lambd C u^k*/
			vx_input->GetInputValue(&vx, gauss);
			vy_input->GetInputValue(&vy, gauss);
			un=bed_normal[0]*vx + bed_normal[1]*vy;
			for(int i=0;i<3;i++) Cuk[i] += - un*rholambda*Jdet*gauss->weight*basis[baselist[i]];
		}

		/*Now update sigmann*/
		Matrix3x3Solve(&deltalambda[0],&Jlambda[0][0],&Cuk[0]);
		delete gauss;
		gauss = element->NewGauss();
		for(int i=0;i<3;i++){
			gauss->GaussNode(P2Enum,baselist[i]);
			sigmann_input->GetInputValue(&sigmann, gauss);
			valueslambda[baselist[i]] = sigmann + deltalambda[i];
		}

		delete gauss;
		xDelete<IssmDouble>(vertexonbase);
		xDelete<IssmDouble>(xyz_list_base);
		xDelete<IssmDouble>(basis);

		element->AddInput(SigmaNNEnum,valueslambda,P2Enum);
	}

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(valueslambda);
	xDelete<IssmDouble>(pressure);
	xDelete<int>(doflist);
}/*}}}*/
void           UzawaPressureAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
