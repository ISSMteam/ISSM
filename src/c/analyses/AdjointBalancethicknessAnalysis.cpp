#include "./AdjointBalancethicknessAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/DatasetInput.h"

/*Model processor*/
void AdjointBalancethicknessAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void AdjointBalancethicknessAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void AdjointBalancethicknessAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
int  AdjointBalancethicknessAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void AdjointBalancethicknessAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void AdjointBalancethicknessAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/

/*Finite Element Analysis*/
void           AdjointBalancethicknessAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           AdjointBalancethicknessAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* AdjointBalancethicknessAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* AdjointBalancethicknessAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* AdjointBalancethicknessAnalysis::CreateKMatrix(Element* element){/*{{{*/

	BalancethicknessAnalysis* analysis = new BalancethicknessAnalysis();
	ElementMatrix* Ke = analysis->CreateKMatrix(element);
	delete analysis;

	/*Transpose and return Ke*/
	Ke->Transpose();
	return Ke;
}/*}}}*/
ElementVector* AdjointBalancethicknessAnalysis::CreatePVector(Element* element){/*{{{*/

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
	int         num_responses,i;
	IssmDouble  dH[2];
	IssmDouble  vx,vy,vel,Jdet;
	IssmDouble  thickness,thicknessobs,weight;
	int        *responses = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe     = basalelement->NewElementVector(SSAApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	basalelement->FindParam(&responses,NULL,InversionCostFunctionsEnum);
	Input* thickness_input    = basalelement->GetInput(ThicknessEnum);                          _assert_(thickness_input);
	Input* thicknessobs_input = basalelement->GetInput(InversionThicknessObsEnum);              _assert_(thicknessobs_input);
	DatasetInput* weights_input      = basalelement->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* vx_input           = basalelement->GetInput(VxEnum);                                 _assert_(vx_input);
	Input* vy_input           = basalelement->GetInput(VyEnum);                                 _assert_(vy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		thickness_input->GetInputValue(&thickness, gauss);
		thickness_input->GetInputDerivativeValue(&dH[0],xyz_list,gauss);
		thicknessobs_input->GetInputValue(&thicknessobs, gauss);

		/*Loop over all requested responses*/
		for(int resp=0;resp<num_responses;resp++){
			weights_input->GetInputValue(&weight,gauss,responses[resp]);

			switch(responses[resp]){
				case ThicknessAbsMisfitEnum:
					for(i=0;i<numnodes;i++) pe->values[i]+=(thicknessobs-thickness)*weight*Jdet*gauss->weight*basis[i];
					break;
				case ThicknessAbsGradientEnum:
					for(i=0;i<numnodes;i++) pe->values[i]+= - weight*dH[0]*dbasis[0*numnodes+i]*Jdet*gauss->weight;
					for(i=0;i<numnodes;i++) pe->values[i]+= - weight*dH[1]*dbasis[1*numnodes+i]*Jdet*gauss->weight;
					break;
				case ThicknessAlongGradientEnum:
					vx_input->GetInputValue(&vx,gauss);
					vy_input->GetInputValue(&vy,gauss);
					vel = sqrt(vx*vx+vy*vy);
					vx  = vx/(vel+1.e-9);
					vy  = vy/(vel+1.e-9);
					for(i=0;i<numnodes;i++) pe->values[i]+= - weight*(dH[0]*vx+dH[1]*vy)*(dbasis[0*numnodes+i]*vx+dbasis[1*numnodes+i]*vy)*Jdet*gauss->weight;
					break;
				case ThicknessAcrossGradientEnum:
					vx_input->GetInputValue(&vx,gauss);
					vy_input->GetInputValue(&vy,gauss);
					vel = sqrt(vx*vx+vy*vy);
					vx  = vx/(vel+1.e-9);
					vy  = vy/(vel+1.e-9);
					for(i=0;i<numnodes;i++) pe->values[i]+= - weight*(dH[0]*(-vy)+dH[1]*vx)*(dbasis[0*numnodes+i]*(-vy)+dbasis[1*numnodes+i]*vx)*Jdet*gauss->weight;
					break;
				case ThicknessPositiveEnum:
					if(thickness<0){
						for(i=0;i<numnodes;i++) pe->values[i]+= - weight*2*thickness*Jdet*gauss->weight*basis[i];
					}
					break;
				default:
					_error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
			}
		}
	}

	/*Clean up and return*/
	xDelete<int>(responses);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete gauss;
	return pe;
}/*}}}*/
void           AdjointBalancethicknessAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           AdjointBalancethicknessAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	/*The gradient of the cost function is calculated in 2 parts.
	 *
	 * dJ    \partial J   \partial lambda^T(KU-F)
	 * --  = ---------- + ------------------------
	 * dk    \partial k   \parial k                  
	 *
	 * */

	/*If on water, grad = 0: */
	if(!element->IsIceInElement()) return;

	/*Get list of cost functions*/
	int *responses = NULL;
	int num_responses,resp;
	element->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	element->FindParam(&responses,NULL,InversionCostFunctionsEnum);

	/*Check that control_type is supported*/
	if(control_type!=VxEnum && 
		control_type!=VyEnum && 
		control_type!=BalancethicknessSpcthicknessEnum && 
		control_type!=BalancethicknessThickeningRateEnum){
		_error_("Control "<<EnumToStringx(control_type)<<" not supported");
	}

	/*Deal with first part (partial derivative a J with respect to k)*/
	for(resp=0;resp<num_responses;resp++) switch(responses[resp]){
		case ThicknessAbsMisfitEnum:      /*Nothing, \partial J/\partial k = 0*/ break;
		case ThicknessAbsGradientEnum:    /*Nothing, \partial J/\partial k = 0*/ break;
		case ThicknessAlongGradientEnum:  /*Nothing, \partial J/\partial k = 0*/ break;
		case ThicknessAcrossGradientEnum: /*Nothing, \partial J/\partial k = 0*/ break;
		case ThicknessPositiveEnum:       /*Nothing, \partial J/\partial k = 0*/ break;
		default: _error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
	}

	/*Deal with second term*/
	switch(control_type){
		case BalancethicknessSpcthicknessEnum:   GradientJDirichlet(element,gradient,control_interp,control_index); break;
		case BalancethicknessThickeningRateEnum: GradientJDhDt(element,gradient,control_interp,control_index); break;
		case VxEnum:                             GradientJVx(  element,gradient,control_interp,control_index); break;
		case VyEnum:                             GradientJVy(  element,gradient,control_interp,control_index); break;
		default: _error_("control type not supported yet: " << EnumToStringx(control_type));
	}

	/*Clean up and return*/
	xDelete<int>(responses);

}/*}}}*/
void           AdjointBalancethicknessAnalysis::GradientJDirichlet(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	IssmDouble* lambda        = xNew<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GradientIndexing(&vertexpidlist[0],control_index);
	element->GetInputListOnVertices(lambda,AdjointEnum);

	BalancethicknessAnalysis* analysis = new BalancethicknessAnalysis();
	ElementMatrix* Ke = analysis->CreateKMatrix(element);
	delete analysis;

	/*Transpose and return Ke*/
	Ke->Transpose();
	_assert_(Ke->nrows == numvertices);

	for(int i=0;i<numvertices;i++){
		for(int j=0;j<numvertices;j++){
			ge[i] += Ke->values[i*Ke->nrows + j] * lambda[j];
		}
		//ge[i]=-lambda[i];
		_assert_(!xIsNan<IssmDouble>(ge[i]));
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,INS_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(ge);
	xDelete<IssmDouble>(lambda);
	xDelete<int>(vertexpidlist);
	delete Ke;
}/*}}}*/
void           AdjointBalancethicknessAnalysis::GradientJDhDt(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	IssmDouble* lambda        = xNew<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GradientIndexing(&vertexpidlist[0],control_index);
	element->GetInputListOnVertices(lambda,AdjointEnum);
	for(int i=0;i<numvertices;i++){
		ge[i]=-lambda[i];
		_assert_(!xIsNan<IssmDouble>(ge[i]));
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,INS_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(ge);
	xDelete<IssmDouble>(lambda);
	xDelete<int>(vertexpidlist);
}/*}}}*/
void           AdjointBalancethicknessAnalysis::GradientJVx(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	IssmDouble thickness,Jdet,Dlambda[3],dp[3];
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* thickness_input = element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* adjoint_input   = element->GetInput(AdjointEnum);   _assert_(adjoint_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		adjoint_input->GetInputDerivativeValue(&Dlambda[0],xyz_list,gauss);
		thickness_input->GetInputValue(&thickness, gauss);
		thickness_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dD): */
		for(int i=0;i<numvertices;i++){
			ge[i]+=thickness*Dlambda[0]*Jdet*gauss->weight*basis[i];
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
}/*}}}*/
void           AdjointBalancethicknessAnalysis::GradientJVy(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	IssmDouble thickness,Jdet,Dlambda[3],dp[3];
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* thickness_input = element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* adjoint_input   = element->GetInput(AdjointEnum);   _assert_(adjoint_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		adjoint_input->GetInputDerivativeValue(&Dlambda[0],xyz_list,gauss);
		thickness_input->GetInputValue(&thickness, gauss);
		thickness_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dvy): */
		for(int i=0;i<numvertices;i++){
			ge[i]+=thickness*Dlambda[1]*Jdet*gauss->weight*basis[i];
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
}/*}}}*/
void           AdjointBalancethicknessAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,AdjointEnum);
			break;
		case Domain3DEnum:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,AdjointEnum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           AdjointBalancethicknessAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
