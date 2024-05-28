#include "./AdjointBalancethickness2Analysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/DatasetInput.h"

/*Model processor*/
void AdjointBalancethickness2Analysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void AdjointBalancethickness2Analysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void AdjointBalancethickness2Analysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
int  AdjointBalancethickness2Analysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void AdjointBalancethickness2Analysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void AdjointBalancethickness2Analysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/

/*Finite Element Analysis*/
void           AdjointBalancethickness2Analysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           AdjointBalancethickness2Analysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* AdjointBalancethickness2Analysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* AdjointBalancethickness2Analysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* AdjointBalancethickness2Analysis::CreateKMatrix(Element* element){/*{{{*/

	Balancethickness2Analysis* analysis = new Balancethickness2Analysis();
	ElementMatrix* Ke = analysis->CreateKMatrix(element);
	delete analysis;

	return Ke;
}/*}}}*/
ElementVector* AdjointBalancethickness2Analysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries */
	int         num_responses,i;
	IssmDouble  vx,vy,vel,Jdet;
	IssmDouble  surface,surfaceobs,weight;
	int        *responses = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe     = element->NewElementVector(SSAApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	element->FindParam(&responses,NULL,InversionCostFunctionsEnum);
	Input* surface_input      = element->GetInput(SurfaceEnum);                          _assert_(surface_input);
	Input* surfaceobs_input   = element->GetInput(InversionSurfaceObsEnum);              _assert_(surfaceobs_input);
	DatasetInput* weights_input      = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* vx_input           = element->GetInput(VxEnum);                                 _assert_(vx_input);
	Input* vy_input           = element->GetInput(VyEnum);                                 _assert_(vy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		surface_input->GetInputValue(&surface, gauss);
		surfaceobs_input->GetInputValue(&surfaceobs, gauss);

		/*Loop over all requested responses*/
		for(int resp=0;resp<num_responses;resp++){
			weights_input->GetInputValue(&weight,gauss,responses[resp]);

			switch(responses[resp]){
				case SurfaceAbsMisfitEnum:
					for(i=0;i<numnodes;i++) pe->values[i]+=(surfaceobs-surface)*weight*Jdet*gauss->weight*basis[i];
					break;
				case OmegaAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case EtaDiffEnum:
					/*Nothing in P vector*/
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
	delete gauss;
	return pe;

}/*}}}*/
void           AdjointBalancethickness2Analysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           AdjointBalancethickness2Analysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
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

	/*Deal with first part (partial derivative a J with respect to k)*/
	for(resp=0;resp<num_responses;resp++) switch(responses[resp]){
		case SurfaceAbsMisfitEnum:
			/*Nothing, \partial J/\partial k = 0*/
			break;
		case OmegaAbsGradientEnum: GradientJOmegaGradient(element,gradient,control_index); break;
		case EtaDiffEnum: GradientJEtaDiff(element,gradient,control_index); break;
		default: _error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
	}

	/*Deal with second term*/
	switch(control_type){
		case BalancethicknessOmegaEnum:           GradientJOmega(element,gradient,control_index); break;
		case BalancethicknessThickeningRateEnum:  GradientJdHdt( element,gradient,control_index); break;
		default: _error_("control type not supported yet: " << EnumToStringx(control_type));
	}

	/*Clean up and return*/
	xDelete<int>(responses);

}/*}}}*/
void           AdjointBalancethickness2Analysis::GradientJdHdt(Element* element,Vector<IssmDouble>* gradient,int control_index){/*{{{*/

	/*Intermediaries*/
	IssmDouble lambda,Jdet; 
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
	Input* adjoint_input = element->GetInput(AdjointEnum);            _assert_(adjoint_input);

	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);

		adjoint_input->GetInputValue(&lambda,gauss);

		/*Build gradient vector (actually -dJ/da): */
		for(int i=0;i<numvertices;i++){
			ge[i]+= -Jdet*gauss->weight*basis[i]*lambda;
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(ge);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<int>(vertexpidlist);
	delete gauss;
}/*}}}*/
void           AdjointBalancethickness2Analysis::GradientJOmega(Element* element,Vector<IssmDouble>* gradient,int control_index){/*{{{*/

	/*Intermediaries*/
	int         n=3;
	IssmDouble  dlambda[2],ds[2],slopex,slopey,slope,omega,Jdet,velobs;
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
	Input* adjoint_input = element->GetInput(AdjointEnum);            _assert_(adjoint_input);
	Input* omega_input   = element->GetInput(BalancethicknessOmegaEnum); _assert_(omega_input);
	Input* surface_input = element->GetInput(SurfaceEnum); _assert_(surface_input); 
	Input* surfaceslopex_input = element->GetInput(SurfaceSlopeXEnum); _assert_(surfaceslopex_input); 
	Input* surfaceslopey_input = element->GetInput(SurfaceSlopeYEnum); _assert_(surfaceslopey_input); 
	Input* velobs_input        = element->GetInput(InversionVelObsEnum); _assert_(velobs_input); 

	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);

		omega_input->GetInputValue(&omega,gauss);
		adjoint_input->GetInputDerivativeValue(&dlambda[0],xyz_list,gauss);
		surface_input->GetInputDerivativeValue(&ds[0],xyz_list,gauss);
		surfaceslopex_input->GetInputValue(&slopex,gauss);
		surfaceslopey_input->GetInputValue(&slopey,gauss);
		velobs_input->GetInputValue(&velobs,gauss);

		slope = sqrt(slopex*slopex + slopey*slopey);
		//if(slope<1.e-5) slope = 1.e-5;

		/*Build gradient vector (actually -dJ/da): */
		for(int i=0;i<numvertices;i++){
			ge[i]+= - Jdet*gauss->weight*basis[i]*velobs/slope*(ds[0]*dlambda[0] + ds[1]*dlambda[1]);
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(ge);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<int>(vertexpidlist);
	delete gauss;
}/*}}}*/
void           AdjointBalancethickness2Analysis::GradientJOmegaGradient(Element* element,Vector<IssmDouble>* gradient,int control_index){/*{{{*/

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble dk[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* dbasis        = xNew<IssmDouble>(2*numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* omega_input = element->GetInput(BalancethicknessOmegaEnum); _assert_(omega_input);
	DatasetInput* weights_input         = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1Derivatives(dbasis,xyz_list,gauss);
		weights_input->GetInputValue(&weight,gauss,OmegaAbsGradientEnum);

		/*Build alpha_complement_list: */
		omega_input->GetInputDerivativeValue(&dk[0],xyz_list,gauss);

		/*Build gradient vector (actually -dJ/ddrag): */
		for(int i=0;i<numvertices;i++){
			ge[i]+=-weight*Jdet*gauss->weight*2*(dk[0]*dk[0] + dk[1]*dk[1])*(dbasis[0*numvertices+i]*dk[0]+dbasis[1*numvertices+i]*dk[1]);
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;

}/*}}}*/
void           AdjointBalancethickness2Analysis::GradientJEtaDiff(Element* element,Vector<IssmDouble>* gradient,int control_index){/*{{{*/

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble omega,omega0; 
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
	Input* omega_input = element->GetInput(BalancethicknessOmegaEnum);   _assert_(omega_input);
	Input* omega0_input = element->GetInput(BalancethicknessOmega0Enum); _assert_(omega0_input);
	DatasetInput* weights_input = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);
		weights_input->GetInputValue(&weight,gauss,EtaDiffEnum);

		/*Build alpha_complement_list: */
		omega_input->GetInputValue(&omega,gauss);
		omega0_input->GetInputValue(&omega0,gauss);

		/*Build gradient vector (actually -dJ/ddrag): */
		for(int i=0;i<numvertices;i++){
			ge[i]+=-weight*Jdet*gauss->weight*(omega - omega0)*basis[i];
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
void           AdjointBalancethickness2Analysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	element->InputUpdateFromSolutionOneDof(solution,AdjointEnum);
}/*}}}*/
void           AdjointBalancethickness2Analysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
