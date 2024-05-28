#include "./Balancethickness2Analysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void Balancethickness2Analysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	int finiteelement = P1Enum;
	IoModelToConstraintsx(constraints,iomodel,"md.balancethickness.spcthickness",Balancethickness2AnalysisEnum,finiteelement);

}/*}}}*/
void Balancethickness2Analysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

}/*}}}*/
void Balancethickness2Analysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int finiteelement = P1Enum;
	::CreateNodes(nodes,iomodel,Balancethickness2AnalysisEnum,finiteelement);
}/*}}}*/
int  Balancethickness2Analysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void Balancethickness2Analysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Finite element type*/
	int finiteelement = P1Enum;

	/*Load variables in element*/
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

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);

			counter++;
		}
	}

}/*}}}*/
void Balancethickness2Analysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           Balancethickness2Analysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           Balancethickness2Analysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* Balancethickness2Analysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* Balancethickness2Analysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* Balancethickness2Analysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  yts = 365*24*3600.;
	IssmDouble  Jdet,vx,vy,vel;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input = element->GetInput(VxEnum); _assert_(vx_input); 
	Input* vy_input = element->GetInput(VyEnum); _assert_(vy_input);

	/*Get element characteristic length*/
	IssmDouble h = element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);

		/*make sure are diffusivisty is large enough*/
		vel = sqrt(vx*vx+vy*vy);
		if(sqrt(vx*vx+vy*vy)==0.){
			vx = 0.1/yts;
			vy = 0.1/yts;
			vel = sqrt(vx*vx+vy*vy);
		}
		else if(vel<30./yts){
			vx = 0.;//vx/vel*0.01;
			vy = 0.;//vy/vel*0.01;
			vel = sqrt(vx*vx+vy*vy);
			vel = 30./yts*500000.;
		}

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += gauss->weight*Jdet*(
							(vx*dbasis[0*numnodes+i] + vy*dbasis[1*numnodes+i])*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j])
				+ vel/500000.*(dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j]));
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* Balancethickness2Analysis::CreatePVector(Element* element){/*{{{*/

	return NULL;
	/*Intermediaries */
	IssmDouble  dhdt[2],mb[2],ms[2],Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* ms_input   = element->GetInput(SmbMassBalanceEnum);                _assert_(ms_input);
	Input* mb_input   = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);       _assert_(mb_input);
	Input* dhdt_input = element->GetInput(BalancethicknessThickeningRateEnum);            _assert_(dhdt_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		ms_input->GetInputDerivativeValue(&ms[0],xyz_list,gauss);

		ms_input->GetInputDerivativeValue(&ms[0],xyz_list,gauss);
		mb_input->GetInputDerivativeValue(&mb[0],xyz_list,gauss);
		dhdt_input->GetInputDerivativeValue(&dhdt[0],xyz_list,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=0*Jdet*gauss->weight*(
					(ms[0]+ms[1]-mb[0]-mb[1]-dhdt[0]-dhdt[1])*basis[i]
					);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           Balancethickness2Analysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
		element->GetSolutionFromInputsOneDof(solution,ThicknessEnum);
}/*}}}*/
void           Balancethickness2Analysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           Balancethickness2Analysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

			element->InputUpdateFromSolutionOneDof(solution,ThicknessEnum);

}/*}}}*/
void           Balancethickness2Analysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
