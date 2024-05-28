#include "./SmoothAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void SmoothAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
}/*}}}*/
void SmoothAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
}/*}}}*/
void SmoothAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	::CreateNodes(nodes,iomodel,SmoothAnalysisEnum,P1Enum);

}/*}}}*/
int  SmoothAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void SmoothAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}
}/*}}}*/
void SmoothAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           SmoothAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           SmoothAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* SmoothAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* SmoothAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* SmoothAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Intermediaries */
	int         domaintype;
	IssmDouble  Jdet,thickness,l;
	IssmDouble *xyz_list = NULL;

	/*Check dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->FindParam(&l,SmoothThicknessMultiplierEnum); _assert_(l>0.);
	element->GetVerticesCoordinates(&xyz_list);
	Input* thickness_input = element->GetInput(ThicknessEnum); _assert_(thickness_input);

	/* Start looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		if(thickness<50.) thickness=50.;

		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += gauss->weight*Jdet*(
							basis[i]*basis[j]
							+(l*thickness)*(l*thickness)*(dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j])
							);
			}
		}
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	return Ke;
}/*}}}*/
ElementVector* SmoothAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Get basal element*/
	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries */
	int         input_enum;
	IssmDouble  Jdet,value;
	IssmDouble *xyz_list  = NULL;
	Input     *input = NULL;

	/*SPECIFICS: Driving stress for balance velocities*/
	Input*      H_input = NULL, *surface_input = NULL, *vx_input = NULL, *vy_input = NULL;
	IssmDouble  taud_x,norms,normv,vx,vy;
	IssmDouble  rho_ice,gravity,slope[2],thickness;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&input_enum,InputToSmoothEnum);

	switch(input_enum){
		case DrivingStressXEnum:
		case DrivingStressYEnum:{
			rho_ice       = element->FindParam(MaterialsRhoIceEnum);
			gravity       = element->FindParam(ConstantsGEnum);
			H_input       = element->GetInput(ThicknessEnum); _assert_(H_input);
			surface_input = element->GetInput(SurfaceEnum);   _assert_(surface_input);
			vx_input      = element->GetInput(VxEnum);
			vy_input      = element->GetInput(VyEnum);
			}
			break;
		case SurfaceSlopeXEnum:
		case SurfaceSlopeYEnum:{
			surface_input = element->GetInput(SurfaceEnum);   _assert_(surface_input);
			}
			break;
		default: input = element->GetInput(input_enum);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		switch(input_enum){
			case DrivingStressXEnum: 
			case DrivingStressYEnum:{
				H_input->GetInputValue(&thickness,gauss);
				surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
				if(vx_input && vy_input){
					vx_input->GetInputValue(&vx,gauss);
					vy_input->GetInputValue(&vy,gauss);
					norms = sqrt(slope[0]*slope[0]+slope[1]*slope[1]+1.e-10);
					normv = sqrt(vx*vx + vy*vy);
					if(normv>15./(365.*24.*3600.)){
						slope[0] = -vx/normv*norms;
						slope[1] = -vy/normv*norms;
					}
				}
				if(input_enum==DrivingStressXEnum)
				 value = rho_ice*gravity*thickness*slope[0];
				else
				 value = rho_ice*gravity*thickness*slope[1];
			}
			break;
			case SurfaceSlopeXEnum: 
				surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
				value = slope[0];
				break;
			case SurfaceSlopeYEnum:
				surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
				value = slope[1];
				break;
			default:
				input->GetInputValue(&value,gauss);
		}

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*value*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           SmoothAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           SmoothAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           SmoothAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	int inputenum,domaintype,elementtype;

	element->FindParam(&inputenum,InputToSmoothEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,inputenum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           SmoothAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
