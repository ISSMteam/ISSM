#include "./DepthAverageAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void DepthAverageAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
}/*}}}*/
void DepthAverageAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
}/*}}}*/
void DepthAverageAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	::CreateNodes(nodes,iomodel,DepthAverageAnalysisEnum,P1Enum);

}/*}}}*/
int  DepthAverageAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void DepthAverageAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	if(iomodel->domaintype==Domain2DverticalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
	}
}/*}}}*/
void DepthAverageAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           DepthAverageAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           DepthAverageAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* DepthAverageAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* DepthAverageAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* DepthAverageAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	int         dim;
	IssmDouble  Jdet,D,dt=1.e+9;
	IssmDouble *xyz_list = NULL;

	/*Get dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		/*vertical diffusion*/
		D=gauss->weight*Jdet*dt;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D*dbasis[(dim-1)*numnodes+i]*dbasis[(dim-1)*numnodes+j];
			}
		}

		/*Next value (Transient)*/
		D=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D*basis[j]*basis[i];
	} 

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* DepthAverageAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	int         input_enum;
	IssmDouble  Jdet,scalar,value;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe     = element->NewElementVector();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&input_enum,InputToDepthaverageInEnum);
	Input* input = element->GetInput(input_enum); _assert_(input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/* Build transient now */
		input->GetInputValue(&value, gauss);
		scalar=value*Jdet*gauss->weight;
		for(int i=0;i<numnodes;i++) pe->values[i]+=scalar*basis[i];

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;
}/*}}}*/
void           DepthAverageAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           DepthAverageAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           DepthAverageAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int inputenum;
	element->FindParam(&inputenum,InputToDepthaverageOutEnum);
	element->InputUpdateFromSolutionOneDof(solution,inputenum);
}/*}}}*/
void           DepthAverageAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
