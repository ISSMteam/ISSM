#include "./ExtrudeFromBaseAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void ExtrudeFromBaseAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
}/*}}}*/
void ExtrudeFromBaseAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
}/*}}}*/
void ExtrudeFromBaseAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	::CreateNodes(nodes,iomodel,ExtrudeFromBaseAnalysisEnum,P1Enum);

}/*}}}*/
int  ExtrudeFromBaseAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void ExtrudeFromBaseAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

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
void ExtrudeFromBaseAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           ExtrudeFromBaseAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           ExtrudeFromBaseAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* ExtrudeFromBaseAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* ExtrudeFromBaseAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* ExtrudeFromBaseAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  Jdet,D;
	IssmDouble *xyz_list = NULL;

	/*Get dimension*/
	int dim;
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		D=gauss->weight*Jdet;
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += gauss->weight*Jdet*(dbasis[(dim-1)*numnodes+i]*dbasis[(dim-1)*numnodes+j]);
			}
		}
	} 

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* ExtrudeFromBaseAnalysis::CreatePVector(Element* element){/*{{{*/
	return NULL;
}/*}}}*/
void           ExtrudeFromBaseAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           ExtrudeFromBaseAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           ExtrudeFromBaseAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int inputenum;
	element->FindParam(&inputenum,InputToExtrudeEnum);
	element->InputUpdateFromSolutionOneDof(solution,inputenum);
}/*}}}*/
void           ExtrudeFromBaseAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
