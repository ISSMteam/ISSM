#include "./L2ProjectionEPLAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void L2ProjectionEPLAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*No constraints*/
}/*}}}*/
void L2ProjectionEPLAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads*/
}/*}}}*/
void L2ProjectionEPLAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	/*Now, do we really want DC?*/
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	bool isefficientlayer;
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(1,"md.mesh.vertexonbase");
	}
	else if(iomodel->domaintype==Domain2DverticalEnum){
		iomodel->FetchData(1,"md.mesh.vertexonbase");
	}
	::CreateNodes(nodes,iomodel,L2ProjectionEPLAnalysisEnum,P1Enum);
	iomodel->DeleteData(1,"md.mesh.vertexonbase");
}/*}}}*/
int  L2ProjectionEPLAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void L2ProjectionEPLAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	bool   isefficientlayer;
	int    hydrology_model;

	/*Now, do we really want DC?*/
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.initialization.epl_head",EplHeadSubstepEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}/*}}}*/
void L2ProjectionEPLAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           L2ProjectionEPLAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           L2ProjectionEPLAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* L2ProjectionEPLAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* L2ProjectionEPLAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* L2ProjectionEPLAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
	bool     active_element;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);

	/* Check that all nodes are active, else return empty matrix */
	if(!active_element){
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		return NULL;
	}

	/*Intermediaries */
	IssmDouble  D,Jdet;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke    = basalelement->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		D=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D*basis[j]*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* L2ProjectionEPLAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
	bool     active_element;
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

	basalelement->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);

	/*Check that all nodes are active, else return empty matrix*/
	if(!active_element) {
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		return NULL;
	}

	/*Intermediaries */
	int         input_enum,index;
	IssmDouble  Jdet,slopes[2];
	Input     *input     = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&input_enum,InputToL2ProjectEnum);
	switch(input_enum){
		case EplHeadSlopeXEnum: input = basalelement->GetInput(EplHeadSubstepEnum); index = 0; _assert_(input); break;
		case EplHeadSlopeYEnum: input = basalelement->GetInput(EplHeadSubstepEnum); index = 1; _assert_(input); break;
		default: _error_("not implemented");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		input->GetInputDerivativeValue(&slopes[0],xyz_list,gauss);
		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*slopes[index]*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void           L2ProjectionEPLAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           L2ProjectionEPLAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           L2ProjectionEPLAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	int inputenum,domaintype;

	element->FindParam(&inputenum,InputToL2ProjectEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,inputenum);
			break;
		case Domain2DverticalEnum:
			element->InputUpdateFromSolutionOneDof(solution,inputenum);
			break;
		case Domain3DEnum:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,inputenum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           L2ProjectionEPLAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
