#include "./L2ProjectionBaseAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void L2ProjectionBaseAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*No constraints*/
}/*}}}*/
void L2ProjectionBaseAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads*/
}/*}}}*/
void L2ProjectionBaseAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	}
	else if(iomodel->domaintype==Domain2DverticalEnum){
		iomodel->FetchData(1,"md.mesh.vertexonbase");
	}
	::CreateNodes(nodes,iomodel,L2ProjectionBaseAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  L2ProjectionBaseAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void L2ProjectionBaseAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum & iomodel->domaintype!=Domain3DsurfaceEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}/*}}}*/
void L2ProjectionBaseAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           L2ProjectionBaseAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           L2ProjectionBaseAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* L2ProjectionBaseAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* L2ProjectionBaseAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* L2ProjectionBaseAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
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

	/*Intermediaries */
	IssmDouble  D,Jdet;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke    = basalelement->NewElementMatrix();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		D=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D*basis[i]*basis[j];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* L2ProjectionBaseAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
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

	/*Intermediaries */
	int         input_enum;
	IssmDouble  Jdet,value,slopes[2];
	Input     *input     = NULL;
	Input     *input2    = NULL;
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
		case SurfaceSlopeXEnum: input2 = basalelement->GetInput(SurfaceEnum); _assert_(input2); break;
		case SurfaceSlopeYEnum: input2 = basalelement->GetInput(SurfaceEnum); _assert_(input2); break;
		case BedSlopeXEnum:     input2 = basalelement->GetInput(BaseEnum);     _assert_(input2); break;
		case BedSlopeYEnum:     input2 = basalelement->GetInput(BaseEnum);     _assert_(input2); break;
		case BaseSlopeXEnum:    input2 = basalelement->GetInput(BaseEnum);    _assert_(input2); break;
		case BaseSlopeYEnum:    input2 = basalelement->GetInput(BaseEnum);    _assert_(input2); break;
		case HydrologyGapHeightXEnum:    input2 = basalelement->GetInput(HydrologyGapHeightEnum);  _assert_(input2); break;
		case HydrologyGapHeightXXEnum:   input2 = basalelement->GetInput(HydrologyGapHeightXEnum); _assert_(input2); break;
		case HydrologyGapHeightYEnum:    input2 = basalelement->GetInput(HydrologyGapHeightEnum);  _assert_(input2); break;
		case HydrologyGapHeightYYEnum:   input2 = basalelement->GetInput(HydrologyGapHeightYEnum); _assert_(input2); break;
		case LevelsetfunctionSlopeXEnum: input2 = basalelement->GetInput(MaskIceLevelsetEnum);     _assert_(input2); break;
		case LevelsetfunctionSlopeYEnum: input2 = basalelement->GetInput(MaskIceLevelsetEnum);     _assert_(input2); break;
		default: input = element->GetInput(input_enum);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		if(input2) input2->GetInputDerivativeValue(&slopes[0],xyz_list,gauss);
		switch(input_enum){
			case SurfaceSlopeXEnum: case BedSlopeXEnum: case BaseSlopeXEnum: case LevelsetfunctionSlopeXEnum: 
			case HydrologyGapHeightXEnum: case HydrologyGapHeightXXEnum:
				value = slopes[0];
				break;
			case SurfaceSlopeYEnum: case BedSlopeYEnum: case BaseSlopeYEnum: case LevelsetfunctionSlopeYEnum: 
			case HydrologyGapHeightYEnum: case HydrologyGapHeightYYEnum:
				value = slopes[1];
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
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
void           L2ProjectionBaseAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           L2ProjectionBaseAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           L2ProjectionBaseAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int inputenum,domaintype,elementtype;

	element->FindParam(&inputenum,InputToL2ProjectEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&elementtype,MeshElementtypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,inputenum);
			break;
		case Domain2DverticalEnum:
			element->InputUpdateFromSolutionOneDof(solution,inputenum);
			break;
		case Domain3DEnum:
			if(elementtype==TetraEnum)
			 element->InputUpdateFromSolutionOneDof(solution,inputenum);
			else
			 element->InputUpdateFromSolutionOneDofCollapsed(solution,inputenum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           L2ProjectionBaseAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
