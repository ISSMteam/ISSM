#include "./GLheightadvectionAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void GLheightadvectionAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*No constraints for now*/
}/*}}}*/
void GLheightadvectionAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads*/
}/*}}}*/
void GLheightadvectionAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*First fetch data: */
	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,GLheightadvectionAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  GLheightadvectionAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void GLheightadvectionAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonboundary",MeshVertexonboundaryEnum);
}/*}}}*/
void GLheightadvectionAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           GLheightadvectionAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           GLheightadvectionAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* GLheightadvectionAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* GLheightadvectionAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* GLheightadvectionAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Spawn basal element */
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/* Check if ice in element */
	if(!basalelement->IsIceInElement()) return NULL;

	/*Intermediaries */
	const IssmPDouble yts = 365*24*3600.;
	int        domaintype,dim;
	IssmDouble Jdet,D_scalar,onboundary;
	IssmDouble vel,vx,vy;
	IssmDouble* xyz_list      = NULL;
	Input* vx_input           = NULL;
	Input* vy_input           = NULL;
	Input* bc_input           = NULL;

	/*Get problem dimension*/
	basalelement->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);
	IssmDouble     D[2][2] = {0.};

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	vx_input=basalelement->GetInput(VxEnum); _assert_(vx_input);
	vy_input=basalelement->GetInput(VyEnum); _assert_(vy_input);
	bc_input=basalelement->GetInput(MeshVertexonboundaryEnum); _assert_(bc_input);

	IssmDouble h = basalelement->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		/*Get velocity*/
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);

		if(false){
			/*Streamline diffusion*/
			vel = sqrt(vx*vx+vy*vy);
			if(vel<10./yts){
				vx = 0.; vy = 0.;
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
		else{
			D_scalar=gauss->weight*Jdet;

			bc_input->GetInputValue(&onboundary,gauss);
			if(onboundary>0.){
				/*We do not want to advect garbage, make sure only diffusion is applied on boundary*/
				vx = 0.; vy = 0.;
			}

			/*Diffusion */
			if(sqrt(vx*vx+vy*vy)<1000./31536000.){
				IssmPDouble kappa = -10.;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += D_scalar*kappa*(dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]);
					}
				}
			}

			/*Advection: */
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += (D_scalar*(vx*dbasis[0*numnodes+j]*basis[i] + vy*dbasis[1*numnodes+j]*basis[i]))*1e-2;
				}
			}

			/*Artificial diffusivity*/
			vel=sqrt(vx*vx + vy*vy)+1.e-14;
			D[0][0]=D_scalar*h/(2.*vel)*fabs(vx*vx);  D[0][1]=D_scalar*h/(2.*vel)*fabs(vx*vy);
			D[1][0]=D_scalar*h/(2.*vel)*fabs(vy*vx);  D[1][1]=D_scalar*h/(2.*vel)*fabs(vy*vy);
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += (
								dbasis[0*numnodes+i] *(D[0][0]*dbasis[0*numnodes+j] + D[0][1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1][0]*dbasis[0*numnodes+j] + D[1][1]*dbasis[1*numnodes+j]) 
								);
				}
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;

}/*}}}*/
ElementVector* GLheightadvectionAnalysis::CreatePVector(Element* element){/*{{{*/
	return NULL;
}/*}}}*/
void           GLheightadvectionAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           GLheightadvectionAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           GLheightadvectionAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,GroundinglineHeightEnum);
			break;
		case Domain3DEnum:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,GroundinglineHeightEnum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           GLheightadvectionAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Deal with ocean constraint*/
	SetActiveNodesLSMx(femmodel);

	/*Constrain all nodes that are grounded and unconstrain the ones that float*/
	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);
		int         numnodes  = element->GetNumberOfNodes();
		IssmDouble *mask      = xNew<IssmDouble>(numnodes);
		IssmDouble *bed       = xNew<IssmDouble>(numnodes);
		IssmDouble *ls_active = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&mask[0],MaskOceanLevelsetEnum);
		element->GetInputListOnNodes(&bed[0],BaseEnum);
		element->GetInputListOnNodes(&ls_active[0],IceMaskNodeActivationEnum);

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]<0. && ls_active[in]==1.){
				node->Activate();
			}
			else{
				node->Deactivate();
				node->ApplyConstraint(0,bed[in]);
			}
		}
		xDelete<IssmDouble>(mask);
		xDelete<IssmDouble>(bed);
		xDelete<IssmDouble>(ls_active);
	}

	return;
}/*}}}*/
