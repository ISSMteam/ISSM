#include "./StressbalanceSIAAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

/*Model processing*/
void StressbalanceSIAAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	bool       isSIA,isSSA,isL1L2,isMOLHO,isHO,isFS,iscoupling;

	/*Fetch parameters: */
	iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
	iomodel->FindConstant(&isSSA,"md.flowequation.isSSA");
	iomodel->FindConstant(&isL1L2,"md.flowequation.isL1L2");
	iomodel->FindConstant(&isMOLHO,"md.flowequation.isMOLHO");
	iomodel->FindConstant(&isHO,"md.flowequation.isHO");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");

	/*Now, is the flag isSIA on? otherwise, do nothing: */
	if (!isSIA) return;

	/*Do we have coupling*/
	if((isSIA?1.:0.) + (isSSA?1.:0.) + (isL1L2?1.:0.) + (isMOLHO?1.:0.) + (isHO?1.:0.) + (isFS?1.:0.) >1.)
	 iscoupling = true;
	else
	 iscoupling = false;

	/*If no coupling, call Regular IoModelToConstraintsx, else, OLD stuff, keep for now*/
	if(!iscoupling){
		IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvx",StressbalanceSIAAnalysisEnum,P1Enum,0);
		IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvy",StressbalanceSIAAnalysisEnum,P1Enum,1);
	}
	else{
		/*Fetch data: */
		iomodel->FetchData(3,"md.stressbalance.spcvx","md.stressbalance.spcvy","md.flowequation.vertex_equation");

		/*Initialize conunter*/
		int count = 0;

		/*vx and vy are spc'd if we are not on nodeonSIA: */
		for(int i=0;i<iomodel->numberofvertices;i++){
			/*keep only this partition's nodes:*/
			if((iomodel->my_vertices[i])){
				if (IoCodeToEnumVertexEquation(reCast<int>(iomodel->Data("md.flowequation.vertex_equation")[i]))!=SIAApproximationEnum){

					constraints->AddObject(new SpcStatic(count+1,i+1,0,0,StressbalanceSIAAnalysisEnum));
					count++;

					constraints->AddObject(new SpcStatic(count+1,i+1,1,0,StressbalanceSIAAnalysisEnum));
					count++;
				}
				else{
					if (!xIsNan<IssmDouble>(iomodel->Data("md.stressbalance.spcvx")[i])){
						constraints->AddObject(new SpcStatic(count+1,i+1,0,iomodel->Data("md.stressbalance.spcvx")[i],StressbalanceSIAAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
					}

					if (!xIsNan<IssmDouble>(iomodel->Data("md.stressbalance.spcvy")[i])){
						constraints->AddObject(new SpcStatic(count+1,i+1,1,iomodel->Data("md.stressbalance.spcvy")[i],StressbalanceSIAAnalysisEnum)); //add count'th spc, on node i+1, setting dof 2 to vy
						count++;
					}
				}
			}
		}
	}

	/*Free data: */
	iomodel->DeleteData(3,"md.stressbalance.spcvx","md.stressbalance.spcvy","md.flowequation.vertex_equation");

}/*}}}*/
void StressbalanceSIAAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads*/

}/*}}}*/
void StressbalanceSIAAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Intermediaries*/
	bool  isSIA;

	/*Fetch parameters: */
	iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");

	/*Now, is the flag isSIA on? otherwise, do nothing: */
	if(!isSIA) return;

	/*First create nodes*/
	int    lid=0;
	iomodel->FetchData(4,"md.flowequation.borderSSA","md.flowequation.borderFS","md.flowequation.vertex_equation","md.stressbalance.referential");
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	}

	::CreateNodes(nodes,iomodel,StressbalanceSIAAnalysisEnum,P1Enum,SIAApproximationEnum);
	for(Object* & object: nodes->objects){
		Node* node=xDynamicCast<Node*>(object);
		int   sid = node->Sid();
		int approximation=IoCodeToEnumVertexEquation(reCast<int>(iomodel->Data("md.flowequation.vertex_equation")[sid]));
		if(approximation!=SIAApproximationEnum) node->Deactivate();
	}
	iomodel->DeleteData(6,"md.mesh.vertexonbase","md.mesh.vertexonsurface","md.flowequation.borderSSA","md.flowequation.borderFS","md.flowequation.vertex_equation","md.stressbalance.referential");

}/*}}}*/
int  StressbalanceSIAAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 2;
}/*}}}*/
void StressbalanceSIAAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	bool   isSIA;
	bool   ismovingfront;
	iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");

	/*Now, is the flag SIA on? otherwise, do nothing: */
	if (!isSIA)return;

	iomodel->FetchData(1,"md.flowequation.element_equation");

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			/*Need to know the type of approximation for this element*/
			if(iomodel->Data("md.flowequation.element_equation")){
				inputs->SetInput(ApproximationEnum,counter,IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[i])));
			}
			counter++;
		}
	}

	/*Free data: */
	iomodel->DeleteData(1,"md.flowequation.element_equation");

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	if(ismovingfront){
		if(iomodel->domaintype!=Domain2DhorizontalEnum)
			iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum); // required for updating active nodes
	}

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);

}/*}}}*/
void StressbalanceSIAAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*Friction*/
	FrictionUpdateParameters(parameters, iomodel);

}/*}}}*/

/*Finite Element Analysis*/
void           StressbalanceSIAAnalysis::Core(FemModel* femmodel){/*{{{*/

		if(VerboseSolution()) _printf0_("   computing SIA velocities\n");
		femmodel->SetCurrentConfiguration(StressbalanceSIAAnalysisEnum);
		solutionsequence_linear(femmodel);
}/*}}}*/
void           StressbalanceSIAAnalysis::PreCore(FemModel* femmodel){/*{{{*/
_error_("not implemented");
}/*}}}*/
ElementVector* StressbalanceSIAAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* StressbalanceSIAAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* StressbalanceSIAAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			return CreateKMatrix2D(element);
		case Domain3DEnum:
			return CreateKMatrix3D(element);
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
ElementMatrix* StressbalanceSIAAnalysis::CreateKMatrix2D(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble connectivity;

	/*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize Element vector*/
	ElementMatrix* Ke=element->NewElementMatrix();
	for(int iv=0;iv<numvertices;iv++){
		connectivity=(IssmDouble)element->VertexConnectivity(iv);
		Ke->values[(2*iv+0)*2*numvertices+(2*iv+0)]=1./connectivity;
		Ke->values[(2*iv+1)*2*numvertices+(2*iv+1)]=1./connectivity;
	}

	/*Clean up and return*/
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceSIAAnalysis::CreateKMatrix3D(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         i0,i1,j0,j1,nodeup,nodedown,numsegments;
	IssmDouble  slope[2],connectivity[2],one0,one1;
	int        *pairindices = NULL;

	/*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();
	int numdof      = 2*numvertices;

	/*Initialize Element vector*/
	ElementMatrix* Ke=element->NewElementMatrix();

	element->VerticalSegmentIndices(&pairindices,&numsegments);
	for(int is=0;is<numsegments;is++){
		nodedown = pairindices[is*2+0];
		nodeup   = pairindices[is*2+1];
		connectivity[0]=(IssmDouble)element->VertexConnectivity(nodedown);
		connectivity[1]=(IssmDouble)element->VertexConnectivity(nodeup);
		one0=1./connectivity[0];
		one1=1./connectivity[1];

		/*2 dofs of first node*/
		i0=2*nodedown;  i1=2*nodedown+1;
		/*2 dofs of second node*/
		j0=2*nodeup;    j1=2*nodeup+1;

		/*Create matrix for these two nodes*/
		if(element->IsOnBase() && element->IsOnSurface()){
			Ke->values[i0*numdof+i0] = +one0;
			Ke->values[i1*numdof+i1] = +one0;
			Ke->values[j0*numdof+i0] = -one1;
			Ke->values[j0*numdof+j0] = +one1;
			Ke->values[j1*numdof+i1] = -one1;
			Ke->values[j1*numdof+j1] = +one1;
		}
		else if(element->IsOnBase()){
			Ke->values[i0*numdof+i0] = one0;
			Ke->values[i1*numdof+i1] = one0;
			Ke->values[j0*numdof+i0] = -2.*one1;
			Ke->values[j0*numdof+j0] = +2.*one1;
			Ke->values[j1*numdof+i1] = -2.*one1;
			Ke->values[j1*numdof+j1] = +2.*one1;
		}
		else if(element->IsOnSurface()){
			Ke->values[j0*numdof+i0] = -one1;
			Ke->values[j0*numdof+j0] = +one1;
			Ke->values[j1*numdof+i1] = -one1;
			Ke->values[j1*numdof+j1] = +one1;
		}
		else{ //node is on two horizontal layers and beams include the values only once, so the have to use half of the connectivity
			Ke->values[j0*numdof+i0] = -2.*one1;
			Ke->values[j0*numdof+j0] = +2.*one1;
			Ke->values[j1*numdof+i1] = -2.*one1;
			Ke->values[j1*numdof+j1] = +2.*one1;
		}
	}

	/*Clean up and return*/
	xDelete<int>(pairindices);
	return Ke;
}/*}}}*/
ElementVector* StressbalanceSIAAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			return CreatePVector2D(element);
		case Domain3DEnum:
			return CreatePVector3D(element);
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
ElementVector* StressbalanceSIAAnalysis::CreatePVector2D(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int        frictionlaw = 1;
	IssmDouble ub,vb,slope2,drag,thickness,surface,connectivity;
	IssmDouble slope[2];

	/*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize Element vector*/
	ElementVector* pe=element->NewElementVector();

	/*Retrieve all inputs and parameters*/
	IssmDouble  rho_ice    = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  gravity    = element->FindParam(ConstantsGEnum);
	IssmDouble  B,n;
	Input* B_input         = element->GetInput(MaterialsRheologyBbarEnum);_assert_(B_input);
	Input* n_input         = element->GetInput(MaterialsRheologyNEnum);   _assert_(n_input);
	Input* slopex_input    = element->GetInput(SurfaceSlopeXEnum);        _assert_(slopex_input);
	Input* slopey_input    = element->GetInput(SurfaceSlopeYEnum);        _assert_(slopey_input);
	Input* thickness_input = element->GetInput(ThicknessEnum);            _assert_(thickness_input);
	Input* surface_input   = element->GetInput(SurfaceEnum);              _assert_(surface_input);
	Input* drag_input      = NULL;
	if(frictionlaw!=5 && frictionlaw!=1){
		drag_input = element->GetInput(FrictionCoefficientEnum);  _assert_(drag_input);
	}

	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		connectivity=(IssmDouble)element->VertexConnectivity(iv);

		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		surface_input->GetInputValue(&surface,gauss);
		slopex_input->GetInputValue(&slope[0],gauss);
		slopey_input->GetInputValue(&slope[1],gauss);
		slope2=slope[0]*slope[0]+slope[1]*slope[1];

		switch(frictionlaw){
			case 1:
				/*Payne 1995 (m = 1, B = 5e-3/yts = 1.58e-10  m/s/Pa*/
				ub=-1.58*1.e-10*rho_ice*gravity*thickness*slope[0];
				vb=-1.58*1.e-10*rho_ice*gravity*thickness*slope[1];
				break;
			case 2:
				/*Ritz et al. 1996*/
				drag_input->GetInputValue(&drag,gauss);
				ub=drag*(rho_ice*gravity*thickness)*(rho_ice*gravity*thickness)*slope[0]/sqrt(slope2);
				vb=drag*(rho_ice*gravity*thickness)*(rho_ice*gravity*thickness)*slope[1]/sqrt(slope2);
				break;
			case 3:
				/*Rutt et al. 2009*/
				drag_input->GetInputValue(&drag,gauss);
				ub=-drag*rho_ice*gravity*thickness*slope[0];
				vb=-drag*rho_ice*gravity*thickness*slope[1];
				break;
			case 4:
				/*Henning Akesson*/
				drag = -4e-15 * surface + 8.6e-12;
				ub=-drag*rho_ice*gravity*thickness*slope[0];
				vb=-drag*rho_ice*gravity*thickness*slope[1];
				break;
			case 6:
				/*No sliding*/
				ub=0.;
				vb=0.;
				break;
			default:
				_error_("Not supported yet");
		}

		pe->values[2*iv+0]=(ub-2.*pow(rho_ice*gravity,n)*pow(slope2,((n-1.)/2.))*pow(thickness,n+1.)/(pow(B,n)*(n+2))*slope[0])/connectivity;
		pe->values[2*iv+1]=(vb-2.*pow(rho_ice*gravity,n)*pow(slope2,((n-1.)/2.))*pow(thickness,n+1.)/(pow(B,n)*(n+2))*slope[1])/connectivity;
	}

	/*Clean up and return*/
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* StressbalanceSIAAnalysis::CreatePVector3D(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         frictionlaw = 1;
	int         nodeup,nodedown,numsegments;
	IssmDouble  ub,vb,slope2,drag,surface,thickness,constant_part,z,Jdet;
	IssmDouble  slope[2],connectivity[2],xyz_list_line[2][3];
	IssmDouble *xyz_list = NULL;
	int        *pairindices = NULL;

	/*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize Element vector*/
	ElementVector* pe=element->NewElementVector();

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble  rho_ice    = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  gravity    = element->FindParam(ConstantsGEnum);
	IssmDouble B,n;
	Input* B_input         = element->GetInput(MaterialsRheologyBEnum);   _assert_(B_input);
	Input* n_input         = element->GetInput(MaterialsRheologyNEnum);   _assert_(n_input);
	Input* surface_input   = element->GetInput(SurfaceEnum);              _assert_(surface_input);
	Input* slopex_input    = element->GetInput(SurfaceSlopeXEnum);        _assert_(slopex_input);
	Input* slopey_input    = element->GetInput(SurfaceSlopeYEnum);        _assert_(slopey_input);
	Input* thickness_input = element->GetInput(ThicknessEnum);            _assert_(thickness_input);
	Input* drag_input      = NULL;
	Friction* friction     = NULL;
	if(frictionlaw!=5 && frictionlaw!=1){
		drag_input = element->GetInput(FrictionCoefficientEnum);  _assert_(drag_input);
	}
	else if(frictionlaw==5){
		friction=new Friction(element,3);
	}

	/*Get Vertical segment indices*/
	element->VerticalSegmentIndices(&pairindices,&numsegments);
	for(int is=0;is<numsegments;is++){
		nodedown = pairindices[is*2+0];
		nodeup   = pairindices[is*2+1];
		connectivity[0]=(IssmDouble)element->VertexConnectivity(nodedown);
		connectivity[1]=(IssmDouble)element->VertexConnectivity(nodeup);
		for(int i=0;i<3;i++){
			xyz_list_line[0][i]=xyz_list[nodedown*3+i];
			xyz_list_line[1][i]=xyz_list[nodeup*3+i];
		}

		Gauss* gauss=element->NewGaussLine(nodedown,nodeup,3);
		while(gauss->next()){

			B_input->GetInputValue(&B,gauss);
			n_input->GetInputValue(&n,gauss);
			slopex_input->GetInputValue(&slope[0],gauss);
			slopey_input->GetInputValue(&slope[1],gauss);
			surface_input->GetInputValue(&surface,gauss);
			thickness_input->GetInputValue(&thickness,gauss);

			slope2=slope[0]*slope[0]+slope[1]*slope[1];
			constant_part=-2.*pow(rho_ice*gravity,n)*pow(slope2,((n-1.)/2.));

			z = element->GetZcoord(xyz_list,gauss);
			element->JacobianDeterminantLine(&Jdet,&xyz_list_line[0][0],gauss);

			if(element->IsOnSurface()){
				pe->values[2*nodeup+0]+=constant_part*pow((surface-z)/B,n)*slope[0]*Jdet*gauss->weight/connectivity[1];
				pe->values[2*nodeup+1]+=constant_part*pow((surface-z)/B,n)*slope[1]*Jdet*gauss->weight/connectivity[1];
			}
			else{/*connectivity is too large, should take only half on it*/
				pe->values[2*nodeup+0]+=constant_part*pow((surface-z)/B,n)*slope[0]*Jdet*gauss->weight*2./connectivity[1];
				pe->values[2*nodeup+1]+=constant_part*pow((surface-z)/B,n)*slope[1]*Jdet*gauss->weight*2./connectivity[1];
			}
		}

		/*Deal with basal velocities*/
		if(element->IsOnBase()){
			/*Make sure to now push the gauss point to the base*/
			delete gauss; gauss=element->NewGauss();
			gauss->GaussVertex(nodedown);
			switch(frictionlaw){
				case 1:
					/*Payne 1995 (m = 1, B = 5e-3/yts = 1.58e-10  m/s/Pa*/
					ub=-1.58*1.e-10*rho_ice*gravity*thickness*slope[0];
					vb=-1.58*1.e-10*rho_ice*gravity*thickness*slope[1];
					break;
				case 2:
					/*Ritz et al. 1996*/
					drag_input->GetInputValue(&drag,gauss);
					ub=drag*(rho_ice*gravity*thickness)*(rho_ice*gravity*thickness)*slope[0]/sqrt(slope2);
					vb=drag*(rho_ice*gravity*thickness)*(rho_ice*gravity*thickness)*slope[1]/sqrt(slope2);
					break;
				case 3:
					/*Rutt et al. 2009*/
					drag_input->GetInputValue(&drag,gauss);
					ub=-drag*rho_ice*gravity*thickness*slope[0];
					vb=-drag*rho_ice*gravity*thickness*slope[1];
					break;
				case 4:
					/*Henning Akesson*/
					drag_input->GetInputValue(&drag,gauss);
					drag = -4e-15 * surface + 8.6e-12;
					ub=-drag*rho_ice*gravity*thickness*slope[0];
					vb=-drag*rho_ice*gravity*thickness*slope[1];
					break;
				case 5: /*Weertman temp for Kevin*/{
					friction->GetAlpha2WeertmanTemp(&drag,gauss);
					ub = -1./drag * rho_ice*gravity*thickness*slope[0];
					vb = -1./drag * rho_ice*gravity*thickness*slope[1];
					}
					break;
				default:
					_error_("Not supported yet");
			}

			pe->values[2*nodedown+0]+=ub/connectivity[0];
			pe->values[2*nodedown+1]+=vb/connectivity[0];
		}
		delete gauss;
	}

	/*Clean up and return*/
	xDelete<int>(pairindices);
	xDelete<IssmDouble>(xyz_list);
	if(frictionlaw==5) delete friction;
	return pe;

}/*}}}*/
void           StressbalanceSIAAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	IssmDouble vx,vy;
	int       *doflist = NULL;

	/*Fetch number of nodes and initialize values*/
	int         numnodes = element->GetNumberOfNodes();
	int         numdof   = numnodes*2;
	IssmDouble* values   = xNew<IssmDouble>(numdof);

	/*Get dof list and inputs */
	element->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);

	/*Ok, we have the velocities in inputs, fill in solution */
	Gauss* gauss=element->NewGauss();
	for(int i=0;i<numnodes;i++){
		gauss->GaussVertex(i);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		values[i*2+0]=vx;
		values[i*2+1]=vy;
	}

	/*Add value to global vector*/
	solution->SetValues(numdof,doflist,values,INS_VAL);

	/*Free resources:*/
	delete gauss;
	xDelete<int>(doflist);
	xDelete<IssmDouble>(values);
}/*}}}*/
void           StressbalanceSIAAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           StressbalanceSIAAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int         i,domaintype;
	IssmDouble  rho_ice,g;
	int*        doflist=NULL;
	IssmDouble* xyz_list=NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*2;

	/*Fetch dof list and allocate solution vectors*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values    = xNew<IssmDouble>(numdof);
	IssmDouble* vx        = xNew<IssmDouble>(numdof);
	IssmDouble* vy        = xNew<IssmDouble>(numdof);
	IssmDouble* vz        = xNew<IssmDouble>(numdof);
	IssmDouble* vel       = xNew<IssmDouble>(numdof);
	IssmDouble* pressure  = xNew<IssmDouble>(numdof);
	IssmDouble* thickness = xNew<IssmDouble>(numdof);
	IssmDouble* surface   = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	element->TransformSolutionCoord(&values[0],XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		vx[i]=values[i*2+0];
		vy[i]=values[i*2+1];

		/*Check solution*/
		if(xIsNan<IssmDouble>(vx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i])) _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vy[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vy[i])) _error_("Inf found in solution vector");
	}

	/*Get Vz and compute vel*/
	element->GetInputListOnNodes(&vz[0],VzEnum,0.);
	for(i=0;i<numnodes;i++) vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

	/*For pressure: we have not computed pressure in this analysis, for this element. We are in 2D,
	 *so the pressure is just the pressure at the bedrock: */
	rho_ice  = element->FindParam(MaterialsRhoIceEnum);
	g        = element->FindParam(ConstantsGEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->GetInputListOnNodes(&thickness[0],ThicknessEnum);
			for(i=0;i<numnodes;i++) pressure[i]=rho_ice*g*thickness[i];
			break;
		case Domain3DEnum:
			element->GetVerticesCoordinates(&xyz_list);
			element->GetInputListOnNodes(&surface[0],SurfaceEnum);
			for(i=0;i<numnodes;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+2]);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Add vx and vy as inputs to the tria element: */
	element->AddInput(VxEnum,vx,P1Enum);
	element->AddInput(VyEnum,vy,P1Enum);
	element->AddInput(VelEnum,vel,P1Enum);
	element->AddInput(PressureEnum,pressure,P1Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
}/*}}}*/
void           StressbalanceSIAAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/
