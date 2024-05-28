#include "./FreeSurfaceTopAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void FreeSurfaceTopAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Use spcthickness for now*/
	IssmDouble* spcthickness = NULL;
	int         M,N;
	iomodel->FetchData(&spcthickness,&M,&N,"md.masstransport.spcthickness");
	if(M!=iomodel->numberofvertices || N!=1) _error_("Size of constraints not supported yet");

	/*Check if there is any NaN*/
	bool isconstraints = false;
	for(int i=0;i<M;i++) if(!xIsNan<IssmDouble>(spcthickness[i])) isconstraints = true;
	if(!isconstraints){
		iomodel->DeleteData(spcthickness,"md.masstransport.spcthickness");
		return;
	}

	_printf0_("   WARNING: using md.geometry to constrain free surface solver\n");

	/*Use spcthickness for now*/
	IssmDouble* surface= NULL;
	iomodel->FetchData(&surface,&M,&N,"md.geometry.surface");
	if(M!=iomodel->numberofvertices || N!=1) _error_("Size of constraints not supported yet");
	for(int i=0;i<M;i++) if(xIsNan<IssmDouble>(spcthickness[i])) surface[i] = NAN;

	/*Create Constraints based on this new vector*/
	IoModelToConstraintsx(constraints,iomodel,surface,M, N, FreeSurfaceTopAnalysisEnum,P1Enum,0);

	/*Cleanup and return*/
	iomodel->DeleteData(spcthickness,"md.masstransport.spcthickness");
	iomodel->DeleteData(surface,"md.geometry.surface");
}/*}}}*/
void FreeSurfaceTopAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int penpair_ids[2];
	int count=0;
	int numvertex_pairing;

	/*Create Penpair for vertex_pairing: */
	IssmDouble *vertex_pairing=NULL;
	IssmDouble *nodeonsurface=NULL;
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(&nodeonsurface,NULL,NULL,"md.mesh.vertexonsurface");
	for(int i=0;i<numvertex_pairing;i++){

		if(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+0])-1]){

			/*In debugging mode, check that the second node is in the same cpu*/
			_assert_(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+1])-1]);

			/*Skip if one of the two is not on the bed*/
			if(iomodel->domaintype!=Domain2DhorizontalEnum){
				if(!(reCast<bool>(nodeonsurface[reCast<int>(vertex_pairing[2*i+0])-1])) || !(reCast<bool>(nodeonsurface[reCast<int>(vertex_pairing[2*i+1])-1]))) continue;
			}

			/*Get node ids*/
			penpair_ids[0]=reCast<int>(vertex_pairing[2*i+0]);
			penpair_ids[1]=reCast<int>(vertex_pairing[2*i+1]);

			/*Create Load*/
			loads->AddObject(new Penpair( count+1, &penpair_ids[0]));
			count++;
		}
	}

	/*Free resources: */
	iomodel->DeleteData(vertex_pairing,"md.masstransport.vertex_pairing");
	iomodel->DeleteData(nodeonsurface,"md.mesh.vertexonsurface");
}/*}}}*/
void FreeSurfaceTopAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,FreeSurfaceTopAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  FreeSurfaceTopAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void FreeSurfaceTopAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Now, is the model 3d? otherwise, do nothing: */
	if (iomodel->domaintype==Domain2DhorizontalEnum)return;

	int smb_model;
	int finiteelement = P1Enum;

	/*Fetch data needed: */
	iomodel->FindConstant(&smb_model,"md.smb.model");

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
	}
	if(iomodel->domaindim==3){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum);
	}
	switch(smb_model){
		case SMBforcingEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.mass_balance",SmbMassBalanceEnum,0.);
			break;
		default:
			/*Nothing for now*/
			;
	}
}/*}}}*/
void FreeSurfaceTopAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           FreeSurfaceTopAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           FreeSurfaceTopAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* FreeSurfaceTopAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* FreeSurfaceTopAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* FreeSurfaceTopAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int         domaintype,dim,stabilization;
	Element*    topelement = NULL;
	IssmDouble *xyz_list  = NULL;
	IssmDouble  Jdet,D_scalar,dt,h;
	IssmDouble  vel,vx,vy,tau;

	/*Get top element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			topelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = topelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = topelement->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNew<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	topelement->GetVerticesCoordinates(&xyz_list);
	topelement->FindParam(&dt,TimesteppingTimeStepEnum);
	topelement->FindParam(&stabilization,MasstransportStabilizationEnum);
	Input* vx_input=topelement->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=NULL;
	if(dim>1){vy_input = topelement->GetInput(VyEnum); _assert_(vy_input);}
	h = topelement->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	while(gauss->next()){

		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		topelement->NodalFunctions(basis,gauss);
		topelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		vx_input->GetInputValue(&vx,gauss);
		if(dim==2) vy_input->GetInputValue(&vy,gauss);

		/*Transient term*/
		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j];

		/*Advection terms*/
		D_scalar=dt*gauss->weight*Jdet;
		for(int i=0;i<dim*dim;i++) D[i]=0.;
		D[0] = D_scalar*vx;
		if(dim==2) D[1*dim+1]=D_scalar*vy;

		if(dim==1){
			/*\phi_i v\cdot\nabla\phi_j*/
			for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[i]*vx*dbasis[0*numnodes+j];
		}
		else{
			_assert_(dim==2);
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					/*\phi_i v\cdot\nabla\phi_j*/
					Ke->values[i*numnodes+j] += D_scalar*basis[i]*(vx*dbasis[0*numnodes+j] + vy*dbasis[1*numnodes+j]);
				}
			}
		}

		if(stabilization==1){
			/*artifical diffusion*/
			if(dim==1){
				vx_input->GetInputAverage(&vx);
				D[0]=h/2.*fabs(vx);
			}
			else{
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);

				D[0*dim+0]=h/2.0*fabs(vx);
				D[1*dim+1]=h/2.0*fabs(vy);
			}
		}
		else if(stabilization==2){
			/*Streamline upwinding*/
			if(dim==1){
				vel=fabs(vx)+1.e-8;
				D[0] = h/(2.*vel)*vx*vx;
			}
			else{
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				D[0*dim+0]=h/(2*vel)*vx*vx;
				D[1*dim+0]=h/(2*vel)*vy*vx;
				D[0*dim+1]=h/(2*vel)*vx*vy;
				D[1*dim+1]=h/(2*vel)*vy*vy;
			}
		}
		else if(stabilization==5){
			/*SUPG*/
			if(dim==1){
				vx_input->GetInputAverage(&vx);
				tau=h/(2.*fabs(vx)+1e-10);
			}
			else{
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);
				tau=1*h/(2.*pow(vx*vx+vy*vy,0.5)+1e-10);
			}
		}
		if(stabilization==1 || stabilization==2){
			for(int i=0;i<dim*dim;i++) D[i]=D_scalar*D[i];
			if(dim==2){
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += (
								dbasis[0*numnodes+i] *(D[0*dim+0]*dbasis[0*numnodes+j] + D[0*dim+1]*dbasis[1*numnodes+j]) +
								dbasis[1*numnodes+i] *(D[1*dim+0]*dbasis[0*numnodes+j] + D[1*dim+1]*dbasis[1*numnodes+j]) 
								);
					}
				}
			}
			else{
				for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += dbasis[0*numnodes+i]*D[0]*dbasis[0*numnodes+j];
			}
		}
		else if(stabilization==5){
			D_scalar=gauss->weight*Jdet*dt;
			if(dim==2){
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=tau*D_scalar*
							(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i])*
							(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j]);
					}
				}
			}
			else{
				for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j]+=tau*D_scalar*(vx*dbasis[0*numnodes+i])*(vx*dbasis[0*numnodes+j]);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(D);
	delete gauss;
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	return Ke;
}/*}}}*/
ElementVector* FreeSurfaceTopAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	int         domaintype,dim,stabilization;
	IssmDouble  Jdet,dt;
	IssmDouble  ms,surface,vx,vy,vz,tau;
	Element*    topelement = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Get top element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			topelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = topelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = topelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	topelement->GetVerticesCoordinates(&xyz_list);
	topelement->FindParam(&dt,TimesteppingTimeStepEnum);
	Input *ms_input      = topelement->GetInput(SmbMassBalanceEnum); _assert_(ms_input);
	Input *surface_input = topelement->GetInput(SurfaceEnum);        _assert_(surface_input);
	Input *vz_input      = NULL;
	Input *vx_input      = NULL;
	Input *vy_input      = NULL;
	switch(dim){
		case 1: 
			vx_input=topelement->GetInput(VxEnum); _assert_(vx_input);
			vz_input = topelement->GetInput(VyEnum)	; _assert_(vz_input); 
			break;
		case 2: 
			vx_input=topelement->GetInput(VxEnum); _assert_(vx_input);
			vy_input = topelement->GetInput(VyEnum); _assert_(vy_input);
			vz_input = topelement->GetInput(VzEnum); _assert_(vz_input); 
			break;
		default: _error_("not implemented");
	}
	IssmDouble h = topelement->CharacteristicLength();

	/*Initialize mb_correction to 0, do not forget!:*/
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	while(gauss->next()){

		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		topelement->NodalFunctions(basis,gauss);

		ms_input->GetInputValue(&ms,gauss);
		vz_input->GetInputValue(&vz,gauss);
		surface_input->GetInputValue(&surface,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(surface + dt*ms + dt*vz)*basis[i];
	}

	if(stabilization==5){
		/*SUPG*/
		topelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		if(dim==1){
			vx_input->GetInputAverage(&vx);
			tau=h/(2.*fabs(vx)+1e-10);
			for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(dt*ms+dt*vz)*tau*(vx*dbasis[0*numnodes+i]);
		}
		else{
			vx_input->GetInputAverage(&vx);
			vy_input->GetInputAverage(&vy);
			tau=h/(2.*pow(vx*vx+vy*vy,0.5)+1e-10);
			for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(dt*ms+dt*vz)*tau*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	return pe;

}/*}}}*/
void           FreeSurfaceTopAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           FreeSurfaceTopAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           FreeSurfaceTopAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	element->InputUpdateFromSolutionOneDof(solution,SurfaceEnum);

	/*Now, we need to do some "processing"*/
	int numvertices = element->GetNumberOfVertices();
	int        migration_style;

	IssmDouble* surface = xNew<IssmDouble>(numvertices);
	IssmDouble* newsurface = xNew<IssmDouble>(numvertices);
	IssmDouble* thickness = xNew<IssmDouble>(numvertices);
	IssmDouble* base = xNew<IssmDouble>(numvertices);
	IssmDouble* bed = xNew<IssmDouble>(numvertices);
	IssmDouble* phi = xNew<IssmDouble>(numvertices);
	IssmDouble* sealevel = xNew<IssmDouble>(numvertices);

	IssmDouble minthickness = element->FindParam(MasstransportMinThicknessEnum);
	IssmDouble rho_ice      = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water    = element->FindParam(MaterialsRhoSeawaterEnum);

	bool isgroundingline;
	element->FindParam(&isgroundingline,TransientIsgroundinglineEnum);
	element->FindParam(&migration_style,GroundinglineMigrationEnum);
	if(isgroundingline) element->GetInputListOnVertices(&bed[0],BedEnum);

	element->GetInputListOnVertices(&base[0],BaseEnum);
	element->GetInputListOnVertices(&surface[0],SurfaceEnum);
	element->GetInputListOnVertices(&phi[0],MaskOceanLevelsetEnum);
	element->GetInputListOnVertices(&sealevel[0],SealevelEnum);

	for(int i=0;i<numvertices;i++){
		newsurface[i]=surface[i];
		thickness[i]=surface[i]-base[i];
		/*Check solution*/
		if(xIsNan<IssmDouble>(thickness[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(thickness[i])) _error_("Inf found in solution vector");

		/* check for thickness<minthickness */
		if(thickness[i]<minthickness){ 
			thickness[i]=minthickness;
			if(phi[i]>0.){
				if(base[i]<=bed[i]) base[i] = bed[i];
				newsurface[i] = base[i]+minthickness;
			}else{
				// assume floatation condition
				newsurface[i] = (1.-rho_ice/rho_water)*minthickness;
				base[i] = -rho_ice/rho_water*minthickness;
			}
		}

		/* update thickness */
		thickness[i]=newsurface[i]-base[i];
		/* some checks */
		if(thickness[i]<0.) _error_("thickness<0");
		if(newsurface[i]<base[i]) _error_("surface<base");
	}

	/* update inputs */
	element->AddInput(BaseEnum,base,element->GetElementType());
	element->AddInput(SurfaceEnum,newsurface,element->GetElementType());
	element->AddInput(ThicknessEnum,thickness,element->GetElementType());

	/* Free resources */
	xDelete<IssmDouble>(newsurface);
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(base);
	xDelete<IssmDouble>(bed);
	xDelete<IssmDouble>(phi);
	xDelete<IssmDouble>(sealevel);
}/*}}}*/
void           FreeSurfaceTopAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
