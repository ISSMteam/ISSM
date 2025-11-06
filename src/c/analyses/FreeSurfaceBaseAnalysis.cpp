#include "./FreeSurfaceBaseAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void FreeSurfaceBaseAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

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

	_printf0_("   WARNING: using md.geometry to constrain free base solver\n");

	/*Use spcthickness for now*/
	IssmDouble* base= NULL;
	iomodel->FetchData(&base,&M,&N,"md.geometry.base");
	if(M!=iomodel->numberofvertices || N!=1) _error_("Size of constraints not supported yet");
	for(int i=0;i<M;i++) if(xIsNan<IssmDouble>(spcthickness[i])) base[i] = NAN;

	/*Create Constraints based on this new vector*/
	IoModelToConstraintsx(constraints,iomodel,base,M, N, FreeSurfaceBaseAnalysisEnum, P1Enum, 0);

	/*Cleanup and return*/
	iomodel->DeleteData(spcthickness,"md.masstransport.spcthickness");
	iomodel->DeleteData(base,"md.geometry.base");
}/*}}}*/
void FreeSurfaceBaseAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int penpair_ids[2];
	int count=0;
	int numvertex_pairing;

	/*Create Penpair for vertex_pairing: */
	IssmDouble *vertex_pairing=NULL;
	IssmDouble *nodeonbase=NULL;
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(&nodeonbase,NULL,NULL,"md.mesh.vertexonbase");
	for(int i=0;i<numvertex_pairing;i++){

		if(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+0])-1]){

			/*In debugging mode, check that the second node is in the same cpu*/
			_assert_(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+1])-1]);

			/*Skip if one of the two is not on the bed*/
			if(iomodel->domaintype!=Domain2DhorizontalEnum){
				if(!(reCast<bool>(nodeonbase[reCast<int>(vertex_pairing[2*i+0])-1])) || !(reCast<bool>(nodeonbase[reCast<int>(vertex_pairing[2*i+1])-1]))) continue;
			}

			/*Get node ids*/
			penpair_ids[0]=reCast<int>(vertex_pairing[2*i+0]);
			penpair_ids[1]=reCast<int>(vertex_pairing[2*i+1]);

			/*Create Load*/
			loads->AddObject(new Penpair(count+1, &penpair_ids[0]));
			count++;
		}
	}

	/*Free resources: */
	iomodel->DeleteData(vertex_pairing,"md.masstransport.vertex_pairing");
	iomodel->DeleteData(nodeonbase,"md.mesh.vertexonbase");
}/*}}}*/
void FreeSurfaceBaseAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,FreeSurfaceBaseAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  FreeSurfaceBaseAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void FreeSurfaceBaseAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Now, is the model 3d? otherwise, do nothing: */
	if (iomodel->domaintype==Domain2DhorizontalEnum)return;

	/*Finite element type*/
	int finiteelement = P1Enum;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum,0.);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.groundingline.intrusion_distance",GroundinglineIntrusionDistanceEnum);

	if(iomodel->domaindim==3){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum);
	}
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}

	/*Get what we need for ocean-induced basal melting*/
	bool isstochastic;
   int basalforcing_model;
   int melt_parameterization;
   iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
   iomodel->FindConstant(&isstochastic,"md.stochasticforcing.isstochasticforcing");
   iomodel->FindConstant(&melt_parameterization,"md.frontalforcings.parameterization");
	switch(basalforcing_model){
		case FloatingMeltRateEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.floatingice_melting_rate",BasalforcingsFloatingiceMeltingRateEnum);
			if(isstochastic){
				iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.floatingice_melting_rate",BaselineBasalforcingsFloatingiceMeltingRateEnum);
			}
			break;
		case LinearFloatingMeltRateEnum:
			break;
		case MismipFloatingMeltRateEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.meltrate_factor",BasalforcingsMeltrateFactorEnum);
			break;
		case MantlePlumeGeothermalFluxEnum:
			break;
		case SpatialLinearFloatingMeltRateEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_melting_rate",BasalforcingsSpatialDeepwaterMeltingRateEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_elevation",BasalforcingsSpatialDeepwaterElevationEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.upperwater_melting_rate",BasalforcingsSpatialUpperwaterMeltingRateEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.upperwater_elevation",BasalforcingsSpatialUpperwaterElevationEnum);
			if(isstochastic){
				iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_melting_rate",BaselineBasalforcingsSpatialDeepwaterMeltingRateEnum);
			}
			break;
		case BasalforcingsPicoEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basin_id",BasalforcingsPicoBasinIdEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.overturning_coeff",BasalforcingsPicoOverturningCoeffEnum);
			break;
		case BasalforcingsIsmip6Enum:
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basin_id",BasalforcingsIsmip6BasinIdEnum);
			break;
		case BeckmannGoosseFloatingMeltRateEnum:
			bool isthermalforcing;
         iomodel->FindConstant(&isthermalforcing,"md.basalforcings.isthermalforcing");
			iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.meltrate_factor",BasalforcingsMeltrateFactorEnum);
         if(isthermalforcing==0){
            iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_salinity",BasalforcingsOceanSalinityEnum);
            iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_temp",BasalforcingsOceanTempEnum);
         }
         else if(melt_parameterization!=FrontalForcingsRignotarmaEnum){
            iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_thermalforcing",ThermalForcingEnum);
         }
			break;
		default:
			_error_("Basal forcing model "<<EnumToStringx(basalforcing_model)<<" not supported yet");
	}
}/*}}}*/
void FreeSurfaceBaseAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           FreeSurfaceBaseAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           FreeSurfaceBaseAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* FreeSurfaceBaseAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* FreeSurfaceBaseAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* FreeSurfaceBaseAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int         domaintype,dim,stabilization;
	Element*    basalelement = NULL;
	IssmDouble *xyz_list  = NULL;
	IssmDouble  Jdet,D_scalar,dt,h;
	IssmDouble  vel,vx,vy,tau;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNew<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&stabilization,MasstransportStabilizationEnum);
	Input* vx_input=basalelement->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=NULL;
	if(dim>1){vy_input = basalelement->GetInput(VyEnum); _assert_(vy_input);}
	h = basalelement->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		vx_input->GetInputValue(&vx,gauss);
		if(dim==2) vy_input->GetInputValue(&vy,gauss);

		/*Transient term*/
		D_scalar=gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D_scalar*basis[i]*basis[j];

		/*Advection terms*/
		D_scalar=dt*gauss->weight*Jdet;
		for(int i=0;i<dim*dim;i++) D[i]=0.;

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
			/*SSA*/
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
			D_scalar=gauss->weight*Jdet*dt*tau;
			if(dim==2){
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=D_scalar*
							(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i])*
							(vx*dbasis[0*numnodes+j]+vy*dbasis[1*numnodes+j]);
					}
				}
			}
			else{
				for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j]+=D_scalar*(vx*dbasis[0*numnodes+i])*(vx*dbasis[0*numnodes+j]);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(D);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* FreeSurfaceBaseAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	int         domaintype,dim,stabilization;
	IssmDouble  Jdet,dt,intrusiondist_avg,factor;
	IssmDouble  gmb,fmb,mb,bed,vx,vy,vz,tau,gldistance;
	Element*    basalelement = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);

	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();
	int         melt_style,point1;
	IssmDouble  fraction1,fraction2;
	bool        mainlyfloating;

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);
	IssmDouble  gllevelset,phi=1.;

	/*Retrieve all inputs and parameters*/
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&melt_style,GroundinglineMeltInterpolationEnum);
	//basalelement->FindParam(&intrusiondist,GroundinglineIntrusionDistanceEnum);

	Input* groundedice_input   = basalelement->GetInput(MaskOceanLevelsetEnum);              _assert_(groundedice_input);
	Input* gmb_input           = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum);  _assert_(gmb_input);
	Input* fmb_input           = basalelement->GetInput(BasalforcingsFloatingiceMeltingRateEnum);  _assert_(fmb_input);
	Input* base_input          = basalelement->GetInput(BaseEnum);                                 _assert_(base_input);
	Input* gllevelset_input = basalelement->GetInput(MaskOceanLevelsetEnum);              _assert_(gllevelset_input);
	Input* vz_input = NULL;
	Input* vx_input = NULL;
	Input* vy_input = NULL;
	//Input* gldistance_input = basalelement->GetInput(DistanceToGroundinglineEnum); _assert_(gldistance_input); 
	Input* intrusiondist_input = basalelement->GetInput(GroundinglineIntrusionDistanceEnum); _assert_(intrusiondist_input);
	
	switch(dim){
		case 1: 
			vx_input=basalelement->GetInput(VxEnum); _assert_(vx_input);
			vz_input=basalelement->GetInput(VyEnum) ; _assert_(vz_input); 
			break;
		case 2: 
			vx_input=basalelement->GetInput(VxEnum); _assert_(vx_input);
			vy_input=basalelement->GetInput(VyEnum); _assert_(vy_input);
			vz_input=basalelement->GetInput(VzEnum); _assert_(vz_input); 
			break;
		default: _error_("not implemented");
	}
	IssmDouble h = basalelement->CharacteristicLength();

	/*Recover portion of element that is grounded*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	phi=basalelement->GetGroundedPortion(xyz_list);
	Gauss*      gauss     = NULL;
	if(melt_style==SubelementMelt2Enum){
		basalelement->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,MaskOceanLevelsetEnum,0);
		gauss = basalelement->NewGauss(point1,fraction1,fraction2,3);
	}
	if(melt_style==IntrusionMeltEnum){
		/* Calculate here the average intrusion distance value over the element to pass to GetGroundedPart*/
		intrusiondist_input->GetInputAverage(&intrusiondist_avg);
		basalelement->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,DistanceToGroundinglineEnum,intrusiondist_avg);
		gauss = basalelement->NewGauss(point1,fraction1,fraction2,3);
	}
	else{
		gauss = basalelement->NewGauss(3);   
	}

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		vz_input->GetInputValue(&vz,gauss);  
		gmb_input->GetInputValue(&gmb,gauss);
		fmb_input->GetInputValue(&fmb,gauss);
		base_input->GetInputValue(&bed,gauss);
		groundedice_input->GetInputValue(&phi,gauss);
		gllevelset_input->GetInputValue(&gllevelset,gauss);
		if(melt_style==SubelementMelt1Enum){ 
			//if (phi>0.999999999) mb=gmb;
			//else mb=(1-phi)*fmb+phi*gmb; // phi is the fraction of grounded ice so (1-phi) is floating
			if(phi>0) mb=gmb;
			else mb=fmb;
		}
		else if(melt_style==SubelementMelt2Enum){
			if(gllevelset>0.) mb=gmb;
			else mb=fmb;
		}
		else if(melt_style==NoMeltOnPartiallyFloatingEnum){
			if (phi<0.00000001) mb=fmb;  
			else mb=gmb;
		}
		else if(melt_style==FullMeltOnPartiallyFloatingEnum){
			if (phi<0.99999999) mb=fmb;  
			else mb=gmb;
		}
		else if(melt_style==IntrusionMeltEnum){
			Input* gldistance_input = basalelement->GetInput(DistanceToGroundinglineEnum); _assert_(gldistance_input); 
			gldistance_input->GetInputValue(&gldistance,gauss);

			if(intrusiondist_avg==0){
				if(gllevelset>0.) mb=gmb;
				else mb=fmb;
			}
			else if(gldistance>intrusiondist_avg) {
				mb=gmb;
			}
			else if(gldistance<=intrusiondist_avg && gldistance>0) {
				mb=fmb*(1-gldistance/intrusiondist_avg); 
			}
			else{
				mb=fmb;
			}
		}
		else{
			_error_("melt interpolation "<<EnumToStringx(melt_style)<<" not implemented yet");
		}

		factor = Jdet*gauss->weight*(bed+dt*(mb) + dt*vz);
		for(int i=0;i<numnodes;i++) pe->values[i]+=factor*basis[i];

		if(stabilization==5){
			/*SUPG*/
			basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			if(dim==1){
				vx_input->GetInputAverage(&vx);
				tau=h/(2.*fabs(vx)+1e-10);
				factor = Jdet*gauss->weight*(dt*mb+dt*vz)*tau;
				for(int i=0;i<numnodes;i++) pe->values[i]+=factor*(vx*dbasis[0*numnodes+i]);
			}
			else{ 
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);
				tau=1*h/(2.*pow(vx*vx+vy*vy,0.5)+1e-10);
				factor = Jdet*gauss->weight*(bed*0.+dt*mb+dt*vz)*tau;
				for(int i=0;i<numnodes;i++) pe->values[i]+=factor*(vx*dbasis[0*numnodes+i]+vy*dbasis[1*numnodes+i]);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;

}/*}}}*/
void           FreeSurfaceBaseAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           FreeSurfaceBaseAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           FreeSurfaceBaseAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	element->InputUpdateFromSolutionOneDof(solution,BaseEnum);
}/*}}}*/
void           FreeSurfaceBaseAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Intermediary*/
	IssmDouble phi,isonbase,base;

	for(Object* & object : femmodel->elements->objects){

		Element* element=xDynamicCast<Element*>(object);
		if(!element->IsOnBase()) continue;

		int             numnodes = element->GetNumberOfNodes();
		Input* groundedice_input = element->GetInput(MaskOceanLevelsetEnum);  _assert_(groundedice_input);
		Input* onbase_input       = element->GetInput(MeshVertexonbaseEnum);          _assert_(onbase_input);
		Input* base_input        = element->GetInput(BaseEnum);                     _assert_(base_input);

		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numnodes;iv++){
			gauss->GaussNode(element->GetElementType(),iv);
			onbase_input->GetInputValue(&isonbase,gauss);
			if(isonbase==1.){
				groundedice_input->GetInputValue(&phi,gauss);
				if(phi>=0.){
					base_input->GetInputValue(&base,gauss);
					element->nodes[iv]->ApplyConstraint(0,base);
				}
				else{
					element->nodes[iv]->DofInFSet(0);
				}
			}
		}
		delete gauss;
	}
}/*}}}*/
