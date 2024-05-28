#include "./StressbalanceVerticalAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

//#define INWOOVZ 1 //For Inwoo Park 0: incompressible assumption(IA), 1: internal deformation(ID), 2: IA+ID

/*Model processing*/
void StressbalanceVerticalAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Intermediary*/
	bool        isSIA,isSSA,isL1L2,isMOLHO,isHO,isFS,iscoupling;
	int         Mz,Nz;
	IssmDouble *spcvz = NULL;

	/*return if not 3d mesh*/
	if(iomodel->domaintype!=Domain3DEnum) return;

	/*Fetch parameters: */
	iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
	iomodel->FindConstant(&isSSA,"md.flowequation.isSSA");
	iomodel->FindConstant(&isL1L2,"md.flowequation.isL1L2");
	iomodel->FindConstant(&isMOLHO,"md.flowequation.isMOLHO");
	iomodel->FindConstant(&isHO,"md.flowequation.isHO");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");

	/*Do we have coupling*/
	if((isSIA?1.:0.) + (isSSA?1.:0.) + (isL1L2?1.:0.) + (isMOLHO?1.:0.) + (isHO?1.:0.) + (isFS?1.:0.) >1.)
	 iscoupling = true;
	else
	 iscoupling = false;

	/*If no coupling, call Regular IoModelToConstraintsx, else, use P1 elements only*/
	if(!iscoupling){
		IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvz",StressbalanceVerticalAnalysisEnum,P1Enum,0);
	}
	else{
		/*Fetch data: */
		iomodel->FetchData(1,"md.flowequation.borderFS");
		/*Fetch Spc*/
		iomodel->FetchData(&spcvz,&Mz,&Nz,"md.stressbalance.spcvz");
		if(Nz>1) _error_("not supported yet (needs to be coded)");

		/*Initialize counter*/
		int count=0;

		/*Create spcs from x,y,z, as well as the spc values on those spcs: */
		for(int i=0;i<iomodel->numberofvertices;i++){

			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){

				if (reCast<int,IssmDouble>(iomodel->Data("md.flowequation.borderFS")[i])){
					constraints->AddObject(new SpcStatic(count+1,i+1,0,0,StressbalanceVerticalAnalysisEnum)); //spc to zero as vertical velocity is done in Horiz for FS
					count++;
				}
				else if (!xIsNan<IssmDouble>(spcvz[i])){
					constraints->AddObject(new SpcStatic(count+1,i+1,0,
									spcvz[i],StressbalanceVerticalAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
					count++;

				}
			} 
		}

		/*Free data: */
		iomodel->DeleteData(1,"md.flowequation.borderFS");
		iomodel->DeleteData(spcvz,"md.stressbalance.spcvz");
	}

}/*}}}*/
void StressbalanceVerticalAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads*/

}/*}}}*/
void StressbalanceVerticalAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*return if not 3d mesh*/
	if(iomodel->domaintype!=Domain3DEnum) return;

	iomodel->FetchData(3,"md.mesh.vertexonbase","md.mesh.vertexonsurface","md.flowequation.vertex_equation");
	::CreateNodes(nodes,iomodel,StressbalanceVerticalAnalysisEnum,P1Enum);
	iomodel->DeleteData(3,"md.mesh.vertexonbase","md.mesh.vertexonsurface","md.flowequation.vertex_equation");
}/*}}}*/
int  StressbalanceVerticalAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void StressbalanceVerticalAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*return if not 3d mesh*/
	if(iomodel->domaintype!=Domain3DEnum) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	#ifdef INWOOVZ
	iomodel->FetchDataToInput(inputs,elements,"md.smb.mass_balance",SmbMassBalanceEnum);
	#endif

	/*Add basal forcings to compute melt rate*/
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
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum,0.);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum,0.);
}/*}}}*/
void StressbalanceVerticalAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/* specific parameters*/

}/*}}}*/

/*Finite Element Analysis*/
void           StressbalanceVerticalAnalysis::Core(FemModel* femmodel){/*{{{*/

		if(VerboseSolution()) _printf0_("   computing vertical velocities\n");
		femmodel->SetCurrentConfiguration(StressbalanceVerticalAnalysisEnum);
		solutionsequence_linear(femmodel);
}/*}}}*/
void           StressbalanceVerticalAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* StressbalanceVerticalAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* StressbalanceVerticalAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* StressbalanceVerticalAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	bool hack = false;

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixVolume(element);
	ElementMatrix* Ke2=NULL;
	if(hack) Ke2=CreateKMatrixBase(element);
	else Ke2=CreateKMatrixSurface(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;

}/*}}}*/
ElementMatrix* StressbalanceVerticalAnalysis::CreateKMatrixBase(Element* element){/*{{{*/

	if(!element->IsOnBase()) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet,normal[3];
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke    = element->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGaussBase(2);
	element->NormalBase(&normal[0],xyz_list);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		D = -gauss->weight*Jdet*normal[2];

		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D*basis[i]*basis[j];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceVerticalAnalysis::CreateKMatrixSurface(Element* element){/*{{{*/

	if(!element->IsOnSurface()) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet,normal[3];
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke    = element->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesTop(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGaussTop(2);
	element->NormalTop(&normal[0],xyz_list);
	while(gauss->next()){

		element->JacobianDeterminantTop(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		D = -gauss->weight*Jdet*normal[2];

		for(int i=0;i<numnodes;i++) for(int j=0;j<numnodes;j++) Ke->values[i*numnodes+j] += D*basis[i]*basis[j];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceVerticalAnalysis::CreateKMatrixVolume(Element* element){/*{{{*/

	/*Intermediaries*/
	IssmDouble  Jdet;
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(3*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += gauss->weight*Jdet*(
							basis[j]*dbasis[2*numnodes+i]
							);
			}
		}
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(basis);
	return Ke;

}/*}}}*/
ElementVector* StressbalanceVerticalAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	bool hack = false;

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorVolume(element);
	ElementVector* pe2=NULL;
	if(hack) pe2=CreatePVectorSurface(element);
	else     pe2=CreatePVectorBase(element);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
ElementVector* StressbalanceVerticalAnalysis::CreatePVectorBase(Element* element){/*{{{*/

	/*Intermediaries */
	int         approximation;
	IssmDouble *xyz_list      = NULL;
	IssmDouble *xyz_list_base = NULL;
	IssmDouble  Jdet,slope[3];
	IssmDouble  vx,vy,vz=0.,dbdx,dbdy;
	IssmDouble  gmb,fmb,phi,basalmeltingvalue;

	if(!element->IsOnBase()) return NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GetInputValue(&approximation,ApproximationEnum);
	Input* base_input=element->GetInput(BaseEnum);                                               _assert_(base_input);
	Input* groundedice_input=element->GetInput(MaskOceanLevelsetEnum);                     _assert_(groundedice_input);
	Input* groundedice_melting_input=element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(groundedice_melting_input);
	Input* floatingice_melting_input=element->GetInput(BasalforcingsFloatingiceMeltingRateEnum); _assert_(floatingice_melting_input);
	Input* vx_input=element->GetInput(VxEnum);                                                   _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);                                                   _assert_(vy_input);
	Input* vzFS_input=NULL;
	if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
		vzFS_input=element->GetInput(VzFSEnum);       _assert_(vzFS_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(2);
	while(gauss->next()){

		groundedice_melting_input->GetInputValue(&gmb,gauss);
		floatingice_melting_input->GetInputValue(&fmb,gauss);
		groundedice_input->GetInputValue(&phi,gauss);
		base_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
			vzFS_input->GetInputValue(&vz,gauss);
		}
		dbdx=slope[0];
		dbdy=slope[1];
		if(phi>0.) basalmeltingvalue=gmb;
		else basalmeltingvalue=fmb;

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=-Jdet*gauss->weight*(vx*dbdx+vy*dbdy-vz-basalmeltingvalue)*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;
}/*}}}*/
ElementVector* StressbalanceVerticalAnalysis::CreatePVectorSurface(Element* element){/*{{{*/

	/*Intermediaries */
	int         approximation;
	IssmDouble *xyz_list      = NULL;
	IssmDouble *xyz_list_surface= NULL;
	IssmDouble  Jdet,slope[3];
	IssmDouble  vx,vy,vz=0.,dsdx,dsdy;
	IssmDouble  smb,smbvalue;

	if(!element->IsOnSurface()) return NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetVerticesCoordinatesTop(&xyz_list_surface);
	element->GetInputValue(&approximation,ApproximationEnum);
	Input* surface_input    =element->GetInput(SurfaceEnum);   _assert_(surface_input);
	Input* smb_input=element->GetInput(SmbMassBalanceEnum);    _assert_(smb_input);
	Input* vx_input=element->GetInput(VxEnum);                 _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);                 _assert_(vy_input);
	Input* vzFS_input=NULL;
	if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
		vzFS_input=element->GetInput(VzFSEnum); _assert_(vzFS_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussTop(2);
	while(gauss->next()){

		smb_input->GetInputValue(&smb,gauss);
		surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
			vzFS_input->GetInputValue(&vz,gauss);
		}
		dsdx=slope[0];
		dsdy=slope[1];

		element->JacobianDeterminantTop(&Jdet,xyz_list_surface,gauss);
		element->NodalFunctions(basis,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=-Jdet*gauss->weight*(vx*dsdx+vy*dsdy-vz+smb)*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_surface);
	return pe;
}/*}}}*/
ElementVector* StressbalanceVerticalAnalysis::CreatePVectorVolume(Element* element){/*{{{*/

	/*Intermediaries*/
	int         approximation;
	IssmDouble  Jdet,dudx,dvdy,dwdz;
	IssmDouble  du[3],dv[3],dw[3];
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and basis functions*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetInputValue(&approximation,ApproximationEnum);
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);
	Input* vzFS_input=NULL;
	if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
		vzFS_input=element->GetInput(VzFSEnum); _assert_(vzFS_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		vx_input->GetInputDerivativeValue(&du[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dv[0],xyz_list,gauss);
		if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
			vzFS_input->GetInputDerivativeValue(&dw[0],xyz_list,gauss);
			dwdz=dw[2];
		}
		else dwdz=0;
		dudx=du[0];
		dvdy=dv[1];

		for(int i=0;i<numnodes;i++) pe->values[i] += (dudx+dvdy+dwdz)*Jdet*gauss->weight*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	return pe;
}/*}}}*/
void           StressbalanceVerticalAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,VzEnum);
}/*}}}*/
void           StressbalanceVerticalAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           StressbalanceVerticalAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int          numnodes = element->GetNumberOfNodes();
	int          numdof=numnodes*1;

	int          i;
	int          approximation;
	int*         doflist  = NULL;
	IssmDouble*  xyz_list = NULL;
	IssmDouble   rho_ice,g;

	/*Get the approximation and do nothing if the element in FS or None*/
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation==FSApproximationEnum || approximation==NoneApproximationEnum){
		return;
	}

	/*Get dof list and vertices coordinates: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble*  values    = xNew<IssmDouble>(numdof);
	IssmDouble*  vx        = xNew<IssmDouble>(numnodes);
	IssmDouble*  vy        = xNew<IssmDouble>(numnodes);
	IssmDouble*  vz        = xNew<IssmDouble>(numnodes);
	IssmDouble*  vzSSA     = xNew<IssmDouble>(numnodes);
	IssmDouble*  vzHO      = xNew<IssmDouble>(numnodes);
	IssmDouble*  vzFS      = xNew<IssmDouble>(numnodes);
	IssmDouble*  vel       = xNew<IssmDouble>(numnodes);
	IssmDouble*  pressure  = xNew<IssmDouble>(numnodes);
	IssmDouble*  surface   = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector vz: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];
	for(i=0;i<numdof;i++){
		vz[i]=values[i*1+0];

		/*Check solution*/
		if(xIsNan<IssmDouble>(vz[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vz[i])) _error_("Inf found in solution vector");
	}

	/*Get Vx and Vy*/
	element->GetInputListOnNodes(&vx[0],VxEnum,0.0); //default is 0
	element->GetInputListOnNodes(&vy[0],VyEnum,0.0); //default is 0

	/*Do some modifications if we actually have a HOFS or SSAFS element*/
	if(approximation==HOFSApproximationEnum){
		Input* vzFS_input=element->GetInput(VzFSEnum);
		if (vzFS_input){
			if (vzFS_input->ObjectEnum()!=PentaInputEnum) _error_("Cannot compute Vel as VzFS is of type " << EnumToStringx(vzFS_input->ObjectEnum()));
			element->GetInputListOnNodes(&vzFS[0],VzFSEnum,0.);
		}
		else _error_("Cannot compute Vz as VzFS in not present in HOFS element");
		for(i=0;i<numnodes;i++){
			vzHO[i]=vz[i];
			vz[i]=vzHO[i]+vzFS[i];
		}
	}
	else if(approximation==SSAFSApproximationEnum){
		Input* vzFS_input=element->GetInput(VzFSEnum);
		if (vzFS_input){
			if (vzFS_input->ObjectEnum()!=PentaInputEnum) _error_("Cannot compute Vel as VzFS is of type " << EnumToStringx(vzFS_input->ObjectEnum()));
			element->GetInputListOnNodes(&vzFS[0],VzFSEnum,0.);
		}
		else _error_("Cannot compute Vz as VzFS in not present in SSAFS element");
		for(i=0;i<numnodes;i++){
			vzSSA[i]=vz[i];
			vz[i]=vzSSA[i]+vzFS[i];
		}
	}

	#ifdef INWOOVZ
	IssmDouble*  thickness = xNew<IssmDouble>(numnodes);
	IssmDouble*  base	     = xNew<IssmDouble>(numnodes);
	IssmDouble*  smb       = xNew<IssmDouble>(numnodes);

	IssmDouble rheology_n=element->material->GetN();
	IssmDouble isvelbc   = 10./(365.25*24*3600); /*10 m/yr*/
	/*Set analytical vertical velocity field at slow flow region.*/
	if(INWOOVZ == 1){
		element->GetInputListOnNodes(&thickness[0],ThicknessEnum,0.);
		element->GetInputListOnNodes(&base[0],BaseEnum,0.);
		element->GetInputListOnNodes(&smb[0],SmbMassBalanceEnum,0.);
		for(i=0;i<numnodes;i++){
			IssmDouble eta = (base[i]+thickness[i]-xyz_list[i*3+2])/thickness[i];
			vz[i] = -(pow(eta,rheology_n+2.0)-1-(rheology_n+2.0)*(eta-1.))/(rheology_n+1.)*smb[i];
		}	
	}
	else if(INWOOVZ== 2) {
		element->GetInputListOnNodes(&thickness[0],ThicknessEnum,0.);
		element->GetInputListOnNodes(&base[0],BaseEnum,0.);
		element->GetInputListOnNodes(&smb[0],SmbMassBalanceEnum,0.);
		for(i=0;i<numnodes;i++){
			if ( vx[i]*vx[i]+vy[i]*vy[i] < isvelbc*isvelbc ){
				IssmDouble eta = (base[i]+thickness[i]-xyz_list[i*3+2])/thickness[i];
				vz[i] = -(pow(eta,rheology_n+2.0)-1-(rheology_n+2.0)*(eta-1.))/(rheology_n+1.)*smb[i];
			}
		}	
	}
	else if(INWOOVZ==0){
		/*nothing to add*/
	}
	else{
		_error_("not supported yet");
	}

	/*Cleanup*/
	xDelete<IssmDouble>(smb);
	xDelete<IssmDouble>(base);
	xDelete<IssmDouble>(thickness);
	#endif

	/*Now Compute vel*/
	for(i=0;i<numnodes;i++) vel[i]=pow( pow(vx[i],2.0) + pow(vy[i],2.0) + pow(vz[i],2.0) , 0.5);

	/*For pressure: we have not computed pressure in this analysis, for this element. We are in 3D, 
	 *so the pressure is just the pressure at the z elevation: except it this is a HOFS element */
	if(approximation!=HOFSApproximationEnum &&  approximation!=SSAFSApproximationEnum){
		rho_ice = element->FindParam(MaterialsRhoIceEnum);
		g       = element->FindParam(ConstantsGEnum);
		element->GetInputListOnNodes(&surface[0],SurfaceEnum,0.);
		for(i=0;i<numnodes;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+2]);
	}
	if(approximation!=HOFSApproximationEnum && approximation!=SSAFSApproximationEnum){
		element->AddInput(PressureEnum,pressure,element->GetElementType());
	}
	else if(approximation==HOFSApproximationEnum){
		element->AddInput(VzHOEnum,vzHO,P1Enum);
	}
	else if(approximation==SSAFSApproximationEnum){
		element->AddInput(VzSSAEnum,vzSSA,P1Enum);
	}
	element->AddInput(VzEnum,vz,P1Enum);
	element->AddInput(VelEnum,vel,P1Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vzSSA);
	xDelete<IssmDouble>(vzHO);
	xDelete<IssmDouble>(vzFS);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
}/*}}}*/
void           StressbalanceVerticalAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	SetActiveNodesLSMx(femmodel);
	return;
}/*}}}*/
