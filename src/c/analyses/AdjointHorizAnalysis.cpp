#include <float.h>
#include "./AdjointHorizAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../classes/Inputs/DatasetInput.h"

/*Model processing*/
void AdjointHorizAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void AdjointHorizAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void AdjointHorizAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
int  AdjointHorizAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void AdjointHorizAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void AdjointHorizAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/

/*Finite Element Analysis*/
void AdjointHorizAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void AdjointHorizAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* AdjointHorizAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateKMatrix(Element* element){/*{{{*/
	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case SSAApproximationEnum: 
			return CreateKMatrixSSA(element);
		case L1L2ApproximationEnum: 
			return CreateKMatrixL1L2(element);
		case HOApproximationEnum: 
			return CreateKMatrixHO(element);
		case FSApproximationEnum: 
			return CreateKMatrixFS(element);
      case MOLHOApproximationEnum:
		// a more accurate option, but integrate in the vertical direction numerically.
		//	return CreateKMatrixMOLHOVerticalIntergrated(element);
			return CreateKMatrixMOLHO(element);
		case NoneApproximationEnum:
			return NULL;
		default:
			_error_("Approximation "<<EnumToStringx(approximation)<<" not supported");
	}
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateKMatrixFS(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	bool        incomplete_adjoint;
	int         dim,epssize;
	IssmDouble  Jdet,mu_prime,factor;
	IssmDouble  eps1dotdphii,eps1dotdphij,eps2dotdphii,eps2dotdphij,eps3dotdphii,eps3dotdphij;
	IssmDouble  eps1[3],eps2[3],eps3[3],epsilon[5];/* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);
	if(dim==2) epssize = 3;
	else       epssize = 6;

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Initialize Jacobian with regular FS (first part of the Gateau derivative)*/
	element->FindParam(&incomplete_adjoint,InversionIncompleteAdjointEnum);
	StressbalanceAnalysis* analysis = new StressbalanceAnalysis();
	ElementMatrix* Ke=analysis->CreateKMatrix(element);
	delete analysis;
	if(incomplete_adjoint) return Ke;

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input = element->GetInput(VxEnum);_assert_(vx_input);
	Input* vy_input = element->GetInput(VyEnum);_assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){
		vz_input = element->GetInput(VzEnum);
	}
	else{
		_error_("Not implemented yet");
	}

	/*Allocate dbasis*/
	IssmDouble* dbasis = xNew<IssmDouble>(dim*vnumnodes);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		element->StrainRateHO(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		element->material->ViscosityFSDerivativeEpsSquare(&mu_prime,&epsilon[0],gauss);
		eps1[0]=epsilon[0];   eps2[0]=epsilon[2];   eps3[0]=epsilon[3];
		eps1[1]=epsilon[2];   eps2[1]=epsilon[1];   eps3[1]=epsilon[4];
		eps1[2]=epsilon[3];   eps2[2]=epsilon[4];   eps3[2]= -epsilon[0] -epsilon[1];

		factor = gauss->weight*Jdet*2*mu_prime;
		for(int i=0;i<vnumnodes;i++){
			for(int j=0;j<vnumnodes;j++){
				eps1dotdphii=eps1[0]*dbasis[0*vnumnodes+i]+eps1[1]*dbasis[1*vnumnodes+i]+eps1[2]*dbasis[2*vnumnodes+i];
				eps1dotdphij=eps1[0]*dbasis[0*vnumnodes+j]+eps1[1]*dbasis[1*vnumnodes+j]+eps1[2]*dbasis[2*vnumnodes+j];
				eps2dotdphii=eps2[0]*dbasis[0*vnumnodes+i]+eps2[1]*dbasis[1*vnumnodes+i]+eps2[2]*dbasis[2*vnumnodes+i];
				eps2dotdphij=eps2[0]*dbasis[0*vnumnodes+j]+eps2[1]*dbasis[1*vnumnodes+j]+eps2[2]*dbasis[2*vnumnodes+j];
				eps3dotdphii=eps3[0]*dbasis[0*vnumnodes+i]+eps3[1]*dbasis[1*vnumnodes+i]+eps3[2]*dbasis[2*vnumnodes+i];
				eps3dotdphij=eps3[0]*dbasis[0*vnumnodes+j]+eps3[1]*dbasis[1*vnumnodes+j]+eps3[2]*dbasis[2*vnumnodes+j];

				Ke->values[numdof*(4*i+0)+4*j+0]+=factor*eps1dotdphij*eps1dotdphii;
				Ke->values[numdof*(4*i+0)+4*j+1]+=factor*eps2dotdphij*eps1dotdphii;
				Ke->values[numdof*(4*i+0)+4*j+2]+=factor*eps3dotdphij*eps1dotdphii;

				Ke->values[numdof*(4*i+1)+4*j+0]+=factor*eps1dotdphij*eps2dotdphii;
				Ke->values[numdof*(4*i+1)+4*j+1]+=factor*eps2dotdphij*eps2dotdphii;
				Ke->values[numdof*(4*i+1)+4*j+2]+=factor*eps3dotdphij*eps2dotdphii;

				Ke->values[numdof*(4*i+2)+4*j+0]+=factor*eps1dotdphij*eps3dotdphii;
				Ke->values[numdof*(4*i+2)+4*j+1]+=factor*eps2dotdphij*eps3dotdphii;
				Ke->values[numdof*(4*i+2)+4*j+2]+=factor*eps3dotdphij*eps3dotdphii;
			}
		}
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,XYZEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(xyz_list);
	return Ke;
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateKMatrixHO(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	bool        incomplete_adjoint;
	IssmDouble  Jdet,mu_prime,factor;
	IssmDouble  eps1dotdphii,eps1dotdphij,eps2dotdphii,eps2dotdphij;
	IssmDouble  eps1[3],eps2[3],epsilon[5];/* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Jacobian with regular HO (first part of the Gateau derivative)*/
	element->FindParam(&incomplete_adjoint,InversionIncompleteAdjointEnum);
	StressbalanceAnalysis* analysis = new StressbalanceAnalysis();
	ElementMatrix* Ke=analysis->CreateKMatrix(element);
	delete analysis;
	if(incomplete_adjoint) return Ke;

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input = element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = element->GetInput(VyEnum); _assert_(vy_input);

	/*Allocate dbasis*/
	IssmDouble* dbasis = xNew<IssmDouble>(3*numnodes);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		element->StrainRateHO(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		element->material->ViscosityHODerivativeEpsSquare(&mu_prime,&epsilon[0],gauss);
		eps1[0]=2.*epsilon[0]+epsilon[1];   eps2[0]=epsilon[2];
		eps1[1]=epsilon[2];                 eps2[1]=epsilon[0]+2.*epsilon[1];
		eps1[2]=epsilon[3];                 eps2[2]=epsilon[4];

		factor = gauss->weight*Jdet*2*mu_prime;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				eps1dotdphii=eps1[0]*dbasis[0*numnodes+i]+eps1[1]*dbasis[1*numnodes+i]+eps1[2]*dbasis[2*numnodes+i];
				eps1dotdphij=eps1[0]*dbasis[0*numnodes+j]+eps1[1]*dbasis[1*numnodes+j]+eps1[2]*dbasis[2*numnodes+j];
				eps2dotdphii=eps2[0]*dbasis[0*numnodes+i]+eps2[1]*dbasis[1*numnodes+i]+eps2[2]*dbasis[2*numnodes+i];
				eps2dotdphij=eps2[0]*dbasis[0*numnodes+j]+eps2[1]*dbasis[1*numnodes+j]+eps2[2]*dbasis[2*numnodes+j];

				Ke->values[2*numnodes*(2*i+0)+2*j+0]+=factor*eps1dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+0)+2*j+1]+=factor*eps2dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+0]+=factor*eps1dotdphij*eps2dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+1]+=factor*eps2dotdphij*eps2dotdphii;
			}
		}
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(xyz_list);
	return Ke;
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateKMatrixMOLHO(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble  Jdet,mu_prime,n,thickness,mu,effmu,factor;
	IssmDouble *xyz_list = NULL;
	IssmDouble  viscosity[9]; //9 mu for different integrand
   int			domaintype;
	int			dim=2;

	IssmDouble  eb1i,eb1j,esh1i,esh1j,eb2i,eb2j,esh2i,esh2j;
	IssmDouble  epsilon[5],epsilonbase[5],epsilonshear[5];/* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble  e1b[2], e2b[2], e1sh[2], e2sh[2];
	IssmDouble  vxshear, vyshear;

   /*Get basal element*/
	Element* basalelement = NULL;
   element->FindParam(&domaintype,DomainTypeEnum);
   switch(domaintype){
      case Domain2DhorizontalEnum:
         basalelement = element;
         break;
      case Domain3DEnum: case Domain2DverticalEnum:
         _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
         break;
      default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
   }

	/*Initialize Jacobian with regular HO (first part of the Gateau derivative)*/
	bool        incomplete_adjoint;
	element->FindParam(&incomplete_adjoint,InversionIncompleteAdjointEnum);
	StressbalanceAnalysis* analysis = new StressbalanceAnalysis();
	ElementMatrix* Ke=analysis->CreateKMatrix(basalelement);
	delete analysis;
	if(incomplete_adjoint) return Ke;

	/*Second part of the Gateau derivative if requested*/

	/*Fetch number of nodes and dof for this finite element*/
   int numnodes = basalelement->GetNumberOfNodes();
   IssmDouble* dbasis = xNew<IssmDouble>(2*numnodes); // like SSA
   IssmDouble* basis  = xNew<IssmDouble>(numnodes); // like SSA

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
   Input* vx_input       = element->GetInput(VxEnum);        _assert_(vx_input); //vertically integrated vx
   Input* vy_input       = element->GetInput(VyEnum);        _assert_(vy_input); //vertically integrated vy
   Input* vxbase_input   = element->GetInput(VxBaseEnum);    _assert_(vxbase_input);
   Input* vybase_input   = element->GetInput(VyBaseEnum);    _assert_(vybase_input);
   Input* vxshear_input  = element->GetInput(VxShearEnum);   _assert_(vxshear_input);
   Input* vyshear_input  = element->GetInput(VyShearEnum);   _assert_(vyshear_input);
   Input* thickness_input= element->GetInput(ThicknessEnum); _assert_(thickness_input);
   Input* n_input        = element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);

	/* Start  looping on the number of gaussian points: */
   Gauss* gauss      = element->NewGauss(5);
   Gauss* gauss_base = basalelement->NewGauss();
   while(gauss->next()){
      gauss->SynchronizeGaussBase(gauss_base);

      element->JacobianDeterminant(&Jdet,xyz_list,gauss_base);
      basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss_base);
      basalelement->NodalFunctions(basis, gauss_base);

      thickness_input->GetInputValue(&thickness, gauss);
      n_input->GetInputValue(&n,gauss);

      vxshear_input->GetInputValue(&vxshear,gauss);
      vyshear_input->GetInputValue(&vyshear,gauss);

		element->material->ViscosityMOLHOAdjoint(&viscosity[0],dim,xyz_list,gauss,vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input);

		effmu = 2.0*(1-n)/2.0/n;

		element->StrainRateHO(&epsilonbase[0],xyz_list,gauss,vxbase_input, vybase_input);
		element->StrainRateHO(&epsilonshear[0],xyz_list,gauss,vxshear_input, vyshear_input);

		e1b[0] = 2*epsilonbase[0]+epsilonbase[1];		e1b[1] = epsilonbase[2];
		e2b[1] = epsilonbase[0]+2*epsilonbase[1];		e2b[0] = epsilonbase[2];

		e1sh[0] = 2*epsilonshear[0]+epsilonshear[1];		e1sh[1] = epsilonshear[2];
		e2sh[1] = epsilonshear[0]+2*epsilonshear[1];		e2sh[0] = epsilonshear[2];

		factor = gauss->weight*Jdet*effmu;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				eb1i = e1b[0]*dbasis[0*numnodes+i]+e1b[1]*dbasis[1*numnodes+i];
				eb1j = e1b[0]*dbasis[0*numnodes+j]+e1b[1]*dbasis[1*numnodes+j];
				eb2i = e2b[0]*dbasis[0*numnodes+i]+e2b[1]*dbasis[1*numnodes+i];
				eb2j = e2b[0]*dbasis[0*numnodes+j]+e2b[1]*dbasis[1*numnodes+j];
            esh1i = e1sh[0]*dbasis[0*numnodes+i]+e1sh[1]*dbasis[1*numnodes+i];
            esh1j = e1sh[0]*dbasis[0*numnodes+j]+e1sh[1]*dbasis[1*numnodes+j];
            esh2i = e2sh[0]*dbasis[0*numnodes+i]+e2sh[1]*dbasis[1*numnodes+i];
            esh2j = e2sh[0]*dbasis[0*numnodes+j]+e2sh[1]*dbasis[1*numnodes+j];

				Ke->values[4*numnodes*(4*i+0)+4*j+0]+=factor*(viscosity[0]*eb1j*eb1i+viscosity[1]*(esh1j*eb1i+eb1j*esh1i)+viscosity[2]*esh1j*esh1i);
				Ke->values[4*numnodes*(4*i+1)+4*j+0]+=factor*(viscosity[1]*eb1j*eb1i+viscosity[2]*(esh1j*eb1i+eb1j*esh1i)+viscosity[4]*esh1j*esh1i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*eb1j*basis[i]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*esh1j*basis[i]);
				Ke->values[4*numnodes*(4*i+2)+4*j+0]+=factor*(viscosity[0]*eb1j*eb2i+viscosity[1]*(esh1j*eb2i+eb1j*esh2i)+viscosity[2]*esh1j*esh2i);
				Ke->values[4*numnodes*(4*i+3)+4*j+0]+=factor*(viscosity[1]*eb1j*eb2i+viscosity[2]*(esh1j*eb2i+eb1j*esh2i)+viscosity[4]*esh1j*esh2i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*eb1j*basis[i]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*esh1j*basis[i]);

				Ke->values[4*numnodes*(4*i+0)+4*j+1]+=factor*(viscosity[1]*eb1j*eb1i+viscosity[2]*(esh1j*eb1i+eb1j*esh1i)+viscosity[4]*esh1j*esh1i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*eb1i*basis[j]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*esh1i*basis[j]);
				Ke->values[4*numnodes*(4*i+1)+4*j+1]+=factor*(viscosity[2]*eb1j*eb1i+viscosity[4]*(esh1j*eb1i+eb1j*esh1i)+viscosity[5]*esh1j*esh1i+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*(eb1i*basis[j]+eb1j*basis[i])+viscosity[7]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*(esh1i*basis[j]+esh1j*basis[i])+viscosity[8]*(n+1)*(n+1)*(n+1)*(n+1)/4.0/thickness/thickness/thickness/thickness*vxshear*vxshear*basis[j]*basis[i]);
				Ke->values[4*numnodes*(4*i+2)+4*j+1]+=factor*(viscosity[1]*eb1j*eb2i+viscosity[2]*(esh1j*eb2i+eb1j*esh2i)+viscosity[4]*esh1j*esh2i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*eb2i*basis[j]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*esh2i*basis[j]);
				Ke->values[4*numnodes*(4*i+3)+4*j+1]+=factor*(viscosity[2]*eb1j*eb2i+viscosity[4]*(esh1j*eb2i+eb1j*esh2i)+viscosity[5]*esh1j*esh2i+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*(vxshear*eb2i*basis[j]+vyshear*eb1j*basis[i])+viscosity[7]*(n+1)*(n+1)/2.0/thickness/thickness*(vxshear*esh2i*basis[j]+vyshear*esh1j*basis[i])+viscosity[8]*(n+1)*(n+1)*(n+1)*(n+1)/4.0/thickness/thickness/thickness/thickness*vxshear*vyshear*basis[j]*basis[i]);

				Ke->values[4*numnodes*(4*i+0)+4*j+2]+=factor*(viscosity[0]*eb2j*eb1i+viscosity[1]*(esh2j*eb1i+eb2j*esh1i)+viscosity[2]*esh2j*esh1i);
				Ke->values[4*numnodes*(4*i+1)+4*j+2]+=factor*(viscosity[1]*eb2j*eb1i+viscosity[2]*(esh2j*eb1i+eb2j*esh1i)+viscosity[4]*esh2j*esh1i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*eb2j*basis[i]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vxshear*esh2j*basis[i]);
				Ke->values[4*numnodes*(4*i+2)+4*j+2]+=factor*(viscosity[0]*eb2j*eb2i+viscosity[1]*(esh2j*eb2i+eb2j*esh2i)+viscosity[2]*esh2j*esh2i);
				Ke->values[4*numnodes*(4*i+3)+4*j+2]+=factor*(viscosity[1]*eb2j*eb2i+viscosity[2]*(esh2j*eb2i+eb2j*esh2i)+viscosity[4]*esh2j*esh2i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*eb2j*basis[i]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*esh2j*basis[i]);

				Ke->values[4*numnodes*(4*i+0)+4*j+3]+=factor*(viscosity[1]*eb2j*eb1i+viscosity[2]*(esh2j*eb1i+eb2j*esh1i)+viscosity[4]*esh2j*esh1i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*eb1i*basis[j]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*esh1i*basis[j]);
				Ke->values[4*numnodes*(4*i+1)+4*j+3]+=factor*(viscosity[2]*eb2j*eb1i+viscosity[4]*(esh2j*eb1i+eb2j*esh1i)+viscosity[5]*esh2j*esh1i+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*(vyshear*eb1i*basis[j]+vxshear*eb2j*basis[i])+viscosity[7]*(n+1)*(n+1)/2.0/thickness/thickness*(vxshear*esh2j*basis[i]+vyshear*esh1i*basis[j])+viscosity[8]*(n+1)*(n+1)*(n+1)*(n+1)/4.0/thickness/thickness/thickness/thickness*vxshear*vxshear*basis[j]*basis[i]);
				Ke->values[4*numnodes*(4*i+2)+4*j+3]+=factor*(viscosity[1]*eb2j*eb2i+viscosity[2]*(esh2j*eb2i+eb2j*esh2i)+viscosity[4]*esh2j*esh2i+viscosity[3]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*eb2i*basis[j]+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*esh2i*basis[j]);
				Ke->values[4*numnodes*(4*i+3)+4*j+3]+=factor*(viscosity[2]*eb2j*eb2i+viscosity[4]*(esh2j*eb2i+eb2j*esh2i)+viscosity[5]*esh2j*esh2i+viscosity[6]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*(eb2i*basis[j]+eb1j*basis[i])+viscosity[7]*(n+1)*(n+1)/2.0/thickness/thickness*vyshear*(esh2i*basis[j]+esh1j*basis[i])+viscosity[8]*(n+1)*(n+1)*(n+1)*(n+1)/4.0/thickness/thickness/thickness/thickness*vyshear*vyshear*basis[j]*basis[i]);
			}
		}
	}

	/*Transform Coordinate System*/
	/*element->TransformStiffnessMatrixCoord(Ke,XYEnum);*/

   /*Clean up and return*/
   delete gauss;
   delete gauss_base;
   if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
   xDelete<IssmDouble>(xyz_list);
   xDelete<IssmDouble>(dbasis);
   xDelete<IssmDouble>(basis);
   return Ke;
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateKMatrixMOLHOVerticalIntergrated(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble  Jdet,mu_prime,n,thickness,mu,effmu;
	IssmDouble *xyz_list = NULL;
	IssmDouble	zeta, epsilon_eff;
   int			domaintype;
	int			dim=2;

	IssmDouble  e1phi1i, e1phi1j, e2phi1i, e2phi1j, e1phi2i, e1phi2j, e2phi2i, e2phi2j;
	IssmDouble  epsilon[5];/* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble	e1[3],e2[3];

   Element* basalelement = NULL;

   /*Get basal element*/
   element->FindParam(&domaintype,DomainTypeEnum);
   switch(domaintype){
      case Domain2DhorizontalEnum:
         basalelement = element;
         break;
      case Domain3DEnum: case Domain2DverticalEnum:
         _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
         break;
      default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
   }

	/*Initialize Jacobian with regular HO (first part of the Gateau derivative)*/
	bool incomplete_adjoint;
	element->FindParam(&incomplete_adjoint,InversionIncompleteAdjointEnum);
	StressbalanceAnalysis* analysis = new StressbalanceAnalysis();
	ElementMatrix* Ke=analysis->CreateKMatrix(basalelement);
	delete analysis;
	if(incomplete_adjoint) return Ke;

	/*Second part of the Gateau derivative if requested*/

	/*Fetch number of nodes and dof for this finite element*/
   int numnodes = basalelement->GetNumberOfNodes();
   IssmDouble* dbasis = xNew<IssmDouble>(2*numnodes); // like SSA
   IssmDouble* basis  = xNew<IssmDouble>(numnodes); // like SSA

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
   Input* vxbase_input   = element->GetInput(VxBaseEnum);    _assert_(vxbase_input);
   Input* vybase_input   = element->GetInput(VyBaseEnum);    _assert_(vybase_input);
   Input* vxshear_input  = element->GetInput(VxShearEnum);   _assert_(vxshear_input);
   Input* vyshear_input  = element->GetInput(VyShearEnum);   _assert_(vyshear_input);
   Input* thickness_input= element->GetInput(ThicknessEnum); _assert_(thickness_input);
   Input* n_input        = element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);

	/* Start  looping on the number of gaussian points: */
   Gauss* gauss      = element->NewGauss(10);
   Gauss* gauss_base = basalelement->NewGauss();

	GaussSeg* gauss_seg=new GaussSeg(5);

   while(gauss->next()){
      gauss->SynchronizeGaussBase(gauss_base);

      element->JacobianDeterminant(&Jdet,xyz_list,gauss_base);
      basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss_base);
      basalelement->NodalFunctions(basis, gauss_base);

      thickness_input->GetInputValue(&thickness, gauss);
      n_input->GetInputValue(&n,gauss);

		effmu = 2*(1-n)/2.0/n;

		/* Get the integration in the vertical direction */
		gauss_seg->Reset();
		while(gauss_seg->next()){
			/*Compute zeta for gauss_seg point (0=surface, 1=base)*/
	      zeta=0.5*(gauss_seg->coord1+1);

		   /* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy (for a given zeta)*/
			element->StrainRateMOLHO(&epsilon[0],xyz_list,gauss,
                  vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input,zeta);
			epsilon_eff=sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[2]*epsilon[2]
                    +  epsilon[3]*epsilon[3] + epsilon[4]*epsilon[4] + epsilon[0]*epsilon[1]);

			/*Get viscosity at zeta*/
			element->material->GetViscosity(&mu,epsilon_eff,gauss);
			/*the adjoint viscosity with zeta dependent term*/
			mu = mu /epsilon_eff/epsilon_eff;

			e1[0] = 2.0*epsilon[0]+epsilon[1];
			e1[1] = epsilon[2];
			e1[2] = epsilon[3];

			e2[0] = epsilon[2];
			e2[1] = epsilon[0]+2.0*epsilon[1];
			e2[2] = epsilon[4];

			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					e1phi1i = e1[0]*dbasis[0*numnodes+i]+e1[1]*dbasis[1*numnodes+i];
					e2phi1i = e2[0]*dbasis[0*numnodes+i]+e2[1]*dbasis[1*numnodes+i];
					e1phi1j = e1[0]*dbasis[0*numnodes+j]+e1[1]*dbasis[1*numnodes+j];
					e2phi1j = e2[0]*dbasis[0*numnodes+j]+e2[1]*dbasis[1*numnodes+j];

					e1phi2i = (1-pow(zeta, n+1))*(e1[0]*dbasis[0*numnodes+i]+e1[1]*dbasis[1*numnodes+i])+(n+1)/thickness*pow(zeta, n)*e1[2]*basis[i];
					e2phi2i = (1-pow(zeta, n+1))*(e2[0]*dbasis[0*numnodes+i]+e2[1]*dbasis[1*numnodes+i])+(n+1)/thickness*pow(zeta, n)*e2[2]*basis[i];
					e1phi2j = (1-pow(zeta, n+1))*(e1[0]*dbasis[0*numnodes+j]+e1[1]*dbasis[1*numnodes+j])+(n+1)/thickness*pow(zeta, n)*e1[2]*basis[j];
					e2phi2j = (1-pow(zeta, n+1))*(e2[0]*dbasis[0*numnodes+j]+e2[1]*dbasis[1*numnodes+j])+(n+1)/thickness*pow(zeta, n)*e2[2]*basis[j];

					Ke->values[4*numnodes*(4*i+0)+4*j+0]+=gauss->weight*Jdet*effmu*mu*e1phi1i*e1phi1j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+1)+4*j+0]+=gauss->weight*Jdet*effmu*mu*e1phi2i*e1phi1j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+2)+4*j+0]+=gauss->weight*Jdet*effmu*mu*e2phi1i*e1phi1j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+3)+4*j+0]+=gauss->weight*Jdet*effmu*mu*e2phi2i*e1phi1j*thickness*0.5*gauss_seg->weight;

					Ke->values[4*numnodes*(4*i+0)+4*j+1]+=gauss->weight*Jdet*effmu*mu*e1phi1i*e1phi2j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+1)+4*j+1]+=gauss->weight*Jdet*effmu*mu*e1phi2i*e1phi2j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+2)+4*j+1]+=gauss->weight*Jdet*effmu*mu*e2phi1i*e1phi2j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+3)+4*j+1]+=gauss->weight*Jdet*effmu*mu*e2phi2i*e1phi2j*thickness*0.5*gauss_seg->weight;

					Ke->values[4*numnodes*(4*i+0)+4*j+2]+=gauss->weight*Jdet*effmu*mu*e1phi1i*e2phi1j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+1)+4*j+2]+=gauss->weight*Jdet*effmu*mu*e1phi2i*e2phi1j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+2)+4*j+2]+=gauss->weight*Jdet*effmu*mu*e2phi1i*e2phi1j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+3)+4*j+2]+=gauss->weight*Jdet*effmu*mu*e2phi2i*e2phi1j*thickness*0.5*gauss_seg->weight;

					Ke->values[4*numnodes*(4*i+0)+4*j+3]+=gauss->weight*Jdet*effmu*mu*e1phi1i*e2phi2j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+1)+4*j+3]+=gauss->weight*Jdet*effmu*mu*e1phi2i*e2phi2j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+2)+4*j+3]+=gauss->weight*Jdet*effmu*mu*e2phi1i*e2phi2j*thickness*0.5*gauss_seg->weight;
					Ke->values[4*numnodes*(4*i+3)+4*j+3]+=gauss->weight*Jdet*effmu*mu*e2phi2i*e2phi2j*thickness*0.5*gauss_seg->weight;
				}
			}
		}
	}
	/*Transform Coordinate System*/
	/*element->TransformStiffnessMatrixCoord(Ke,XYEnum);*/

   /*Clean up and return*/
   delete gauss;
   delete gauss_base;
	delete gauss_seg;
   if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
   xDelete<IssmDouble>(xyz_list);
   xDelete<IssmDouble>(dbasis);
   xDelete<IssmDouble>(basis);
   return Ke;
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateKMatrixL1L2(Element* element){/*{{{*/

	/*Intermediaries*/
	bool incomplete_adjoint;

	/*Initialize Jacobian with regular L1L2 (first part of the Gateau derivative)*/
	StressbalanceAnalysis* analysis = new StressbalanceAnalysis();
	ElementMatrix* Ke=analysis->CreateKMatrix(element);
	delete analysis;

	/*return*/
	element->FindParam(&incomplete_adjoint,InversionIncompleteAdjointEnum);
	if(!incomplete_adjoint){
		_error_("Exact adjoint not supported yet for L1L2 model");
	}
	return Ke;
}/*}}}*/
ElementMatrix* AdjointHorizAnalysis::CreateKMatrixSSA(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement(true);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries */
	bool        incomplete_adjoint;
	IssmDouble  Jdet,thickness,mu_prime;
	IssmDouble  eps1dotdphii,eps1dotdphij,eps2dotdphii,eps2dotdphij;
	IssmDouble  eps1[2],eps2[2],epsilon[3];
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Jacobian with regular SSA (first part of the Gateau derivative)*/
	basalelement->FindParam(&incomplete_adjoint,InversionIncompleteAdjointEnum);
	StressbalanceAnalysis* analysis = new StressbalanceAnalysis();
	ElementMatrix* Ke=analysis->CreateKMatrix(element);
	delete analysis;
	if(incomplete_adjoint){
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		return Ke;
	}

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	Input* vx_input        = basalelement->GetInput(VxEnum);       _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyEnum);       _assert_(vy_input);
	Input* thickness_input = basalelement->GetInput(ThicknessEnum); _assert_(thickness_input);

	/*Allocate dbasis*/
	IssmDouble* dbasis = xNew<IssmDouble>(2*numnodes);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		thickness_input->GetInputValue(&thickness, gauss);
		basalelement->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		basalelement->material->ViscositySSADerivativeEpsSquare(&mu_prime,&epsilon[0],gauss);
		eps1[0]=2.*epsilon[0]+epsilon[1]; eps2[0]=epsilon[2];
		eps1[1]=epsilon[2];               eps2[1]=epsilon[0]+2*epsilon[1];

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				eps1dotdphii=eps1[0]*dbasis[0*numnodes+i]+eps1[1]*dbasis[1*numnodes+i];
				eps1dotdphij=eps1[0]*dbasis[0*numnodes+j]+eps1[1]*dbasis[1*numnodes+j];
				eps2dotdphii=eps2[0]*dbasis[0*numnodes+i]+eps2[1]*dbasis[1*numnodes+i];
				eps2dotdphij=eps2[0]*dbasis[0*numnodes+j]+eps2[1]*dbasis[1*numnodes+j];

				Ke->values[2*numnodes*(2*i+0)+2*j+0]+=gauss->weight*Jdet*2*mu_prime*thickness*eps1dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+0)+2*j+1]+=gauss->weight*Jdet*2*mu_prime*thickness*eps2dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+0]+=gauss->weight*Jdet*2*mu_prime*thickness*eps1dotdphij*eps2dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+1]+=gauss->weight*Jdet*2*mu_prime*thickness*eps2dotdphij*eps2dotdphii;
			}
		}
	}

	/*Transform Coordinate System*/
	basalelement->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(xyz_list);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* AdjointHorizAnalysis::CreatePVector(Element* element){/*{{{*/

	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case SSAApproximationEnum: 
			return CreatePVectorSSA(element);
		case L1L2ApproximationEnum: 
			return CreatePVectorL1L2(element);
		case HOApproximationEnum: 
			return CreatePVectorHO(element);
		case FSApproximationEnum: 
			return CreatePVectorFS(element);
      case MOLHOApproximationEnum:
			return CreatePVectorMOLHO(element);
		case NoneApproximationEnum:
			return NULL;
		default:
			_error_("Approximation "<<EnumToStringx(approximation)<<" not supported");
	}
}/*}}}*/
ElementVector* AdjointHorizAnalysis::CreatePVectorFS(Element* element){/*{{{*/

	/*Nothing to be done if not on surface*/
	if(!element->IsOnSurface()) return NULL;

	/*Intermediaries */
	int        num_responses,i,domaintype;
	IssmDouble Jdet,obs_velocity_mag,velocity_mag;
	IssmDouble vx,vy,vxobs,vyobs,dux,duy,weight;
	IssmDouble scalex,scaley,scale,S;
	int        *responses    = NULL;
	IssmDouble *xyz_list_top = NULL;

	/* Get domaintype*/
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(domaintype==Domain2DverticalEnum) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize Element vector and vectors*/
	ElementVector* pe     = element->NewElementVector(FSApproximationEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesTop(&xyz_list_top);
	element->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	element->FindParam(&responses,NULL,InversionCostFunctionsEnum);
	DatasetInput* weights_input = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* vx_input      = element->GetInput(VxEnum);             _assert_(vx_input);
	Input* vxobs_input   = element->GetInput(InversionVxObsEnum); _assert_(vxobs_input);
	Input* vy_input    = NULL;
	Input* vyobs_input = NULL;
	if(domaintype!=Domain2DverticalEnum){
		vy_input      = element->GetInput(VyEnum);             _assert_(vy_input);
		vyobs_input   = element->GetInput(InversionVyObsEnum); _assert_(vyobs_input);
	}
	IssmDouble epsvel  = DBL_EPSILON;
	IssmDouble meanvel = 3.170979198376458e-05; /*1000 m/yr*/

	/*Get Surface if required by one response*/
	Input* S_input = NULL;
	for(int resp=0;resp<num_responses;resp++){
		if(responses[resp]==SurfaceAverageVelMisfitEnum){
			S_input = element->GetInput(SurfaceAreaEnum);  _assert_(S_input); break;
		}
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussTop(4);
	while(gauss->next()){

		element->JacobianDeterminantTop(&Jdet,xyz_list_top,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);

		vx_input->GetInputValue(&vx,gauss);
		vxobs_input->GetInputValue(&vxobs,gauss);
		if(domaintype!=Domain2DverticalEnum) {
			vy_input->GetInputValue(&vy,gauss);
			vyobs_input->GetInputValue(&vyobs,gauss);
		}

		/*Loop over all requested responses*/
		for(int resp=0;resp<num_responses;resp++){
			weights_input->GetInputValue(&weight,gauss,responses[resp]);

			switch(responses[resp]){
				case SurfaceAbsVelMisfitEnum:
					/*
					 *      1  [           2              2 ]
					 * J = --- | (u - u   )  +  (v - v   )  |
					 *      2  [       obs            obs   ]
					 *
					 *        dJ
					 * DU = - -- = (u   - u )
					 *        du     obs
					 */
					for(i=0;i<vnumnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux=vxobs-vx;
							pe->values[i*3+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
							duy=vyobs-vy;
							pe->values[i*3+1]+=duy*weight*Jdet*gauss->weight*vbasis[i]; 
						}
						else{
							dux=vxobs-vx;
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
						}
					}
					break;
				case SurfaceRelVelMisfitEnum:
					/*
					 *      1  [     \bar{v}^2             2   \bar{v}^2              2 ]
					 * J = --- | -------------  (u - u   ) + -------------  (v - v   )  |
					 *      2  [  (u   + eps)^2       obs    (v   + eps)^2       obs    ]
					 *              obs                        obs                      
					 *
					 *        dJ     \bar{v}^2
					 * DU = - -- = ------------- (u   - u )
					 *        du   (u   + eps)^2    obs
					 *               obs
					 */
					for(i=0;i<vnumnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							scalex=pow(meanvel/(vxobs+epsvel),2); if(vxobs==0)scalex=0;
							scaley=pow(meanvel/(vyobs+epsvel),2); if(vyobs==0)scaley=0;
							dux=scalex*(vxobs-vx);
							duy=scaley*(vyobs-vy);
							pe->values[i*3+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
							pe->values[i*3+1]+=duy*weight*Jdet*gauss->weight*vbasis[i]; 
						}
						else{
							scalex=pow(meanvel/(vxobs+epsvel),2); if(vxobs==0)scalex=0;
							dux=scalex*(vxobs-vx);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*vbasis[i]; 
						}
					}
					break;
				case SurfaceLogVelMisfitEnum:
					/*
					 *                 [        vel + eps     ] 2
					 * J = 4 \bar{v}^2 | log ( -----------  ) |  
					 *                 [       vel   + eps    ]
					 *                            obs
					 *
					 *        dJ                 2 * log(...)
					 * DU = - -- = - 4 \bar{v}^2 -------------  u
					 *        du                 vel^2 + eps
					 *            
					 */
					for(i=0;i<vnumnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							velocity_mag    =sqrt(vx*vx+vy*vy)+epsvel;
							obs_velocity_mag=sqrt(vxobs*vxobs+vyobs*vyobs)+epsvel;
							scale=-8.*meanvel*meanvel/(velocity_mag*velocity_mag)*log(velocity_mag/obs_velocity_mag);
							dux=scale*vx;
							duy=scale*vy;
							pe->values[i*3+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
							pe->values[i*3+1]+=duy*weight*Jdet*gauss->weight*vbasis[i]; 
						}
						else{
							velocity_mag    =fabs(vx)+epsvel;
							obs_velocity_mag=fabs(vxobs)+epsvel;
							scale=-8.*meanvel*meanvel/(velocity_mag*velocity_mag)*log(velocity_mag/obs_velocity_mag);
							dux=scale*vx;
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
						}
					}
					break;
				case SurfaceAverageVelMisfitEnum:
					/*
					 *      1                    2              2
					 * J = ---  sqrt(  (u - u   )  +  (v - v   )  )
					 *      S                obs            obs
					 *
					 *        dJ      1       1 
					 * DU = - -- = - --- ----------- * 2 (u - u   )
					 *        du      S  2 sqrt(...)           obs
					 */
					S_input->GetInputValue(&S,gauss);
					for(i=0;i<vnumnodes;i++){
						if (domaintype!=Domain2DverticalEnum){
							scale=1./(S*sqrt(pow(vx-vxobs,2)+pow(vy-vyobs,2))+epsvel);
							dux=scale*(vxobs-vx);
							duy=scale*(vyobs-vy);
							pe->values[i*3+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
							pe->values[i*3+1]+=duy*weight*Jdet*gauss->weight*vbasis[i]; 
						}
						else{
							scale=1./(S*fabs(vx-vxobs)+epsvel);
							dux=scale*(vxobs-vx);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
						}
					}
					break;
				case SurfaceLogVxVyMisfitEnum:
					/*
					 *      1            [        |u| + eps     2          |v| + eps     2  ]
					 * J = --- \bar{v}^2 | log ( -----------  )   +  log ( -----------  )   |  
					 *      2            [       |u    |+ eps              |v    |+ eps     ]
					 *                              obs                       obs
					 *        dJ                              1      u                             1
					 * DU = - -- = - \bar{v}^2 log(u...) --------- ----  ~ - \bar{v}^2 log(u...) ------
					 *        du                         |u| + eps  |u|                           u + eps
					 */
					for(i=0;i<vnumnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux = - meanvel*meanvel * log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)) / (vx+epsvel);
							duy = - meanvel*meanvel * log((fabs(vy)+epsvel)/(fabs(vyobs)+epsvel)) / (vy+epsvel);
							pe->values[i*3+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
							pe->values[i*3+1]+=duy*weight*Jdet*gauss->weight*vbasis[i]; 
						}
						else{
							dux = - meanvel*meanvel * log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)) / (vx+epsvel);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*vbasis[i]; 
						}
					}
					break;
				case DragCoefficientAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAlongGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAcrossGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBbarAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBInitialguessMisfitEnum:
					/*Nothing in P vector*/
					break;
				default:
					_error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
			}
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,cs_list);

	/*Clean up and return*/
	xDelete<int>(cs_list);
	xDelete<int>(responses);
	xDelete<IssmDouble>(xyz_list_top);
	xDelete<IssmDouble>(vbasis);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* AdjointHorizAnalysis::CreatePVectorL1L2(Element* element){/*{{{*/

	/*Same as SSA*/
	return this->CreatePVectorSSA(element);
}/*}}}*/
ElementVector* AdjointHorizAnalysis::CreatePVectorHO(Element* element){/*{{{*/

	/*Nothing to be done if not on surface*/
	if(!element->IsOnSurface()) return NULL;

	/*Intermediaries */
	int        num_responses,i,domaintype;
	IssmDouble Jdet,obs_velocity_mag,velocity_mag;
	IssmDouble vx,vy,vxobs,vyobs,dux,duy,weight;
	IssmDouble scalex,scaley,scale,S;
	int        *responses    = NULL;
	IssmDouble *xyz_list_top = NULL;

	/* Get domaintype*/
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe    = element->NewElementVector(HOApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesTop(&xyz_list_top);
	element->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	element->FindParam(&responses,NULL,InversionCostFunctionsEnum);
	DatasetInput* weights_input = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* vx_input      = element->GetInput(VxEnum);                                 _assert_(vx_input);
	Input* vxobs_input   = element->GetInput(InversionVxObsEnum);                     _assert_(vxobs_input);
	Input* vy_input=NULL;
	Input* vyobs_input=NULL;
	if(domaintype!=Domain2DverticalEnum){
		vy_input      = element->GetInput(VyEnum);                                 _assert_(vy_input);
		vyobs_input   = element->GetInput(InversionVyObsEnum);                     _assert_(vyobs_input);
	}
	IssmDouble epsvel  = DBL_EPSILON;
	IssmDouble meanvel = 3.170979198376458e-05; /*1000 m/yr*/

	/*Get Surface if required by one response*/
	Input* S_input = NULL;
	for(int resp=0;resp<num_responses;resp++){
		if(responses[resp]==SurfaceAverageVelMisfitEnum){
			S_input = element->GetInput(SurfaceAreaEnum);  _assert_(S_input); break;
		}
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussTop(4);
	while(gauss->next()){

		element->JacobianDeterminantTop(&Jdet,xyz_list_top,gauss);
		element->NodalFunctions(basis, gauss);

		vx_input->GetInputValue(&vx,gauss);
		vxobs_input->GetInputValue(&vxobs,gauss);
		if(domaintype!=Domain2DverticalEnum){
			vy_input->GetInputValue(&vy,gauss);
			vyobs_input->GetInputValue(&vyobs,gauss);
		}
		/*Loop over all requested responses*/
		for(int resp=0;resp<num_responses;resp++){
			weights_input->GetInputValue(&weight,gauss,responses[resp]);

			switch(responses[resp]){
				case SurfaceAbsVelMisfitEnum:
					/*
					 *      1  [           2              2 ]
					 * J = --- | (u - u   )  +  (v - v   )  |
					 *      2  [       obs            obs   ]
					 *
					 *        dJ
					 * DU = - -- = (u   - u )
					 *        du     obs
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux=vxobs-vx;
							duy=vyobs-vy;
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{ 
							dux=vxobs-vx;
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceRelVelMisfitEnum:
					/*
					 *      1  [     \bar{v}^2             2   \bar{v}^2              2 ]
					 * J = --- | -------------  (u - u   ) + -------------  (v - v   )  |
					 *      2  [  (u   + eps)^2       obs    (v   + eps)^2       obs    ]
					 *              obs                        obs                      
					 *
					 *        dJ     \bar{v}^2
					 * DU = - -- = ------------- (u   - u )
					 *        du   (u   + eps)^2    obs
					 *               obs
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							scalex=pow(meanvel/(vxobs+epsvel),2); if(vxobs==0)scalex=0;
							scaley=pow(meanvel/(vyobs+epsvel),2); if(vyobs==0)scaley=0;
							dux=scalex*(vxobs-vx);
							duy=scaley*(vyobs-vy);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							scalex=pow(meanvel/(vxobs+epsvel),2); if(vxobs==0)scalex=0;
							dux=scalex*(vxobs-vx);
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceLogVelMisfitEnum:
					/*
					 *                 [        vel + eps     ] 2
					 * J = 4 \bar{v}^2 | log ( -----------  ) |  
					 *                 [       vel   + eps    ]
					 *                            obs
					 *
					 *        dJ                 2 * log(...)
					 * DU = - -- = - 4 \bar{v}^2 -------------  u
					 *        du                 vel^2 + eps
					 *            
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							velocity_mag    =sqrt(pow(vx,   2)+pow(vy,   2))+epsvel;
							obs_velocity_mag=sqrt(pow(vxobs,2)+pow(vyobs,2))+epsvel;
							scale=-8*pow(meanvel,2)/pow(velocity_mag,2)*log(velocity_mag/obs_velocity_mag);
							dux=scale*vx;
							duy=scale*vy;
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							velocity_mag    =fabs(vx)+epsvel;
							obs_velocity_mag=fabs(vxobs)+epsvel;
							scale=-8*pow(meanvel,2)/pow(velocity_mag,2)*log(velocity_mag/obs_velocity_mag);
							dux=scale*vx;
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceAverageVelMisfitEnum:
					/*
					 *      1                    2              2
					 * J = ---  sqrt(  (u - u   )  +  (v - v   )  )
					 *      S                obs            obs
					 *
					 *        dJ      1       1 
					 * DU = - -- = - --- ----------- * 2 (u - u   )
					 *        du      S  2 sqrt(...)           obs
					 */
					S_input->GetInputValue(&S,gauss);
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							scale=1./(S*sqrt(pow(vx-vxobs,2)+pow(vy-vyobs,2))+epsvel);
							dux=scale*(vxobs-vx);
							duy=scale*(vyobs-vy);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							scale=1./(S*2*fabs(vx-vxobs)+epsvel);
							dux=scale*(vxobs-vx);
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceLogVxVyMisfitEnum:
					/*
					 *      1            [        |u| + eps     2          |v| + eps     2  ]
					 * J = --- \bar{v}^2 | log ( -----------  )   +  log ( -----------  )   |  
					 *      2            [       |u    |+ eps              |v    |+ eps     ]
					 *                              obs                       obs
					 *        dJ                              1      u                             1
					 * DU = - -- = - \bar{v}^2 log(u...) --------- ----  ~ - \bar{v}^2 log(u...) ------
					 *        du                         |u| + eps  |u|                           u + eps
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux = - meanvel*meanvel * log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)) / (vx+epsvel);
							duy = - meanvel*meanvel * log((fabs(vy)+epsvel)/(fabs(vyobs)+epsvel)) / (vy+epsvel);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							dux = - meanvel*meanvel * log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)) / (vx+epsvel);
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case DragCoefficientAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAlongGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAcrossGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBbarAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBInitialguessMisfitEnum:
					/*Nothing in P vector*/
					break;
				default:
					_error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
			}
		}
	}

	/*Transform coordinate system*/
	if(domaintype!=Domain2DverticalEnum) element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<int>(responses);
	xDelete<IssmDouble>(xyz_list_top);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;

}/*}}}*/
ElementVector* AdjointHorizAnalysis::CreatePVectorMOLHO(Element* element){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

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
			basalelement = element->SpawnBasalElement(true);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries */
	int         num_responses,i;
	IssmDouble  Jdet,obs_velocity_mag,velocity_mag;
	IssmDouble  vx,vy,vxobs,vyobs,dux,duy,weight;
	IssmDouble scalex,scaley,scale,S,thickness,n;
	int        *responses = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe    = basalelement->NewElementVector(MOLHOApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	basalelement->FindParam(&responses,NULL,InversionCostFunctionsEnum);
	DatasetInput* weights_input = basalelement->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* vx_input      = basalelement->GetInput(VxSurfaceEnum);                                        _assert_(vx_input);
	Input* vxobs_input   = basalelement->GetInput(InversionVxObsEnum);                                   _assert_(vxobs_input);
   Input* thickness_input=basalelement->GetInput(ThicknessEnum);													  _assert_(thickness_input);
   Input*     n_input        =element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
	Input* vy_input=NULL;
	Input* vyobs_input=NULL;
	if(domaintype!=Domain2DverticalEnum){
		vy_input      = basalelement->GetInput(VySurfaceEnum);       _assert_(vy_input);
		vyobs_input   = basalelement->GetInput(InversionVyObsEnum);  _assert_(vyobs_input);
	}
	IssmDouble epsvel  = DBL_EPSILON;
	IssmDouble meanvel = 3.170979198376458e-05; /*1000 m/yr*/

	/*Get Surface if required by one response*/
	Input* S_input = NULL;
	for(int resp=0;resp<num_responses;resp++){
		if(responses[resp]==SurfaceAverageVelMisfitEnum){
			S_input = element->GetInput(SurfaceAreaEnum);  _assert_(S_input); break;
		}
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis, gauss);

		thickness_input->GetInputValue(&thickness,gauss);
      n_input->GetInputValue(&n,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vxobs_input->GetInputValue(&vxobs,gauss);
		if(domaintype!=Domain2DverticalEnum){
			vy_input->GetInputValue(&vy,gauss);
			vyobs_input->GetInputValue(&vyobs,gauss);
		}
		/*Loop over all requested responses*/
		for(int resp=0;resp<num_responses;resp++){
			weights_input->GetInputValue(&weight,gauss,responses[resp]);

			switch(responses[resp]){
				case SurfaceAbsVelMisfitEnum:
					/*
					 *      1  [           2              2 ]
					 * J = --- | (u - u   )  +  (v - v   )  |
					 *      2  [       obs            obs   ]
					 *
					 *        dJ
					 * DU = - -- = (u   - u )
					 *        du     obs
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux=vxobs-vx;
							duy=vyobs-vy;
							pe->values[i*4+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+1]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+2]+=duy*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+3]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else {
							_error_("2D vertical is not implemented for MOLHO");
						}
					}
					break;
				case SurfaceRelVelMisfitEnum:
					/*
					 *      1  [     \bar{v}^2             2   \bar{v}^2              2 ]
					 * J = --- | -------------  (u - u   ) + -------------  (v - v   )  |
					 *      2  [  (u   + eps)^2       obs    (v   + eps)^2       obs    ]
					 *              obs                        obs                      
					 *
					 *        dJ     \bar{v}^2
					 * DU = - -- = ------------- (u   - u )
					 *        du   (u   + eps)^2    obs
					 *               obs
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							scalex=pow(meanvel/(vxobs+epsvel),2); if(vxobs==0)scalex=0;
							scaley=pow(meanvel/(vyobs+epsvel),2); if(vyobs==0)scaley=0;
							dux=scalex*(vxobs-vx);
							duy=scaley*(vyobs-vy);
							pe->values[i*4+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+1]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+2]+=duy*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+3]+=duy*weight*Jdet*gauss->weight*basis[i];
						}
						else{
							_error_("2D vertical is not implemented for MOLHO");
						}
					}
					break;
				case SurfaceLogVelMisfitEnum:
					/*
					 *                 [        vel + eps     ] 2
					 * J = 4 \bar{v}^2 | log ( -----------  ) |  
					 *                 [       vel   + eps    ]
					 *                            obs
					 *
					 *        dJ                 2 * log(...)
					 * DU = - -- = - 4 \bar{v}^2 -------------  u
					 *        du                 vel^2 + eps
					 *            
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							velocity_mag    =sqrt(pow(vx,   2)+pow(vy,   2))+epsvel;
							obs_velocity_mag=sqrt(pow(vxobs,2)+pow(vyobs,2))+epsvel;
							scale=-8*pow(meanvel,2)/pow(velocity_mag,2)*log(velocity_mag/obs_velocity_mag);
							dux=scale*vx;
							duy=scale*vy;
							pe->values[i*4+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+1]+=dux*weight*Jdet*gauss->weight*basis[i];
							pe->values[i*4+2]+=duy*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+3]+=duy*weight*Jdet*gauss->weight*basis[i];
						}
						else{
							_error_("2D vertical is not implemented for MOLHO");
						}
					}
					break;
				case SurfaceAverageVelMisfitEnum:
					/*
					 *      1                    2              2
					 * J = ---  sqrt(  (u - u   )  +  (v - v   )  )
					 *      S                obs            obs
					 *
					 *        dJ      1       1 
					 * DU = - -- = - --- ----------- * 2 (u - u   )
					 *        du      S  2 sqrt(...)           obs
					 */
					S_input->GetInputValue(&S,gauss);
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							scale=1./(S*sqrt(pow(vx-vxobs,2)+pow(vy-vyobs,2))+epsvel);
							dux=scale*(vxobs-vx);
							duy=scale*(vyobs-vy);
							pe->values[i*4+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+1]+=dux*weight*Jdet*gauss->weight*basis[i];
							pe->values[i*4+2]+=duy*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+3]+=duy*weight*Jdet*gauss->weight*basis[i];
						}
						else{
							_error_("2D vertical is not implemented for MOLHO");
						}
					}
					break;
				case SurfaceLogVxVyMisfitEnum:
					/*
					 *      1            [        |u| + eps     2          |v| + eps     2  ]
					 * J = --- \bar{v}^2 | log ( -----------  )   +  log ( -----------  )   |  
					 *      2            [       |u    |+ eps              |v    |+ eps     ]
					 *                              obs                       obs
					 *        dJ                              1      u                             1
					 * DU = - -- = - \bar{v}^2 log(u...) --------- ----  ~ - \bar{v}^2 log(u...) ------
					 *        du                         |u| + eps  |u|                           u + eps
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux = - meanvel*meanvel * log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)) / (vx+epsvel);
							duy = - meanvel*meanvel * log((fabs(vy)+epsvel)/(fabs(vyobs)+epsvel)) / (vy+epsvel);
							pe->values[i*4+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+1]+=dux*weight*Jdet*gauss->weight*basis[i];
							pe->values[i*4+2]+=duy*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*4+3]+=duy*weight*Jdet*gauss->weight*basis[i];
						}
						else{
							_error_("2D vertical is not implemented for MOLHO");
						}
					}
					break;
				case DragCoefficientAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAlongGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAcrossGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBbarAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBInitialguessMisfitEnum:
					/*Nothing in P vector*/
					break;
				default:
					_error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
			}
		}
	}

	/*Transform coordinate system*/
	//if(domaintype!=Domain2DverticalEnum)	basalelement->TransformLoadVectorCoord(pe,XYEnum);
	/*Clean up and return*/
	xDelete<int>(responses);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* AdjointHorizAnalysis::CreatePVectorSSA(Element* element){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

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
			basalelement = element->SpawnBasalElement(true);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries */
	int         num_responses,i;
	IssmDouble  Jdet,obs_velocity_mag,velocity_mag;
	IssmDouble  vx,vy,vxobs,vyobs,dux,duy,weight;
	IssmDouble scalex,scaley,scale,S;
	int        *responses = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe    = basalelement->NewElementVector(SSAApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	basalelement->FindParam(&responses,NULL,InversionCostFunctionsEnum);
	DatasetInput* weights_input = basalelement->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* vx_input      = basalelement->GetInput(VxEnum);                                               _assert_(vx_input);
	Input* vxobs_input   = basalelement->GetInput(InversionVxObsEnum);                                   _assert_(vxobs_input);
	Input* vy_input=NULL;
	Input* vyobs_input=NULL;
	if(domaintype!=Domain2DverticalEnum){
		vy_input      = basalelement->GetInput(VyEnum);              _assert_(vy_input);
		vyobs_input   = basalelement->GetInput(InversionVyObsEnum);  _assert_(vyobs_input);
	}
	IssmDouble epsvel  = DBL_EPSILON;
	IssmDouble meanvel = 3.170979198376458e-05; /*1000 m/yr*/

	/*Get Surface if required by one response*/
	Input* S_input = NULL;
	for(int resp=0;resp<num_responses;resp++){
		if(responses[resp]==SurfaceAverageVelMisfitEnum){
			S_input = element->GetInput(SurfaceAreaEnum);  _assert_(S_input); break;
		}
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis, gauss);

		vx_input->GetInputValue(&vx,gauss);
		vxobs_input->GetInputValue(&vxobs,gauss);
		if(domaintype!=Domain2DverticalEnum){
			vy_input->GetInputValue(&vy,gauss);
			vyobs_input->GetInputValue(&vyobs,gauss);
		}
		/*Loop over all requested responses*/
		for(int resp=0;resp<num_responses;resp++){
			weights_input->GetInputValue(&weight,gauss,responses[resp]);

			switch(responses[resp]){
				case SurfaceAbsVelMisfitEnum:
					/*
					 *      1  [           2              2 ]
					 * J = --- | (u - u   )  +  (v - v   )  |
					 *      2  [       obs            obs   ]
					 *
					 *        dJ
					 * DU = - -- = (u   - u )
					 *        du     obs
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux=vxobs-vx;
							duy=vyobs-vy;
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else {
							dux=vxobs-vx;
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceRelVelMisfitEnum:
					/*
					 *      1  [     \bar{v}^2             2   \bar{v}^2              2 ]
					 * J = --- | -------------  (u - u   ) + -------------  (v - v   )  |
					 *      2  [  (u   + eps)^2       obs    (v   + eps)^2       obs    ]
					 *              obs                        obs                      
					 *
					 *        dJ     \bar{v}^2
					 * DU = - -- = ------------- (u   - u )
					 *        du   (u   + eps)^2    obs
					 *               obs
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							scalex=pow(meanvel/(vxobs+epsvel),2); if(vxobs==0)scalex=0;
							scaley=pow(meanvel/(vyobs+epsvel),2); if(vyobs==0)scaley=0;
							dux=scalex*(vxobs-vx);
							duy=scaley*(vyobs-vy);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							scalex=pow(meanvel/(vxobs+epsvel),2); if(vxobs==0)scalex=0;
							dux=scalex*(vxobs-vx);
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceLogVelMisfitEnum:
					/*
					 *                 [        vel + eps     ] 2
					 * J = 4 \bar{v}^2 | log ( -----------  ) |  
					 *                 [       vel   + eps    ]
					 *                            obs
					 *
					 *        dJ                 2 * log(...)
					 * DU = - -- = - 4 \bar{v}^2 -------------  u
					 *        du                 vel^2 + eps
					 *            
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							velocity_mag    =sqrt(pow(vx,   2)+pow(vy,   2))+epsvel;
							obs_velocity_mag=sqrt(pow(vxobs,2)+pow(vyobs,2))+epsvel;
							scale=-8*pow(meanvel,2)/pow(velocity_mag,2)*log(velocity_mag/obs_velocity_mag);
							dux=scale*vx;
							duy=scale*vy;
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							velocity_mag    =fabs(vx)+epsvel;
							obs_velocity_mag=fabs(vxobs)+epsvel;
							scale=-8*pow(meanvel,2)/pow(velocity_mag,2)*log(velocity_mag/obs_velocity_mag);
							dux=scale*vx;
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceAverageVelMisfitEnum:
					/*
					 *      1                    2              2
					 * J = ---  sqrt(  (u - u   )  +  (v - v   )  )
					 *      S                obs            obs
					 *
					 *        dJ      1       1 
					 * DU = - -- = - --- ----------- * 2 (u - u   )
					 *        du      S  2 sqrt(...)           obs
					 */
					S_input->GetInputValue(&S,gauss);
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							scale=1./(S*sqrt(pow(vx-vxobs,2)+pow(vy-vyobs,2))+epsvel);
							dux=scale*(vxobs-vx);
							duy=scale*(vyobs-vy);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							scale=1./(S*fabs(vx-vxobs)+epsvel);
							dux=scale*(vxobs-vx);
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case SurfaceLogVxVyMisfitEnum:
					/*
					 *      1            [        |u| + eps     2          |v| + eps     2  ]
					 * J = --- \bar{v}^2 | log ( -----------  )   +  log ( -----------  )   |  
					 *      2            [       |u    |+ eps              |v    |+ eps     ]
					 *                              obs                       obs
					 *        dJ                              1      u                             1
					 * DU = - -- = - \bar{v}^2 log(u...) --------- ----  ~ - \bar{v}^2 log(u...) ------
					 *        du                         |u| + eps  |u|                           u + eps
					 */
					for(i=0;i<numnodes;i++){
						if(domaintype!=Domain2DverticalEnum){
							dux = - meanvel*meanvel * log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)) / (vx+epsvel);
							duy = - meanvel*meanvel * log((fabs(vy)+epsvel)/(fabs(vyobs)+epsvel)) / (vy+epsvel);
							pe->values[i*2+0]+=dux*weight*Jdet*gauss->weight*basis[i]; 
							pe->values[i*2+1]+=duy*weight*Jdet*gauss->weight*basis[i]; 
						}
						else{
							dux = - meanvel*meanvel * log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)) / (vx+epsvel);
							pe->values[i]+=dux*weight*Jdet*gauss->weight*basis[i]; 
						}
					}
					break;
				case DragCoefficientAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAlongGradientEnum:
					/*Nothing in P vector*/
					break;
				case ThicknessAcrossGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBbarAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBAbsGradientEnum:
					/*Nothing in P vector*/
					break;
				case RheologyBInitialguessMisfitEnum:
					/*Nothing in P vector*/
					break;
				default:
					_error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
			}
		}
	}

	/*Transform coordinate system*/
	if(domaintype!=Domain2DverticalEnum)	basalelement->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<int>(responses);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete gauss;
	return pe;
}/*}}}*/
void           AdjointHorizAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           AdjointHorizAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	/*The gradient of the cost function is calculated in 2 parts.
	 *
	 * dJ    \partial J   \partial lambda^T(KU-F)
	 * --  = ---------- + ------------------------
	 * dk    \partial k   \parial k                  
	 *
	 * */

	/*If on water, grad = 0: */
	if(!element->IsIceInElement()) return;

	/*Get Approximation*/
	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);

	/*Get list of cost functions*/
	int *responses = NULL;
	int num_responses,resp;
	element->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	element->FindParam(&responses,NULL,InversionCostFunctionsEnum);

	/*Check that control_type is supported*/
	if(control_type!=MaterialsRheologyBbarEnum && 
		control_type!=FrictionCoefficientEnum   &&
		control_type!=FrictionCEnum             &&
		control_type!=FrictionAsEnum            && 
		control_type!=DamageDbarEnum            &&
		control_type!=MaterialsRheologyBEnum){
		_error_("Control "<<EnumToStringx(control_type)<<" not supported");
	}

	/*Deal with first part (partial derivative a J with respect to k)*/
	for(resp=0;resp<num_responses;resp++) switch(responses[resp]){
		case SurfaceAbsVelMisfitEnum:     /*Nothing, \partial J/\partial k = 0*/ break;
		case SurfaceRelVelMisfitEnum:     /*Nothing, \partial J/\partial k = 0*/ break;
		case SurfaceLogVelMisfitEnum:     /*Nothing, \partial J/\partial k = 0*/ break;
		case SurfaceLogVxVyMisfitEnum:    /*Nothing, \partial J/\partial k = 0*/ break;
		case SurfaceAverageVelMisfitEnum: /*Nothing, \partial J/\partial k = 0*/ break;
		case DragCoefficientAbsGradientEnum:   GradientJDragGradient(element,gradient,control_interp,control_index); break;
		case RheologyBbarAbsGradientEnum:      GradientJBbarGradient(element,gradient,control_interp,control_index); break;
		case RheologyBAbsGradientEnum:         GradientJBGradient(element,gradient,control_interp,control_index);    break;
		case RheologyBInitialguessMisfitEnum:  GradientJBinitial(element,gradient,control_interp,control_index);     break;
		default: _error_("response " << EnumToStringx(responses[resp]) << " not supported yet");
	}

	/*Deal with second term*/
	switch(control_type){
		case FrictionCoefficientEnum:
		case FrictionCEnum:
			switch(approximation){
				case SSAApproximationEnum: GradientJDragSSA(element,gradient,control_interp,control_index); break;
				case L1L2ApproximationEnum:GradientJDragL1L2(element,gradient,control_interp,control_index); break;
				case HOApproximationEnum:  GradientJDragHO( element,gradient,control_interp,control_index); break;
				case MOLHOApproximationEnum:  GradientJDragMOLHO( element,gradient,control_interp,control_index); break;
				case FSApproximationEnum:  GradientJDragFS( element,gradient,control_interp,control_index); break;
				case NoneApproximationEnum: /*Gradient is 0*/                    break;
				default: _error_("approximation " << EnumToStringx(approximation) << " not supported yet");
			}
			break;
		case FrictionAsEnum:
			switch(approximation){
				case SSAApproximationEnum: GradientJDragHydroSSA(element,gradient,control_interp,control_index); break;
				case L1L2ApproximationEnum:GradientJDragHydroL1L2(element,gradient,control_interp,control_index); break;
				case HOApproximationEnum:  GradientJDragHydroHO( element,gradient,control_interp,control_index); break;
				case FSApproximationEnum:  GradientJDragHydroFS( element,gradient,control_interp,control_index); break;
				case NoneApproximationEnum: /*Gradient is 0*/                    break;
				default: _error_("approximation " << EnumToStringx(approximation) << " not supported yet");
			}
			break;
		case MaterialsRheologyBbarEnum:
			switch(approximation){
				case SSAApproximationEnum: GradientJBbarSSA(element,gradient,control_interp,control_index); break;
				case L1L2ApproximationEnum:GradientJBbarL1L2(element,gradient,control_interp,control_index); break;
				case HOApproximationEnum:  GradientJBbarHO( element,gradient,control_interp,control_index); break;
				case MOLHOApproximationEnum:  GradientJBbarMOLHO( element,gradient,control_interp,control_index); break;
				case FSApproximationEnum:  GradientJBbarFS( element,gradient,control_interp,control_index); break;
				case NoneApproximationEnum: /*Gradient is 0*/                    break;
				default: _error_("approximation " << EnumToStringx(approximation) << " not supported yet");
			}
			break;
		case MaterialsRheologyBEnum:
			switch(approximation){
				case SSAApproximationEnum: GradientJBSSA(element,gradient,control_interp,control_index); break;
				case HOApproximationEnum:  GradientJBHO( element,gradient,control_interp,control_index); break;
			//	case MOLHOApproximationEnum:  GradientJBMOLHO( element,gradient,control_interp,control_index); break;
				case FSApproximationEnum:  GradientJBFS( element,gradient,control_interp,control_index); break;
				case NoneApproximationEnum: /*Gradient is 0*/                    break;
				default: _error_("approximation " << EnumToStringx(approximation) << " not supported yet");
			}
			break;
		case DamageDbarEnum:
			switch(approximation){
				case SSAApproximationEnum: GradientJDSSA(element,gradient,control_interp,control_index); break;
				case NoneApproximationEnum: /*Gradient is 0*/                 break;
				default: _error_("approximation " << EnumToStringx(approximation) << " not supported yet");
			}
			break;
		default: _error_("control type not supported yet: " << EnumToStringx(control_type));
	}

	/*Clean up and return*/
	xDelete<int>(responses);

}/*}}}*/
void           AdjointHorizAnalysis::GradientJBbarFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/
	/*WARNING: We use SSA as an estimate for now*/
	this->GradientJBbarSSA(element,gradient,control_interp,control_index);
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBbarGradient(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	if(control_interp==P0Enum) _error_("cannot require regularization for P0 controls");

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement(true);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble dk[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* dbasis        = xNew<IssmDouble>(2*numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* rheologyb_input = basalelement->GetInput(MaterialsRheologyBbarEnum);              _assert_(rheologyb_input);
	DatasetInput* weights_input   = basalelement->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1Derivatives(dbasis,xyz_list,gauss);
		weights_input->GetInputValue(&weight,gauss,RheologyBbarAbsGradientEnum);

		/*Build alpha_complement_list: */
		rheologyb_input->GetInputDerivativeValue(&dk[0],xyz_list,gauss);

		/*Build gradje_g_gaussian vector (actually -dJ/ddrag): */
		for(int i=0;i<numvertices;i++){
			if(domaintype!=Domain2DverticalEnum){
				ge[i]+=-weight*Jdet*gauss->weight*(dbasis[0*numvertices+i]*dk[0]+dbasis[1*numvertices+i]*dk[1]);
			}
			else{
				ge[i]+=-weight*Jdet*gauss->weight*dbasis[0*numvertices+i]*dk[0];
			}
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBbarL1L2(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	   /*Same as SSA*/
	   return this->GradientJBbarSSA(element,gradient,control_interp,control_index);
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBbarHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*WARNING: We use SSA as an estimate for now*/
	this->GradientJBbarSSA(element,gradient,control_interp,control_index);
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBbarMOLHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

   if(control_interp!=P1Enum) _error_("not implemented yet...");
	/*Intermediaries*/
	int      domaintype,dim;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement(true);
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble thickness,dmudB,n;
	IssmDouble dvx[3],dvy[3],dadjbx[3],dadjby[3],dadjshx[3],dadjshy[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);
	IssmDouble  zeta;
   IssmDouble  e1[3],e2[3], phishx[3], phishy[3];
   IssmDouble  epsilon[5];/* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble  adjshx, adjshy;

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* thickness_input	= basalelement->GetInput(ThicknessEnum);					_assert_(thickness_input);
   Input* n_input				= basalelement->GetInput(MaterialsRheologyNEnum);		_assert_(n_input);
   Input* vxbase_input		= basalelement->GetInput(VxBaseEnum);						_assert_(vxbase_input);
   Input* vybase_input		= basalelement->GetInput(VyBaseEnum);						_assert_(vybase_input);
   Input* vxshear_input		= basalelement->GetInput(VxShearEnum);						_assert_(vxshear_input);
   Input* vyshear_input		= basalelement->GetInput(VyShearEnum);						_assert_(vyshear_input);

	Input* adjointbx_input	= basalelement->GetInput(AdjointxBaseEnum);				_assert_(adjointbx_input);
	Input* adjointby_input	= basalelement->GetInput(AdjointyBaseEnum);				_assert_(adjointby_input);
	Input* adjointshx_input	= basalelement->GetInput(AdjointxShearEnum);				_assert_(adjointshx_input);
	Input* adjointshy_input	= basalelement->GetInput(AdjointyShearEnum);				_assert_(adjointshy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(5);
   GaussSeg* gauss_seg=new GaussSeg(2);

	while(gauss->next()){

		thickness_input->GetInputValue(&thickness,gauss);
	   n_input->GetInputValue(&n,gauss);
		adjointbx_input->GetInputDerivativeValue(&dadjbx[0],xyz_list,gauss);
		adjointby_input->GetInputDerivativeValue(&dadjby[0],xyz_list,gauss);
		adjointshx_input->GetInputDerivativeValue(&dadjshx[0],xyz_list,gauss);
		adjointshy_input->GetInputDerivativeValue(&dadjshy[0],xyz_list,gauss);
		adjointshx_input->GetInputValue(&adjshx, gauss);
		adjointshy_input->GetInputValue(&adjshy, gauss);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1(basis,gauss);

	   /* Get the integration in the vertical direction */
      gauss_seg->Reset();
      while(gauss_seg->next()){
			zeta=0.5*(gauss_seg->coord1+1);

			basalelement->StrainRateMOLHO(&epsilon[0],xyz_list,gauss,
                  vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input,zeta);

			basalelement->dViscositydBMOLHO(&dmudB,dim,xyz_list,gauss,vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input,zeta);

			e1[0] = 2.0*epsilon[0]+epsilon[1];
         e1[1] = epsilon[2];
         e1[2] = epsilon[3];

         e2[0] = epsilon[2];
         e2[1] = epsilon[0]+2.0*epsilon[1];
         e2[2] = epsilon[4];

			phishx[0] = dadjbx[0] + (1-pow(zeta, n+1))*dadjshx[0];
			phishx[1] = dadjbx[1] + (1-pow(zeta, n+1))*dadjshx[1];
			phishx[2] = (n+1)/thickness*pow(zeta, n)*adjshx;
			phishy[0] = dadjby[0] + (1-pow(zeta, n+1))*dadjshy[0];
			phishy[1] = dadjby[1] + (1-pow(zeta, n+1))*dadjshy[1];
			phishy[2] = (n+1)/thickness*pow(zeta, n)*adjshy;

			/*Build gradient vector (actually -dJ/dB): */
			for(int i=0;i<numvertices;i++){
				for (int j=0;j<3;j++){
					ge[i]+=(-dmudB)*2*(e1[j]*phishx[j]+e2[j]*phishy[j])*Jdet*gauss->weight*basis[i]*thickness*0.5*gauss_seg->weight;
					_assert_(!xIsNan<IssmDouble>(ge[i]));
				}
			}
		}
	}
	if(control_interp==P1Enum){
		gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);
	}
	else{
		_error_("not supported");
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete gauss_seg;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBbarSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	int      domaintype,dim;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement(true);
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble thickness,dmudB;
	IssmDouble dvx[3],dvy[3],dadjx[3],dadjy[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* thickness_input = basalelement->GetInput(ThicknessEnum);             _assert_(thickness_input);
	Input* vx_input        = basalelement->GetInput(VxEnum);                    _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyEnum);                    _assert_(vy_input);
	Input* adjointx_input  = basalelement->GetInput(AdjointxEnum);              _assert_(adjointx_input);
	Input* adjointy_input  = basalelement->GetInput(AdjointyEnum);              _assert_(adjointy_input);
	Input* rheologyb_input = basalelement->GetInput(MaterialsRheologyBbarEnum); _assert_(rheologyb_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		adjointx_input->GetInputDerivativeValue(&dadjx[0],xyz_list,gauss);
		adjointy_input->GetInputDerivativeValue(&dadjy[0],xyz_list,gauss);

		basalelement->dViscositydBSSA(&dmudB,dim,xyz_list,gauss,vx_input,vy_input);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dB): */
		if(control_interp==P1Enum){
			for(int i=0;i<numvertices;i++){
				ge[i]+=-dmudB*thickness*(
							(2*dvx[0]+dvy[1])*2*dadjx[0]+(dvx[1]+dvy[0])*(dadjx[1]+dadjy[0])+(2*dvy[1]+dvx[0])*2*dadjy[1]
							)*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
		else if(control_interp==P0Enum){
			ge[0]+=-dmudB*thickness*(
						(2*dvx[0]+dvy[1])*2*dadjx[0]+(dvx[1]+dvy[0])*(dadjx[1]+dadjy[0])+(2*dvy[1]+dvx[0])*2*dadjy[1]
						)*Jdet*gauss->weight;
			_assert_(!xIsNan<IssmDouble>(ge[0]));
		}
		else{
			_error_("not supported");
		}
	}
	if(control_interp==P1Enum){
		gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);
	}
	else if(control_interp==P0Enum){
		gradient->SetValue(vertexpidlist[0],ge[0],ADD_VAL);
	}
	else{
		_error_("not supported");
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/
	/*WARNING: We use HO as an estimate for now*/
	this->GradientJBHO(element,gradient,control_interp,control_index);
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBGradient(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	int      domaintype;

	if(control_interp==P0Enum) _error_("cannot require regularization for P0 controls");

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			break;
		case Domain2DverticalEnum:
			break;
		case Domain3DEnum:
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble dk[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* dbasis        = xNew<IssmDouble>(3*numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* rheology_input = element->GetInput(MaterialsRheologyBEnum);              _assert_(rheology_input);
	DatasetInput* weights_input   = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1Derivatives(dbasis,xyz_list,gauss);
		weights_input->GetInputValue(&weight,gauss,RheologyBAbsGradientEnum);

		/*Build alpha_complement_list: */
		rheology_input->GetInputDerivativeValue(&dk[0],xyz_list,gauss);

		/*Build gradje_g_gaussian vector (actually -dJ/ddrag): */
		for(int i=0;i<numvertices;i++){
			if(domaintype!=Domain2DverticalEnum){
				ge[i]+=-weight*Jdet*gauss->weight*(dbasis[0*numvertices+i]*dk[0]+dbasis[1*numvertices+i]*dk[1]);
			}
			else{
				ge[i]+=-weight*Jdet*gauss->weight*(dbasis[0*numvertices+i]*dk[0]+dbasis[1*numvertices+i]*dk[1]+dbasis[2*numvertices+i]*dk[2]);
			}
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBinitial(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			break;
		case Domain2DverticalEnum:
			break;
		case Domain3DEnum:
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble B,B0; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis        = xNew<IssmDouble>(numvertices);
	IssmDouble* ge           = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* rheology_input  = element->GetInput(MaterialsRheologyBbarEnum);              _assert_(rheology_input);
	Input* rheology0_input = element->GetInput(RheologyBInitialguessEnum);              _assert_(rheology0_input);
	DatasetInput* weights_input   = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);
		weights_input->GetInputValue(&weight,gauss,RheologyBInitialguessMisfitEnum);

		/*Build alpha_complement_list: */
		rheology_input->GetInputValue(&B,gauss);
		rheology0_input->GetInputValue(&B0,gauss);

		/*Build gradje_g_gaussian vector (actually -dJ/ddrag): */
		for(int i=0;i<numvertices;i++){
			ge[i]+=-weight*Jdet*gauss->weight*basis[i]*(B-B0);
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/
	/*Intermediaries*/
	int      domaintype,dim;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Get domaintype*/
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble thickness,dmudB;
	IssmDouble dvx[3],dvy[3],dadjx[3],dadjy[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinates(&xyz_list);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* thickness_input = element->GetInput(ThicknessEnum);             _assert_(thickness_input);
	Input* vx_input        = element->GetInput(VxEnum);                    _assert_(vx_input);
	Input* vy_input        = NULL;
	Input* adjointx_input  = element->GetInput(AdjointxEnum);              _assert_(adjointx_input);
	Input* adjointy_input  = NULL;
	Input* rheologyb_input = element->GetInput(MaterialsRheologyBEnum); _assert_(rheologyb_input);
	if(domaintype!=Domain2DverticalEnum){
		vy_input        = element->GetInput(VyEnum);                   _assert_(vy_input);
		adjointy_input  = element->GetInput(AdjointyEnum);             _assert_(adjointy_input);
	}
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	while(gauss->next()){

		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		adjointx_input->GetInputDerivativeValue(&dadjx[0],xyz_list,gauss);
		dim=2;
		if(domaintype!=Domain2DverticalEnum){
			adjointy_input->GetInputDerivativeValue(&dadjy[0],xyz_list, gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			dim=3;
		}

		element->dViscositydBHO(&dmudB,dim,xyz_list,gauss,vx_input,vy_input);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dB): */
		for(int i=0;i<numvertices;i++){
			if(domaintype!=Domain2DverticalEnum){
				ge[i]+=-dmudB*thickness*(
							(2*dvx[0]+dvy[1])*2*dadjx[0]+(dvx[1]+dvy[0])*(dadjx[1]+dadjy[0])+(2*dvy[1]+dvx[0])*2*dadjy[1]
							)*Jdet*gauss->weight*basis[i];
			}
			else{
				ge[i]+=-dmudB*thickness*4*dvx[0]*dadjx[0]*Jdet*gauss->weight*basis[i];
			}
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBMOLHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/
	_error_("not implemented yet...");
}/*}}}*/
void           AdjointHorizAnalysis::GradientJBSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	int      domaintype,dim;
	Element* basalelement;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble thickness,dmudB;
	IssmDouble dvx[3],dvy[3],dadjx[3],dadjy[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* thickness_input = basalelement->GetInput(ThicknessEnum);             _assert_(thickness_input);
	Input* vx_input        = basalelement->GetInput(VxEnum);                    _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyEnum);                    _assert_(vy_input);
	Input* adjointx_input  = basalelement->GetInput(AdjointxEnum);              _assert_(adjointx_input);
	Input* adjointy_input  = basalelement->GetInput(AdjointyEnum);              _assert_(adjointy_input);
	Input* rheologyb_input = basalelement->GetInput(MaterialsRheologyBEnum); _assert_(rheologyb_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		adjointx_input->GetInputDerivativeValue(&dadjx[0],xyz_list,gauss);
		adjointy_input->GetInputDerivativeValue(&dadjy[0],xyz_list,gauss);

		basalelement->dViscositydBSSA(&dmudB,dim,xyz_list,gauss,vx_input,vy_input);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dB): */
		for(int i=0;i<numvertices;i++){
			ge[i]+=-dmudB*thickness*(
						(2*dvx[0]+dvy[1])*2*dadjx[0]+(dvx[1]+dvy[0])*(dadjx[1]+dadjy[0])+(2*dvy[1]+dvx[0])*2*dadjy[1]
						)*Jdet*gauss->weight*basis[i];
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragGradient(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating (gradient is 0)*/
	if(element->IsAllFloating()) return;
	if(control_interp==P0Enum) _error_("cannot require regularization for P0 controls");

	/*Intermediaries*/
	int      domaintype,dim;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble dk[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* dbasis        = xNew<IssmDouble>(2*numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);

	/* get the friction law: if 2-Weertman, 11-Schoof, use a special name for the coefficient*/
	int frictionlaw; element->FindParam(&frictionlaw, FrictionLawEnum);
	Input* dragcoefficient_input;
	switch(frictionlaw) {
		case 2:
		case 11:
		case 13:
		case 14:
			dragcoefficient_input = basalelement->GetInput(FrictionCEnum); _assert_(dragcoefficient_input);
			break;
		default:
			dragcoefficient_input = basalelement->GetInput(FrictionCoefficientEnum); _assert_(dragcoefficient_input);
	}

	DatasetInput* weights_input         = basalelement->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1Derivatives(dbasis,xyz_list,gauss);
		weights_input->GetInputValue(&weight,gauss,DragCoefficientAbsGradientEnum);

		/*Build alpha_complement_list: */
		dragcoefficient_input->GetInputDerivativeValue(&dk[0],xyz_list,gauss);

		/*Build gradient vector (actually -dJ/ddrag): */
		for(int i=0;i<numvertices;i++){
			if(dim==2){
				ge[i]+=-weight*Jdet*gauss->weight*(dbasis[0*numvertices+i]*dk[0]+dbasis[1*numvertices+i]*dk[1]);
			}
			else{
				ge[i]+=-weight*Jdet*gauss->weight*dbasis[0*numvertices+i]*dk[0];
			}
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};

}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating or not on bed (gradient is 0)*/
	if(element->IsAllFloating()) return;
	if(!element->IsOnBase()) return;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Intermediaries*/
	int        domaintype,dim;
	IssmDouble Jdet,weight;
	IssmDouble drag,dalpha2dk,normal[3];
	IssmDouble vx,vy,vz,lambda,mu,xi;
	IssmDouble *xyz_list_base= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/* get domaintype */
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Build friction element, needed later: */
	if(domaintype!=Domain2DverticalEnum) dim=3;
	else dim=2;
	Friction* friction=new Friction(element,dim);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* vx_input        = element->GetInput(VxEnum);                   _assert_(vx_input);
	Input* vy_input        = element->GetInput(VyEnum);                   _assert_(vy_input);
	Input* adjointx_input  = element->GetInput(AdjointxEnum);             _assert_(adjointx_input);
	Input* adjointy_input  = element->GetInput(AdjointyEnum);             _assert_(adjointy_input);
	Input* vz_input        = NULL;
	Input* adjointz_input  = NULL;
	if(domaintype!=Domain2DverticalEnum){
		vz_input        = element->GetInput(VzEnum);                   _assert_(vy_input);
		adjointz_input  = element->GetInput(AdjointzEnum);             _assert_(adjointz_input);
	}
	/* get the friction law: 1- Budd, 11-Schoof*/
	int frictionlaw; element->FindParam(&frictionlaw, FrictionLawEnum);
	Input* dragcoeff_input = NULL;
	switch(frictionlaw) {
		case 1:
			dragcoeff_input = element->GetInput(FrictionCoefficientEnum); _assert_(dragcoeff_input);
			break;
		case 2:
		case 11:
		case 14:
			dragcoeff_input = element->GetInput(FrictionCEnum); _assert_(dragcoeff_input);
			break;
		default:
			_error_("Friction law "<< frictionlaw <<" not supported in the inversion.");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	while(gauss->next()){

		adjointx_input->GetInputValue(&lambda, gauss);
		adjointy_input->GetInputValue(&mu, gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		if(domaintype!=Domain2DverticalEnum){
			adjointz_input->GetInputValue(&xi    ,gauss);
			vz_input->GetInputValue(&vz,gauss);
		}
		dragcoeff_input->GetInputValue(&drag, gauss);

		friction->GetAlphaComplement(&dalpha2dk,gauss);
		element->NormalBase(&normal[0],xyz_list_base);

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dk): */
		if(domaintype!=Domain2DverticalEnum){
			for(int i=0;i<numvertices;i++){
				ge[i]+=(
							-lambda*(2*drag*dalpha2dk*(vx - vz*normal[0]*normal[2]))
							-mu    *(2*drag*dalpha2dk*(vy - vz*normal[1]*normal[2]))
							-xi    *(2*drag*dalpha2dk*(-vx*normal[0]*normal[2]-vy*normal[1]*normal[2]))
						 )*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
		else{
			for(int i=0;i<numvertices;i++){
				ge[i]+=(
							-lambda*2*drag*dalpha2dk*vx
							-mu    *2*drag*dalpha2dk*vy
						 )*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete friction;
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragL1L2(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Same as SSA*/
	return this->GradientJDragSSA(element,gradient,control_interp,control_index);
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating or not on bed (gradient is 0)*/
	if(element->IsAllFloating()) return;
	if(!element->IsOnBase()) return;

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble drag,dalpha2dk;
	IssmDouble vx,vy,lambda,mu;
	IssmDouble *xyz_list_base= NULL;

	int      domaintype,dim;
	element->FindParam(&domaintype,DomainTypeEnum);
	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Build friction element, needed later: */
	if(domaintype!=Domain2DverticalEnum) dim=3;
	else dim=2;
	Friction* friction=new Friction(element,dim);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* vx_input        = element->GetInput(VxEnum);                   _assert_(vx_input);
	Input* vy_input        = NULL;
	Input* adjointx_input  = element->GetInput(AdjointxEnum);             _assert_(adjointx_input);
	Input* adjointy_input  = NULL;
	if(domaintype!=Domain2DverticalEnum){
		vy_input        = element->GetInput(VyEnum);                   _assert_(vy_input);
		adjointy_input  = element->GetInput(AdjointyEnum);             _assert_(adjointy_input);
	}
	/* get the friction law: 1- Budd, 11-Schoof*/
	int frictionlaw; element->FindParam(&frictionlaw, FrictionLawEnum);
	Input* dragcoeff_input = NULL;
	switch(frictionlaw) {
		case 1:
			dragcoeff_input = element->GetInput(FrictionCoefficientEnum); _assert_(dragcoeff_input);
			break;
		case 2:
		case 11:
		case 14:
			dragcoeff_input = element->GetInput(FrictionCEnum); _assert_(dragcoeff_input);
			break;
		default:
			_error_("Friction law "<< frictionlaw <<" not supported in the inversion.");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	while(gauss->next()){

		adjointx_input->GetInputValue(&lambda, gauss);
		vx_input->GetInputValue(&vx,gauss);
		if(domaintype!=Domain2DverticalEnum){
			adjointy_input->GetInputValue(&mu, gauss);
			vy_input->GetInputValue(&vy,gauss);
		}
		dragcoeff_input->GetInputValue(&drag, gauss);

		friction->GetAlphaComplement(&dalpha2dk,gauss);

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dD): */
		if(control_interp==P1Enum){
			for(int i=0;i<numvertices;i++){
				if(domaintype!=Domain2DverticalEnum) ge[i]+=-2.*drag*dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight*basis[i];
				else ge[i]+=-2.*drag*dalpha2dk*(lambda*vx)*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
		else if(control_interp==P0Enum){
			if(domaintype!=Domain2DverticalEnum) ge[0]+=-2.*drag*dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight;
			else ge[0]+=-2.*drag*dalpha2dk*(lambda*vx)*Jdet*gauss->weight;
			_assert_(!xIsNan<IssmDouble>(ge[0]));
		}
		else{
			_error_("not supported");
		}
	}
	if(control_interp==P1Enum){
		gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);
	}
	else if(control_interp==P0Enum){
		gradient->SetValue(vertexpidlist[0],ge[0],ADD_VAL);
	}
	else{
		_error_("not supported");
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete friction;
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragMOLHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating (gradient is 0)*/
	if(element->IsAllFloating()) return;

	/*Intermediaries*/
	int      domaintype,dim;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble drag,dalpha2dk;
	IssmDouble vx,vy,lambda,mu;
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(basalelement,dim);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* vx_input        = basalelement->GetInput(VxBaseEnum);                   _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyBaseEnum);                   _assert_(vy_input);
	Input* adjointx_input  = basalelement->GetInput(AdjointxBaseEnum);             _assert_(adjointx_input);
	Input* adjointy_input  = basalelement->GetInput(AdjointyBaseEnum);             _assert_(adjointy_input);

	/* get the friction law: 1- Budd, 11-Schoof*/
	int frictionlaw;element->FindParam(&frictionlaw, FrictionLawEnum);
	Input* dragcoeff_input = NULL;
	switch(frictionlaw) {
		case 1:
			dragcoeff_input = basalelement->GetInput(FrictionCoefficientEnum); _assert_(dragcoeff_input);
			break;
		case 2:
		case 11:
		case 14:
			dragcoeff_input = basalelement->GetInput(FrictionCEnum); _assert_(dragcoeff_input);
			break;
		default:
			_error_("Friction law "<< frictionlaw <<" not supported in the inversion.");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		adjointx_input->GetInputValue(&lambda, gauss);
		adjointy_input->GetInputValue(&mu, gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		dragcoeff_input->GetInputValue(&drag, gauss);

		friction->GetAlphaComplement(&dalpha2dk,gauss);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dD): */
		if(control_interp==P1Enum){
			for(int i=0;i<numvertices;i++){
				ge[i]+=-2.*drag*dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
		else if(control_interp==P0Enum){
			ge[0]+=-2.*drag*dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight;
			_assert_(!xIsNan<IssmDouble>(ge[0]));
		}
		else{
			_error_("not supported");
		}
	}
	if(control_interp==P1Enum){
		gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);
	}
	else if(control_interp==P0Enum){
		gradient->SetValue(vertexpidlist[0],ge[0],ADD_VAL);
	}
	else{
		_error_("not supported");
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete friction;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating (gradient is 0)*/
	if(element->IsAllFloating()) return;

	/*Intermediaries*/
	int      domaintype,dim;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble drag,dalpha2dk;
	IssmDouble vx,vy,lambda,mu;
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(basalelement,dim);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* vx_input        = basalelement->GetInput(VxEnum);                   _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyEnum);                   _assert_(vy_input);
	Input* adjointx_input  = basalelement->GetInput(AdjointxEnum);             _assert_(adjointx_input);
	Input* adjointy_input  = basalelement->GetInput(AdjointyEnum);             _assert_(adjointy_input);

	/* get the friction law: 1- Budd, 11-Schoof*/
	int frictionlaw;element->FindParam(&frictionlaw, FrictionLawEnum);
	Input* dragcoeff_input = NULL;
	switch(frictionlaw) {
		case 1:
			dragcoeff_input = basalelement->GetInput(FrictionCoefficientEnum); _assert_(dragcoeff_input);
			break;
		case 2:
		case 11:
		case 13:
		case 14:
			dragcoeff_input = basalelement->GetInput(FrictionCEnum); _assert_(dragcoeff_input);
			break;
		default:
			_error_("Friction law "<< frictionlaw <<" not supported in the inversion.");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		adjointx_input->GetInputValue(&lambda, gauss);
		adjointy_input->GetInputValue(&mu, gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		dragcoeff_input->GetInputValue(&drag, gauss);

		friction->GetAlphaComplement(&dalpha2dk,gauss);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dD): */
		if(control_interp==P1Enum){
			for(int i=0;i<numvertices;i++){
				ge[i]+=-2.*drag*dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
		else if(control_interp==P0Enum){
			ge[0]+=-2.*drag*dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight;
			_assert_(!xIsNan<IssmDouble>(ge[0]));
		}
		else{
			_error_("not supported");
		}
	}
	if(control_interp==P1Enum){
		gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);
	}
	else if(control_interp==P0Enum){
		gradient->SetValue(vertexpidlist[0],ge[0],ADD_VAL);
	}
	else{
		_error_("not supported");
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete friction;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragHydroFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating or not on bed (gradient is 0)*/
	if(element->IsAllFloating()) return;
	if(!element->IsOnBase()) return;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Intermediaries*/
	int        domaintype,dim;
	IssmDouble Jdet,weight;
	IssmDouble drag,dalpha2dk,normal[3];
	IssmDouble vx,vy,vz,lambda,mu,xi;
	IssmDouble *xyz_list_base= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/* get domaintype */
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Build friction element, needed later: */
	if(domaintype!=Domain2DverticalEnum) dim=3;
	else dim=2;
	Friction* friction=new Friction(element,dim);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* vx_input        = element->GetInput(VxEnum);                   _assert_(vx_input);
	Input* vy_input        = element->GetInput(VyEnum);                   _assert_(vy_input);
	Input* adjointx_input  = element->GetInput(AdjointxEnum);             _assert_(adjointx_input);
	Input* adjointy_input  = element->GetInput(AdjointyEnum);             _assert_(adjointy_input);
	Input* vz_input        = NULL;
	Input* adjointz_input  = NULL;
	if(domaintype!=Domain2DverticalEnum){
		vz_input        = element->GetInput(VzEnum);                   _assert_(vy_input);
		adjointz_input  = element->GetInput(AdjointzEnum);             _assert_(adjointz_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	while(gauss->next()){

		adjointx_input->GetInputValue(&lambda, gauss);
		adjointy_input->GetInputValue(&mu, gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		if(domaintype!=Domain2DverticalEnum){
			adjointz_input->GetInputValue(&xi    ,gauss);
			vz_input->GetInputValue(&vz,gauss);
		}

		friction->GetAlphaComplement(&dalpha2dk,gauss);
		element->NormalBase(&normal[0],xyz_list_base);

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dk): */
		if(domaintype!=Domain2DverticalEnum){
			for(int i=0;i<numvertices;i++){
				ge[i]+=(
							-lambda*(dalpha2dk*(vx - vz*normal[0]*normal[2]))
							-mu    *(dalpha2dk*(vy - vz*normal[1]*normal[2]))
							-xi    *(dalpha2dk*(-vx*normal[0]*normal[2]-vy*normal[1]*normal[2]))
						 )*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
		else{
			for(int i=0;i<numvertices;i++){
				ge[i]+=(
							-lambda*dalpha2dk*vx
							-mu    *dalpha2dk*vy
						 )*Jdet*gauss->weight*basis[i];
				_assert_(!xIsNan<IssmDouble>(ge[i]));
			}
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete friction;
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragHydroL1L2(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Same as SSA*/
	return this->GradientJDragSSA(element,gradient,control_interp,control_index);
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragHydroHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating or not on bed (gradient is 0)*/
	if(element->IsAllFloating()) return;
	if(!element->IsOnBase()) return;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble drag,dalpha2dk;
	IssmDouble vx,vy,lambda,mu;
	IssmDouble *xyz_list_base= NULL;

	int      domaintype,dim;
	element->FindParam(&domaintype,DomainTypeEnum);
	/*Fetch number of vertices for this finite element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Build friction element, needed later: */
	if(domaintype!=Domain2DverticalEnum) dim=3;
	else dim=2;
	Friction* friction=new Friction(element,dim);

	/*Retrieve all inputs we will be needing: */
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GradientIndexing(&vertexpidlist[0],control_index);
	Input* vx_input        = element->GetInput(VxEnum);                   _assert_(vx_input);
	Input* vy_input        = NULL;
	Input* adjointx_input  = element->GetInput(AdjointxEnum);             _assert_(adjointx_input);
	Input* adjointy_input  = NULL;
	if(domaintype!=Domain2DverticalEnum){
		vy_input        = element->GetInput(VyEnum);                   _assert_(vy_input);
		adjointy_input  = element->GetInput(AdjointyEnum);             _assert_(adjointy_input);
	}
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	while(gauss->next()){

		adjointx_input->GetInputValue(&lambda, gauss);
		vx_input->GetInputValue(&vx,gauss);
		if(domaintype!=Domain2DverticalEnum){
			adjointy_input->GetInputValue(&mu, gauss);
			vy_input->GetInputValue(&vy,gauss);
		}

		friction->GetAlphaComplement(&dalpha2dk,gauss);

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsP1(basis,gauss);

		/*Build gradient vector (actually -dJ/dD): */
		for(int i=0;i<numvertices;i++){
			if(domaintype!=Domain2DverticalEnum) ge[i]+=-dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight*basis[i];
			else ge[i]+=-dalpha2dk*(lambda*vx)*Jdet*gauss->weight*basis[i];
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete friction;
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDragHydroSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*return if floating (gradient is 0)*/
	if(element->IsAllFloating()) return;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Intermediaries*/
	int      domaintype,dim;

	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble dalpha2dk;
	IssmDouble vx,vy,lambda,mu;
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(basalelement,dim);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* vx_input        = basalelement->GetInput(VxEnum);          _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyEnum);          _assert_(vy_input);
	Input* adjointx_input  = basalelement->GetInput(AdjointxEnum);    _assert_(adjointx_input);
	Input* adjointy_input  = basalelement->GetInput(AdjointyEnum);    _assert_(adjointy_input);

	IssmDouble  q_exp;
	IssmDouble  C_param;
	IssmDouble  As;
	IssmDouble  Neff;
	IssmDouble  n;
	IssmDouble  alpha;
	IssmDouble  Chi,Gamma;
	IssmDouble  vz,vmag;
	IssmDouble  Uder;

	/*Recover parameters: */
	Input* qinput = basalelement->GetInput(FrictionQEnum);
	Input* cinput = basalelement->GetInput(FrictionCEnum);
	Input* Asinput = basalelement->GetInput(FrictionAsEnum);
	Input* nInput =basalelement->GetInput(MaterialsRheologyNEnum);
	Input* Ninput = basalelement->GetInput(FrictionEffectivePressureEnum);	
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		adjointx_input->GetInputValue(&lambda, gauss);
		adjointy_input->GetInputValue(&mu, gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);

		friction->GetAlphaComplement(&dalpha2dk,gauss);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1(basis,gauss);

		/*Dealing with dalpha/du*/
		qinput->GetInputValue(&q_exp,gauss);
		cinput->GetInputValue(&C_param,gauss);
		Asinput->GetInputValue(&As,gauss);
		Ninput->GetInputValue(&Neff,gauss);
		nInput->GetInputValue(&n,gauss);

		if (q_exp==1){
			alpha=1;
		}
		else{
			alpha=(pow(q_exp-1,q_exp-1))/pow(q_exp,q_exp);
		}

		vmag  = sqrt(vx*vx + vy*vy);
		Chi   = vmag/(pow(C_param,n)*pow(Neff,n)*As);
		Gamma = (Chi/(1.+alpha*pow(Chi,q_exp)));

		Uder =Neff*C_param/(vmag*vmag*n) *
			(Gamma-alpha*q_exp*pow(Chi,q_exp-1.)*Gamma*Gamma* pow(Gamma,(1.-n)/n) -
			 n* pow(Gamma,1./n));

		/*Build gradient vector (actually -dJ/dD): */
		for(int i=0;i<numvertices;i++){
			ge[i]+=-dalpha2dk*((lambda*vx+mu*vy))*Jdet*gauss->weight*basis[i];
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	delete friction;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::GradientJDSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index){/*{{{*/

	/*Intermediaries*/
	int      domaintype,dim;
	Element* basalelement;
	if(control_interp!=P1Enum) _error_("not implemented yet...");

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim          = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			dim          = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement(true);
			dim          = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries*/
	IssmDouble Jdet,weight;
	IssmDouble thickness,dmudD;
	IssmDouble dvx[3],dvy[3],dadjx[3],dadjy[3]; 
	IssmDouble *xyz_list= NULL;

	/*Fetch number of vertices for this finite element*/
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize some vectors*/
	IssmDouble* basis         = xNew<IssmDouble>(numvertices);
	IssmDouble* ge            = xNewZeroInit<IssmDouble>(numvertices);
	int*        vertexpidlist = xNew<int>(numvertices);

	/*Retrieve all inputs we will be needing: */
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->GradientIndexing(&vertexpidlist[0],control_index);
	Input* thickness_input = basalelement->GetInput(ThicknessEnum);             _assert_(thickness_input);
	Input* vx_input        = basalelement->GetInput(VxEnum);                    _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyEnum);                    _assert_(vy_input);
	Input* adjointx_input  = basalelement->GetInput(AdjointxEnum);              _assert_(adjointx_input);
	Input* adjointy_input  = basalelement->GetInput(AdjointyEnum);              _assert_(adjointy_input);
	Input* rheologyb_input = basalelement->GetInput(MaterialsRheologyBbarEnum); _assert_(rheologyb_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(4);
	while(gauss->next()){

		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		adjointx_input->GetInputDerivativeValue(&dadjx[0],xyz_list,gauss);
		adjointy_input->GetInputDerivativeValue(&dadjy[0],xyz_list,gauss);

		basalelement->dViscositydDSSA(&dmudD,dim,xyz_list,gauss,vx_input,vy_input);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsP1(basis,gauss);

		for(int i=0;i<numvertices;i++){
			ge[i]+=-dmudD*thickness*(
						(2*dvx[0]+dvy[1])*2*dadjx[0]+(dvx[1]+dvy[0])*(dadjx[1]+dadjy[0])+(2*dvy[1]+dvx[0])*2*dadjy[1]
						)*Jdet*gauss->weight*basis[i];
			_assert_(!xIsNan<IssmDouble>(ge[i]));
		}
	}
	gradient->SetValues(numvertices,vertexpidlist,ge,ADD_VAL);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(ge);
	xDelete<int>(vertexpidlist);
	delete gauss;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           AdjointHorizAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation==FSApproximationEnum || approximation==NoneApproximationEnum){
		InputUpdateFromSolutionFS(solution,element);
	}
	else if (approximation==MOLHOApproximationEnum) {
		InputUpdateFromSolutionMOLHO(solution, element);
	}
	else{
		InputUpdateFromSolutionHoriz(solution,element);
	}
}/*}}}*/
void           AdjointHorizAnalysis::InputUpdateFromSolutionFS(IssmDouble* solution,Element* element){/*{{{*/
	int          i,fe_FS;
	int*         vdoflist=NULL;
	int*         pdoflist=NULL;
	IssmDouble   FSreconditioning;

	int      domaintype,dim;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			dim          = 3;
			break;
		case Domain2DverticalEnum:
			dim          = 2;
			break;
		case Domain3DEnum:
			dim          = 3;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int vnumdof   = vnumnodes*dim;
	int pnumdof   = pnumnodes*1;

	/*Initialize values*/
	IssmDouble* values  = xNew<IssmDouble>(vnumdof+pnumdof);
	IssmDouble* lambdax = xNew<IssmDouble>(vnumnodes);
	IssmDouble* lambday = xNew<IssmDouble>(vnumnodes);
	IssmDouble* lambdaz = xNew<IssmDouble>(vnumnodes);
	IssmDouble* lambdap = xNew<IssmDouble>(pnumnodes);

	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Get dof list: */
	element->GetDofListLocalVelocity(&vdoflist,GsetEnum);
	element->GetDofListLocalPressure(&pdoflist,GsetEnum);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<vnumdof;i++) values[i]        =solution[vdoflist[i]];
	for(i=0;i<pnumdof;i++) values[vnumdof+i]=solution[pdoflist[i]];

	/*Transform solution in Cartesian Space*/
	element->TransformSolutionCoord(values,cs_list);

	/*fill in all arrays: */
	for(i=0;i<vnumnodes;i++){
		lambdax[i] = values[i*dim+0];
		if(xIsNan<IssmDouble>(lambdax[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambdax[i])) _error_("Inf found in solution vector");
		lambday[i] = values[i*dim+1];
		if(xIsNan<IssmDouble>(lambday[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambday[i])) _error_("Inf found in solution vector");
		if(dim==3){
			lambdaz[i] = values[i*dim+2];
			if(xIsNan<IssmDouble>(lambdaz[i])) _error_("NaN found in solution vector");
			if(xIsInf<IssmDouble>(lambdaz[i])) _error_("Inf found in solution vector");
		}
	}
	for(i=0;i<pnumnodes;i++){
		lambdap[i] = values[vnumdof+i];
		if(xIsNan<IssmDouble>(lambdap[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambdap[i])) _error_("Inf found in solution vector");
	}

	/*Recondition pressure and compute vel: */
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	for(i=0;i<pnumnodes;i++) lambdap[i]=lambdap[i]*FSreconditioning;

	/*Add vx and vy as inputs to the tria element: */
	element->AddInput(AdjointxEnum,lambdax,element->VelocityInterpolation());
	element->AddInput(AdjointyEnum,lambday,element->VelocityInterpolation());
	if(domaintype!=Domain2DverticalEnum) element->AddInput(AdjointzEnum,lambdaz,element->VelocityInterpolation());

	element->FindParam(&fe_FS,FlowequationFeFSEnum);
	if(fe_FS!=LATaylorHoodEnum && fe_FS!=LACrouzeixRaviartEnum)	
	 element->AddInput(AdjointpEnum,lambdap,element->PressureInterpolation());	

	/*Free resources:*/
	xDelete<int>(vdoflist);
	xDelete<int>(pdoflist);
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(lambdap);
	xDelete<IssmDouble>(lambdaz);
	xDelete<IssmDouble>(lambday);
	xDelete<IssmDouble>(lambdax);
	xDelete<IssmDouble>(values);
}/*}}}*/
void           AdjointHorizAnalysis::InputUpdateFromSolutionHoriz(IssmDouble* solution,Element* element){/*{{{*/
	int  i;
	int* doflist=NULL;

	int    domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof;
	if(domaintype!=Domain2DverticalEnum)  numdof   = numnodes*2;
	else			                          numdof   = numnodes*1;
	/*Fetch dof list and allocate solution vectors*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values  = xNew<IssmDouble>(numdof);
	IssmDouble* lambdax = xNew<IssmDouble>(numnodes);
	IssmDouble* lambday = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	if(domaintype!=Domain2DverticalEnum)	element->TransformSolutionCoord(&values[0],XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		if(domaintype!=Domain2DverticalEnum){
			lambdax[i]=values[i*2+0];
			lambday[i]=values[i*2+1];
		}
		else {lambdax[i]=values[i];lambday[i]=0;}
		/*Check solution*/
		if(xIsNan<IssmDouble>(lambdax[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambdax[i])) _error_("Inf found in solution vector");
		if(domaintype!=Domain2DverticalEnum && xIsNan<IssmDouble>(lambday[i])) _error_("NaN found in solution vector");
		if(domaintype!=Domain2DverticalEnum && xIsInf<IssmDouble>(lambday[i])) _error_("Inf found in solution vector");
	}

	/*Add vx and vy as inputs to the tria element: */
	element->AddInput(AdjointxEnum,lambdax,element->GetElementType());
	element->AddInput(AdjointyEnum,lambday,element->GetElementType());

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(lambdax);
	xDelete<IssmDouble>(lambday);
	xDelete<int>(doflist);
}/*}}}*/
void           AdjointHorizAnalysis::InputUpdateFromSolutionMOLHO(IssmDouble* solution,Element* element){/*{{{*/
	int  i;
	int* doflist=NULL;

	int    domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if (domaintype!=Domain2DhorizontalEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof = numnodes * 4;

	/*Fetch dof list and allocate solution vectors*/
	element->GetDofListLocal(&doflist,MOLHOApproximationEnum,GsetEnum);
	IssmDouble* values    = xNew<IssmDouble>(numdof);
	IssmDouble* lambdax   = xNew<IssmDouble>(numnodes);
	IssmDouble* lambday   = xNew<IssmDouble>(numnodes);
	IssmDouble* lambdabx  = xNew<IssmDouble>(numnodes);
	IssmDouble* lambdaby  = xNew<IssmDouble>(numnodes);
	IssmDouble* lambdashx = xNew<IssmDouble>(numnodes);
	IssmDouble* lambdashy = xNew<IssmDouble>(numnodes);
	IssmDouble* n         = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	if(domaintype!=Domain2DverticalEnum)	element->TransformSolutionCoord(&values[0],XYEnum);

	element->GetInputListOnNodes(&n[0],MaterialsRheologyNEnum,0.);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		lambdabx[i] =values[i*4+0];
		lambdashx[i]=values[i*4+1];
		/*Check solution*/
		if(xIsNan<IssmDouble>(lambdabx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambdabx[i])) _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(lambdashx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambdashx[i])) _error_("Inf found in solution vector");
		/* adjoint for the surface velocity */
		lambdax[i] = lambdabx[i] + lambdashx[i]*(n[i]+1)/(n[i]+2);

		lambdaby[i] =values[i*4+2];
		lambdashy[i]=values[i*4+3];
		/*Check solution*/
		if(xIsNan<IssmDouble>(lambdaby[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambdaby[i])) _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(lambdashy[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(lambdashy[i])) _error_("Inf found in solution vector");
		/* adjoint for the surface velocity */
		lambday[i] = lambdaby[i] + lambdashy[i]*(n[i]+1)/(n[i]+2);
	}

	/*Add vx and vy as inputs to the tria element: */
	element->AddInput(AdjointxBaseEnum,lambdabx,element->GetElementType());
	element->AddInput(AdjointyBaseEnum,lambdaby,element->GetElementType());
	element->AddInput(AdjointxShearEnum,lambdashx,element->GetElementType());
	element->AddInput(AdjointyShearEnum,lambdashy,element->GetElementType());

	element->AddInput(AdjointxEnum,lambdax,element->GetElementType());
	element->AddInput(AdjointyEnum,lambday,element->GetElementType());

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(lambdax);
	xDelete<IssmDouble>(lambday);
	xDelete<IssmDouble>(lambdabx);
	xDelete<IssmDouble>(lambdaby);
	xDelete<IssmDouble>(lambdashx);
	xDelete<IssmDouble>(lambdashy);
   xDelete<IssmDouble>(n);
	xDelete<int>(doflist);
}/*}}}*/
void           AdjointHorizAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
