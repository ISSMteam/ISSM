/*!\file Matice.c
 * \brief: implementation of the Matice object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Matice.h"
#include "./Materials.h"
#include "../Elements/Element.h"
#include "../Elements/Tria.h"
#include "../Elements/Penta.h"
#include "../Params/Parameters.h"
#include "../Vertex.h"
#include "../Hook.h"
#include "../Node.h"
#include "../IoModel.h"
#include "../../shared/shared.h"

/*Matice constructors and destructor*/
Matice::Matice(){/*{{{*/
	this->helement=NULL;
	this->element=NULL;
	this->isdamaged=false;
	this->isenhanced=false;
	return;
}
/*}}}*/
Matice::Matice(int matice_mid,int index, IoModel* iomodel){/*{{{*/

	 /*Get material type and initialize object*/
   int materialtype;
   iomodel->FindConstant(&materialtype,"md.materials.type");
	this->Init(matice_mid,index,materialtype);

}
/*}}}*/
Matice::Matice(int matice_mid,int index,int materialtype){/*{{{*/

	this->Init(matice_mid,index,materialtype);
	return;
} /*}}}*/
Matice::~Matice(){/*{{{*/
	delete helement;
	return;
}
/*}}}*/
void Matice::Init(int matice_mid,int index,int materialtype){/*{{{*/

	/*Initialize id*/
	this->mid=matice_mid;

	/*Hooks: */
	int matice_eid=index+1;
	this->helement=new Hook(&matice_eid,1);
	this->element=NULL;

	/*Material specific properties*/
	switch(materialtype){
		case MatdamageiceEnum:
			this->isdamaged = true;
			this->isenhanced = false;
			break;
		case MaticeEnum:
			this->isdamaged = false;
			this->isenhanced = false;
			break;
		case MatenhancediceEnum:
			this->isdamaged = false;
			this->isenhanced = true;
			break;
		default:
			_error_("Material type not recognized");
	}

	return;
} /*}}}*/

/*Object virtual functions definitions:*/
Object*   Matice::copy() {/*{{{*/

	/*Output*/
	Matice* matice=NULL;

	/*Initialize output*/
	matice=new Matice();

	/*copy fields: */
	matice->mid=this->mid;
	matice->helement=(Hook*)this->helement->copy();
	matice->element =(Element*)this->helement->delivers();
	matice->isdamaged = this->isdamaged;
	matice->isenhanced = this->isenhanced;

	return matice;
}
/*}}}*/
Material* Matice::copy2(Element* element_in) {/*{{{*/

	/*Output*/
	Matice* matice=NULL;

	/*Initialize output*/
	matice=new Matice();

	/*copy fields: */
	matice->mid=this->mid;
	matice->helement=(Hook*)this->helement->copy();
	matice->element =element_in;
	matice->isdamaged = this->isdamaged;
	matice->isenhanced = this->isenhanced;

	return matice;
}
/*}}}*/
void      Matice::DeepEcho(void){/*{{{*/

	_printf_("Matice:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   isdamaged: " << isdamaged << "\n");
	_printf_("   isenhanced: " << isenhanced << "\n");

	/*helement and element DeepEcho were commented to avoid recursion.*/
	/*Example: element->DeepEcho calls matice->DeepEcho which calls element->DeepEcho etc*/
	_printf_("   helement:\n");
	_printf_("		note: helement not printed to avoid recursion.\n");
	//if(helement) helement->DeepEcho();
	//else _printf_("   helement = NULL\n");

	_printf_("   element:\n");
	_printf_("     note: element not printed to avoid recursion.\n");
	//if(element) element->DeepEcho();
	//else _printf_("   element = NULL\n");
}		
/*}}}*/
void      Matice::Echo(void){/*{{{*/

	_printf_("Matice:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   isdamaged: " << isdamaged << "\n");
	_printf_("   isenhanced: " << isenhanced << "\n");

	/*helement and element Echo were commented to avoid recursion.*/
	/*Example: element->Echo calls matice->Echo which calls element->Echo etc*/
	_printf_("   helement:\n");
	_printf_("     note: helement not printed to avoid recursion.\n");
	//if(helement) helement->Echo();
	//else _printf_("   helement = NULL\n");

	_printf_("   element:\n");
	_printf_("     note: element not printed to avoid recursion.\n");
	//if(element) element->Echo();
	//else _printf_("   element = NULL\n");
}
/*}}}*/
int       Matice::Id(void){ return mid; }/*{{{*/
/*}}}*/
void      Matice::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD)helement=new Hook(); 

	int object_enum = MaticeEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->mid);
	marshallhandle->call(this->isdamaged);
	marshallhandle->call(this->isenhanced);
	this->helement->Marshall(marshallhandle);
	this->element=(Element*)this->helement->delivers();
}/*}}}*/
int       Matice::ObjectEnum(void){/*{{{*/

	return MaticeEnum;

}
/*}}}*/

/*Matice management*/
void  Matice::Configure(Elements* elementsin){/*{{{*/

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	helement->configure((DataSet*)elementsin);
	this->element  = (Element*)helement->delivers();
}
/*}}}*/
IssmDouble Matice::GetA(Gauss* gauss){/*{{{*/
	/*
	 * A = 1/B^n
	 */

	IssmDouble B=this->GetB(gauss);
	IssmDouble n=this->GetN();

	return pow(B,-n);
}
/*}}}*/
IssmDouble Matice::GetAbar(Gauss* gauss){/*{{{*/
	/*
	 * A = 1/B^n
	 */

	IssmDouble B=this->GetBbar(gauss);
	IssmDouble n=this->GetN();

	return pow(B,-n);
}
/*}}}*/
IssmDouble Matice::GetB(Gauss* gauss){/*{{{*/

	_assert_(gauss); 

	/*Output*/
	IssmDouble B;
	Input* B_input = element->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	B_input->GetInputValue(&B,gauss);
	return B;
}
/*}}}*/
IssmDouble Matice::GetBbar(Gauss* gauss){/*{{{*/

	_assert_(gauss); 

	/*Output*/
	IssmDouble Bbar;

	Input* B_input = element->GetInput(MaterialsRheologyBbarEnum); _assert_(B_input);
	B_input->GetInputValue(&Bbar,gauss);
	return Bbar;
}
/*}}}*/
IssmDouble Matice::GetD(Gauss* gauss){/*{{{*/

	_assert_(this->isdamaged);
	/*Output*/
	IssmDouble D;
	if(this->isdamaged){
		Input* D_input = element->GetInput(DamageDEnum); _assert_(D_input);
		D_input->GetInputValue(&D,gauss);
	}
	else{
		_error_("Cannot get DamageD for non damaged ice");
	}
	return D;
}
/*}}}*/
IssmDouble Matice::GetDbar(Gauss* gauss){/*{{{*/

	_assert_(this->isdamaged);
	/*Output*/
	IssmDouble Dbar;
	if(this->isdamaged){
		Input* D_input = element->GetInput(DamageDbarEnum); _assert_(D_input);
		D_input->GetInputValue(&Dbar,gauss);
	}
	else{
		_error_("Cannot get DamageD for non damaged ice");
	}
	return Dbar;
}
/*}}}*/
IssmDouble Matice::GetE(Gauss* gauss){/*{{{*/

	_assert_(this->isenhanced);
	/*Output*/
	IssmDouble E;
	Input* E_input = element->GetInput(MaterialsRheologyEEnum); _assert_(E_input);
	E_input->GetInputValue(&E,gauss);
	return E;
}
/*}}}*/
IssmDouble Matice::GetEbar(Gauss* gauss){/*{{{*/

	_assert_(this->isenhanced);
	/*Output*/
	IssmDouble Ebar;
	Input* E_input = element->GetInput(MaterialsRheologyEbarEnum); _assert_(E_input);
	E_input->GetInputValue(&Ebar,gauss);
	return Ebar;
}
/*}}}*/
IssmDouble Matice::GetN(){/*{{{*/

	/*Output*/
	IssmDouble n;
	Input* n_input = element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
	n_input->GetInputAverage(&n);
	return n;
}
/*}}}*/
bool       Matice::IsDamage(){/*{{{*/

	return this->isdamaged;
}
/*}}}*/
bool       Matice::IsEnhanced(){/*{{{*/

	return this->isenhanced;
}
/*}}}*/
void  Matice::GetViscosity(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){/*{{{*/
	/*From a string tensor and a material object, return viscosity, using Glen's flow law.
								(1-D) B
	  viscosity= -------------------------
						  2 E^[1/n] eps_eff ^[(n-1)/n]

	  where viscosity is the viscosity, B the flow law parameter , eps_eff is the effective strain rate,
	  n the flow law exponent, and E is the enhancement factor.

	  If eps_eff = 0 , it means this is the first time SystemMatrices is being run, and we 
	  return 10^14, initial viscosity.
	  */

	/*output: */
	IssmDouble viscosity;

	/*Intermediary: */
	IssmDouble B,D=0.,E=1.,n;

	/*Get B and n*/
	B=GetB(gauss); _assert_(B>0.);
	n=GetN(); _assert_(n>0.);
	if(this->isdamaged){
		D=GetD(gauss);
		_assert_(D>=0. && D<1.);
	}
	if(this->isenhanced){
		E=GetE(gauss);
		_assert_(E>0.);
	}

	if (n==1.){
		/*Linear Viscous behavior (Newtonian fluid) viscosity=B/2E: */
		viscosity=(1.-D)*B/(2.*E);
	}
	else{

		/*if no strain rate, return maximum viscosity*/
		if(eps_eff==0.){
			viscosity = 1.e+14/2.;
			//viscosity = B;
			//viscosity=2.5*pow(10.,17);
		}

		else{
			viscosity=(1.-D)*B/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
		}
	}

	/*Checks in debugging mode*/
	if(viscosity<=0) _error_("Negative viscosity");

	/*Return: */
	*pviscosity=viscosity;
}
/*}}}*/
void  Matice::GetViscosityBar(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){/*{{{*/
	/*From a string tensor and a material object, return viscosity, using Glen's flow law.
								(1-D) B
	  viscosity= -------------------------
						  2 E^[1/n] eps_eff ^[(n-1)/n]

	  where B the flow law parameter, eps_eff is the effective strain rate, n the flow law exponent,
	  and E is the enhancement factor.

	  If eps_eff = 0 , it means this is the first time SystemMatrices is being run, and we 
	  return 10^14, initial viscosity.
	  */

	/*output: */
	IssmDouble viscosity;

	/*Intermediary: */
	IssmDouble B,D=0.,E=1.,n;

	/*Get B and n*/
	B=GetBbar(gauss); _assert_(B>0.);
	n=GetN();    _assert_(n>0.);
	if(this->isdamaged){
		D=GetDbar(gauss);
		_assert_(D>=0. && D<1.);
	}
	if(this->isenhanced){
		E=GetEbar(gauss);
		_assert_(E>0.);
	}

	if (n==1.){
		/*Linear Viscous behavior (Newtonian fluid) viscosity=B/2E: */
		viscosity=(1.-D)*B/(2.*E);
	}
	else{
		/*if strain rate is 0., it is probably our initial guess, use an average
		 * viscosity instead of a large one*/
		if(eps_eff==0.) viscosity = 1.e+14/2.;
		else{
			/*if no strain rate, return maximum viscosity*/
			//if(eps_eff<1.e-6) eps_eff = 1e-6;
			viscosity=(1.-D)*B/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
		}
	}

	/*Checks in debugging mode*/
	if(viscosity<=0) _error_("Negative viscosity");

	/*Return: */
	*pviscosity=viscosity;
}
/*}}}*/
void  Matice::GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	/*Return viscosity accounting for steady state power law creep [Thomas and SSA, 1982]: 
	 *
	 *  										                (1-D)
	 * viscosity= -------------------------------------------------------------------
	 *  				  2[ exx^2+eyy^2+exx*eyy+exy^2+exz^2+eyz^2 ]^[(n-1)/2n]
	 *
	 * If epsilon is NULL, it means this is the first time Gradjb is being run, and we 
	 * return mu20, initial viscosity.
	 */

	/*output: */
	IssmDouble viscosity_complement;

	/*input strain rate: */
	IssmDouble exx,eyy,exy;

	/*Intermediary value A and exponent e: */
	IssmDouble A,e;
	IssmDouble D=0.,n;

	/*Get D and n*/
	if(this->isdamaged){
		D=GetDbar(gauss); /* GetD()? */
		_assert_(D>=0. && D<1.);
	}
	n=GetN();

	if(epsilon){
		exx=*(epsilon+0);
		eyy=*(epsilon+1);
		exy=*(epsilon+2);

		/*Build viscosity: mu2=(1-D)/(2*A^e) */
		A=pow(exx,2)+pow(eyy,2)+pow(exy,2)+exx*eyy;
		if(A==0){
			/*Maximum viscosity_complement for 0 shear areas: */
			viscosity_complement=2.25*pow(10.,17);
		}
		else{
			e=(n-1)/(2*n);

			viscosity_complement=(1-D)/(2*pow(A,e));
		}
	}
	else{
		viscosity_complement=4.5*pow(10.,17);
	}

	/*Checks in debugging mode*/
	_assert_(D>=0 && D<1);
	_assert_(n>0);
	_assert_(viscosity_complement>0);

	/*Return: */
	*pviscosity_complement=viscosity_complement;
}
/*}}}*/
void  Matice::GetViscosityDComplement(IssmDouble* pviscosity_complement, IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	/*Return viscosity derivative for control method d(mu)/dD: 
	 *
	 *  										               B 
	 * dviscosity= - -------------------------------------------------------------------
	 *  				  2[ exx^2+eyy^2+exx*eyy+exy^2+exz^2+eyz^2 ]^[(n-1)/2n]
	 *
	 * If epsilon is NULL, it means this is the first time Gradjb is being run, and we 
	 * return mu20, initial viscosity.
	 */

	/*output: */
	IssmDouble viscosity_complement;

	/*input strain rate: */
	IssmDouble exx,eyy,exy;

	/*Intermediary value A and exponent e: */
	IssmDouble A,e;
	IssmDouble B,n;

	/*Get B and n*/
	B=GetBbar(gauss);
	n=GetN();

	if(epsilon){
		exx=*(epsilon+0);
		eyy=*(epsilon+1);
		exy=*(epsilon+2);

		/*Build viscosity: mu2=B/(2*A^e) */
		A=pow(exx,2)+pow(eyy,2)+pow(exy,2)+exx*eyy;
		if(A==0){
			/*Maximum viscosity_complement for 0 shear areas: */
			viscosity_complement=- 2.25*pow(10.,17);
		}
		else{
			e=(n-1)/(2*n);

			viscosity_complement=- B/(2*pow(A,e));
		}
	}
	else{
		viscosity_complement=- 4.5*pow(10.,17);
	}

	/*Checks in debugging mode*/
	_assert_(B>0);
	_assert_(n>0);
	_assert_(viscosity_complement<0);

	/*Return: */
	*pviscosity_complement=viscosity_complement;
}
/*}}}*/
void  Matice::GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* epsilon,Gauss* gauss){/*{{{*/

	/*output: */
	IssmDouble mu_prime;
	IssmDouble mu,n,eff2;

	/*input strain rate: */
	IssmDouble exx,eyy,exy,exz,eyz;

	if((epsilon[0]==0) && (epsilon[1]==0) && (epsilon[2]==0) && 
				(epsilon[3]==0) && (epsilon[4]==0)){
		mu_prime=0.5*pow(10.,14);
	}
	else{

		/*Retrieve strain rate components: */
		exx=epsilon[0];
		eyy=epsilon[1];
		exy=epsilon[2];
		exz=epsilon[3];
		eyz=epsilon[4];
		eff2 = exx*exx + eyy*eyy + exx*eyy + exy*exy + exz*exz + eyz*eyz;

		GetViscosity(&mu,sqrt(eff2),gauss);
		n=GetN();
		mu_prime=(1.-n)/(2.*n) * mu/eff2;
	}

	/*Assign output pointers:*/
	*pmu_prime=mu_prime;
}
/*}}}*/
void  Matice::GetViscosity_B(IssmDouble* pdmudB,IssmDouble eps_eff,Gauss* gauss){/*{{{*/

	/*output: */
	IssmDouble dmudB;

	/*Intermediary: */
	IssmDouble D=0.,E=1.,n;

	/*Get B and n*/
	n=GetN(); _assert_(n>0.);
	if(this->isdamaged){
		D=GetD(gauss);
		_assert_(D>=0. && D<1.);
	}
	if(this->isenhanced){
		E=GetE(gauss);
		_assert_(E>0.);
	}

	if(n==1.){
		/*Linear Viscous behavior (Newtonian fluid) dmudB=B/2E: */
		dmudB=(1.-D)/(2.*E);
	}
	else{
		if(eps_eff==0.) dmudB = 0.;
		else            dmudB = (1.-D)/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
	}

	/*Return: */
	*pdmudB=dmudB;
}
/*}}}*/
void  Matice::GetViscosity_D(IssmDouble* pdmudD,IssmDouble eps_eff,Gauss* gauss){/*{{{*/

	/*output: */
	IssmDouble dmudD;

	/*Intermediary: */
	IssmDouble n,B,E=1.;

	/*Get B and n*/
	n=GetN(); _assert_(n>0.);
	B=GetBbar(gauss);
	_assert_(this->isdamaged);
	if(this->isenhanced){
		E=GetE(gauss);
		_assert_(E>0.);
	}

	if(n==1.){
		/*Linear Viscous behavior (Newtonian fluid) dmudB=B/2E: */
		dmudD=-B/(2.*E);
	}
	else{
		if(eps_eff==0.) dmudD = 0.;
		else            dmudD = -B/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
	}

	/*Return: */
	*pdmudD=dmudD;
}
/*}}}*/
void  Matice::GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* epsilon,Gauss* gauss){/*{{{*/

	/*output: */
	IssmDouble mu_prime;
	IssmDouble mu,n,eff2;

	/*input strain rate: */
	IssmDouble exx,eyy,exy;

	if((epsilon[0]==0) && (epsilon[1]==0) && (epsilon[2]==0)){
		mu_prime=0.5*pow(10.,14);
	}
	else{
		/*Retrive strain rate components: */
		exx=epsilon[0];
		eyy=epsilon[1];
		exy=epsilon[2];
		eff2 = exx*exx + eyy*eyy + exx*eyy + exy*exy ;

		GetViscosityBar(&mu,sqrt(eff2),gauss);
		n=GetN();
		mu_prime=(1.-n)/(2.*n)*mu/eff2;
	}

	/*Assign output pointers:*/
	*pmu_prime=mu_prime;
}
/*}}}*/
void  Matice::ResetHooks(){/*{{{*/

	this->element=NULL;

	/*Get Element type*/
	this->helement->reset();

}
/*}}}*/
void  Matice::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Matice::ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	this->GetViscosityDerivativeEpsSquare(pmu_prime,epsilon,gauss);
}/*}}}*/
void  Matice::ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	this->GetViscosityDerivativeEpsSquare(pmu_prime,epsilon,gauss);
}/*}}}*/
void  Matice::ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	this->GetViscosity2dDerivativeEpsSquare(pmu_prime,epsilon,gauss);
}/*}}}*/

void  Matice::ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/
	/*The effective strain rate is defined in Paterson 3d Ed p 91 eq 9,
	 * and Cuffey p 303 eq 8.18:
	 *
	 *  2 eps_eff^2 = eps_xx^2 + eps_yy^2 + eps_zz^2 + 2(eps_xy^2 + eps_xz^2 + eps_yz^2)
	 *
	 *  or
	 *
	 *  eps_eff = 1/sqrt(2) sqrt( \sum_ij eps_ij^2 )
	 *
	 *          = 1/sqrt(2) ||eps||_F
	 *
	 *  where ||.||_F is the Frobenius norm */

	/*Intermediaries*/
	IssmDouble viscosity;
	IssmDouble epsilon3d[6]; /* epsilon=[exx,eyy,ezz,exy,exz,eyz];*/
	IssmDouble epsilon2d[3]; /* epsilon=[exx,eyy,exy];            */
	IssmDouble eps_eff;
	IssmDouble eps0=1.e-27;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		element->StrainRateFS(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[3]*epsilon3d[3] +  epsilon3d[4]*epsilon3d[4] + epsilon3d[5]*epsilon3d[5] + epsilon3d[0]*epsilon3d[1]+eps0*eps0);
	}
	else{
		/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
		element->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + 2.*epsilon2d[2]*epsilon2d[2]);
	}

	/*Get viscosity*/
	this->GetViscosity(&viscosity,eps_eff,gauss);

	/*Assign output pointer*/
	*pviscosity=viscosity;
}
/*}}}*/
void  Matice::ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble viscosity;
	IssmDouble epsilon3d[5];/* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble epsilon2d[2];/* epsilon=[exx,exy];            */
	IssmDouble eps_eff;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		element->StrainRateHO(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[2]*epsilon3d[2] +  epsilon3d[3]*epsilon3d[3] + epsilon3d[4]*epsilon3d[4] + epsilon3d[0]*epsilon3d[1]);
	}
	else{
		/* eps_eff^2 = 1/2 (2*exx^2 + 2*exy^2 ) (since eps_zz = - eps_xx)*/
		element->StrainRateHO2dvertical(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(2*epsilon2d[0]*epsilon2d[0] + 2*epsilon2d[1]*epsilon2d[1]);
	}

	/*Get viscosity*/
	this->GetViscosity(&viscosity,eps_eff,gauss);
	_assert_(!xIsNan<IssmDouble>(viscosity));

	/*Assign output pointer*/
	*pviscosity=viscosity;
}/*}}}*/
void  Matice::ViscosityMOLHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vxbase_input,Input* vybase_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble epsilon[5];	/* epsilon=[exx,eyy,exy,exz,eyz]; */
	IssmDouble epsilon_eff;
	IssmDouble zeta,H,n;
	IssmDouble f[4],F[4];
	IssmDouble mubar[4];
	IssmDouble mu;
	int order=5; 

	for(int i=0;i<4;++i) mubar[i]=0;

	GaussSeg* gauss_seg=new GaussSeg(order);
	//IssmDouble eps_eff_averaged=0;
	while(gauss_seg->next()){

		/*Compute zeta for gauss_seg point (0=surface, 1=base)*/
		zeta=0.5*(gauss_seg->coord1+1);	

		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy (for a given zeta)*/
		element->StrainRateMOLHO(&epsilon[0],xyz_list,gauss,
						vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input,zeta);
		epsilon_eff=sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[2]*epsilon[2] 
						  +  epsilon[3]*epsilon[3] + epsilon[4]*epsilon[4] + epsilon[0]*epsilon[1]);

		/*Get viscosity at zeta */
//?? need to use Bar for the current inversion
//    this->GetViscosity(&mu,epsilon_eff,gauss);
		this->GetViscosityBar(&mu, epsilon_eff,gauss);
		thickness_input->GetInputValue(&H, gauss);
      n_input->GetInputValue(&n,gauss);

		/*Compute fi and Fi at zeta*/
		f[0]=1;
		f[1]=(1-pow(zeta,n+1));
		f[2]=(1-pow(zeta,n+1))*(1-pow(zeta,n+1));
		f[3]=((n+1)/H)*pow(zeta,n) * ((n+1)/H)*pow(zeta,n);

		F[0]=H;
		F[1]=H*(n+1)/(n+2);
		F[2]=2*H*(n+1)*(n+1)/( (2*n+3)*(n+2) );
		F[3]=(n+1)*(n+1)/( H*(2*n+1) );

		/*Sum the viscosity*/
		mubar[0]+=(H/(2*F[0]))*gauss_seg->weight*mu*f[0];
		mubar[1]+=(H/(2*F[1]))*gauss_seg->weight*mu*f[1];
		mubar[2]+=(H/(2*F[2]))*gauss_seg->weight*mu*f[2];
		mubar[3]+=(H/(2*F[3]))*gauss_seg->weight*mu*f[3];

	}//while

	/*Assign output pointer*/
	pviscosity[0]=mubar[0];
	pviscosity[1]=mubar[1];
	pviscosity[2]=mubar[2];
	pviscosity[3]=mubar[3];

	/*Clean up*/
	delete gauss_seg;
}/*}}}*/
void  Matice::ViscosityMOLHOAdjoint(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vxbase_input,Input* vybase_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input){/*{{{*/

	/* To compute the additional 5 terms in the viscosity appear in the adjoint equation*/
	/*Intermediaries*/
	IssmDouble epsilon[5];	/* epsilon=[exx,eyy,exy,exz,eyz]; */
	IssmDouble epsilon_eff;
	IssmDouble zeta,H,n;
	IssmDouble f[9],F[9];
	IssmDouble mubar[9];
	IssmDouble mu;
	int order=10; 

	for(int i=0;i<9;++i) mubar[i]=0;

	GaussSeg* gauss_seg=new GaussSeg(order);
	//IssmDouble eps_eff_averaged=0;
	while(gauss_seg->next()){
		
		/*Compute zeta for gauss_seg point (0=surface, 1=base)*/
		zeta=0.5*(gauss_seg->coord1+1);	

		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy (for a given zeta)*/
		element->StrainRateMOLHO(&epsilon[0],xyz_list,gauss,
						vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input,zeta);
		epsilon_eff=sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[2]*epsilon[2] 
						  +  epsilon[3]*epsilon[3] + epsilon[4]*epsilon[4] + epsilon[0]*epsilon[1]);

		this->GetViscosity(&mu,epsilon_eff,gauss);
		/*the adjoint viscosity with zeta dependent term*/
		mu = mu /epsilon_eff/epsilon_eff;

		thickness_input->GetInputValue(&H, gauss);
      n_input->GetInputValue(&n,gauss);

		/*Compute fi and Fi at zeta*/
		f[0]=1;
		f[1]=(1-pow(zeta,n+1));
		f[2]=(1-pow(zeta,n+1))*(1-pow(zeta,n+1));
		f[3]=pow(zeta,2*n); // NOTE: this is different from the forward formulation
		f[4]=(1-pow(zeta,n+1))*(1-pow(zeta,n+1))*(1-pow(zeta,n+1));
		f[5]=(1-pow(zeta,n+1))*(1-pow(zeta,n+1))*(1-pow(zeta,n+1))*(1-pow(zeta,n+1));
		f[6]=(1-pow(zeta,n+1))*pow(zeta,2*n);
		f[7]=(1-pow(zeta,n+1))*(1-pow(zeta,n+1))*pow(zeta,2*n);
		f[8]=pow(zeta,4*n);

	
		F[0]=H;
		F[1]=H*(n+1)/(n+2);
		F[2]=2*H*(n+1)*(n+1)/( (2*n+3)*(n+2) );
		F[3]=H/(2*n+1);
		F[4]=6*H*(n+1)*(n+1)*(n+1)/( (n+2)*(2*n+3)*(3*n+4) );
		F[5]=24*H*(n+1)*(n+1)*(n+1)*(n+1)/( (n+2)*(2*n+3)*(3*n+4)*(4*n+5) );
		F[6]=H*(n+1)/( (2*n+1)*(3*n+2) );
		F[7]=2*H*(n+1)*(n+1)/( (2*n+1)*(3*n+2)*(4*n+3) );
		F[8]=H/(4*n+1);

		/*Sum the viscosity*/
		for(int i=0;i<9;i++) {
			mubar[i]+=(H/2)*gauss_seg->weight*mu*f[i];
		}
	}//while

	/*Assign output pointer*/
	for(int i=0;i<9;i++) pviscosity[i]=mubar[i];

	/*Clean up*/
	delete gauss_seg;
}/*}}}*/
void  Matice::ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surface_input){/*{{{*/
	/*Compute the L1L2 viscosity
	 *
	 *      1
	 * mu = - A^-1 (sigma'_e)^(1-n)
	 *      2
	 *
	 * sigma'_e^2 = |sigma'_//|^2 + |sigma'_perp|^2 (see Perego 2012 eq. 17,18)
	 *
	 * L1L2 assumptions:
	 *
	 * (1) |eps_b|_// = A (|sigma'_//|^2 + |sigma'_perp|^2)^((n-1)/2) |sigma'_//|
	 * (2) |sigma'_perp|^2 = |rho g (s-z) grad(s)|^2
	 *
	 * Assuming that n = 3, we have a polynom of degree 3 to solve (the only unkown is X=|sigma'_//|)
	 *
	 * A X^3 + A |rho g (s-z) grad(s)|^2 X - |eps_b|_// = 0     */

	IssmDouble z,s,viscosity,p,q,delta;
	IssmDouble tau_perp,tau_par,eps_b,A;
	IssmDouble epsilon[5];   /*exx eyy exy exz eyz*/
	IssmDouble slope[3];

	/*Check that both inputs have been found*/
	if (!vx_input || !vy_input || !surface_input) _error_("Input missing");

	/*Get tau_perp*/
	surface_input->GetInputValue(&s,gauss);
	surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
	z=this->element->GetZcoord(xyz_list,gauss);
	tau_perp = element->FindParam(MaterialsRhoIceEnum) * element->FindParam(ConstantsGEnum) * fabs(s-z)*sqrt(slope[0]*slope[0]+slope[1]*slope[1]);

	/* Get eps_b*/
	element->StrainRateHO(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
	eps_b = sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[0]*epsilon[1] + epsilon[2]*epsilon[2]);
	if(eps_b==0.){
		*pviscosity = 2.5e+17;
		return;
	}

	/*Get A*/
	_assert_(this->GetN()==3.0);
	A=this->GetA(gauss);

	/*Solve for tau_perp (http://fr.wikipedia.org/wiki/Méthode_de_Cardan)*/
	p     = tau_perp *tau_perp;
	q     = - eps_b/A;
	delta = q *q + p*p*p*4./27.;
	_assert_(delta>0);
	tau_par = pow(0.5*(-q+sqrt(delta)),1./3.) - pow(0.5*(q+sqrt(delta)),1./3.);

	/*Viscosity*/
	viscosity = 1./(2.*A*(tau_par*tau_par + tau_perp*tau_perp));
	_assert_(!xIsNan(viscosity));
	_assert_(viscosity > 0.);

	/*Assign output pointer*/
	*pviscosity = viscosity;
}/*}}}*/
void  Matice::ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble viscosity;
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble epsilon1d;   /* epsilon=[exx];    */
	IssmDouble eps_eff;

	if(dim==2){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exx*eyy*/
		element->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + epsilon2d[2]*epsilon2d[2] + epsilon2d[0]*epsilon2d[1]);
	}
	else{
		/* eps_eff^2 = exx^2*/
		element->StrainRateSSA1d(&epsilon1d,xyz_list,gauss,vx_input);
		eps_eff = fabs(epsilon1d);
	}

	/*Get viscosity*/
	this->GetViscosityBar(&viscosity,eps_eff,gauss);

	/*Assign output pointer*/
	*pviscosity=viscosity;
}/*}}}*/
