/*!\file Matestar.c
 * \brief: implementation of the Matestar object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Matestar.h"
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

/*Matestar constructors and destructor*/
Matestar::Matestar(){/*{{{*/
	this->helement=NULL;
	this->element=NULL;
	return;
}
/*}}}*/
Matestar::Matestar(int matestar_mid,int index, IoModel* iomodel){/*{{{*/

	/*Intermediaries:*/
	int    matestar_eid;

	/*Initialize id*/
	this->mid=matestar_mid;

	/*Hooks: */
	matestar_eid=index+1;
	this->helement=new Hook(&matestar_eid,1);
	this->element=NULL;

	return;
}
/*}}}*/
Matestar::~Matestar(){/*{{{*/
	delete helement;
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object*   Matestar::copy() {/*{{{*/

	/*Output*/
	Matestar* matestar=NULL;

	/*Initialize output*/
	matestar=new Matestar();

	/*copy fields: */
	matestar->mid=this->mid;
	matestar->helement=(Hook*)this->helement->copy();
	matestar->element =(Element*)this->helement->delivers();

	return matestar;
}
/*}}}*/
Material* Matestar::copy2(Element* element_in) {/*{{{*/

	/*Output*/
	Matestar* matestar=NULL;

	/*Initialize output*/
	matestar=new Matestar();

	/*copy fields: */
	matestar->mid=this->mid;
	matestar->helement=(Hook*)this->helement->copy();
	matestar->element =element_in;

	return matestar;
}
/*}}}*/
void      Matestar::DeepEcho(void){/*{{{*/

	_printf_("Matestar:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   element:\n");
	helement->Echo();
}		
/*}}}*/
void      Matestar::Echo(void){/*{{{*/

	_printf_("Matestar:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   element:\n");
	helement->Echo();
}
/*}}}*/
int       Matestar::Id(void){ return mid; }/*{{{*/
/*}}}*/
void      Matestar::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD)helement=new Hook(); 

	int object_enum = MatestarEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->mid);
	this->helement->Marshall(marshallhandle);
	this->element=(Element*)this->helement->delivers();

}
/*}}}*/
int       Matestar::ObjectEnum(void){/*{{{*/

	return MatestarEnum;

}
/*}}}*/

/*Matestar management*/
void  Matestar::Configure(Elements* elementsin){/*{{{*/

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	helement->configure((DataSet*)elementsin);
	this->element  = (Element*)helement->delivers();
}
/*}}}*/
IssmDouble Matestar::GetA(Gauss* gauss){/*{{{*/
	/*
	 * A = 1/B^n
	 */

	IssmDouble B=this->GetB(gauss);
	IssmDouble n=this->GetN();

	return pow(B,-n);
}
/*}}}*/
IssmDouble Matestar::GetAbar(Gauss* gauss){/*{{{*/
	/*
	 * A = 1/B^n
	 */

	IssmDouble B=this->GetBbar(gauss);
	IssmDouble n=this->GetN();

	return pow(B,-n);
}
/*}}}*/
IssmDouble Matestar::GetB(Gauss* gauss){/*{{{*/

	/*Output*/
	IssmDouble B;

	Input* B_input = element->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	B_input->GetInputValue(&B,gauss);
	return B;
}
/*}}}*/
IssmDouble Matestar::GetBbar(Gauss* gauss){/*{{{*/

	/*Output*/
	IssmDouble Bbar;

	Input* B_input = element->GetInput(MaterialsRheologyBbarEnum); _assert_(B_input);
	B_input->GetInputValue(&Bbar,gauss);
	return Bbar;
}
/*}}}*/
IssmDouble Matestar::GetD(Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}
/*}}}*/
IssmDouble Matestar::GetDbar(Gauss* gauss){/*{{{*/

	_error_("not implemented yet");
}
/*}}}*/
IssmDouble Matestar::GetEc(Gauss* gauss){/*{{{*/

	/*Output*/
	IssmDouble Ec;

	Input* Ec_input = element->GetInput(MaterialsRheologyEcEnum); _assert_(Ec_input);
	Ec_input->GetInputValue(&Ec,gauss);
	return Ec;
}
/*}}}*/
IssmDouble Matestar::GetEcbar(Gauss* gauss){/*{{{*/

	/*Output*/
	IssmDouble Ecbar;

	Input* Ecbar_input = element->GetInput(MaterialsRheologyEcbarEnum); _assert_(Ecbar_input);
	Ecbar_input->GetInputValue(&Ecbar,gauss);
	return Ecbar;
}
/*}}}*/
IssmDouble Matestar::GetEs(Gauss* gauss){/*{{{*/

	/*Output*/
	IssmDouble Es;

	Input* Es_input = element->GetInput(MaterialsRheologyEsEnum); _assert_(Es_input);
	Es_input->GetInputValue(&Es,gauss);
	return Es;
}
/*}}}*/
IssmDouble Matestar::GetEsbar(Gauss* gauss){/*{{{*/

	/*Output*/
	IssmDouble Esbar;

	Input* Esbar_input = element->GetInput(MaterialsRheologyEsbarEnum); _assert_(Esbar_input);
	Esbar_input->GetInputValue(&Esbar,gauss);
	return Esbar;
}
/*}}}*/
IssmDouble Matestar::GetN(){/*{{{*/

	/*Output*/
	IssmDouble n=3.0;
	return n;
}
/*}}}*/
void  Matestar::GetViscosity(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}
/*}}}*/
void  Matestar::GetViscosityBar(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}
/*}}}*/
void  Matestar::GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}
/*}}}*/
void  Matestar::GetViscosityDComplement(IssmDouble* pviscosity_complement, IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}
/*}}}*/
void  Matestar::GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}
/*}}}*/
IssmDouble Matestar::GetViscosityGeneral(IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz,IssmDouble eps_eff,bool isdepthaveraged,Gauss* gauss){/*{{{*/

	/*output: */
	IssmDouble viscosity;

	/*Intermediaries*/
	IssmDouble epsprime_norm;
	IssmDouble lambdas;
	IssmDouble vmag,dvmag[3];
	IssmDouble B,Ec,Es,E,n;

	/*Calculate velocity magnitude and its derivative*/
	vmag = sqrt(vx*vx+vy*vy+vz*vz);
	if(vmag<1e-12){
		dvmag[0]=0;
		dvmag[1]=0;
		dvmag[2]=0;
	}
	else{
		dvmag[0]=1./(2*sqrt(vmag))*(2*vx*dvx[0]+2*vy*dvy[0]+2*vz*dvz[0]);
		dvmag[1]=1./(2*sqrt(vmag))*(2*vx*dvx[1]+2*vy*dvy[1]+2*vz*dvz[1]);
		dvmag[2]=1./(2*sqrt(vmag))*(2*vx*dvx[2]+2*vy*dvy[2]+2*vz*dvz[2]);
	}

	EstarStrainrateQuantities(&epsprime_norm,vx,vy,vz,vmag,dvx,dvy,dvz,&dvmag[0]);
	lambdas=EstarLambdaS(eps_eff,epsprime_norm);

	/*Get B and enhancement*/
	n=GetN(); _assert_(n>0.);
	if (isdepthaveraged==0.){
		B=GetB(gauss);   _assert_(B>0.);
		Ec=GetEc(gauss); _assert_(Ec>=0.); Es=GetEs(gauss); _assert_(Es>=0.);
	}
	else{
		B=GetBbar(gauss);   _assert_(B>0.);
		Ec=GetEcbar(gauss); _assert_(Ec>=0.);
		Es=GetEsbar(gauss); _assert_(Es>=0.);
	}

	/*Get total enhancement factor E(lambdas)*/
	E = Ec + (Es-Ec)*lambdas*lambdas; _assert_(E>0.);

	/*Compute viscosity*/
	/*if no strain rate, return maximum viscosity*/
	if(eps_eff==0.){
		viscosity = 1.e+14/2.;
		//viscosity = B;
		//viscosity=2.5*pow(10.,17);
	}
	else{
		viscosity = B/(2.*pow(E,1./n)*pow(eps_eff,2./n));
	}

   /*Checks in debugging mode*/
	if(viscosity<=0) _error_("Negative viscosity");

	/*Assign output pointer*/
	return viscosity;
}
/*}}}*/
IssmDouble Matestar::GetViscosity_BGeneral(IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz,IssmDouble eps_eff,bool isdepthaveraged,Gauss* gauss){/*{{{*/

	/*Intermediaries*/
	IssmDouble dmudB;
	IssmDouble epsprime_norm;
	IssmDouble lambdas;
	IssmDouble vmag,dvmag[3];
	IssmDouble Ec,Es,E;

	/*Calculate velocity magnitude and its derivative*/
	vmag = sqrt(vx*vx+vy*vy+vz*vz);
	if(vmag<1e-12){
		dvmag[0]=0;
		dvmag[1]=0;
		dvmag[2]=0;
	}
	else{
		dvmag[0]=1./(2*sqrt(vmag))*(2*vx*dvx[0]+2*vy*dvy[0]+2*vz*dvz[0]);
		dvmag[1]=1./(2*sqrt(vmag))*(2*vx*dvx[1]+2*vy*dvy[1]+2*vz*dvz[1]);
		dvmag[2]=1./(2*sqrt(vmag))*(2*vx*dvx[2]+2*vy*dvy[2]+2*vz*dvz[2]);
	}

	EstarStrainrateQuantities(&epsprime_norm,vx,vy,vz,vmag,dvx,dvy,dvz,&dvmag[0]);
	lambdas=EstarLambdaS(eps_eff,epsprime_norm);

	/*Get enhancement*/
	if (isdepthaveraged==0.){
		Ec=GetEc(gauss); _assert_(Ec>=0.);
		Es=GetEs(gauss); _assert_(Es>=0.);
	}
	else{
		Ec=GetEcbar(gauss); _assert_(Ec>=0.);
		Es=GetEsbar(gauss); _assert_(Es>=0.);
	}

	/*Get total enhancement factor E(lambdas)*/
	E = Ec + (Es-Ec)*lambdas*lambdas; _assert_(E>0.);

	/*Compute dmudB*/
	if(eps_eff==0.) dmudB = 0.;
	else            dmudB = 1./(2.*pow(E,1./3.)*pow(eps_eff,2./3.));

	/*Assign output*/
	return dmudB;

}
/*}}}*/
void  Matestar::GetViscosity_B(IssmDouble* pdmudB,IssmDouble eps_eff,Gauss* gauss){/*{{{*/
	 _error_("not implemented yet");
}
/*}}}*/
void  Matestar::GetViscosity_D(IssmDouble* pdmudD,IssmDouble eps_eff,Gauss* gauss){/*{{{*/
	 _error_("not implemented yet");
}
/*}}}*/
void  Matestar::GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}
/*}}}*/
bool Matestar::IsDamage(){/*{{{*/

	_error_("not implemented yet");
}
/*}}}*/
void  Matestar::ResetHooks(){/*{{{*/

	this->element=NULL;

	/*Get Element type*/
	this->helement->reset();

}
/*}}}*/
void  Matestar::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Matestar::ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	this->GetViscosityDerivativeEpsSquare(pmu_prime,epsilon,gauss);
}/*}}}*/
void  Matestar::ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void  Matestar::ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/

void  Matestar::ViscosityBFS(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input,IssmDouble eps_eff){/*{{{*/

	/*Intermediaries*/
	IssmDouble vx,vy,vz;
	IssmDouble dvx[3],dvy[3],dvz[3];
	bool isdepthaveraged=0.;

	/*Get velocity derivatives in all directions*/
	_assert_(dim>1);
	_assert_(vx_input);
	vx_input->GetInputValue(&vx,gauss);
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	_assert_(vy_input);
	vy_input->GetInputValue(&vy,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	if(dim==3){
		_assert_(vz_input);
		vz_input->GetInputValue(&vz,gauss);
		vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
	}
	else{
		vz = 0.;
		dvz[0] = 0.; dvz[1] = 0.; dvz[2] = 0.;
	}

	/*Compute dmudB*/
	*pdmudB=GetViscosity_BGeneral(vx,vy,vz,&dvx[0],&dvy[0],&dvz[0],eps_eff,isdepthaveraged,gauss);
}
/*}}}*/
void  Matestar::ViscosityBHO(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble eps_eff){/*{{{*/

	/*Intermediaries*/
	IssmDouble vx,vy,vz;
	IssmDouble dvx[3],dvy[3],dvz[3];
	bool isdepthaveraged=0.;

	/*Get velocity derivatives in all directions*/
	_assert_(dim==2 || dim==3);
	_assert_(vx_input);
	vx_input->GetInputValue(&vx,gauss);
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	if(dim==3){
		_assert_(vy_input);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	}
	else{
		dvx[2] = 0.;
		vy = 0.;
		dvy[0] = 0.; dvy[1] = 0.; dvy[2] = 0.;
	}
	vz = 0.;
	dvz[0] = 0.; dvz[1] = 0.; dvz[2] = -dvx[0]-dvy[1];

	/*Compute viscosity*/
	*pdmudB=GetViscosity_BGeneral(vx,vy,vz,&dvx[0],&dvy[0],&dvz[0],eps_eff,isdepthaveraged,gauss);
}/*}}}*/
void  Matestar::ViscosityBSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble eps_eff){/*{{{*/
	/*Intermediaries*/
	IssmDouble vx,vy,vz;
	IssmDouble dvx[3],dvy[3],dvz[3];
	bool isdepthaveraged=1.;

	/*Get velocity derivatives in all directions*/
	_assert_(dim==1 || dim==2);
	_assert_(vx_input);
	vx_input->GetInputValue(&vx,gauss);
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	if(dim==2){
		_assert_(vy_input);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	}
	else{
		dvx[1] = 0.;
		dvx[2] = 0.;
		vy = 0.;
		dvy[0] = 0.; dvy[1] = 0.; dvy[2] = 0.;
	}
	dvx[2] = 0.;
	dvy[2] = 0.;
	vz = 0.;
	dvz[0] = 0.; dvz[1] = 0.; dvz[2] = -dvx[0]-dvy[1];

	/*Compute viscosity*/
	*pdmudB=GetViscosity_BGeneral(vx,vy,vz,&dvx[0],&dvy[0],&dvz[0],eps_eff,isdepthaveraged,gauss);
}/*}}}*/
void  Matestar::ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble vx,vy,vz;
	IssmDouble dvx[3],dvy[3],dvz[3];
	IssmDouble epsilon3d[6]; /* epsilon=[exx,eyy,ezz,exy,exz,eyz];*/
	IssmDouble epsilon2d[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble eps_eff,eps0=1.e-27;
	bool isdepthaveraged=0.;

	/*Get velocity derivatives in all directions*/
	_assert_(dim>1);
	_assert_(vx_input);
	vx_input->GetInputValue(&vx,gauss);
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	_assert_(vy_input);
	vy_input->GetInputValue(&vy,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	if(dim==3){
		_assert_(vz_input);
		vz_input->GetInputValue(&vz,gauss);
		vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
	}
	else{
		vz = 0.;
		dvz[0] = 0.; dvz[1] = 0.; dvz[2] = 0.;
	}

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

	/*Compute viscosity*/
	*pviscosity=GetViscosityGeneral(vx,vy,vz,&dvx[0],&dvy[0],&dvz[0],eps_eff,isdepthaveraged,gauss);
}
/*}}}*/
void  Matestar::ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble vx,vy,vz;
	IssmDouble dvx[3],dvy[3],dvz[3];
	IssmDouble epsilon3d[5]; /* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble epsilon2d[5]; /* epsilon=[exx,exy];*/
	IssmDouble eps_eff;
	bool isdepthaveraged=0.;

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

	/*Get velocity derivatives in all directions*/
	_assert_(dim==2 || dim==3);
	_assert_(vx_input);
	vx_input->GetInputValue(&vx,gauss);
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	if(dim==3){
		_assert_(vy_input);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	}
	else{
		dvx[2] = 0.;
		vy = 0.;
		dvy[0] = 0.; dvy[1] = 0.; dvy[2] = 0.;
	}
	vz = 0.;
	dvz[0] = 0.; dvz[1] = 0.; dvz[2] = -dvx[0]-dvy[1];

	/*Compute viscosity*/
	*pviscosity=GetViscosityGeneral(vx,vy,vz,&dvx[0],&dvy[0],&dvz[0],eps_eff,isdepthaveraged,gauss);
}/*}}}*/
void  Matestar::ViscosityMOLHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void  Matestar::ViscosityMOLHOAdjoint(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void  Matestar::ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surface_input){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void  Matestar::ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble vx,vy,vz;
	IssmDouble dvx[3],dvy[3],dvz[3];
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy]; */
	IssmDouble epsilon1d;   /* epsilon=[exx];         */
	IssmDouble eps_eff;
	bool isdepthaveraged=1.;

	/*Get velocity derivatives in all directions*/
	_assert_(dim==1 || dim==2);
	_assert_(vx_input);
	vx_input->GetInputValue(&vx,gauss);
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	if(dim==2){
		_assert_(vy_input);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	}
	else{
		dvx[1] = 0.;
		dvx[2] = 0.;
		vy = 0.;
		dvy[0] = 0.; dvy[1] = 0.; dvy[2] = 0.;
	}
	dvx[2] = 0.;
	dvy[2] = 0.;
	vz = 0.;
	dvz[0] = 0.; dvz[1] = 0.; dvz[2] = -dvx[0]-dvy[1];

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

	/*Compute viscosity*/
	*pviscosity=GetViscosityGeneral(vx,vy,vz,&dvx[0],&dvy[0],&dvz[0],eps_eff,isdepthaveraged,gauss);
}/*}}}*/
