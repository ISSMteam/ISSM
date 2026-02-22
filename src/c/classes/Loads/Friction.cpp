/*!\file Friction.c
 * \brief: implementation of the Friction object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
#include "../../modules/InputUpdateFromConstantx/InputUpdateFromConstantx.h"
/*}}}*/

/*Constructors/destructors*/
Friction::Friction(){/*{{{*/
	this->element=NULL;
	this->law        = 0;
	this->linearize  = 0;
	this->apply_dim  = 1.;
	this->domaintype = -1;
	this->vx_input=NULL;
	this->vy_input=NULL;
	this->vz_input=NULL;
	this->alpha2_list=NULL;
	this->alpha2_complement_list=NULL;
}
/*}}}*/
Friction::Friction(Element* element_in){/*{{{*/
	/* Determine the dimension according to the domain type automatically. 
	 * There are exceptions, e.g. HO, which needs the user to specify the dimension used in Friciton.*/

	/*Intermediaries*/
	int linearization_type;

	this->element=element_in;
	this->linearize  = 0;

	/* Load necessary parameters */
	element_in->FindParam(&this->law,FrictionLawEnum);
	element_in->FindParam(&this->domaintype,DomainTypeEnum);

	/* Load VxBase and VyBase for this special case */
	switch(this->domaintype){
		case Domain2DhorizontalEnum: 
			this->apply_dim = 2.;
			this->vx_input = element_in->GetInput(VxBaseEnum);	_assert_(this->vx_input); 
			this->vy_input = element_in->GetInput(VyBaseEnum);	_assert_(this->vy_input);
			this->vz_input = NULL;
			break;
      case Domain2DverticalEnum:
			this->apply_dim = 2.;
			this->vx_input = element_in->GetInput(VxEnum);	_assert_(this->vx_input);
			this->vy_input = element_in->GetInput(VyEnum);	_assert_(this->vy_input);
			this->vz_input = NULL;
			break;
      case Domain3DEnum:           
			this->apply_dim = 3.;
			this->vx_input = element_in->GetInput(VxEnum);	_assert_(this->vx_input);
			this->vy_input = element_in->GetInput(VyEnum);	_assert_(this->vy_input);
			this->vz_input = element_in->GetInput(VzEnum);	_assert_(this->vz_input);
			break;
      default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	if(this->law==1 || this->law==2){
		element_in->FindParam(&linearization_type,FrictionLinearizeEnum);
		if(linearization_type==0){
			/*Don't do anything*/
		}
		else if(linearization_type==1){
			int numvertices = this->element->GetNumberOfVertices();
			this->alpha2_list            = xNew<IssmDouble>(numvertices);
			this->alpha2_complement_list = xNew<IssmDouble>(numvertices);
			Gauss* gauss=this->element->NewGauss();
			for(int iv=0;iv<numvertices;iv++){
				gauss->GaussVertex(iv);
				this->GetAlpha2(&this->alpha2_list[iv], gauss);
				this->GetAlphaComplement(&this->alpha2_complement_list[iv], gauss);
			}
			this->linearize = linearization_type; /*Change back, we are now all set!*/
			delete gauss;
		}
		else if(linearization_type==2){
			this->alpha2_list            = xNew<IssmDouble>(1);
			this->alpha2_complement_list = xNew<IssmDouble>(1);
			Gauss* gauss=element->NewGauss(1); gauss->GaussPoint(0);
			this->GetAlpha2(&this->alpha2_list[0], gauss);
			this->GetAlphaComplement(&this->alpha2_complement_list[0], gauss);
			this->linearize = linearization_type; /*Change back, we are now all set!*/
			delete gauss;
		}
		else{
			_error_("not supported yet");
		}
	}
}
/*}}}*/
Friction::Friction(Element* element_in,int dim) : Friction(element_in) {/*{{{*/
	this->apply_dim = reCast<IssmPDouble>(dim);
}
/*}}}*/
Friction::Friction(Element* element_in,IssmPDouble dim) : Friction(element_in) {/*{{{*/
	this->apply_dim = dim;
}
/*}}}*/
Friction::~Friction(){/*{{{*/
	if(this->linearize){
		xDelete<IssmDouble>(this->alpha2_list);
		xDelete<IssmDouble>(this->alpha2_complement_list);
	}
}
/*}}}*/

/*methods: */
void Friction::Echo(void){/*{{{*/
	_printf_("Friction:\n");
	_printf_("   Domain type: " << this->domaintype<< "\n");
}
/*}}}*/
void Friction::GetAlphaComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	if(this->linearize==0){
		switch(this->law){
			case 1:
				GetAlphaBuddComplement(palpha_complement,gauss);
				break;
			case 2:
				GetAlphaWeertmanComplement(palpha_complement, gauss);
				break;
			case 3:
				GetAlphaHydroComplement(palpha_complement,gauss);
				break;
			case 4:
				GetAlphaTempComplement(palpha_complement,gauss);
				break;
			case 11:
				GetAlphaSchoofComplement(palpha_complement,gauss);
				break;
			case 13:
				GetAlphaCoulomb2Complement(palpha_complement,gauss);
				break;
			case 14:
				GetAlphaRegCoulombComplement(palpha_complement,gauss);
				break;
			default:
				_error_("not supported");
		}
	}
	else if(this->linearize==1){
		this->element->ValueP1OnGauss(palpha_complement, this->alpha2_complement_list, gauss);
	}
	else if(this->linearize==2){
		*palpha_complement = this->alpha2_complement_list[0];
	}
	else{
		_error_("not supported yet");
	}

	/*Checks*/
	_assert_(!xIsNan<IssmDouble>(*palpha_complement));
	_assert_(!xIsInf<IssmDouble>(*palpha_complement));

}/*}}}*/
void Friction::GetAlphaHydroComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/*diverse: */
	IssmDouble  q_exp;
	IssmDouble  C_param;
	IssmDouble  As;
	IssmDouble  n;
	IssmDouble  alpha;
	IssmDouble  Chi,Gamma;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&q_exp,gauss,FrictionQEnum);
	element->GetInputValue(&C_param,gauss,FrictionCEnum);
	element->GetInputValue(&As,gauss,FrictionAsEnum);
	element->GetInputValue(&n,gauss,MaterialsRheologyNEnum);

	/*Get effective pressure and velocity magnitude*/
	IssmDouble Neff = EffectivePressure(gauss);
	IssmDouble vmag = VelMag(gauss);

	if (q_exp==1){
		alpha=1;
	}
	else{
		alpha=(pow(q_exp-1,q_exp-1))/pow(q_exp,q_exp);
	}
	Chi   = vmag/(pow(C_param,n)*pow(Neff,n)*As);
	Gamma = (Chi/(1.+alpha*pow(Chi,q_exp)));
	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0.) alpha_complement=0.;
	else	if(Neff==0.) alpha_complement=0.;
	else	alpha_complement=-(C_param*Neff/(n*vmag)) *
					pow(Gamma,((1.-n)/n)) *
					(Gamma/As - (alpha*q_exp*pow(Chi,q_exp-1.)* Gamma * Gamma/As));

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}
/*}}}*/
void Friction::GetAlphaTempComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_viscous * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  f,T,pressure,Tpmp,gamma;
	IssmDouble  alpha_complement;

	/*Get viscous part*/
	this->GetAlphaBuddComplement(&alpha_complement,gauss);

	/*Get pressure melting point (Tpmp) for local pressure and get current temperature*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);

	/*Compute scaling parameter*/
	element->parameters->FindParam(&gamma,FrictionGammaEnum);
	alpha_complement = alpha_complement/ (exp((T-Tpmp)/gamma)+1e-3);

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}/*}}}*/
void Friction::GetAlphaBuddComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/* FrictionGetAlpha2 computes alpha2= drag^2 * Neff ^r * vel ^(s-1), with Neff=rho_ice*g*thickness+rho_ice*g*base, r=q/p and s=1/p.
	 * FrictionGetAlphaComplement is used in control methods on drag, and it computes:
	 * alpha_complement= Neff ^r * vel ^(s-1)*/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  drag_p,drag_q;
	IssmDouble  drag_coefficient;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&drag_p,gauss,FrictionPEnum);
	element->GetInputValue(&drag_q,gauss,FrictionQEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);

	/*compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	/*Get effective pressure*/
	IssmDouble Neff = EffectivePressure(gauss);
	IssmDouble vmag = VelMag(gauss);

	if(s==1.){
		/*This is to make AD happy and avoid 0^0*/
		alpha_complement=pow(Neff,r);
	}
	else{
		/*Check to prevent dividing by zero if vmag==0*/
		if(vmag==0. && (s-1.)<0.) alpha_complement=0.;
		else alpha_complement=pow(Neff,r)*pow(vmag,(s-1.));
	}

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}
/*}}}*/
void Friction::GetAlphaSchoofComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/* Compute the complement of Schoof's law for inversion
	 * d alpha2                       
	 * -------- = |u_b|^(m-1) *(1+ (C/(Cmax N))^(1/m)|u_b|)^(-m-1)
	 *  dC                           
	*/
	/*diverse: */
	IssmDouble  m,Cmax;
	IssmDouble  C, coeff;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&m,gauss,FrictionMEnum);
	element->GetInputValue(&Cmax,gauss,FrictionCmaxEnum);
	element->GetInputValue(&coeff, gauss,FrictionCEnum);

	C = coeff*coeff;
	/*Get effective pressure*/
	IssmDouble Neff = EffectivePressure(gauss);
	IssmDouble vmag = VelMag(gauss);

	/*Check to prevent dividing by zero if vmag==0*/
	if((vmag==0.) || (Neff == 0.)) {
		alpha_complement=0.;
	}
	else {
		alpha_complement= pow(vmag, m-1.)*pow((1 + pow(C/(Cmax*Neff),1./m)*vmag), -m-1.);
	}
	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}/*}}}*/
void Friction::GetAlphaWeertmanComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/* Compute the complement of Weertman's law for inversion
	 * alpha2 = C^2 * vel^(1/m-1)
	 * alpha_complement = vel^(1/m-1)
	*/
	/*diverse: */
	IssmDouble  m;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&m,gauss,FrictionMEnum);

	/*Get effective pressure*/
	IssmDouble vmag = VelMag(gauss);

	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0. && (1./m-1.)<0.) alpha_complement=0.;
	else alpha_complement= pow(vmag, 1.0/m-1.);

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}/*}}}*/
void Friction::GetAlphaRegCoulombComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/* Compute the complement of regularised Coulombs law for inversion
	 * d alpha2                       
	 * -------- = |u_b|^(1/m-1) * (|u_b|/u_0 + 1)^(-1/m)
	 *  dC                           
	 */

	/*diverse: */
	IssmDouble  m, u0;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&m,gauss,FrictionMEnum);
	element->parameters->FindParam(&u0,FrictionU0Enum);

	/*Get velocity magnitude*/
	IssmDouble ub = VelMag(gauss);

	/*Check to prevent dividing by zero if vmag==0*/
	if(ub==0.) {
		alpha_complement=0.;
	}
	else {
		/*Compute friction complement*/
		alpha_complement= (pow(ub,1./m-1.)) / pow(ub/u0 + 1.,1./m);	
	}

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}/*}}}*/
void Friction::GetAlphaCoulomb2Complement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/* Compute the complement of Cornford's friction law for inversion
	 * d alpha2                       
	 * ------ = (C*N*v^m)/(C^(2/m)*v + (N/2)^(1/m))^m - (C^(2/m - 1)*C^2*N*v*v^m)/(C^(2/m)*v + (N/2)^(1/m))^(m + 1)
	 *  dC                           
	 */

	/*diverse: */
	IssmDouble  m, C;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&C,gauss,FrictionCEnum);
	element->GetInputValue(&m,gauss,FrictionMEnum);

	/*Get effective pressure and velocity magnitude*/
	IssmDouble N = EffectivePressure(gauss);
	IssmDouble v = VelMag(gauss);

	/*Compute alpha*/
	if(v<1e-10){
		alpha_complement = 0.;
	}
	else{
		alpha_complement= pow(0.5*N,1./m+1)* pow(v,m-1.) * pow(v*pow(C,1./m) +pow(0.5*N,1./m) ,-m-1.);
	}

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement/2.;
}/*}}}*/
void Friction::GetAlpha2(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	if(this->linearize==0){
		switch(this->law){
			case 1:
				GetAlpha2Budd(palpha2,gauss);
				break;
			case 2:
				GetAlpha2Weertman(palpha2,gauss);
				break;
			case 3:
				GetAlpha2Hydro(palpha2,gauss);
				break;
			case 4:
				GetAlpha2Temp(palpha2,gauss);
				break;
			case 5:
				GetAlpha2WaterLayer(palpha2,gauss);
				break;
			case 6:
				GetAlpha2WeertmanTemp(palpha2,gauss);
				break;
			case 7:
				GetAlpha2Coulomb(palpha2,gauss);
				break;
			case 8:
				GetAlpha2Shakti(palpha2,gauss);
				break;
			case 9:
				GetAlpha2Josh(palpha2,gauss);
				break;
			case 10:
				GetAlpha2PISM(palpha2,gauss);
				break;
			case 11:
				GetAlpha2Schoof(palpha2,gauss);
				break;
			case 12:
				GetAlpha2Tsai(palpha2,gauss);
				break;
			case 13:
				GetAlpha2Coulomb2(palpha2,gauss);
				break;
			case 14:
				GetAlpha2RegCoulomb(palpha2,gauss);
				break;
			case 15:
				GetAlpha2RegCoulomb2(palpha2,gauss);
				break;
			default:
				_error_("Friction law "<< this->law <<" not supported");
		}
	}
	else if(this->linearize==1){
		this->element->ValueP1OnGauss(palpha2, this->alpha2_list, gauss);
	}
	else if(this->linearize==2){
		*palpha2 = this->alpha2_list[0];
	}
	else{
		_error_("not supported yet");
	}

	/*Checks*/
	_assert_(!xIsNan<IssmDouble>(*palpha2));
	_assert_(!xIsInf<IssmDouble>(*palpha2));
	_assert_(*palpha2>=0);

}/*}}}*/
void Friction::GetAlpha2Coulomb(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	  alpha2= drag^2 * Neff ^r * | vel | ^(s-1), with Neff=rho_ice*g*thickness+rho_ice*g*base, r=q/p and s=1/p
	  alpha2= min(drag^2 * Neff ^r * | vel | ^(s-1), drag_coulomb^2 * Neff*/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  drag_p, drag_q;
	IssmDouble  drag_coefficient,drag_coefficient_coulomb;
	IssmDouble  alpha2,alpha2_coulomb;

	/*Recover parameters: */
	element->GetInputValue(&drag_p,gauss,FrictionPEnum);
	element->GetInputValue(&drag_q,gauss,FrictionQEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	element->GetInputValue(&drag_coefficient_coulomb, gauss,FrictionCoefficientcoulombEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity   = element->FindParam(ConstantsGEnum);

	//compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	/*Get effective pressure*/
	bool ispwHydro,ispwStochastic;
   IssmDouble Neff;
   element->parameters->FindParam(&ispwStochastic,StochasticForcingIsWaterPressureEnum);
   element->parameters->FindParam(&ispwHydro,HydrologyIsWaterPressureArmaEnum);
   if(ispwStochastic || ispwHydro){
      /*Retrieve pre-computed water pressure and compute ice pressure*/
      IssmDouble p_ice,p_water,Neff_limit;
      element->GetInputValue(&p_water,gauss,FrictionWaterPressureEnum);
      element->parameters->FindParam(&Neff_limit,FrictionEffectivePressureLimitEnum);
      p_ice = IcePressure(gauss);
      Neff  = max(Neff_limit*p_ice, p_ice - p_water);
   }
   else{
      /*Compute effective pressure directly*/
      Neff = EffectivePressure(gauss);
   }

	/*Get velocity magnitude*/
	IssmDouble vmag = VelMag(gauss);

	if(s==1.){
		/*This is to make AD happy and avoid 0^0*/
		alpha2=drag_coefficient*drag_coefficient*pow(Neff,r);
	}
	else{
		/*Check to prevent dividing by zero if vmag==0*/
		if(vmag==0. && (s-1.)<0.) alpha2=0.;
		else alpha2=drag_coefficient*drag_coefficient*pow(Neff,r)*pow(vmag,(s-1.));
	}

	if(vmag==0.){
		alpha2_coulomb=0.;
	}
	else{
		//alpha2_coulomb=drag_coefficient_coulomb*drag_coefficient_coulomb*rho_ice*gravity*max(0.,thickness-floatation_thickness)/vmag;
		alpha2_coulomb=drag_coefficient_coulomb*drag_coefficient_coulomb*Neff/vmag;
	}

	if(alpha2_coulomb<alpha2) alpha2=alpha2_coulomb;

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Hydro(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
		Based on Gagliardini 2007, needs a good effective pressure computation
		Not tested so far so use at your own risks
	  alpha2= NeffC[Chi/(1+alpha*Chi^q)]^(1/n)*|vel|^(-1)  with
		 Chi=|vel|/(C^n*Neff^n*As)
		 alpha=(q-1)^(q-1)/q^q */

	/*diverse: */
	IssmDouble  q_exp;
	IssmDouble  C_param;
	IssmDouble  As;
	IssmDouble  n;
	IssmDouble  alpha;
	IssmDouble  Chi,Gamma;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&q_exp,gauss,FrictionQEnum);
	element->GetInputValue(&C_param,gauss,FrictionCEnum);
	element->GetInputValue(&As,gauss,FrictionAsEnum);
	element->GetInputValue(&n,gauss,MaterialsRheologyNEnum);

	/*Get effective pressure*/
	IssmDouble Neff = EffectivePressure(gauss);
	IssmDouble vmag = VelMag(gauss);

	//compute alpha and Chi coefficients: */
	if (q_exp==1){
		alpha=1;
	}
	else{
		alpha=(pow(q_exp-1,q_exp-1))/pow(q_exp,q_exp);
	}
	Chi=vmag/(pow(C_param,n)*pow(Neff,n)*As);
	Gamma=(Chi/(1. + alpha * pow(Chi,q_exp)));
	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0.) alpha2=0.;
	else	if (Neff==0) alpha2=0.0;
	else	alpha2=Neff * C_param * pow(Gamma,1./n) * pow(vmag,-1);

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Shakti(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/* FrictionGetAlpha2 computes alpha2= drag^2 * Neff, with Neff=rho_ice*g*thickness+rho_ice*g*(head-base)*/

	/*diverse: */
	IssmDouble  pressure_ice,pressure_water;
	IssmDouble  Neff;
	IssmDouble  drag_coefficient;
	IssmDouble  base,thickness,head,sealevel;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&base, gauss,BaseEnum);
	element->GetInputValue(&head, gauss,HydrologyHeadEnum);
	element->GetInputValue(&sealevel, gauss,SealevelEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	IssmDouble rho_water   = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble rho_ice     = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->FindParam(ConstantsGEnum);

	//From base and thickness, compute effective pressure when drag is viscous:
	pressure_ice   = rho_ice*gravity*thickness;
	pressure_water = rho_water*gravity*(head-base+sealevel);
	Neff=pressure_ice-pressure_water;
	if(Neff<0.) Neff=0.;

	alpha2=drag_coefficient*drag_coefficient*Neff;

	/*Assign output pointers:*/
	*palpha2=alpha2;
}
/*}}}*/
void Friction::GetAlpha2Temp(IssmDouble* palpha2, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_viscous * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  f,T,pressure,Tpmp,gamma;
	IssmDouble  alpha2;

	/*Get viscous part*/
	this->GetAlpha2Budd(&alpha2,gauss);

	/*Get pressure melting point (Tpmp) for local pressure and get current temperature*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);

	/*Compute scaling parameter*/
	element->parameters->FindParam(&gamma,FrictionGammaEnum);
	alpha2 = alpha2 / (exp((T-Tpmp)/gamma)+1e-3);

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Josh(IssmDouble* palpha2, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_viscous * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  T,Tpmp,deltaT,deltaTref,pressure,diff,drag_coefficient;
	IssmDouble  alpha2,time,gamma,ref,alp_new,alphascaled,max_coefficient;
	const IssmDouble yts = 365*24*3600.;

	/*Get viscous part*/
	this->GetAlpha2Budd(&alpha2,gauss);

	/*Get delta Refs*/
	element->GetInputValue(&deltaTref,gauss,FrictionPressureAdjustedTemperatureEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	/*New*/
	/*element->GetInputValue(&deltaTrefsfc,gauss,FrictionSurfaceTemperatureEnum);
	 *    element->GetInputValue(&Tpdd,gauss,TemperaturePDDEnum);
	 *       */

	element->parameters->FindParam(&max_coefficient,FrictionMaxCoefficientEnum);

	/*Compute delta T*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);
	deltaT = T-Tpmp;

	/*Compute gamma*/
	element->parameters->FindParam(&time,TimeEnum);
	element->parameters->FindParam(&gamma,FrictionGammaEnum);

	ref = exp(deltaTref/gamma);
	alp_new = ref/exp(deltaT/gamma);

	alphascaled = sqrt(alp_new)*drag_coefficient;
	if (alphascaled > max_coefficient) alp_new = (max_coefficient/drag_coefficient)*(max_coefficient/drag_coefficient);

	alp_new=alp_new*alpha2;

	/*Assign output pointers:*/
	*palpha2=alp_new;
}/*}}}*/
void Friction::GetAlpha2Budd(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	  alpha2= drag^2 * Neff ^r * | vel | ^(s-1), with Neff=rho_ice*g*thickness+rho_ice*g*base, r=q/p and s=1/p**/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  drag_p, drag_q;
	IssmDouble  drag_coefficient;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&drag_p,gauss,FrictionPEnum);
	element->GetInputValue(&drag_q,gauss,FrictionQEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);

	/*compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	/*Get effective pressure and basal velocity*/
	IssmDouble vmag = VelMag(gauss);

	bool ispwHydro,ispwStochastic;
   IssmDouble Neff;
   element->parameters->FindParam(&ispwStochastic,StochasticForcingIsWaterPressureEnum);
   element->parameters->FindParam(&ispwHydro,HydrologyIsWaterPressureArmaEnum);
   if(ispwStochastic || ispwHydro){
      /*Retrieve pre-computed water pressure and compute ice pressure*/
      IssmDouble p_ice,p_water,Neff_limit;
      element->GetInputValue(&p_water,gauss,FrictionWaterPressureEnum);
		element->parameters->FindParam(&Neff_limit,FrictionEffectivePressureLimitEnum);
      p_ice = IcePressure(gauss);
      Neff  = max(Neff_limit*p_ice, p_ice - p_water);
   }	
	else{
		/*Compute effective pressure directly*/
		Neff = EffectivePressure(gauss);
	}

	/*Check to prevent dividing by zero if vmag==0*/
	if(s==1.){
		/*This is to make AD happy and avoid 0^0*/
		alpha2=drag_coefficient*drag_coefficient*pow(Neff,r);
	}
	else{
		if(vmag==0. && (s-1.)<0.) alpha2=0.;
		else alpha2=drag_coefficient*drag_coefficient*pow(Neff,r)*pow(vmag,(s-1.));
	}

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2WaterLayer(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	  alpha2= drag^2 * Neff ^r * | vel | ^(s-1), with Neff=rho_ice*g*thickness+rho_ice*g*base, r=q/p and s=1/p**/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  drag_p, drag_q;
	IssmDouble  Neff,F;
	IssmDouble  thickness,base,sealevel;
	IssmDouble  drag_coefficient,water_layer;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->parameters->FindParam(&F,FrictionFEnum);
	element->GetInputValue(&drag_p,gauss,FrictionPEnum);
	element->GetInputValue(&drag_q,gauss,FrictionQEnum);
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&base, gauss,BaseEnum);
	element->GetInputValue(&sealevel, gauss,SealevelEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	element->GetInputValue(&water_layer, gauss,FrictionWaterLayerEnum);
	IssmDouble rho_water   = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->FindParam(ConstantsGEnum);

	//compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	//From base and thickness, compute effective pressure when drag is viscous:
	if(base>0) base=0;
	if(water_layer==0) Neff=gravity*rho_ice*thickness+gravity*rho_water*(base-sealevel);
	else if(water_layer>0) Neff=gravity*rho_ice*thickness*F;
	else _error_("negative water layer thickness");
	if(Neff<0) Neff=0;

	IssmDouble vmag = VelMag(gauss);

	if(s==1.){
		/*This is to make AD happy and avoid 0^0*/
		alpha2=drag_coefficient*drag_coefficient*pow(Neff,r);
	}
	else{
		/*Check to prevent dividing by zero if vmag==0*/
		if(vmag==0. && (s-1.)<0.) alpha2=0.;
		else alpha2=drag_coefficient*drag_coefficient*pow(Neff,r)*pow(vmag,(s-1.));
	}

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Weertman(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient alpha2= C^2 |v|^(1/m-1) */

	/*diverse: */
	IssmDouble  C,m;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&C,gauss,FrictionCEnum);
	element->GetInputValue(&m,gauss,FrictionMEnum);

	/*Get velocity magnitude*/
	IssmDouble vmag = VelMag(gauss);

	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0. && (1./m-1.)<0.) alpha2=0.;
	else alpha2=C*C*pow(vmag,(1./m-1.));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2WeertmanTemp(IssmDouble* palpha2, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_weertman * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  f,T,pressure,Tpmp,gamma;
	IssmDouble  alpha2;

	/*Get viscous part*/
	this->GetAlpha2Weertman(&alpha2,gauss);

	/*Get pressure melting point (Tpmp) for local pressure and get current temperature*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);

	/*Compute scaling parameter*/
	element->parameters->FindParam(&gamma,FrictionGammaEnum);
	alpha2 = alpha2 / exp((T-Tpmp)/gamma);

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2PISM(IssmDouble* palpha2, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction using a pseudoplastic friction law,
	 * computing the basal shear stress as
	 *
	 * alpha2 = tau_c (u_b/(abs(u_b)^(1-q)*u_0^q))
	 *
	 * The yield stress tau_c is a function of the effective pressure N
	 * using a Mohr-Coloumb criterion, so that
	 * tau_c = tan(phi)*N,
	 * where phi is the till friction angle, representing sediment strength
	 *
	 * The effective pressure is given by Eq. (5) in Aschwanden et al. 2016:
	 *
	 * N = delta * P0 * 10^((e_0/Cc)(1-(W/Wmax)))
	 *
	 * W is calculated by a non-conserving hydrology model in HydrologyPismAnalysis.cpp
	 *
	 * see Aschwanden et al. 2016 and Bueler and Brown, 2009 for more details
	 */

	/*compute ice overburden pressure P0*/
	IssmDouble thickness,base,P0;
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&base, gauss,BaseEnum);
	//element->GetInputValue(&sealevel, gauss,SealevelEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity   = element->FindParam(ConstantsGEnum);
	P0 = gravity*rho_ice*thickness;

	/*Compute effective pressure*/
	IssmDouble  N,delta,W,Wmax,e0,Cc;
	element->parameters->FindParam(&delta,FrictionDeltaEnum);
	element->parameters->FindParam(&e0,FrictionVoidRatioEnum);
	element->GetInputValue(&Cc,gauss,FrictionSedimentCompressibilityCoefficientEnum);
	element->GetInputValue(&W,gauss,WatercolumnEnum);
	element->GetInputValue(&Wmax,gauss,HydrologyWatercolumnMaxEnum);

 	/*Check that water column height is within 0 and upper bound, correct if needed*/
 	if(W>Wmax) W=Wmax;
 	if(W<0)    W=0.;

	N = delta*P0*pow(10.,(e0/Cc)*(1.-W/Wmax));

	/*Get till friction angles, defined by user [deg]*/
	IssmDouble phi;
	element->GetInputValue(&phi,gauss,FrictionTillFrictionAngleEnum);

	/*Convert till friction angle from user-defined deg to rad, which Matlab uses*/
	phi = phi*PI/180.;

	/*Compute yield stress following a Mohr-Colomb criterion*/
	IssmDouble tau_c = N*tan(phi);

	/*Compute basal speed*/
	IssmDouble ub;
	element->GetInputValue(&ub,gauss,VelEnum);

	/*now compute alpha^2*/
	IssmDouble u0,q;
	element->parameters->FindParam(&u0,FrictionThresholdSpeedEnum);
	element->parameters->FindParam(&q,FrictionPseudoplasticityExponentEnum);
	IssmDouble alpha2 = tau_c/(pow(ub+1.e-10,1.-q)*pow(u0,q));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Schoof(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	 *
	 *               C^2 |u_b|^(m-1)
	 * alpha2= __________________________
	 *          (1+(C^2/(Cmax Neff))^1/m |u_b| )^m
	 *
	 * */

	/*diverse: */
	IssmDouble  C,Cmax,m,alpha2;

	/*Recover parameters: */
	element->GetInputValue(&Cmax,gauss,FrictionCmaxEnum);
	element->GetInputValue(&C,gauss,FrictionCEnum);
	element->GetInputValue(&m,gauss,FrictionMEnum);

	/*Get effective pressure*/
	bool ispwStochastic;
	IssmDouble Neff;
	element->parameters->FindParam(&ispwStochastic,StochasticForcingIsWaterPressureEnum);
	if(ispwStochastic){
		/*Retrieve stochastic water pressure and compute ice pressure*/
		IssmDouble p_ice,p_water,Neff_limit;
		element->GetInputValue(&p_water,gauss,FrictionWaterPressureEnum);
		element->parameters->FindParam(&Neff_limit,FrictionEffectivePressureLimitEnum);
		p_ice = IcePressure(gauss);
		Neff  = max(Neff_limit*p_ice, p_ice - p_water);
	}	
	else{
		/*Compute effective pressure directly*/
		Neff = EffectivePressure(gauss);
	}

	/*Get velocity magnitude*/
	IssmDouble ub = VelMag(gauss);

	/*Compute alpha^2*/
	if((ub<1e-10) ||(Neff==0.0)){
		alpha2 = 0.;
	}
	else{
		alpha2= (C*C*pow(ub,m-1.)) / pow(1.+  pow(C*C/(Cmax*Neff),1./m)*ub,m);
	}

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Tsai(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	 *
	 * alpha2= min(C |ub|^m , f N ) / |ub|
	 *
	 * */

	/*diverse: */
	IssmDouble  C,f,m,alpha2;

	/*Recover parameters: */
	element->GetInputValue(&f,gauss,FrictionfEnum);
	element->GetInputValue(&C,gauss,FrictionCEnum);
	element->GetInputValue(&m,gauss,FrictionMEnum);

	/*Get effective pressure and velocity magnitude*/
	IssmDouble N  = EffectivePressure(gauss);
	IssmDouble ub = VelMag(gauss);

	/*Compute alpha^2*/
	if(ub<1e-10){
		alpha2 = 0.;
	}
	else{
		alpha2= C*pow(ub,m);

		if(alpha2>f*N) alpha2 = f*N;

		alpha2 = alpha2/ub;
	}

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Coulomb2(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	 *
	 *               C^2 |u_b|^(m-1) * (.5*N)
	 * alpha2= ___________________________________
	 *          (C^(2/m) |u_b| + (0.5*N)^(1/m) )^m
	 *
	 * */

	/*diverse: */
	IssmDouble  C,m,alpha2;

	/*Recover parameters: */
	element->GetInputValue(&C,gauss,FrictionCEnum);
	element->GetInputValue(&m,gauss,FrictionMEnum);

	/*Get effective pressure and velocity magnitude*/
	IssmDouble N  = EffectivePressure(gauss);
	IssmDouble ub = VelMag(gauss);

	/*Compute alpha^2*/
	if(ub<1e-10){
		alpha2 = 0.;
	}
	else{
		alpha2= (pow(C,2)*pow(ub,m-1.)*(0.5*N)) / pow(pow(C,2./m)*ub + pow(0.5*N,1./m),m);
	}

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2RegCoulomb(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	 *
	 *               C |u_b|^(1/m-1)
	 * alpha2= __________________________
	 *          (|u_b|/u0 + 1 )^(1/m)
	 *
	 * */

	/*diverse: */
	IssmDouble  C,coeff,u0,m,alpha2;

	/*Recover parameters: */
	element->GetInputValue(&coeff,gauss,FrictionCEnum);
	element->GetInputValue(&m,gauss,FrictionMEnum);
	element->parameters->FindParam(&u0,FrictionU0Enum);

	/* scale C for a better inversion */
	C = coeff*coeff;

	/*Get velocity magnitude*/
	IssmDouble ub = VelMag(gauss);

	/*Check to prevent dividing by zero if vmag==0*/
	if(ub==0.) {
		alpha2=0.;
	}
	else {
		/*Compute alpha^2*/
		alpha2= (C*pow(ub,1./m-1.)) / pow(ub/u0 + 1.,1./m);
	}

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2RegCoulomb2(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient
	 *
	 *               C N |u_b|^(1/m-1)
	 * alpha2= __________________________
	 *          (|u_b| + K N^m )^(1/m)
	 *
	 * */

	/*diverse: */
	IssmDouble  C,K,m,alpha2;

	/*Recover parameters: */
	element->GetInputValue(&C,gauss,FrictionCEnum);
	element->GetInputValue(&m,gauss,FrictionMEnum);
	element->GetInputValue(&K,gauss,FrictionKEnum);

	/*Get velocity magnitude*/
	IssmDouble ub = VelMag(gauss);
	IssmDouble Neff = EffectivePressure(gauss);

	/*Check to prevent dividing by zero if vmag==0*/
	if(ub==0. && (m-1.)<0) {
		alpha2=0.;
	}
	else {
		/*Compute alpha^2*/
		alpha2= (C*pow(ub,1./m-1.)) * Neff / pow((ub+pow(K*Neff,m)),1./m);
	}

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
IssmDouble Friction::EffectivePressure(Gauss* gauss){/*{{{*/
	/*Get effective pressure as a function of  flag */

	/*diverse: */
	int         coupled_flag;
	IssmDouble  thickness,base,sealevel;
	IssmDouble  p_ice,p_water;
	IssmDouble  Neff,Neff_limit;

	/*Recover parameters: */
	element->parameters->FindParam(&coupled_flag,FrictionCouplingEnum);
	element->parameters->FindParam(&Neff_limit,FrictionEffectivePressureLimitEnum);

	/*Compute ice pressure*/
	p_ice = IcePressure(gauss);

	/*From base and thickness, compute effective pressure when drag is viscous, or get Neff from forcing:*/
	switch(coupled_flag){
		case 0:{
			element->GetInputValue(&base, gauss,BaseEnum);
			element->GetInputValue(&sealevel, gauss,SealevelEnum);
			IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);
			IssmDouble gravity   = element->FindParam(ConstantsGEnum);
			p_water = rho_water*gravity*(sealevel-base);
			Neff = p_ice - p_water;
		}
			break;
		case 1:{
			p_water = 0.;
			Neff = p_ice - p_water;
		}
			break;
		case 2:{
			element->GetInputValue(&base, gauss,BaseEnum);
			element->GetInputValue(&sealevel, gauss,SealevelEnum);
			IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);
			IssmDouble gravity   = element->FindParam(ConstantsGEnum);
			p_water = max(0.,rho_water*gravity*(sealevel-base));
			Neff = p_ice - p_water;
		}
			break;
		case 3:{
			element->GetInputValue(&Neff,gauss,FrictionEffectivePressureEnum);
		}
			break;
		case 4:{
			element->GetInputValue(&Neff,gauss,EffectivePressureEnum);
		}
			break;
		default:
			_error_("not supported");
	}

	/*Make sure Neff is positive*/
	if(Neff<Neff_limit*p_ice) Neff=Neff_limit*p_ice;

	/*Return effective pressure*/
	return Neff;

}/*}}}*/
IssmDouble Friction::IcePressure(Gauss* gauss){/*{{{*/
	/*Get ice pressure*/

	IssmDouble  thickness,p_ice;
	/*Recover Inputs and Parameters*/
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	IssmDouble rho_ice = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity = element->FindParam(ConstantsGEnum);

	/*Compute*/
	p_ice = gravity*rho_ice*thickness;

	/*Return ice pressure*/
	return p_ice;

}/*}}}*/
IssmDouble Friction::SubglacialWaterPressure(Gauss* gauss){/*{{{*/
	/*Get water pressure as a function of  flag */

	int         coupled_flag;
	IssmDouble  base,sealevel,p_water;

	/*Recover parameters: */
	element->parameters->FindParam(&coupled_flag,FrictionCouplingEnum);

	switch(coupled_flag){
		case 0:{
			element->GetInputValue(&base, gauss,BaseEnum);
			element->GetInputValue(&sealevel, gauss,SealevelEnum);
			IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);
			IssmDouble gravity   = element->FindParam(ConstantsGEnum);
			p_water = rho_water*gravity*(sealevel-base);
		}
			break;
		case 1:{
			p_water = 0.;
		}
			break;
		case 2:{
			element->GetInputValue(&base, gauss,BaseEnum);
			element->GetInputValue(&sealevel, gauss,SealevelEnum);
			IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);
			IssmDouble gravity   = element->FindParam(ConstantsGEnum);
			p_water = max(0.,rho_water*gravity*(sealevel-base));
		}
			break;
		case 3:{
			_error_("water pressure not computed for coupling==3 in friction law");
		}
			break;
		case 4:{
			_error_("water pressure not computed for coupling==4 in friction law");
		}
			break;
		default:
			_error_("not supported");
	}

	/*Return water pressure*/
	return p_water;

}/*}}}*/
IssmDouble Friction::VelMag(Gauss* gauss){/*{{{*/
	/*Get the velocity magnitude as a function of flag */

	/*diverse*/
	IssmDouble vx,vy,vz,vmag;

	this->vx_input->GetInputValue(&vx, gauss);
	this->vy_input->GetInputValue(&vy, gauss);

	if ((this->vz_input == NULL) || (this->apply_dim<3.)) vz = 0.0;
	else this->vz_input->GetInputValue(&vz, gauss);

	if (this->apply_dim<2.) vy = 0.0;

	vmag = sqrt(vx*vx+vy*vy+vz*vz);
	return vmag;
}/*}}}*/
void Friction::GetBasalSlidingSpeeds(IssmDouble* pvx, Gauss* gauss){/*{{{*/

	this->vx_input->GetInputValue(pvx, gauss);
	/*Checks*/
	_assert_(!xIsNan<IssmDouble>(*pvx));
	_assert_(!xIsInf<IssmDouble>(*pvx));
}/*}}}*/
void Friction::GetBasalSlidingSpeeds(IssmDouble* pvx, IssmDouble* pvy, Gauss* gauss){/*{{{*/

	this->vx_input->GetInputValue(pvx, gauss);
	this->vy_input->GetInputValue(pvy, gauss);
	/*Checks*/
	_assert_(!xIsNan<IssmDouble>(*pvx));
	_assert_(!xIsInf<IssmDouble>(*pvx));
	_assert_(!xIsNan<IssmDouble>(*pvy));
	_assert_(!xIsInf<IssmDouble>(*pvy));
}/*}}}*/
void Friction::GetBasalSlidingSpeeds(IssmDouble* pvx, IssmDouble* pvy, IssmDouble* pvz, Gauss* gauss){/*{{{*/

	this->vx_input->GetInputValue(pvx, gauss);
	this->vy_input->GetInputValue(pvy, gauss);
	this->vz_input->GetInputValue(pvz, gauss);
	/*Checks*/
	_assert_(!xIsNan<IssmDouble>(*pvx));
	_assert_(!xIsInf<IssmDouble>(*pvx));
	_assert_(!xIsNan<IssmDouble>(*pvy));
	_assert_(!xIsInf<IssmDouble>(*pvy));
	_assert_(!xIsNan<IssmDouble>(*pvz));
	_assert_(!xIsInf<IssmDouble>(*pvz));
}/*}}}*/

/*IO*/
void FrictionUpdateInputs(Elements* elements,Inputs* inputs,IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int    frictionlaw;
	int    frictioncoupling;

	/*Friction law variables*/
	iomodel->FindConstant(&frictionlaw,"md.friction.law");
	switch(frictionlaw){
		case 1:
			iomodel->FindConstant(&frictioncoupling,"md.friction.coupling");
			iomodel->FetchDataToInput(inputs,elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.q",FrictionQEnum);
			if(frictioncoupling==3){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",FrictionEffectivePressureEnum);}
			else if(frictioncoupling==4){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",EffectivePressureEnum);
			}
			break;
		case 2:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.m",FrictionMEnum);
			break;
		case 3:
			iomodel->FindConstant(&frictioncoupling,"md.friction.coupling");
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.As",FrictionAsEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.q",FrictionQEnum);
			if(frictioncoupling==3){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",FrictionEffectivePressureEnum);}
			else if(frictioncoupling==4){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",EffectivePressureEnum);
			}
			break;
		case 4:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.q",FrictionQEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.temperature",TemperatureEnum);
			iomodel->FindConstant(&frictioncoupling,"md.friction.coupling");
			break;
		case 5:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.q",FrictionQEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.water_layer",FrictionWaterLayerEnum);
			break;
		case 6:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.m",FrictionMEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.temperature",TemperatureEnum);
			break;
		case 7:
			iomodel->FindConstant(&frictioncoupling,"md.friction.coupling");
			iomodel->FetchDataToInput(inputs,elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.coefficientcoulomb",FrictionCoefficientcoulombEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.q",FrictionQEnum);
			if(frictioncoupling==3){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",FrictionEffectivePressureEnum);}
			else if(frictioncoupling==4){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",EffectivePressureEnum);

			}
			break;
		case 8:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.coefficient",FrictionCoefficientEnum);
			break;
		case 9:
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.temperature",TemperatureEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.pressure_adjusted_temperature",FrictionPressureAdjustedTemperatureEnum);
			InputUpdateFromConstantx(inputs,elements,1.,FrictionPEnum);
			InputUpdateFromConstantx(inputs,elements,1.,FrictionQEnum);
			break;
		case 10:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.till_friction_angle",FrictionTillFrictionAngleEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.sediment_compressibility_coefficient",FrictionSedimentCompressibilityCoefficientEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.hydrology.watercolumn_max",HydrologyWatercolumnMaxEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",WatercolumnEnum,0.);
			break;
		case 11:
			iomodel->FindConstant(&frictioncoupling,"md.friction.coupling");
			iomodel->FetchDataToInput(inputs,elements,"md.friction.m",FrictionMEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.Cmax",FrictionCmaxEnum);
			if(frictioncoupling==3){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",FrictionEffectivePressureEnum);}
			else if(frictioncoupling==4){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",EffectivePressureEnum);
			}
			break;
		case 12:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.m",FrictionMEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.f",FrictionfEnum);
			break;
		case 13:
			iomodel->FindConstant(&frictioncoupling,"md.friction.coupling");
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.m",FrictionMEnum);
			if(frictioncoupling==3){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",FrictionEffectivePressureEnum);}
			else if(frictioncoupling==4){
				iomodel->FetchDataToInput(inputs,elements,"md.friction.effective_pressure",EffectivePressureEnum);
			}
			break;
		case 14:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.m",FrictionMEnum);
			break;
		case 15:
			iomodel->FetchDataToInput(inputs,elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.m",FrictionMEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.friction.K",FrictionKEnum);
			break;
		default:
			_error_("friction law "<< frictionlaw <<" not supported");
	}

}/*}}}*/
void FrictionUpdateParameters(Parameters* parameters,IoModel* iomodel){/*{{{*/

	parameters->AddObject(iomodel->CopyConstantObject("md.friction.law",FrictionLawEnum));

	/*Set default linearize parameter to 0 for now*/
	parameters->AddObject(new IntParam(FrictionLinearizeEnum,0));

	int frictionlaw;
	iomodel->FindConstant(&frictionlaw,"md.friction.law");
	switch(frictionlaw){
		case 1:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.linearize",FrictionLinearizeEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.coupling",FrictionCouplingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		case 2:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.linearize",FrictionLinearizeEnum));
			break;
		case 3:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.coupling",FrictionCouplingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		case 4:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.gamma",FrictionGammaEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.coupling",FrictionCouplingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		case 5:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.f",FrictionFEnum));
			break;
		case 6:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.gamma",FrictionGammaEnum));
			break;
		case 7:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.coupling",FrictionCouplingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		case 8:
			break;
		case 9:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.gamma",FrictionGammaEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.coefficient_max",FrictionMaxCoefficientEnum));
			parameters->AddObject(new IntParam(FrictionCouplingEnum,2));/*comment this line to use effective pressure from Beuler and Pelt (2015)*/
			break;
		case 10:
			parameters->AddObject(new IntParam(FrictionCouplingEnum,2)); /*comment this line to use effective pressure from Beuler and Pelt (2015)*/
			parameters->AddObject(new DoubleParam(FrictionEffectivePressureLimitEnum,0.));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.pseudoplasticity_exponent",FrictionPseudoplasticityExponentEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.threshold_speed",FrictionThresholdSpeedEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.delta",FrictionDeltaEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.void_ratio",FrictionVoidRatioEnum));
			break;
		case 11:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.coupling",FrictionCouplingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		case 12:
			parameters->AddObject(new IntParam(FrictionCouplingEnum,2));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		case 13:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.coupling",FrictionCouplingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		case 14:
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.u0",FrictionU0Enum));
			break;
		case 15:
			parameters->AddObject(new IntParam(FrictionCouplingEnum,2));
			parameters->AddObject(iomodel->CopyConstantObject("md.friction.effective_pressure_limit",FrictionEffectivePressureLimitEnum));
			break;
		default: _error_("Friction law "<<frictionlaw<<" not implemented yet");
	}

}/*}}}*/
