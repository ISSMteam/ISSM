/* \file LliboutryDuval.cpp
 * \brief figure out B of ice for a certain temperature and water fraction or enthalpy
 */

#include <math.h>
#include "../Numerics/types.h"
#include "../Exceptions/exceptions.h"

/* get ice stiffness B from enthalpy, pressure and flow law exponent*/
IssmDouble LliboutryDuval(IssmDouble enthalpy, IssmDouble pressure, IssmDouble n, IssmDouble betaCC, IssmDouble referencetemperature, IssmDouble heatcapacity, IssmDouble latentheat){
  /*Use Lliboutry & Duval's 1985 parameterization for the rheology: 
	* see also: Grewe/Blatter 2009, Aschwanden et al. 2012
	*
	* ISSM uses enthalpy/temperature values that are not corrected for pressure.
   *
   *  A(H,p) = A0 exp(-Q/RT(H,p)), if H < H_s(p)
   *         = A0 exp(-Q/RTpmp) (1+181.25w(H,p)), if H_s(p) \le H < H_l(p)
   *  
   *  T(H,p) = Tref + H/c_i, if H < H_s(p)
   *         = Tpmp , if H_s(p) \le H \le H_l(p)
   *
   *  w(H,p) = 0, if H < H_s(p)
   *         = (H - H_s(p))/L
   *
   *  H_s(p) = c_i (Tpmp - Tref)
   *
   *  Tpmp   = T - betaCC p;
   *
   *  A0 constant of proportionality
   *     = 3.61 * 10^-13   if T*<263.15K
   *     = 1.73 * 10^3     if T*>263.15K
   *  Q  Activation energy for creep
   *     = 6.0  * 10^4     if T*<263.15K
   *     = 13.9 * 10^4     if T*>263.15K
   *  R  Universal gas constant
   *     = 8.314
   *  
   *  Convert A to B :  B = A^(-1/n) */
	/*check feasibility*/
	_assert_(pressure+1.e-4>=0); // deal with pressure instability at ice surface
	_assert_(n>0);
	_assert_(betaCC>=0);
	_assert_(referencetemperature>=0);
	_assert_(heatcapacity>0);
	_assert_(latentheat>0);

	/*Some physical constants*/
	IssmDouble R=8.314; 

	/*Intermediaries*/
	IssmDouble A,B,Tstar,Tpmp,H_sp,waterfraction;

	Tpmp=273.15-betaCC*pressure; //pressure melting point temperature
	H_sp=heatcapacity*(Tpmp-referencetemperature); //pressure melting point enthalpy
	if (enthalpy<H_sp){ //compute homologous temperature and water fraction
		Tstar=referencetemperature+enthalpy/heatcapacity+betaCC*pressure; 
		waterfraction=0.;
	}
	else{
		Tstar=273.15;
		waterfraction=(enthalpy-H_sp)/latentheat;
		if (waterfraction>0.01) waterfraction=0.01; // limit softness of ice
	}

	/*Get A*/
	if(Tstar<=263.15){A=3.61e-13*exp(-6.e+4/(R*Tstar));}
	else{A=1.73e3*exp(-13.9e+4/(R*Tstar));}
	A*=(1.+181.25*waterfraction);

	/*Convert to B*/
	B=pow(A,-1./n);
	return B;
}
