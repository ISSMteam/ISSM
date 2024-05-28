/* \file Arrhenius.cpp
 * \brief figure out B of ice for a certain temperature
 */

#include <math.h>
#include "../Numerics/types.h"
#include "../Exceptions/exceptions.h"

IssmDouble Arrhenius(IssmDouble temperature,IssmDouble depth,IssmDouble n){
	/*Use EISMINT Parameterization for the rheology: Payne2000
	 *
	 *  A(T*) = A0 exp(-Q/RT*)
	 *
	 *  A0 constant of proportionality
	 *     = 3.61 * 10^-13   if T*<263.15K
	 *     = 1.73 * 10^3     if T*>263.15K
	 *  Q  Activation energy for creep
	 *     = 6.0  * 10^4     if T*<263.15K
	 *     = 13.9 * 10^4     if T*>263.15K
	 *  R  Universal gas constant
	 *     = 8.314
	 *  T* Absolute temperature corrected for the dependence of Tpmp on P
	 *     = T - beta (s-z)
	 *
	 *  Convert A to B :  B = A^(-1/n) */

	/*Some physical constants (Payne2000)*/
	IssmDouble beta=8.66*1.e-4;
	IssmDouble R=8.314;

	/*Intermediaries*/
	IssmDouble A,B,Tstar;

	/*convert temperature to absolute temperature*/
	_assert_(depth>0);
	Tstar=temperature-beta*depth;
	_assert_(Tstar>0);

	/*Get A*/
	if(Tstar<263.15){
		A=3.61e-13*exp(  -6.e+4/(R*Tstar));
	}
	else{
		A=1.73e+3 *exp(-13.9e+4/(R*Tstar));
	}

	/*Convert to B*/
	B=pow(A,-1./n);

	return B;
}
