/* \file NyeCO2.cpp
 * \brief figure out B of CO2 ice for a certain temperature
 *		INPUT function B=NyeCO2(temperature)
 *    	where rigidigty (in s^(1/n)Pa) is the flow law paramter in the flow law sigma=B*e(1/n) (Nye, p2000). 
 */

#include "../io/io.h" 
#include <math.h> 
#include "../Numerics/types.h"

IssmDouble NyeCO2(IssmDouble temperature){

	/*Coefficients*/
	const IssmPDouble Rg      = 8.3144598;     /* J mol^-1 K^-1   */ 
	const IssmPDouble A_const = pow(10.,13.0); /* s^-1 MPa        */ 
	const IssmPDouble Q       = 66900.;        /* J mol^-1        */ 
	const IssmPDouble n       = 8.;            /* Glen's exponent */

	/*Arrhenius Law*/
	IssmDouble A = A_const *exp(-Q/(temperature*Rg));  /* s^-1 MPa   */
	IssmDouble B = 1e6*pow(A,-1/n);                    /* s^(1/n) Pa */

	/*Beyond-melting-point cases*/
	if((temperature>200.)&&(temperature<220.)) _printf0_("CO2 ICE - POSSIBLE MELTING. Some temperature values are between 200K and 220K.\n");
	else if(temperature>=220.) _printf0_("CO2 ICE - GUARANTEED MELTING. Some temperature values are beyond 220K.\n");

	/*Return output*/
	return B; 
}
