/* \file NyeH2O.cpp
 * brief figure out B of H2O ice for a certain temperature
 *		INPUT function B=NyeH2O(temperature)
 *    	where rigidigty (in s^(1/n)Pa) is the flow law paramter in the flow law sigma=B*e(1/n) (Nye, p2000). 
 */

#include "../io/io.h"
#include <math.h>
#include "../Numerics/types.h"

IssmDouble NyeH2O(IssmDouble temperature){

	/*Coefficients*/
	const IssmPDouble Rg      = 8.3144598; /* J mol^-1 K^-1 */
	const IssmPDouble A_const = 9.e4;      /*s^-1 MPa       */
	const IssmPDouble Q       = 60000.;    /*J mol^-1       */
	const IssmPDouble n       = 3.;        /*Glen's exponent*/

	/*Arrhenius Law*/
	IssmDouble A = A_const *exp(-Q/(temperature*Rg));  /*s^-1 MPa   */
	IssmDouble B = 1e6*pow(A,-1/n);                    /*s^(1/n) Pa */

	/*Beyond-melting-point case*/
	if(temperature>=273.15) _printf0_("H2O ICE - GUARANTEED MELTING. Some temperature values are beyond 273.15K.\n");

	/*Return output*/
	return B;
}
