/* \file StressIntensityIntegralWeight.cpp
 * \Weight to integrate the stress along the ice flow direction to compute the stress intensity factor : Linear Fracture Mechanics (see "Combining damage and fracture mechanics to model calving",Krug 2014 in the appendix) 
 */

#include <math.h>

#include "../Numerics/types.h"

IssmDouble StressIntensityIntegralWeight(IssmDouble depth, IssmDouble water_depth, IssmDouble thickness){

	/*output: */
	IssmDouble beta;

	/*intermediaries: */
	IssmDouble M1,M2,M3,x,y,d;
	const double pi = 3.141592653589793;
	x    = water_depth/thickness;
	y    = depth;
	d    = water_depth;

	M1   = 0.0719768-1.513476*x-61.1001*pow(x,2)+1554.95*pow(x,3)-14583.8*pow(x,4)+71590.7*pow(x,5)-205384*pow(x,6)+356469*pow(x,7)-368270*pow(x,8)+208233*pow(x,9)-49544*pow(x,10);
	//printf("M1 : %g",M1);
	M2   = 0.246984+6.47583*x+176.456*pow(x,2)-4058.76*pow(x,3)+37303.8*pow(x,4)-181755*pow(x,5)+520551*pow(x,6)-904370*pow(x,7)+936863*pow(x,8)-531940*pow(x,9)+127291*pow(x,10);
	//printf("M2 : %g",M2);
	M3   = 0.529659-22.3235*x+532.074*pow(x,2)-5479.53*pow(x,3)+28592.2*pow(x,4)-81388.6*pow(x,5)+128746*pow(x,6)-106246*pow(x,7)+35780.7*pow(x,8);
	//printf("M3 : %g",M3);

	beta = 2/sqrt(2*pi*(d-y))*(1+M1*sqrt(1-y/d)+M2*(1-y/d)+M3*pow((1-y/d),1.5));

	return beta;
}
