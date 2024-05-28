#include <math.h>
#include "./types.h"
#include "../Exceptions/exceptions.h"

int NewtonSolveDnorm(IssmDouble* pdnorm,IssmDouble c1,IssmDouble c2,IssmDouble c3,IssmDouble n,IssmDouble dnorm){
	/* solve the following equation using Newton-Raphson
	 *
	 *   c1*x^(s-1) + c2*x = c3
	 *
	 *   s = (1+n)/n
	 *
	 *   we solve y = 10^x:
	 */

	/*trivial solution*/
	if(c3==0.){
		*pdnorm = 0.;
		return 0;
	}

	/*Intermediaries*/
	int        counter = 0;
	IssmDouble s = (1.+n)/n;
	IssmDouble y2;
	IssmDouble threshold = 1.e-12;

	/*Initial guess*/
	_assert_(dnorm>0.); 
	IssmDouble y1 = log10(dnorm);

	while(true){

		/*Newton step*/
		y2 = y1 - (c1*pow(pow(10.,y1),s-1.) + c2*pow(10.,y1) - c3)/((s-1)*c1*log(10.)*pow(pow(10.,y1),s-1.) + c2*log(10.)*pow(10.,y1));

		if( fabs(y2-y1)/(fabs(y2))<threshold ){
			break;
		}
		else{
			y1=y2;
			counter++;
		}

		if(counter>50) break;
	}

	/*Avoid extremely large values that indicate non convergence*/
	if(y2>50.) y2 = 50;

	/*Assign output pointer*/
	*pdnorm = pow(10.,y2);
	return 0;
}
