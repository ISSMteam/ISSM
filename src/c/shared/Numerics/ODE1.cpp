#include <math.h>
#include "./types.h"
#include "../Exceptions/exceptions.h"

IssmDouble ODE1(IssmDouble alpha,IssmDouble beta,IssmDouble Si, IssmDouble dt,int method){
	/* solve the following equation:
	 *
	 *   dS/dt = alpha S + beta
	 *
	 *   method 0: Forward Euler (explicit)
	 *   method 1: backward Euler (implicit)
	 *   method 2: Crank Nicolson
	 *
	 *   return S^{i+1} based on  Si, dt, alpha and beta
	 */

	switch(method){
		case 0: return Si*(1.+alpha*dt) + beta*dt;
		case 1: return (Si+beta*dt)/(1.-alpha*dt);
		case 2: return (Si*(1.+alpha*dt/2.) + beta*dt)/(1-alpha*dt/2.);
		default: _error_("not supported yet");
	}
}
