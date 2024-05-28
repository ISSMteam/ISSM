#include <math.h>

#include "./numerics.h"

IssmDouble CBRT(IssmDouble Z){

	IssmDouble ret;

	if (Z> 0.0){
		ret = fabs(pow(fabs(Z),1./3.));
	}
	else if(Z< 0.0){
		ret = - fabs(pow(fabs(Z),1./3.));
	}
	else{
		ret = 0.;
	}
	return ret;
}

int cubic(IssmDouble a,IssmDouble b,IssmDouble c,IssmDouble d, IssmDouble x[3], int* num){
	/* Find the real roots of linear/quadratic and cubic equations:
	 *
	 *   a x^3 + bx^2 + cx + d = 0
	 *
	 *   returns the roots in x
	 *   num is the number of roots */

	/*Some useful constants*/
	const IssmDouble pi    = 3.1415926535897932;
	const IssmDouble third = 1./3.;

	/*Intermediaries*/
	IssmDouble U[3],W,P,Q,delta,phi;

	/* determine the degree of the polynomial */
	if (a != 0.0){
		//cubic problem
		W     = b/a *third;
		P     = pow((c/a *third - pow(W,2)),3);
		Q     = -.5 *(2.0*pow(W,3)-(c*W-d)/a );
		delta = pow(Q,2)+P;
		if ( delta < 0.0 ){
			//three real solutions!
			//Confine the argument of coeffCOS to the interval [-1;1]!
			phi = acos(min(1.0,max(-1.0,Q/sqrt(-P))));
			P   = 2.0*pow((-P),(5.e-1*third));
			for(int i=0;i<3;i++)	U[i] = P*cos((phi+2*((IssmDouble)i)*pi)*third)-W;
			x[0] = min(U[0], min(U[1], U[2]));
			x[1] = max(min(U[0], U[1]),max( min(U[0], U[2]), min(U[1], U[2])));
			x[2] = max(U[0], max(U[1], U[2]));
			*num = 3;
		}
		else{
			// only one real solution!
			delta = sqrt(delta);
			x[0] = CBRT(Q+delta)+CBRT(Q-delta)-W;
			*num=1;
		}
	}
	else if (b != 0.0){
		// quadratic problem
		P     = 0.5*c/b;
		delta = pow(P,2)-d/b;
		if (delta > 0.0){
			// 2 real solutions
			x[0] = -P - sqrt(delta);
			x[1] = -P + sqrt(delta);
			*num = 2;
		}
		else{
			// no real solution
			*num = 0;
		}
	}
	else if (c != 0.0){
		//linear equation
		x[0] = d/c;
		*num = 1;
	}
	else{
		//no equation
		*num = 0;
	}

	/* perform one step of a newton iteration in order to minimize round-off errors */
	for(int i=0;i<*num;i++){
		x[i] = x[i] - (d+x[i]*(c+x[i]*(b+x[i]*a)))/(c+x[i]*(2.0*b+x[i]*3.0*a));
	}

	return 0;
}
