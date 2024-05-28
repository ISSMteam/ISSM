#include <math.h>

#include "../MemOps/MemOps.h"
#include "../Exceptions/exceptions.h"
#include "../Numerics/types.h"
#include "./isnan.h"

void XZvectorsToCoordinateSystem(IssmDouble* T,IssmDouble* xzvectors){

	IssmDouble	x[3],y[3],z[3];
	IssmDouble	x_norm, y_norm, z_norm;

	for(int i=0;i<6;i++){
		if(xIsNan<IssmDouble>(xzvectors[i])){
			/*At least one NaN found: default to Id*/
			T[0*3+0] = 1.0;	T[0*3+1] = 0.0;	T[0*3+2] = 0.0;
			T[1*3+0] = 0.0;	T[1*3+1] = 1.0;	T[1*3+2] = 0.0;
			T[2*3+0] = 0.0;	T[2*3+1] = 0.0;	T[2*3+2] = 1.0;

			return;
		}
	}

	/* get input {x} (vector in local x-z plane): */
	x[0] = xzvectors[0];
	x[1] = xzvectors[1];
	x[2] = xzvectors[2];

	/* get input {z} (local tangent plane normal vector): */
	z[0] = xzvectors[3];
	z[1] = xzvectors[4];
	z[2] = xzvectors[5];

	/* compute {y} = {z} x {x}: */
	y[0] =  x[2]*z[1] - x[1]*z[2];
	y[1] = -x[2]*z[0] + x[0]*z[2];
	y[2] =  x[1]*z[0] - x[0]*z[1];

	/* normalise {x}, {y} and {z} to form unit vectors {i_hat}, {j_hat} and {k_hat};
		store in {x}, {y}, and {z}: */
	x_norm = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	y_norm = sqrt( y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
	z_norm = sqrt( z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);

	x[0] = x[0]/x_norm;		x[1] = x[1]/x_norm;		x[2] = x[2]/x_norm;
	y[0] = y[0]/y_norm;		y[1] = y[1]/y_norm;		y[2] = y[2]/y_norm;
	z[0] = z[0]/z_norm;		z[1] = z[1]/z_norm;		z[2] = z[2]/z_norm;

	/* Tlg columns are just {i_hat}, {j_hat} and {k_hat}, respectively: */
	T[0*3+0] = x[0];	T[0*3+1] = y[0];	T[0*3+2] = z[0];
	T[1*3+0] = x[1];	T[1*3+1] = y[1];	T[1*3+2] = z[1];
	T[2*3+0] = x[2];	T[2*3+1] = y[2];	T[2*3+2] = z[2];
}
