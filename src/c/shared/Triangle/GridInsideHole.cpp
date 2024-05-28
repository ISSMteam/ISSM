/*
 * GridInsideHole.c:
 * from a convex set of points, figure out a point that for sure lies inside the profile.
 */

#include <math.h>

#include "./triangle.h"
#include "../Exp/exp.h"

#undef M_PI
#define M_PI 3.141592653589793238462643

int GridInsideHole(double* px0,double* py0,int n,double* x,double* y){

	double flag=0.0;
	double xA,xB,xC,xD,xE;
	double yA,yB,yC,yD,yE;

	/*Take first and last vertices: */
	xA=x[0];
	yA=y[0];
	xB=x[n-1];
	yB=y[n-1];

	/*Figure out middle of segment [A B]: */
	xC=(xA+xB)/2;
	yC=(yA+yB)/2;

	/*D and E are on each side of segment [A B], on the median line between segment [A  B], 
	 *at an angle of 10 degree (less than the minimum 30 enforced by the quality of the mesh: */
	xD=xC+tan(10./180.*M_PI)*(yC-yA);
	yD=yC+tan(10./180.*M_PI)*(xA-xC);
	xE=xC-tan(10./180.*M_PI)*(yC-yA);
	yE=yC-tan(10./180.*M_PI)*(xA-xC);

	/*Either E or D is inside profile (x,y): */
	IsInPolySerial(&flag,&xD,&yD,1,x,y,n,2);
	/*FIXME: used to be 'flag' and not '!flag', check*/
	if(!flag){
		/*D is inside the poly: */
		*px0=xD;
		*py0=yD;
	}
	else{
		/*E is inside the poly: */
		*px0=xE;
		*py0=yE;
	}
	return 1;
}
