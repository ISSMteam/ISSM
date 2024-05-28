/*!\file Xy2llx.cpp
 */

#include "../../shared/shared.h"
#include "./latlong.h"
#include <math.h>

int Xy2llx(double* lat, double* lon, double* x, double* y, int ncoord, int sgn){
/*  This is a cpp conversion of the following:
%XY2LL - converts xy to lat long
%
%   Converts Polar  Stereographic (X,Y) coordinates for the polar regions to
%   latitude and longitude Stereographic (X,Y) coordinates for the polar
%   regions.
%   Author: Michael P. Schodlok, December 2003 (map2xy.m)
%
%   Usage:
%      [lat,lon] = xy2ll(x,y,sgn);
%      [lat,lon] = xy2ll(x,y,sgn,central_meridian,standard_parallel);
%
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)
*/
	double  central_meridian,standard_parallel;

	Xy2lldef(&central_meridian,&standard_parallel,sgn);

	return(Xy2llx(lat,lon,x,y,ncoord,sgn,central_meridian,standard_parallel));
}

int Xy2llx(double* lat, double* lon, double* x, double* y, int ncoord, int sgn, double central_meridian, double standard_parallel){
/*  This is a cpp conversion of the following:
%XY2LL - converts xy to lat long
%
%   Converts Polar  Stereographic (X,Y) coordinates for the polar regions to
%   latitude and longitude Stereographic (X,Y) coordinates for the polar
%   regions.
%   Author: Michael P. Schodlok, December 2003 (map2xy.m)
%
%   Usage:
%      [lat,lon] = xy2ll(x,y,sgn);
%      [lat,lon] = xy2ll(x,y,sgn,central_meridian,standard_parallel);
%
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)
*/

	int     i,iret=0;
	double  delta,slat;
	double  re,ex2,ex;
	double  sl,rho,cm,T,chi;

	if((sgn!=1) && (sgn!=-1)) _error_("Sign should be either +1 or -1.\n");

	delta = central_meridian;
	slat  = standard_parallel;

	/*  Radius of the earth in meters  */
	re   = 6378.273e+3;
	/*  Eccentricity of the Hughes ellipsoid squared  */
	ex2  = 0.006693883;
	/*  Eccentricity of the Hughes ellipsoid  */
	ex   =  sqrt(ex2);

	/*  loop over all the coordinate pairs  */
	for(i=0; i<ncoord; i++){
		sl = slat*PI/180.;
		cm = cos(sl)/sqrt(1.0-ex2*(pow(sin(sl),2)));
		rho= sqrt(pow(x[i],2) + pow(y[i],2));
		T  = tan((PI/4.0) - (sl/2.0))/pow(((1.0-ex*sin(sl))/(1.0+ex*sin(sl))),(ex/2.0));

		if(fabs(slat-90.) < 1.e-5)
			T =rho*sqrt(pow((1.+ex),(1.+ex))*pow((1.-ex),(1.-ex)))/2./re;
		else
			T =rho*T/(re*cm);

		chi = (PI / 2.0) - 2.0 * atan(T);
		lat[i] = chi + ((ex2 / 2.0) + (5.0 * pow(ex2,2.0) / 24.0) + (pow(ex2,3.0) / 12.0)) *
			   sin(2.0 * chi) + ((7.0 * pow(ex2,2.0) / 48.0) + (29.0 * pow(ex2,3.0) / 240.0)) *
			   sin(4.0 * chi) + (7.0 * pow(ex2,3.0) / 120.0) * sin(6.0 * chi) ;

		lat[i] = (double)sgn * lat[i];
		lon[i] = atan2((double)sgn * x[i],-(double)sgn * y[i]);
		lon[i] = (double)sgn * lon[i];

		if(rho <= 0.1){
			lat[i] = 90. * (double)sgn;
			lon[i] = 0.0;
			iret=1;
		}

		lon[i] = lon[i] * 180. / PI;
		lat[i] = lat[i] * 180. / PI;
		lon[i] = lon[i] - delta; 
	}

	return(iret);
}

void Xy2lldef(double* pdelta, double* pslat, int sgn){
/*  This is a cpp conversion of the following:
%XY2LL - converts xy to lat long
%
%   Converts Polar  Stereographic (X,Y) coordinates for the polar regions to
%   latitude and longitude Stereographic (X,Y) coordinates for the polar
%   regions.
%   Author: Michael P. Schodlok, December 2003 (map2xy.m)
%
%   Usage:
%      [lat,lon] = xy2ll(x,y,sgn);
%      [lat,lon] = xy2ll(x,y,sgn,central_meridian,standard_parallel);
%
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)
*/
	bool    flag=true;

	/*  Get central_meridian and standard_parallel depending on hemisphere  */
	if (sgn == 1) {
		*pdelta= 45;
		*pslat = 70;
		if(flag) _printf0_("Warning: expecting coordinates in polar stereographic (Std Latitude: 70N Meridian: 45).\n");
	}
	else if (sgn == -1) {
		*pdelta= 0;
		*pslat = 71;
		if(flag) _printf0_("Warning: expecting coordinates in polar stereographic (Std Latitude: 71S Meridian: 0).\n");
	}
	else _error_("Sign should be either +1 or -1.\n");

	return;
}
