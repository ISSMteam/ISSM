/*!\file Ll2xyx.cpp
 */

#include "../../shared/shared.h"
#include "./latlong.h"
#include <math.h>

int Ll2xyx(double* x, double* y, double* lat, double* lon, int ncoord, int sgn){
/*  This is a cpp conversion of the following:
%LL2XY - converts lat long to polar stereographic
%
%   Converts from geodetic latitude and longitude to Polar 
%   Stereographic (X,Y) coordinates for the polar regions.
%   Author: Michael P. Schodlok, December 2003 (map2ll)
%
%   Usage:
%      [x,y] = ll2xy(lat,lon,sgn)
%      [x,y] = ll2xy(lat,lon,sgn,central_meridian,standard_parallel)
%
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)
*/
	double  central_meridian,standard_parallel;

	Ll2xydef(&central_meridian,&standard_parallel,sgn);

	return(Ll2xyx(x,y,lat,lon,ncoord,sgn,central_meridian,standard_parallel));
}

int Ll2xyx(double* x, double* y, double* lat, double* lon, int ncoord, int sgn, double central_meridian, double standard_parallel){
/*  This is a cpp conversion of the following:
%LL2XY - converts lat long to polar stereographic
%
%   Converts from geodetic latitude and longitude to Polar 
%   Stereographic (X,Y) coordinates for the polar regions.
%   Author: Michael P. Schodlok, December 2003 (map2ll)
%
%   Usage:
%      [x,y] = ll2xy(lat,lon,sgn)
%      [x,y] = ll2xy(lat,lon,sgn,central_meridian,standard_parallel)
%
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)
*/

	int     i,iret=0;
	double  delta,slat;
	double  re,ex2,ex;
	double  latitude,longitude;
	double  T,rho,sl,tc,mc;

	if((sgn!=1) && (sgn!=-1)) _error_("Sign should be either +1 or -1.\n");

	delta = central_meridian;
	slat  = standard_parallel;

	/*  Radius of the earth in meters  */
	re  = 6378.273*1.e3;
	/*  Eccentricity of the Hughes ellipsoid squared  */
	ex2 = 0.006693883;
	/*  Eccentricity of the Hughes ellipsoid  */
	ex  =  sqrt(ex2);

	/*  loop over all the coordinate pairs  */
	for(i=0; i<ncoord; i++){
		latitude  = fabs(lat[i]) * PI/180.;
		longitude = (lon[i] + delta) * PI/180.;

		/*  compute X and Y in grid coordinates.  */
		T = tan(PI/4.-latitude/2.) / pow(((1.-ex*sin(latitude))/(1.+ex*sin(latitude))),(ex/2.));

		if ((90. - slat) < 1.e-5)
			rho = 2.*re*T/sqrt(pow((1.+ex),(1.+ex))*pow((1.-ex),(1.-ex)));
		else {
			sl  = slat*PI/180.;
			tc  = tan(PI/4.-sl/2.)/pow(((1.-ex*sin(sl))/(1.+ex*sin(sl))),(ex/2.));
			mc  = cos(sl)/sqrt(1.0-ex2*(pow(sin(sl),2)));
			rho = re*mc*T/tc;
		}

		y[i]= -rho*(double)sgn*cos(sgn*longitude);
		x[i]=  rho*(double)sgn*sin(sgn*longitude);

		if (latitude>= PI/2.){
			x[i] = 0.0;
			y[i] = 0.0;
			iret=1;
		}
	}
	return(iret);
}

void Ll2xydef(double* pdelta, double* pslat, int sgn){
/*  This is a cpp conversion of the following:
%LL2XY - converts lat long to polar stereographic
%
%   Converts from geodetic latitude and longitude to Polar 
%   Stereographic (X,Y) coordinates for the polar regions.
%   Author: Michael P. Schodlok, December 2003 (map2ll)
%
%   Usage:
%      [x,y] = ll2xy(lat,lon,sgn)
%      [x,y] = ll2xy(lat,lon,sgn,central_meridian,standard_parallel)
%
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)
*/
	bool    flag=true;

	/*  Get central_meridian and standard_parallel depending on hemisphere  */
	if (sgn ==  1) {
		*pdelta= 45;
		*pslat = 70;
		if(flag) _printf0_("Info: creating coordinates in polar stereographic (Std Latitude: 70N Meridian: 45).\n");
	}
	else if (sgn == -1) {
		*pdelta= 0;
		*pslat = 71;
		if(flag) _printf0_("Info: creating coordinates in polar stereographic (Std Latitude: 71S Meridian: 0).\n");
	}
	else _error_("Sign should be either +1 or -1.\n");

	return;
}
