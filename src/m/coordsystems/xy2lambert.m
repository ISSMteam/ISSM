function [lat,lon] = xy2lambert(x,y,sgn,projection_center_lat,projection_center_lon)  
%XY2LAMBERT - converts xy to lat lon in Lambert Azimuthal
%
%   Converts from Ploar Stereographic (X,Y) coordinates to geodetic 
%   lat lon that are in Lambert Azimuthal (equal area) projection.
%
%   Usage:
%      [lat,lon] = xy2lambert(x,y,sgn)
%      [lat,lon] = xy2lambert(x,y,sgn,projection_center_lat,projection_center_lon)
%
%      - provide lat in [-90,90] and lon in [-180,180].
%
%      - sgn = +1 N hemisphere [default projection center lat = 90 lon=0]
%              -1 S hemisphere [default projection center lat = -90 lon=0]

%Get projection_center_lat and projection_center_lon 
if nargin==5,
	latitude0  = projection_center_lat;
	longitude0 = projection_center_lon;
elseif nargin==3,
	if sgn==1,
		latitude0 = 90; longitude0 = 0;
		disp('Info: creating coordinates in Lambert Azimuthal equal-area (Projection center lat: 90N lon: 0)');
	elseif sgn==-1,
		latitude0 = -90; longitude0 = 0;
		disp('Info: creating coordinates in Lambert Azimuthal equal-area (Projection center lat: 90S lon: 0)');
	else
		error('Sign should be either +1 or -1');
	end
else
	help xy2lambert
	error('bad usage');
end

% Radius of the earth in meters 
a = 6378137.0;
% Eccentricity of the Hughes ellipsoid squared
e = 0.081819191;

% Projection center latitude and longitude in radians 
phi0 = latitude0 * pi/180; 
lam0 = longitude0 * pi/180; 

% Some constants based on phi0 and lam0
% (as in forward calculation)
qp= (1-e^2)*((1/(1-e^2))-((1/(2*e))*log((1-e)/(1+e))));
q0=(1-e^2)*((sin(phi0)/(1-e^2*sin(phi0)*sin(phi0)))-((1/(2*e))*log((1-e*sin(phi0))/(1+e*sin(phi0)))));
Rq=a*sqrt(qp/2);
b0=asin(q0/qp);
D =a*(cos(phi0)/sqrt(1-e^2*sin(phi0)*sin(phi0)))/(Rq*cos(b0));

% Some other (x,y) dependent parameters 
rho=sqrt((x/D)^2+(D*y)^2);
C=2*asin(rho/(2*Rq));
b_prime=asin((cos(C)*sin(b0))+((D*y*sin(C)*cos(b0))/rho));

% Calculation of lat and lon 
dist=sqrt(x^2+y^2);
if(dist<=0.1)
	lat=sgn*90.0;
	lon=0.0;
else
	lat_rad=b_prime+((e^2/3+31*e^4/180+517*e^6/5040)*sin(2*b_prime))+((23*e^4/360+251*e^6/3780)*sin(4*b_prime))+((761*e^6/45360)*sin(6*b_prime));
	lon_rad=lam0+atan(x*sin(C)/(D*rho*cos(b0)*cos(C)-D^2*y*sin(b0)*sin(C)));
	% in degrees 
	lat=lat_rad*180/pi;
	lon=lon_rad*180/pi;
end
