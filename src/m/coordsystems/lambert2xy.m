function [x,y] = lambert2xy(lat,lon,sgn,projection_center_lat,projection_center_lon)  
%LAMBERT2XY - converts lat long from Lambert Azimuthal to Polar Stereographic
%
%   Converts from geodetic latitude and longitude that are 
%   in Lambert Azimuthal (equal area) projections to Polar 
%   Stereographic (X,Y) coordinates for the polar regions.
%
%   Usage:
%      [x,y] = lambert2xy(lat,lon,sgn)
%      [x,y] = lambert2xy(lat,lon,sgn,projection_center_lat,projection_center_lon)
%
%      - provide lat in [-90,90] and lon in [-180,180].

%      - sgn = +1 N hemisphere [default projection center lat = 90 lon=0]
%              -1 S hemisphere [default projection center lat = -90 lon=0]

%Get projection_center_lat and projection_center_lon 
if nargin==5,
	latitude0  = projection_center_lat;
	longitude0 = projection_center_lon;
elseif nargin==3,
	if sgn==1,
		latitude0 = 90; longitude0 = 0;
		disp('Info: creating coordinates in polar stereographic (Projection center lat: 90N lon: 0)');
	elseif sgn==-1,
		latitude0 = -90; longitude0 = 0;
		disp('Info: creating coordinates in polar stereographic (Projection center lat: 90S lon: 0)');
	else
		error('Sign should be either +1 or -1');
	end
else
	help lambert2xy
	error('bad usage');
end

% Radius of the earth in meters 
a = 6378137.0;
% Eccentricity of the Hughes ellipsoid squared
e = 0.081819191;

% Projection center latitude and longitude in radians 
phi0 = latitude0 * pi/180; 
lam0 = longitude0 * pi/180; 

% Some constant based on phi0 and lam0
qp= (1-e^2)*((1/(1-e^2))-((1/(2*e))*log((1-e)/(1+e))));
q0=(1-e^2)*((sin(phi0)/(1-e^2*sin(phi0)*sin(phi0)))-((1/(2*e))*log((1-e*sin(phi0))/(1+e*sin(phi0)))));
Rq=a*sqrt(qp/2);
b0=asin(q0/qp);
D =a*(cos(phi0)/sqrt(1-e^2*sin(phi0)*sin(phi0)))/(Rq*cos(b0));

% Latitude and longitude in radians 
phi = lat*pi/180;
lam = lon*pi/180;

% Some other phi,lam dependent parameters 
q=(1-e^2)*((sin(phi)/(1-e^2*sin(phi)*sin(phi)))-((1/(2*e))*log((1-e*sin(phi))/(1+e*sin(phi)))));
b =asin(q/qp);
B =Rq*sqrt(2/(1+sin(b0)*sin(b)+(cos(b0)*cos(b)*cos(lam-lam0))));

% Calculation of x and y
if(abs(lat)==90)
	x=0.0; 
	y=0.0; 
else
	x=(B*D)*(cos(b)*sin(lam-lam0));
	y=(B/D)*((cos(b0)*sin(b))-(sin(b0)*cos(b)*cos(lam-lam0)));
end
