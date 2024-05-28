function [lat,lon,scale_factor] = xy2ll(x,y,sgn,central_meridian,standard_parallel)
%XY2LL - converts xy to lat long
%
%   Converts Polar  Stereographic (X,Y) coordinates for the polar regions to
%   latitude and longitude Stereographic (X,Y) coordinates for the polar
%   regions.
%   Optional scale factor provides the scaling factor needed to correct projection error
%   in areas and volumes
%   Author: Michael P. Schodlok, December 2003 (map2xy.m)
%
%   Usage:
%      [lat,lon] = xy2ll(x,y,sgn);
%      [lat,lon,scale_factor] = xy2ll(x,y,sgn);
%      [lat,lon] = xy2ll(x,y,sgn,central_meridian,standard_parallel);
%
%      - sgn = Sign of latitude		1 : north latitude (default is mer=45 lat=70)
%                              	   -1 : south latitude (default is mer=0  lat=71)

%Get central_meridian and standard_parallel depending on hemisphere
if nargin==5,
	delta = central_meridian;
	slat  = standard_parallel;
elseif nargin==3
	if sgn == 1,
		delta = 45; slat = 70;
		disp('Warning: expecting coordinates in polar stereographic (Std Latitude: 70ºN Meridian: 45º)');
	elseif sgn==-1,
		delta = 0;  slat = 71;
		disp('Warning: expecting coordinates in polar stereographic (Std Latitude: 71ºS Meridian: 0º)');
	else
		error('Sign should be either 1 or -1');
	end
else
	help xy2ll
	error('bad usage');
end

%Choose ellipsoid
if 0
	%Hughes ellipsoid
	re   = 6378.273*10^3; % Radius of the earth in meters
	ex2 = .006693883;     % Eccentricity of the Hughes ellipsoid squared
else
	%WGS84 ellipsoid
	re = 6378137;         % Radius of the earth in meters
	f  = 1./298.257223563;% Earth flattening
	ex2 = 2*f-f^2;        % Eccentricity squared
end

% Eccentricity 
ex = sqrt(ex2);

sl  = slat*pi/180.;
rho = sqrt(x.^2 + y.^2);
cm = cos(sl) / sqrt(1.0 - ex2 * (sin(sl)^2));
T = tan((pi / 4.0) - (sl / 2.0)) / ((1.0 - ex * sin(sl)) / (1.0 + ex * sin(sl)))^(ex / 2.0);

if  abs(slat-90.) < 1.e-5
	T = rho * sqrt((1. + ex)^(1. + ex) * (1. - ex)^(1. - ex)) / 2. / re;
else
	T = rho * T / (re * cm);
end

chi = (pi / 2.0) - 2.0 * atan(T);
lat = chi + ((ex2 / 2.0) + (5.0 * ex2^2.0 / 24.0) + (ex2^3.0 / 12.0)) * ...
	sin(2 * chi) + ((7.0 * ex2^2.0 / 48.0) + (29.0 * ex2^3 / 240.0)) * ...
	sin(4.0 * chi) + (7.0 * ex2^3.0 / 120.0) * sin(6.0 * chi) ;

lat = sgn * lat;
lon = atan2(sgn * x,-sgn * y);
lon = sgn * lon;

[res1,res2] = find(rho(:) <= 0.1);
if res1
	lat(res1) = pi/2. * sgn;
	lon(res1) = 0.0;
end

lon = lon * 180. / pi;
lat = lat * 180. / pi;
lon = lon - delta;
if nargout==3,
	m=((1+sin(abs(slat)*pi/180))*ones(length(lat),1)./(1+sin(abs(lat)*pi/180)));
	scale_factor=(1./m).^2;
end
