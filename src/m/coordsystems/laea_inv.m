function [lat,lon] = laea_inv(x,y,lat1,lon0,falseeasting,falsenorthing)
% laea_inv performs an inverse Lambert Azimuthal Equal Area projection
% for a simple spherical Earth of radius 6371000 meters. 
% 
%% Syntax 
% 
%  [lat,lon] = laea_inv(x,y,lat1,lon0)
%  [lat,lon] = laea_inv(x,y,lat1,lon0,falseeasting,falsenorthing)
% 
%% Description 
% 
% [lat,lon] = laea_inv(x,y,lat1,lon0) transforms the coordinates x,y (in 
% meters) into geographic coordinates lat,lon. Inputs lat1 and lon0 specify 
% the origin.
%
% [lat,lon] = laea_inv(x,y,lat1,lon0,falseeasting,falsenorthing) also allows
% inclusion of false eastings and northings in meters. 
% 
%% Author Info 
% Function written by Chad A. Greene of NASA Jet Propulsion Laboratory. 
% December 2020. 
% Formulas taken directly from Snyder 1987's classic tome, "Map Projections 
% A Working Manual" starting around page 185. 

%% Parse inputs: 

narginchk(4,6)
assert(isequal(size(lat1),size(lon0),[1 1]),'Error: Inputs lat1 and lon0 must both be scalars.')
assert(abs(lat1)<=90,'lat1 cannot exceed +/-90 degrees.') 
assert(abs(lon0)<=360,'lon0 cannot exceed +/-360 degrees.') 
assert(isequal(size(x),size(y)),'Dimensions of x and y must match.') 

if nargin<5
   falseeasting = 0; 
   falsenorthing = 0; 
end

%% 
% Projection formulas:
% From Snyder 1987, MAP PROJECTIONS-A WORKING MANUAL, page 185: https://pubs.usgs.gov/pp/1395/report.pdf

% Account for false easting, northing: 
x = x - falseeasting; 
y = y - falsenorthing; 

% Define constants: 
R = 6371000; % earth radius (meters) 
rho = hypot(x,y); 
c = 2*asind(rho./(2*R)); 

% Unproject: 
lat = asind(cosd(c) .* sind(lat1) + (y .* sind(c) .* cosd(lat1))./rho); 

switch lon0
   case 90
      lon = lon0 + atand(-x./y);
   case -90
      lon = lon0 + atand(x./y); 
   otherwise
      lon = lon0 + atand((x .* sind(c))./(rho.*cosd(lat1) .*cosd(c) - y.*sind(lat1).*sind(c))); 
end

end