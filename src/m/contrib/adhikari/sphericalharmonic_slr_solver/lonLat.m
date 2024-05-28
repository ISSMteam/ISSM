function lon_lat = lonLat(nSide,nPix)

%---------------------------------------------------------------------
% lonLat :: a function to compute (long, lat) in degrees for GMT 
%---------------------------------------------------------------------
% This code is written as a part of the ISSM-PSL project
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     November 3, 2014
%---------------------------------------------------------------------

thetaLambda = pix2ang(nSide); 
data = cell2mat(thetaLambda); 
data = data*180/pi;  % in degrees 
data(1,:) = -data(1,:)+90;  % lat \in [-90, 90] in GMT 
% 
lon_lat = zeros(2,nPix); 
lon_lat(1,:) = data(2,:); 
lon_lat(2,:) = data(1,:); 

