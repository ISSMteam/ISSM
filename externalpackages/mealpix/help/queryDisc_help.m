%% queryDisc
% Find all pixels whose centers are within a specified disc

%% Syntax
%  pix = queryDisc(nSide,c0,r,'Param1',Value1,...)

%% Input Arguments
%  nSide     HEALPix resolution parameter
%  c0        disc center. May be specified as a pixel number, spherical
%            coordinate location [theta, phi] in rads, or as a cartesian
%            unit vector [x;y;z].   
%  r         disc radius (radians)
% 
%  Param     Value
%  'nest'    pixels are returned using nested indexing (true | {false})

%% Return Arguments
%  pix      list of pixels in disc

%% Example
nSide=32;
% Pixels within pi/10 radians of north pole
pix = queryDisc(nSide,[0;0;1],pi/10);
numel(pix)
nSide2nPix(nSide)*((pi/10)/(4*pi))
% Pixels within pi/10 radians of south pole
pix = queryDisc(nSide,[0;0;-1],pi/10);
numel(pix)
% Pixels within pi/10 radians of 0 long, lat
pix = queryDisc(nSide,[1;0;0],pi/10);
numel(pix)

%% See also
% queryPolygon, queryStrip, queryTriangle

%% Requires
% pix2vec, ang2vec, angDist

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.