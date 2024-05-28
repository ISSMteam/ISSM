%% xy2pix
% Find (nest indexed) MEALPix pixel number corresponding to a pixel index
% on a HEALPix face 

%% Syntax
% ipix = xy2pix(nSide, ix, iy, nFace)

%% Input Arguments
%  nside     HEALPix resolution parameter (power of 2)
%  ix,iy     integer pixel location on face (1 <= ix, iy <= nSide)
%  nFace     HEALPix face  number (1 <= nFace <= 12)
%
% One or both ix and iy may be scalars. If neither is a scalar both must
% have the same size. 

%% Return Arguments
%  ipix      NESTED pixel index for given x,y on the face

%% Description
% At lowest resolution HEALPix divides the sphere into 12 faces, each with
% its own coordinate system. Each face is further subdivided in a regular
% fashion as the resolution parameter is increased. Given a face number
% (from 1 to 12), a resolution parameter, and a location on the face XY2PIX
% identifies the pixel number (in the nested indexing scheme).

%% Example
nSide = 4;
nFace = 3;
ix = [1 2; 1 2]; 
iy = [1 1; 2 2];
nPix = xy2pix(nSide,ix,iy,nFace)

%% See also
% pix2xy

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.