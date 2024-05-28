%% ring2z
% Find HEALPix ring number corresponding to a given z coordinate

%% Syntax
%  z = ring2z(nSide,nRing)

%% Input Arguments
%  nSide    HEALPix resolution parameter
%  nRing    numeric array of ring numbers (in range [1,4*nSide-1])

%% Return Arguments
%  z        size(nRing) numeric array of pixel center z coordinates

%% Example
nRing = reshape(1:15,3,5)
z = ring2z(4,nRing)

%% See also
% z2ring

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.