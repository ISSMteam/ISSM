%% ang2pix
% Convert spherical coordinate location to HEALPix pixel number

%% Syntax
%  pix = ang2pix(nSide, v);
%  pix = ang2pix(nSide, v, nest);

%% Input Arguments
%  nSide      HEALPix resolution parameter
%  v          (theta, phi) pairs in radians as either [2,n] or as cellarray
%             of [2,1] 
%  nest       (optional) nested pixel number scheme (true | {false})

%% Return Arguments
%  pix             array of pix numbers:
%                 v on input     output
%                 iscell(v)      size(v) cell array of pixel numbers
%                 isnumeric(v)   size(v,1) numeric array of pixel numbers
%
% Note that Matlab numbers pixels from 1

%% Example
vTheta = acos(2*rand(3,5)-1);
vPhi = 2*pi*rand(3,5);
vc = arrayfun(@(x,y)([x;y]),vTheta,vPhi,'UniformOutput',false);
va = reshape(cell2mat(vc),[2,3,5]);
nSide = 8;
size(vc)
size(va)
pix = ang2pix(nSide,vc)
pix = ang2pix(nSide,va)

%% See also
% vec2pix, pix2ang, pix2vec

%% Requires
% ring2nest

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.