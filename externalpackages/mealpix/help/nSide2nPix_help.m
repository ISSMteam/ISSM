%% nSide2nPix
% Convert HEALPix resolution parameter to number of pixels on sphere

%% Syntax
%  nPix = nSide2nPix(nSide)

%% Input Arguments
%  nSide    sub-divisions of facet sides. May be an array. 

%% Return Arguments
%  nPix     size(nSide) number of pixels on sphere corresponding to nSide

%% Description
% HEALPix divides the sphere into 12 basic facets. Each facet is further
% divided into nSide sub-divisions for a total of 12*nSide^2 pixels. 

%% Example
nSide = 2.^[1:4];
nPix = nSide2nPix(nSide)

%% See also
% nPix2nSide

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.