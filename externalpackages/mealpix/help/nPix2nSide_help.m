%% nPix2nSide 
% Finds number of HEALPix facet side subdivisions from number of pixels on
% sphere

%% Syntax
%  nSide = nPix2nSide(nPix);

%% Input Arguments
%  nPix     number of pixels on sphere

%% Return Arguments
%  nSide    sub-divisions of HEALPix facet side

%% Example
nPix = 12*[1:4].^2;
nSide = nPix2nSide(nPix)

%% See also
% nSide2nPix

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.