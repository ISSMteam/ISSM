%% nest2ring
% Convert MEALPix pixel numbers from nest to ring indexing

%% Syntax
%  rPix = nest2ring(nSide,nPix)

%% Input Arguments
%  nSide     HEALPix resolution parameters
%  nPix      numeric array of MEALPix nest indexed pixel numbers

%% Return Arguments
%  rPix      size(nPix) numeric array of ring indexed pixel numbers

%% Example
nPix = 1:10;
rPix = nest2ring(2,nPix)
nPix = ring2nest(2,rPix)

%% See also
% ring2nest

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.