%% pix2xy
% Find pixel cartesian face coordinates

%% Syntax
%  xy = pix2xy(nSide,nPix)

%% Input Arguments
%  nSide     HEALPix resolution parameter
%  nPix      nested indexing pixel numbers to find coordinates of

%% Return Arguments
%  xy        size(nPix) cell array of pixel coordinates on face

%% Example
nPix = reshape(1:24,[4,3,2]);
xy = pix2xy(2,nPix)

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.