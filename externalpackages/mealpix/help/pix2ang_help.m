%% pix2ang
% Find spherical coordinate location of HEALPix pixels

%% Syntax
%  tp = pix2ang(nSide,pix,'Param1',Value1,...)

%% Input Arguments
%  nSide    HEALPix resolution parameter (positive integer power of 2)
%  pix      (optional) array of pixel values to find coordinates of
%
%  Param    Values
%  'nest'   pixels are given in nested indexing (true | {false})

%% Return Arguments
%  tp       cell array of [theta; phi] pixel locations (radians)

%% Example
% Return [theta; phi] postion angles for all pix of nSide=4, nest scheme
nSide = 4;
tp = pix2ang(nSide,'nest',true);
size(tp)
tp{1}
% Return position angles for a random assortment of pixels
pix = randi(nSide2nPix(nSide),3,4);
tp = pix2ang(nSide,pix);
size(tp)
tp{1}

%% See also
% See also ang2pix, pix2vec, vec2pix

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.