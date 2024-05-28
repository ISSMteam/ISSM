%% pix2vec
% Convert MEALPix pixel numbers to cartesian location vectors

%% Syntax
%  [x,y,z] = pix2vec(nSide,pix,'Param1',Value1,...)
%  xyz = pix2vec(nSide,pix,'Param1',Value1,...)

%% Input Arguments
%  nSide       HEALPix resolution parameter
%  pix         (optional) numeric array of pixel numbers
%
%  Param       Value
%  'nest'      indexing scheme is nested (true | {false}) 

%% Return Arguments
%  x,y,z       size(pix) size(pix) numeric array of unit vector coordinate
%              components corresponding to pixels pix
%  xyz         size(pix) size(pix) cell array of cartesian unit vectors
%              corresponding to pixels pix

%% Example
nSide = 2^4;
nPix = nSide2nPix(nSide);
pix0 = reshape(1:6,3,2);
pix1 = reshape(nPix+(-5:0),3,2);
pix = [pix0,pix1];
xyz = pix2vec(nSide,pix);
xyzn = pix2vec(nSide,pix,'nest',true);
size(xyz)
[xyz{1},xyz{end}]
[xyzn{1},xyzn{end}]

%% See also
% ang2pix, pix2ang, vec2pix

%% Requires
% nest2ring

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.