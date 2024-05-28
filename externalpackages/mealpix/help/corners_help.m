%% corners
% Find pixel corners

%% Syntax
%  xyz = corners(nSide,nPix,'Param1',Value1,...);

%% Input Arguments
%  nSide    HEALPix resolution parameter
%  nPix     (optional) pixel list. Default all pixels.
%  
%  Param    Value
%  'nest'   nested indexing (true | {false})

%% Return Arguments
%  xyz    size(nPix) cell array of [3,1] cartesian vector pixel corner
%         locations

%% Example:
nSide = 2^2;
% Corners for all 48 nSide=4 pixels in ring-indexed order
xyzR = corners(nSide);
size(xyzR)
xyzR{1} 
xyzR{end-1}

nPix = randi(nSide2nPix(nSide), 3, 4);
% Corners for an assortment of pixels in nested index order
xyzN = corners(nSide,nPix,'nest',true);
% Convert to pixel list to ring indexing
r = xyzR{nest2ring(nSide,nPix(2,3))}
% Corners for same pixels requested via ring indexing
n = xyzN{2,3}
% Compare
all(r == n)

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.