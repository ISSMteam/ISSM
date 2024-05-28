%% inRing 
% Return list of pixels in a ring or ring section

%% Syntax
%  pixList = inRing(nSide, nRing, phi0, dPhi, 'Param1', Value1, ...);

%% Input Arguments
%  nSide       HEALPix resolution parameter
%  nRing       Ring number (1 <= nRing <= 4*nSide-1)
%  phi0, dPhi  (optional) Ring section longitudes [phi-dPhi,phi+dPhi]
%
%  Param       Value
%  'nest'      return pixels with nested indexing (true | {false})
%  'strict'    return only pixels whose center is in longitude range
%              (true | {false})

%% Return Arguments
%  pixList     size(nRing) cell array of pixels in ring defined by nRing,
%              phi0, dPhi

%% Description
%  nRing may be a numeric array, in which case phi0, dPhi and nest may each
%  be either scalar or a size(nRing) array. If strict is false (default)
%  then all pixels intersected by longitude range are included in pixList.

%% Example

% All pixels in ring 2 of nSide = 4 pixelization (ring indexed)
pix = inRing(4,2);
pix{:}

% Same, but nested indexing
pix = inRing(4,2,'nest',true);
pix{:}
nest2ring(4,pix{:})

% Pixels in rings [2,4;8,10;12,14] and the longitude band
% [7*pi/4,9*pi/4] 
nRing = [2,3;5,12;13,14];
pix = inRing(4,nRing,2*pi,pi/4);
pix{1}
pix{end}

% Same, with strict selection
pix = inRing(4,nRing,2*pi,pi/4,'strict',true);
pix{1}
pix{end}


%% See also 
% ringNum

%% Requires 
% ring2z, pix2vec

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.