%% vec2pix
% Convert cartesian direction vectors to MEALPix pixel numbers

%% Syntax
%  nPix = vec2pix(nSide,xyz,'Param1',Value1,...);

%% Input Arguments
%  nSide     HEALPix resolution parameter
%  xyz       cell array of [3,1] cartesian direction vectors
%  
%  Param     Value
%  'nest'    use nested indexing (true | {false})

%% Return Arguments
%  nPix      size(xyz) pixel number array

%% Example
% Find pixels corresponding to a [4,5] set of (unnormalized) cartesian
% vectors:
nSide = 4;
xyz0 = squeeze(num2cell(2*rand(3,4,5)-1,1));
nPix = vec2pix(nSide,xyz0)

% Find the normalized vectors corresponding to the pixels
xyz1 = pix2vec(nSide,nPix);

% Find angular distance from vectors to pixels
d = angDist(xyz0,xyz1)

%% See also
% ang2pix, pix2ang, pix2vec

%% Requires
% ring2nest

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.