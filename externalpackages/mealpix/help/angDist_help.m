%% angDist 
% Compute angular distance between pixel centers or points on the sphere 

%% Syntax
%   rad = angDist(v0,v1,'Param1',Value1,'Param2',Value2,...)

%% Input Arguments
%  v0, v1      (numerical or cell) array of points on sphere
%
%  Param       Value
%  'nSide'     HEALPix resolution parameter (integer power 2)
%  'nest'      nested indexing (true | {false})

%% Return Arguments
%   rad     array of angular distances (radians)

%% Description
% angDist computes the angular distances between points on the sphere.
% Points may be specified as MEALPix pixel numbers, spherical coordinate
% location vectors, or cartesian coordinate vectors. Cartesian coordinate
% input vectors need not be normalized. 
%
% Either or both of v0, v1  may be a scalar; otherwise they must both be
% the same size. When one is a scalar than the distance is computed between
% that point and all the points specified in the complementary argument.
% When pix0 and pix1 are arrays then the distances are returned between the
% corresponding points in the two arrays.
%
% If any location is specified as a pixel number than nSide must be
% specified. The indexing scheme for pixels center locations defaults to
% ring; to specify nested indexing set the 'nest' parameter to true. 

%% v0 or v1 specified as numeric arrays
% If all points in v0 are pixel numbers then v0 may be specified as a
% numeric array (similarly v1, mutatis mundis). then nSide must be specified.
% Pixels numbers may be specified in ring indexing or nested indexing. For
% nested indexing the nest flag must be specified as true. 

%% pix0 or pix1 specified as cell arrays
% Either or both pix0 or pix1 may be specified as a cell array, in which
% case each cell specifies a location on the sphere as either a pixel
% number, a two-component ([theta; phi]) spherical coordinate vector, or a
% three component ([x;y;z]) cartesian vector. nSide must be specified if
% any location is identified as a pixel center. 

%% Example
 % angular distance from pole to equator:
 v0 = {[0;0;1]}; % north pole
 v1 = {[1;0;0]}; % point on equator
 rad = angDist(v0, v1);
 
 % angular distance from pole to pixels 1, 17, 33 for nSide = 4:
 v0 = {[0;0;1]};
 rad = angDist([1, 17, 33], v0, 'nSide', 4)

 % angular distance from pixels 1 to 4, 2 to 5, and 3 to 6 for nSide = 8:
 rad = angDist([1, 2, 3], [4, 5, 6],'nSide',8)

%% See also
% queryDisc

%% Requires
% pix2vec, ang2vec

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.