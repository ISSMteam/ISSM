function rad = angDist(varargin)
%ANGDIST Computes angular distance on the unit sphere
%
%   rad = angDist(v0,v1,'Param1',Value1,'Param2',Value2,...)
%
%   v0, v1      (numerical or cell) array of points on sphere
%
%   Param       Value
%   'nSide'     HEALPix resolution parameter (integer power 2)
%   'nest'      nested indexing (true | {false})
%
%   rad         array of angular distances (radians)
%
% Description
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
% If any location is specified as a pixel number then nSide must be
% specified. The indexing scheme for pixels center locations defaults to
% ring; to specify nested indexing set the 'nest' parameter to true. 
%
% v0 or v1 specified as numeric arrays:
% If all points in v0 are pixel numbers then v0 may be specified as a
% numeric array (similarly v1, mutatis mundis). nSide must be specified.
% Pixels numbers may be specified in ring indexing or nested indexing. For
% nested indexing the nest parameter must be set to true. 
%
% pix0 or pix1 specified as cell arrays:
% Either or both pix0 or pix1 may be specified as a cell array, in which
% case each cell specifies a location on the sphere as either a pixel
% number, a two-component ([theta; phi]) spherical coordinate vector, or a
% three component ([x;y;z]) cartesian vector. nSide must be specified if
% any location is identified as a pixel center. 
%
% Example
%  % angular distance from pole to equator:
%  v0 = {[0;0;1]}; % north pole
%  v1 = {[1;0;0]}; % point on equator
%  rad = angDist(v0, v1);
%  
%  % angular distance from pole to pixels 1, 17, 33 for nSide = 4:
%  v0 = {[0;0;1]};
%  rad = angDist([1, 17, 33], v0, 'nSide', 4);
% 
%  % angular distance from pixels 1 to 4, 2 to 5, and 3 to 6 for nSide = 8:
%  rad = angDist(8, [1, 2, 3], [4, 5, 6]);
%
% Requires pix2vec, ang2vec
%
% See also queryDisc
% 
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: angDist.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

% Copyright 2010-2011 by Lee Samuel Finn. 
% This program file is part of MEALPix. It is licensed under the Apache
% License, Version 2.0 (the  "License"); you may not use MEALPix except in
% compliance with the License. You may obtain a copy of the License at 
% <http://www.apache.org/licenses/LICENSE-2.0> 
%
% Unless required by applicable law or agreed to in writing, MEALPix
% software distributed under the License is distributed on an "AS IS"
% BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
% implied. See the License for the specific language governing permissions
% and limitations under the License.


%% Verify and parse arguments
error(nargoutchk(1,1,nargout,'struct'));

p = inputParser;
p.addRequired('v0',@(x)(iscell(x) || isnumeric(x)));
p.addRequired('v1',@(x)(iscell(x) || isnumeric(x)));
p.addParamValue('nest',false,@islogical);
% nSide scalar non-negative integer
p.addParamValue('nSide',0,...
  @(x)(isscalar(x) && isnumeric(x) && floor(abs(x)) == x));
p.parse(varargin{:});

v0 = p.Results.v0;
v1 = p.Results.v1;
nest = p.Results.nest;
nSide = p.Results.nSide;

if ~(1 == numel(v0) || 1 == numel(v1))
  nV0 = size(v0);
  nV1 = size(v1);
  assert(all(nV0 == nV1),[mfilename ':vDim']);
end

%%
if isnumeric(v0)
  v0 = pix2vec(nSide,v0);
else
  v0 = cellfun(@(x)(v2vec(x,nSide,nest)),v0,'UniformOutput',false);
end
if isnumeric(v1)
  v1 = pix2vec(nSide,v1);
else
  v1 = cellfun(@(x)(v2vec(x,nSide,nest)),v1,'UniformOutput',false);
end

if 1 == numel(v0)
  rad = cellfun(@(x)(f_angdist(v0{1},x)),v1);
elseif 1 == numel(v1)
  rad = cellfun(@(x)(f_angdist(x,v1{1})),v0);
else
  rad = cellfun(@(x,y)(f_angdist(x,y)),v0,v1);
end

return

function  dist = f_angdist(v1, v2)
% use pix_tools
%=======================================================================
% call angdist(v1, v2, dist)
% computes the angular distance dist (in rad) between 2 vectors v1 and v2
% in general dist = acos ( v1 . v2 )
% except if the 2 vectors are almost aligned.
%=======================================================================
% persistent diff r1 r2 sprod vdiff ;

% if isempty(r1), r1=zeros(1,3); end;
% if isempty(r2), r2=zeros(1,3); end;
% if isempty(vdiff), vdiff=zeros(1,3); end;
% if isempty(diff), diff=0; end;
% if isempty(sprod), sprod=0; end;
%=======================================================================
% normalize both vectors
r1 = v1/norm(v1);
r2 = v2/norm(v2);
sprod = dot(r1, r2);
if(sprod > 0.999)
  % almost colinear vectors
  vdiff = r1 - r2;
  % norm of difference
  diff = sqrt(dot(vdiff,vdiff));
  dist = 2*asin(diff/2);
elseif(sprod < -0.999) ;
  % almost anti-colinear vectors
  vdiff = r1 + r2;
  % norm of sum
  diff = sqrt(dot(vdiff,vdiff));
  dist = pi - 2*asin(diff/2);
else
  % other cases
  dist = acos( sprod );
end;
return

function v = v2vec(v,nSide,nest)
% V2VEC pixel, angular locations or vectors to vectors

switch numel(v)
  case 1, % pixel
    v = pix2vec(nSide,v,'nest',nest);
  case 2, % angular coordinate location
    z = cos(v(1));
    s = sin(v(1));
    v = [s*cos(v(2));s*sin(v(2));z];
  case 3,
    % nop
  otherwise,
    msgid = [mfilename ':v2vec'];
    error(msgid,msgid);
end
v = reshape(v,[],1);

return