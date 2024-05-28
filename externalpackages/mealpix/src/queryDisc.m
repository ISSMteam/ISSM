function pix = queryDisc( nSide, c0, r0, varargin )
% QUERYDISC Find all pixels whose centers are within a specified disc
%
%
%  pix = queryDisc(nSide,c0,r,'Param1',Value1,...)
%
%  nSide
%  c0        disc center. May be specified as a pixel number, anuglar
%            position vector [theta, phi] in rads, or as a cartesian unit
%            vector [x;y;z]. 
%  r         disc radius (rads)
%
%  Param     Value
%  'nest'    pixel indexing is nested (true | {false})
%
%  pix       list of pixels in disc
%
% Example
% nSide=8;
% pix = queryDisc(nSide,1,pi/10,'nest');
% numel(pix)
% nSide2nPix(nSide)*((pi/10)/(4 Pi))
%
% See also queryPolygon, queryStrip, queryTriangle
%
% Requires pix2vec, ang2vec, angDist
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: queryDisc.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

%% validate and parse arguments
p = inputParser;
p.addRequired('nSide',@qNSide);
p.addRequired('c0',@(x)(isnumeric(x) && any(numel(x) == 1:3)));
p.addRequired('r0',@(x)(isnumeric(x) && isscalar(x) && x > 0 && x <= pi));
p.addParamValue('nest',false,@(n)(isscalar(n) && islogical(n)));
p.parse(nSide,c0,r0,varargin{:});
nest = p.Results.nest;

%% Setup vIn
switch numel(c0)
  case 1 %pix
    vc = cell2mat(pix2vec(nSide,c0,'nest',nest));
  case 2 %ang
    vc = cell2mat(ang2vec(c0));
  case 3 %vec
    vc = c0/norm(c0);
end

%% Computation
% Find max, min z of disk
% Find corresponding ring numbers
% Find pixels in rings
% Select pixels within ring

th = acos(vc(3));
zMax = cos(max(th-r0,0));
zMin = cos(min(th+r0,pi));

ringMin = z2ring(nSide,zMax);
ringMax = z2ring(nSide,zMin);

rings = min(ringMin,ringMax):max(ringMin,ringMax);
pix = arrayfun(@(r)(cell2mat(inRing(nSide,r))),rings,'UniformOutput',false);
pix = cell2mat(pix);

d = angDist(pix,{vc},'nSide',nSide);
pix = pix(d < r0); % closer than r0
if nest
  pix = ring2nest(nSide,pix);
end
return