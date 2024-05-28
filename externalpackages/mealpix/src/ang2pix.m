function pix = ang2pix(varargin)
% ANG2PIX Convert spherical coordinate location to HEALPix pixel number
%
% pix = ang2pix(nSide, v);
% pix = ang2pix(nSide, v, nest);
%
%  nSide      HEALPix resolution parameter
%  v          (theta, phi) pairs in radians as either size(v,1) == 2 or as
%             cellarray of [2,1] 
%  nest       (optional) nested pixel number scheme (true | {false})
%
%  pix        array of pix numbers:
%             v on input     output
%             iscell(v)      size(v) cell array of pixel numbers
%             isnumeric(v)   numeric array of pixel numbers
%
% Example
% vTheta = acos(2*rand(3,5)-1);
% vPhi = 2*pi*rand(3,5);
% v = arrayfun(@(x,y)([x,y]),vTheta,vPhi,'UniformOutput',false);
% nSide = 8;
% pix = ang2pix(nSide,v);
%
% NB: Matlab indexes pixels from 1
%
% See also vec2pix, pix2ang, pix2vec
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: ang2pix.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

%% Parse and validate arguments
p = inputParser;
p.addRequired('nSide',@qNSide);
p.addRequired('v',@qV);
p.addOptional('nest',false,@islogical)
p.parse(varargin{:});
nSide = p.Results.nSide;
v = p.Results.v;
nest = p.Results.nest;

%% unpack v
if isnumeric(v)
  ns = size(v);
  ns = ns(2:end);
  vt = reshape(v(1,:),ns);
  vp = reshape(v(2,:),ns);
else
  vt = cellfun(@(x)(x(1)),v);
  vp = cellfun(@(x)(x(2)),v);
end

%% Compute pix numbers
pix = f_ang2pix_ring(nSide,vt,vp);

%% If required convert to nest indexing
if nest
  pix = ring2nest(nSide,pix);
end

return

%% qV
function tf = qV(v)
% QV test v for inputParser

if isnumeric(v)
  tf = 2 == size(v,1);
elseif iscell(v)
  nSize = cellfun('prodofsize',v);
  tf = all(2 == nSize(:));
else
  tf = false;
end

return

%% f_ang2pix_ring
function ipix = f_ang2pix_ring(nside, theta, phi)

%=======================================================================
% renders the pixel number ipix (RING scheme) for a pixel which contains
% a point on a sphere at coordinates theta and phi, given the map
% resolution parameter nside
%=======================================================================

%% f90 I/O Check/Init
if(nside<1 || nside>8192)
  msgid = [mfilename ':nSide'];
  error(msgid,'%s nSide (%d) out of range',msgid,nside);
end;
if(theta<0 | theta>pi) %#ok<OR2>
  msgid = [mfilename ':theta'];
  error(msgid,'%s %e is out of range [0,pi)',msgid,theta);
end;

z = cos(theta);
za = abs(z);
% in [0,4)
tt = mod(phi,2*pi)/(pi/2);

ip = zeros(size(z));
ipix = zeros(size(z));

% Logical arrays of pix regions (to remove if statement and vectorize calc)
emask = za <= 2/3;
pmask = ~emask;

%% Equatorial Region Calculation
temp1(emask) = nside*(0.5+tt(emask));
temp2(emask) = nside*0.75*z(emask);

% index of  ascending edge line
jp(emask) = fix(temp1(emask)-temp2(emask));

% index of descending edge line
jm(emask) = fix(temp1(emask)+temp2(emask));

% in {1,2n+1} (ring number counted from z=2/3)
ir = zeros(size(emask));
ir(emask) = fix(nside + 1 + jp(emask) - jm(emask));

% kshift=1 if ir even, 0 otherwise
kshift(emask) = fix(1 - mod(ir(emask),2));

%% What is this all about?
nl4 = 4*nside;
% in {0,4n-1}
ip(emask) = fix(fix((jp(emask)+jm(emask)-nside+kshift(emask)+1)/2));
% if(ip(q) >= nl4)
ip(emask & (ip >= nl4)) = fix(ip(emask & (ip >= nl4)) - nl4);
% end;
ipix(emask) = fix(2*nside*(nside-1) + nl4*(ir(emask)-1) + ip(emask));

%% NorthSouth Polar Cap Calculation
%else;
%tt=mod(tt,1);
tp(pmask) = tt(pmask) - fix(tt(pmask));
tmp(pmask) = nside*sqrt(3*(1-za(pmask)));
% increasing edge line index
jp(pmask) = fix(tp(pmask).*tmp(pmask));
% decreasing edge line index
jm(pmask) = fix((1-tp(pmask)).*tmp(pmask));
% ring number counted from the closest pole
ir = zeros(size(pmask));
ir(pmask) = fix(jp(pmask) + jm(pmask) + 1);
% in {0,4*ir-1}
ip(pmask) = fix(fix(tt(pmask).*ir(pmask)));
newmask1 = (pmask&(ip >= 4*ir));
%if(ip(pmask) >= 4.*ir(pmask))
ip(newmask1) = fix(ip(newmask1) - 4*ir(newmask1));
%end;
newmask2 = (pmask&(z>0));
%if(z(pmask)>0.)
ipix(newmask2) = fix(2*ir(newmask2).*(ir(newmask2)-1) + ip(newmask2));
newmask3 = (pmask&~(z>0));
%else;
nPix = 12*nside^2;
ipix(newmask3) = fix(nPix-2*ir(newmask3).*(ir(newmask3)+1)+ip(newmask3)); 
%end;

ipix = ipix+1; %MEALPix numbering
return