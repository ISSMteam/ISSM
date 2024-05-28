function nPix = vec2pix(nSide,xyz,varargin)
% VEC2PIX Convert cartesian direction vectors to MEALPix pixel numbers
%
%  nPix = vec2pix(nSide,xyz,'Param1',Value1,...);
%
%  nSide     HEALPix resolution parameter
%  xyz       cell array of [3,1] cartesian direction vectors
%
%  Param     Value
%  'nest'    use nested indexing (true | {false})
%
%  nPix      size(xyz) pixel number array
%
% Example
%   xyz = {[1;0;0],[0;1;0],[0;0;1]};
%   nPix = vec2pix(4,xyz);
%
% See also ang2pix, pix2ang, pix2vec
%
% Requires ring2nest
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: vec2pix.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


%% Validate and parse input
p = inputParser;

p.addRequired('nSide', ...
  @(x)(isscalar(x) && isnumeric(x) && x>0 && nextpow2(x) == log2(x)));
p.addRequired('xyz',@iscell);
p.addParamValue('nest',false,@(n)(isscalar(n) && islogical(n)));

p.parse(nSide,xyz,varargin{:});
nest = p.Results.nest;

%% Convert vector(s) to pixel numbers
% Process xyz as a vector then restore shape

% Save xyz shape
sXyz = size(xyz);

% Insure cell vectors are [3,1]
xyz = cellfun(@(v)(reshape(v,[],1)),xyz,'UniformOutput',false);

% Re-express xyz as a column vector
xyz = cell2mat(reshape(xyz,1,[]));

% Find pixels
nPix = f_vec2pix_ring(nSide,xyz(1,:),xyz(2,:),xyz(3,:));
nPix = nPix + 1;
if nest
  nPix = ring2nest(nSide,nPix);
end

% Restore shape
nPix = reshape(nPix,sXyz);

return

function  iPix = f_vec2pix_ring(nSide, x, y, zin)

%=======================================================================
%     renders the pixel number iPix (RING scheme) for a pixel containing
%     a point on a sphere at coordinate vector (=x,y,z), given the map
%     resolution parameter nside
%=======================================================================

dnorm = arrayfun(@(x,y,z)(norm([x,y,z])),x,y,zin);
z = zin./dnorm;

phi = atan2(y,x);
za = abs(z);
% phi(phi < 0.0) = phi(phi < 0.0) + 2*pi;
phi = mod(phi,2*pi);
tt = phi/(pi/2);
nl4 = fix(4*nSide);

%% preallocate result
iPix = zeros(size(x));

%% Equatorial region
% Identify pixels
pixMask = 3*za <= 2;
if ~isempty(pixMask)
  ttPix = tt(pixMask);
  zPix = z(pixMask);
  
  temp1 = nSide*(ttPix + 0.5);
  temp2 = nSide*zPix*0.75;
  
  % find ascending, descending line indices
  jp = fix(temp1-temp2);
  jm = fix(temp1+temp2);
  
  ir = fix(nSide+jp-jm);
  kShift = fix(bitand(ir,1));
  ip = fix((jp+jm-nSide+kShift+1)/2);
  ip(ip>=nl4) = fix(ip(ip>=nl4)-nl4); % can this be replaced by a mod?
  iPix(pixMask) = fix(2*nSide*(nSide-1) + nl4*fix(ir) + ip);
end
%% North, south polar cap
% Identify pixels
pixMask = ~pixMask;
if ~isempty(pixMask)
  ttPix = tt(pixMask);
  tp = ttPix - fix(ttPix);
  tmp = nSide*sqrt(3*(1-za(pixMask)));
  
  % find increasing, decreasing edge line indices
  jp = fix(tp.*tmp);
  jm = fix((1-tp).*tmp);
  
  % ring number counted from the closest pole
  ir = fix(jp+jm+1);
  ip = fix(ttPix.*ir);
  tMask = ip>=4*ir;
  ip(tMask) = fix(ip(tMask)-4*ir(tMask));
  
  % north pole pixels
  nPixMask = pixMask & z > 0;
  tMask = z(pixMask) > 0;
  iPix(nPixMask) = fix(2*ir(tMask).*(ir(tMask)-1) + ip(tMask));
  
  % south pole pixels
  sPixMask = pixMask & z <= 0;
  tMask = ~tMask;
  iPix(sPixMask) = fix(...
    3*nSide*fix(nl4) - 2*ir(tMask).*(ir(tMask)+1)+ip(tMask)...
    );
end

return