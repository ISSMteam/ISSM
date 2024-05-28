function nRing = ringNum(nSide, z, varargin) 
%% RINGNUM Find ring number for given z coordinate
%
% ring = ringNum(nSide, z, shift)
%
% nSide    HEALPix resolution parameter (integer power 2)
% z        numeric array of z coordinates on unit sphere
% shift    (optional) ring number "rounding" mode
%          shift = 0: return ring number whose center is closest to z
%          shift < 0: return ring number immediately north of closest ring
%          shift > 0: return ring number immediately south of closest ring
%
% nRing    size(z) numeric array of ring numbers
%
% Author: Lee Samuel Finn
% Copyright: 2011

% $Id: ringNum.m 5792 2011-02-20 03:21:25Z lsf@GRAVITY.PSU.EDU $

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


%% Parse and validate input
p = inputParser;
p.addRequired('nSide',@qNSide);
p.addRequired('z',@(x)(all(abs(x(:))<=1)));
p.addOptional('shift',0,@isnumeric);
p.parse(nSide,z,varargin{:});
shift = p.Results.shift/2;
  
%% preallocate 
nRing = zeros(size(z));

%% equatorial regime
mask = abs(z) <= 2/3;
nRing(mask) = round(nSide*(2-1.5*z(mask)) + shift);

%% north polar cap
mask = z > 2/3;
nTmp = round(nSide*sqrt(3*(1-z(mask))) + shift);
nTmp(0 == nTmp) = 1;
nRing(mask) = nTmp;

%% south polar cap
mask = z < -2/3;
nTmp = round(nSide*sqrt(3*(1+z(mask))) - shift); % Yes: this is a -shift
nTmp(0 == nTmp) = 1;
nRing(mask) = 4*nSide - nTmp;

return