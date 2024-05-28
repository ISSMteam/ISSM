function nSide = res2nSide(r)
% RES2NSIDE Finds minimum base pixel sub-division for required resolution
%
%  nSide = res2nSide(r)
%
%  r        angular resolution (arc seconds)
%
%  nSide    minimum base pixel sub-divisions corresponding to resolution r
%
% See also nSide2res
% 
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: res2nSide.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


nPix0 = 4*pi./(r*pi/(180*3600)).^2;
nSide0 = realsqrt(nPix0/12);
nSide = pow2(ceil(log2(nSide0)-2*eps));
nSide(nSide<1) = 1;

return