function nRing = z2ring(nSide,z)
%% Z2RING Find HEALPix ring numbers from direction vector z coordinate
%
%  nRing = z2ring(nSide,z)
%
%  nSide      HEALPix resolution parameter
%  z          numeric array of z values
%
%  nRing      size(z) numeric array of ring numbers
%
% Description
% HEALPix pixels centers are arranged in 4*nSide-1 rings of constant z.
% Each ring spans a given range of z. z2ring returns the ring numbers that
% cover the input z values. 
% 
% Example:
% nRing = reshape(1:15,3,5);
% z = ring2z(4,nRing);
% 
% See also ring2z
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: z2ring.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.FunctionName = mfilename;
p.addRequired('nSide',@qNSide);
p.addRequired('z',@(x)(all(abs(x(:))<=1)));
p.parse(nSide,z);

%% Work
% Save input size
% Flatten input
% Find ring numbers
% Reshape ring numbers to reflect input

zSz = size(z);
nRing = reshape(ringNum(nSide, z(:)),zSz);

return