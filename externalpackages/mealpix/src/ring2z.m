function z = ring2z(nSide,nRing)
% RING2Z Find HEALPix ring number corresponding to a given z coordinate
%
%  z = ring2z(nSide,nRing)
%
%  nSide    HEALPix base-pixel sub-divisions
%  nRing    numeric array of ring numbers
%
%  z        size(nRing) numeric array of pixel center z coordinates
%
% Example
% z = ring2z(4,reshape(1:15,3,5));
%
% See also z2ring
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: ring2z.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.addRequired('nRing',@(n)(qNRing(nSide,nRing))); 
p.parse(nSide,nRing);

%% pre-allocate result
z = zeros(size(nRing));

%% north polar cap
mask = nRing < nSide;
z(mask) = 1 - nRing(mask).^2/(3*nSide^2);

%% tropical band
mask = ~mask & nRing < 3*nSide; 
z(mask) = (2*nSide-nRing(mask))*2/(3*nSide);

%% south polar cap
mask = nRing >= 3*nSide;
z(mask) = -1 + (4*nSide-nRing(mask)).^2/(3*nSide^2);

return

function tf = qNRing(nSide,nRing)
% QNRING true if valid ring number

tf = isnumeric(nRing);
if tf
  nRing = nRing(:);
  tf = ...
    all(nRing == fix(abs(nRing))) && ...
    all(nRing > 0) && ...
    all(nRing < 4*nSide);
end

return