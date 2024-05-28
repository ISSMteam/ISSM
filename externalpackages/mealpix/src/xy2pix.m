function ipix = xy2pix(nSide,ix,iy,nFace)
% XY2PIX gives the pixel number ipix (NESTED) of ix, iy and face_num
%
% ipix = xy2pix(nSide, ix, iy, nFace)
%
% nside     HEALPix resolution parameter (power of 2)
% ix,iy     x,y location on the Healpix face (1 <= ix, iy <= nSide)
% nFace     HEALPix face  number (1 <= nFace <= 12)
%
% ipix      MEALPix NESTED pixel index for given x,y on the face
%
% See also pix2xy
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: xy2pix.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.KeepUnmatched = false;

p.addRequired('nSide',@qNSide);
qNXY = @(n)(all(n(:)>0) && all(n(:)<=nSide));
p.addRequired('ix',qNXY);
p.addRequired('iy',qNXY);
p.addRequired('nFace',...
  @(n)(isscalar(n) && isnumeric(n) && ...
  n == fix(n) && 0 < n && n <= 12));
p.parse(nSide,ix,iy,nFace);

%% Initialization
[x2pix y2pix]=mk_xy2pix;

%% Convert to HEALPix indexing
ix = ix - 1;
iy = iy - 1;

scale = 1;
scaleFactor = 16384;
ipf = 0;
ismax = 1;
if nSide > 1048576
  ismax = 3;
end
for k = 1:ismax
  ix_low = mod(ix,128);
  iy_low = mod(iy,128);
  ipf = ipf + (x2pix(ix_low+1) + y2pix(iy_low+1))*scale;
  scale = scale*scaleFactor;
  ix = fix(ix/128);
  iy = fix(iy/128);
end
ipf = ipf + (x2pix(ix+1) + y2pix(iy+1))*scale;

ipix = ipf + nFace*nSide^2;

ipix = ipix+1; % MEALPix numbering

return