function ipring = f_nest2ring(nSide, ipnest)

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

% $Id: f_nest2ring.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

nPix = nSide2nPix(nSide);
if ~(all(ipnest>=0) && all(ipnest<nPix))
  msgid = [mfilename ':index'];
  error(msgid,msgid);
end

% initializations
jrll = [ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 ]';
jpll = [ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 ]';
[pix2x,pix2y] = mk_pix2xy;

npface = nSide^2;  % number of pixels in each face
nl4    = 4*nSide;

%% find the face and the number in the face
% face number in {0,11}
face_num = fix(ipnest/npface);
ipf = mod(ipnest,npface);

%     finds the x,y on the face (starting from the lowest corner)
%     from the pixel number
ix = 0;
iy = 0;
scalemlv = 1;
ismax = 4;

for i = 0:ismax;
  ip_low = mod(ipf,1024);
  ix = fix(ix + scalemlv*pix2x(ip_low+1));
  iy = fix(iy + scalemlv*pix2y(ip_low+1));
  scalemlv = scalemlv*32;
  ipf   = fix(ipf/1024);
end; 
ix = fix(ix + scalemlv*pix2x(ipf+1));
iy = fix(iy + scalemlv*pix2y(ipf+1));

%% transforms to (horizontal, vertical) coordinates
% 'vertical' in {0,2*(nSide-1)}
jrt = fix(ix + iy);
% 'horizontal' in {-nSide+1,nSide-1}
jpt = fix(ix - iy);

%% Find z coordinate on the sphere
% ring number in {1,4*nSide-1}
jr =  fix(jrll(face_num+1)*nSide - jrt - 1);

% initialization
nr = zeros(size(ipnest));
kshift = zeros(size(ipnest));
n_before = zeros(size(ipnest));

% north pole region
mask = jr < nSide;
jrN = jr(mask);
nrN = fix(jrN);
n_before(mask) = fix(2*nrN.*(nrN - 1));
kshift(mask) = 0;
nr(mask) = nrN;
  
% equatorial region (the most frequent)
mask = jr >= nSide & (jr <= 3*nSide);
jrE = jr(mask);
nrE = fix(nSide);
n_before(mask) = fix(2*nrE*(2*jrE - nrE - 1));
kshift(mask) = mod(jrE - nSide, 2);
nr(mask) = nrE;

% south pole region
mask = jr > 3*nSide;
jrS = jr(mask);
nrS = fix(nl4 - jrS);
n_before(mask) = fix(nPix - 2*nrS.*(nrS + 1));
kshift(mask) = 0;
nr(mask) = nrS;

% computes the phi coordinate on the sphere, in [0,2Pi]
% 'phi' number in the ring in {1,4*nr}
jp = fix((jpll(face_num+1).*nr + jpt + 1 + kshift)/2);
maskH = jp>nl4;
maskL = jp<1;
jp(maskH) = fix(jp(maskH) - nl4);
jp(maskL) = fix(jp(maskL) + nl4);

% in {0, nPix-1}
ipring = fix(n_before + jp - 1);

return

