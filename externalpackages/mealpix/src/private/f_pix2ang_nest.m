function [theta, phi] = f_pix2ang_nest(nSide, ipix)

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

% $Id: f_pix2ang_nest.m 6016 2011-04-27 02:52:56Z lsf@GRAVITY.PSU.EDU $

% coordinate of the lowest corner of each face
% in unit of nSide
jrll =[ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 ]';
% in unit of nSide/2
jpll =[ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 ]';
%-----------------------------------------------------------------------
% total number of points
npix = nSide2nPix(nSide);
if ~(all(ipix>=0) && all(ipix<npix))
  msgid = [mfilename ':ndx'];
  error(msgid,msgid);
end

% initialize the array for the pixel number -> (x,y) mapping
[pix2x,pix2y] = mk_pix2xy;

npface = nSide^2;
nl4    = 4*nSide;

% find face and pixel number in the face
% face number in {0,11}
face_num = fix(ipix/npface);
% pixel number in the face {0,npface-1}
ipf = mod(ipix,npface);

fact1 = 1/(3*nSide^2);
fact2 = 2/(3*nSide);
%     finds the x,y on the face (starting from the lowest corner)
%     from the pixel number
ix = 0;
iy = 0;
scalemlv = 1;
ismax = 4;
for i=0: ismax;
  ip_low = mod(ipf,1024);
  ix = fix(ix + scalemlv*pix2x(ip_low+1));
  iy = fix(iy + scalemlv*pix2y(ip_low+1));
  scalemlv = scalemlv*32;
  ipf   = fix(ipf/1024);
end
ix = fix(ix + scalemlv*pix2x(ipf+1));
iy = fix(iy + scalemlv*pix2y(ipf+1));

%% transforms to (horizontal, vertical) coordinates
% 'vertical' in {0,2*(nSide-1)}
jrt = fix(ix + iy);
% 'horizontal' in {-nSide+1,nSide-1}
jpt = fix(ix - iy);

%% compute z coordinate on the sphere
% ring number in {1,4*nSide-1}
jr =  fix(jrll(face_num+1)*nSide - jrt - 1);

z = zeros(size(jr));
nr = zeros(size(jr));
kshift = zeros(size(jr));

% north pole region
mask = jr<nSide;
nrM = fix(jr(mask));
z(mask) = 1 - nrM.^2*fact1;
kshift(mask) = 0;
nr(mask) = nrM;

% equatorial region
mask = ~mask & jr <= 3*nSide;
nr(mask) = fix(nSide);
z(mask)  =(2*nSide-jr(mask))*fact2;
kshift(mask) = mod(jr(mask)-nSide,2);

% south pole region
mask = jr > 3*nSide;
nrM = fix(nl4 - jr(mask));
z(mask) = - 1 + nrM.^2*fact1;
kshift(mask) = 0;
nr(mask) = nrM;

%% Convert z to theta
theta = acos(z);

%% Compute phi in [0,2*pi)
% 'phi' number in the ring in {1,4*nr}
jp =fix(((jpll(face_num+1).*nr) + jpt + 1 + kshift)/2);
if(jp > nl4)
  jp = fix(jp - nl4);
end;
if(jp < 1)
  jp = fix(jp + nl4);
end;
phi = (pi/2)*(jp -(kshift+1)/2)./nr;

return