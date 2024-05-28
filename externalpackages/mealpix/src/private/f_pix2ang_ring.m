function [theta, phi] = f_pix2ang_ring(nSide, iPix)
% F_PIX2ANG_RING convert ring-indexed pixel to angular position vector
%
% [theta, phi] = f_pix2ang_ring(nSide,iPix)
%
% nSide    HEALPix resolution parameter
% iPix     HEALPix ring-indexed pixel number
% 
% Hand conversion of HEALPix distribution f90 code pix2ang_ring
%
% Author: Lee Samuel Finn
% Copyring 2011

% $Id: f_pix2ang_ring.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


npix = nSide2nPix(nSide);
if any(iPix<0) || any(iPix>npix-1)
  msgid = [mfilename ':ndx'];
  error(msgid,msgid);
end

nl2 = 2*nSide;

% points in each polar cap, =0 for nSide =1
nCap = nl2*(nSide-1);

% preallocate results
theta = zeros(size(iPix));
phi = zeros(size(iPix));

%% North Polar cap -------------
mask = iPix < nCap;
% counted from North pole
iPixM = iPix(mask);
iRingM = round(sqrt((iPixM+1)/2));
iPhiM  = fix(iPixM - 2*iRingM.*(iRingM - 1));
% fix round-off error appearing at large nSide
[iRingM, iPhiM]=correct_ring_phi(1, iRingM, iPhiM);
theta(mask) = acos(1-(iRingM/nSide).^2/3);
phi(mask)   = (pi/2)*(iPhiM + 0.5)./iRingM;

%% Equatorial region ------
mask = ~mask & (iPix < npix-nCap);
ipM  = fix(iPix(mask) - nCap);
nl4   = 4*nSide;
% counted from North pole
iRingM = fix(ipM/nl4) + nSide;
iPhiM  = mod(ipM,nl4);
% 0 if iring+nSide is odd, 1/2 otherwise
fodd  = 0.5*mod(iRingM+nSide+1,2);
theta(mask) = acos((nl2 - iRingM)/(1.5*nSide));
phi(mask)   = (pi/2)*(iPhiM + fodd)/nSide;

%% South Polar cap -----------------------------------
mask = npix-nCap <= iPix;
ipM   = fix(npix - iPix(mask));
% counted from South pole
iRingM = round(sqrt(ipM/2));
iPhiM  = fix(2*iRingM.*(iRingM + 1) - ipM);
% fix round-off error appearing at large nSide
[iRingM, iPhiM]=correct_ring_phi(-1, iRingM, iPhiM);
theta(mask) = acos((iRingM/nSide).^2/3 - 1);
phi(mask)   = (pi/2)*(iPhiM + 0.5)./iRingM;

return