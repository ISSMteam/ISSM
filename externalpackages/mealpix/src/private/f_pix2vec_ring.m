function varargout = f_pix2vec_ring(nSide, hPix)
% F_PIX2VEC_RING conversion of f90 PIX2VEC_RING HEALPix code
%
% [vector, vertex] = f_pix2vec_ring(nSide,hPix)
%
% nSide    (integer) sub-divisions of HEALPix facets
% hPix     (vector integers) HEALPix pixel numbers 
%          (0 <= hPix < 12*nSide^2)
%
% vector   cell(size(hPix)) of [3,1] unit vectors to pixel centers
% vertex   (optional) cell(size(hPix)) of [3,4] unit vectors to pixel
%          vertices (corners)
% 
%     renders vector (x,y,z) coordinates of the nominal pixel center
%     for the pixel number hPix (RING scheme)
%     given the map resolution parameter nSide
%     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
%     in the order N,W,S,E
%
% Hand rewrite of f90 code
%
% Author: Lee Samuel Finn
% Copyright 2010

% $Id: f_pix2vec_ring.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

%% I/O Check/Init
error(nargoutchk(1,2,nargout,'struct'));
error(nargchk(2,2,nargin,'struct'));

nPix = nSide2nPix(nSide);
if ~(...
    all(floor(hPix(:)) == abs(hPix(:))) && all(hPix(:) < nPix))
  msgid = [mfilename ':hPix'];
  error(msgid,msgid);
end

% in {1, nPix}
hPix1 = fix(hPix + 1);
nl2 = fix(2*nSide);
nl4 = fix(4*nSide);
% points in each polar cap, = 0 for nSide = 1
nCap = fix(2*nSide*(nSide-1));
fact1 = 1.5*nSide;
fact2 = 3.0*nSide^2;

% Calculate corners (vertices)?
do_vertex = (nargout == 2);

% Associate pixels with regions (north cap, equatorial, south cap)
ncMask = hPix1 <= nCap;
eqMask = hPix1 <= nl2*(5*nSide+1) & ~ncMask;
scMask = ~(eqMask | ncMask);

% Pre-allocate 
z   = zeros(size(hPix));
phi = zeros(size(hPix));
hDelta_Phi = zeros(size(hPix));
z_nv = zeros(size(hPix));
z_sv = zeros(size(hPix));
phi_nv = zeros(size(hPix));
phi_sv = zeros(size(hPix));

% North polar cap
if any(ncMask(:))
  hip = hPix1(ncMask)/2;
  fihip = fix(hip);
  iRing = fix(sqrt(hip-sqrt(fihip)))+1;
  iPhi = hPix1(ncMask) - 2*iRing.*(iRing-1);
  
  z(ncMask) = 1-iRing.^2/fact2;
  phi(ncMask) = (iPhi-0.5).*pi./(2*iRing);
  
  if (do_vertex)
    hDelta_Phi(ncMask) = pi./(4*iRing);     % half pixel width
    z_nv(ncMask) = 1-(iRing-1).^2/fact2;
    z_sv(ncMask) = 1-(iRing+1).^2/fact2;
    iPhi_Mod = mod(iPhi-1,iRing);    % in {0..iRing-1}
    iPhi_Rat = fix((iPhi-1)./iRing);      % in {0,1,2,3}
    phi_nv_m = phi_nv(ncMask);
    ndx = (1 < iRing);
    phi_nv_m(ndx) = (iPhi_Rat(ndx)+iPhi_Mod(ndx)./(iRing(ndx)-1))*pi/2;
    phi_nv(ncMask) = phi_nv_m;
    phi_sv(ncMask) = (iPhi_Rat+(iPhi_Mod+1)./(iRing+1))*pi/2;
  end
end

% Equatorial region
if any(eqMask(:))
  ip = hPix1(eqMask) - nCap -1;
  iRing = fix(ip/nl4) + nSide;      % counted from north pole
  iPhi = mod(ip,nl4)+1;
  
  fOdd = 0.5*(1+mod(iRing+nSide,2)); % 1 (1/2) if iRing+nSide is odd (even)
  z(eqMask) = (nl2-iRing)/fact1;
  phi(eqMask) = (iPhi-fOdd)*pi/(2*nSide);
  
  if (do_vertex)
    hDelta_Phi(eqMask) = pi./(4*nSide);     % half pixel width
    phi_nv_m = phi(eqMask);
    phi_sv_m = phi(eqMask);
    z_nv_m = (nl2-iRing+1)/fact1;
    z_sv_m = (nl2-iRing-1)/fact1;
    nNdx = (iRing == nSide);   % northern transition
    sNdx = (iRing == 3*nSide); % southern transition
    ndx = sNdx|nNdx;           % transition
    z_nv_m(nNdx) =  1-(nSide-1)^2/fact2;
    z_sv_m(sNdx) = -1+(nSide-1)^2/fact2;
    iPhi_Mod(ndx) = mod(iPhi(ndx)-1,nSide);   % in {0..nSide-1}
    iPhi_Rat(ndx) = fix((iPhi(ndx)-1)/nSide);      % in {0..3} 
    if nSide > 1
      phi_nv_m(nNdx) = (iPhi_Rat(nNdx)+iPhi_Mod(nNdx)/(nSide-1))*pi/2;
      phi_sv_m(sNdx) = (iPhi_Rat(sNdx)+iPhi_Mod(sNdx)/(nSide-1))*pi/2;
    end
    phi_nv(eqMask) = phi_nv_m;
    phi_sv(eqMask) = phi_sv_m;
    z_nv(eqMask) = z_nv_m;
    z_sv(eqMask) = z_sv_m;
  end
end

% South polar cap
if any(scMask(:))
  ip = nPix - hPix1(scMask) + 1;
  hip = ip/2;
  fihip = fix(hip);
  iRing = fix(sqrt(hip-sqrt(fihip)))+1;    % counted from south pole
  iPhi = 4*iRing+1-(ip-2*iRing.*(iRing-1));
  
  z(scMask) = -1+iRing.^2/fact2;
  phi(scMask) = (iPhi-0.5).*pi./(2*iRing);
  
  if (do_vertex)
    hDelta_Phi(scMask) = pi./(4*iRing);
    z_nv(scMask) = -1+(iRing+1).^2/fact2;
    z_sv(scMask) = -1+(iRing-1).^2/fact2;
    iPhi_Mod = mod(iPhi-1,iRing);  % in {0..iRing-1}
    iPhi_Rat = fix((iPhi-1)./iRing);    % in {0..3}
    phi_nv(scMask) = (iPhi_Rat+(iPhi_Mod+1)./(iRing+1))*pi/2;
    ndx = iRing > 1; 
    phi_sv_m = phi_sv(scMask);
    phi_sv_m(ndx) = (iPhi_Rat(ndx)+iPhi_Mod(ndx)./(iRing(ndx)-1))*pi/2;
    phi_sv(scMask) = phi_sv_m;
  end
end

sth = sqrt(1-z).*sqrt(1+z);

%% Vertices
if (do_vertex)
  v = zeros(horzcat([3,4],size(hPix)));
  
  % West vertex
  phi_wv = phi - hDelta_Phi;
  v(1,2,:) = sth(:).*cos(phi_wv(:));
  v(2,2,:) = sth(:).*sin(phi_wv(:));
  v(3,2,:) = z(:);
  
  % East vertex
  phi_ev = phi + hDelta_Phi;
  v(1,4,:) = sth(:).*cos(phi_ev(:));
  v(2,4,:) = sth(:).*sin(phi_ev(:));
  v(3,4,:) = z(:);
  
  % North vertex
  sth_nv = sqrt(1-z_nv).*sqrt(1+z_nv);
  v(1,1,:) = sth_nv(:).*cos(phi_nv(:));
  v(2,1,:) = sth_nv(:).*sin(phi_nv(:));
  v(3,1,:) = z_nv(:);
  
  % South vertex
  sth_sv = sqrt(1-z_sv).*sqrt(1+z_sv);
  v(1,3,:) = sth_sv(:).*cos(phi_sv(:));
  v(2,3,:) = sth_sv(:).*sin(phi_sv(:));
  v(3,3,:) = z_sv(:);
  
  varargout{2} = v;
end

%% Pixel centers
sth = reshape(sth,[],1);
phi = reshape(phi,[],1);
vector = zeros(horzcat(3,size(hPix)));
vector(1,:) = sth.*cos(phi);
vector(2,:) = sth.*sin(phi);
vector(3,:) = z(:);

varargout{1} = vector;

return
