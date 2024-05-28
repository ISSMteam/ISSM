function alm = pix2alm(varargin)
% PIX2ALM Find spherical harmonic decomposition of function on sphere
%
%  alm = pix2alm(v)
%  alm = pix2alm(v, lmax)
%
%  v       array of pixel values
%  lMax    (optional) max order of harmonic to calculate
%
% nPix = numel(v), with nPix = 12*nSide^2 for nSide a power of 2. 
% lMax defaults to 2*floor(nSide/3)
%
%  alm     Spherical harmonic expansion coefficients
%
% alm(j) is coefficient corresponding to order L, M with j = (L+1)^2+M-L
%
% Example
% % estimates |a_l^m| of dummy data
% pix = ylm(16,4,2) + ylm(16,7,2) + ylm(16,16,-4) + .1.*rand(1,12*16*16);
% alm = pix2alm(16,pix,100);
% plot(abs(alm));
%
% Requires ylm
% See also alm2pix
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: pix2alm.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.addRequired('v',@qPix);
hLMax = @(v)(2*floor(nPix2nSide(numel(v))/3));
p.addOptional('lMax', hLMax(varargin{1}), @qLMax);

p.parse(varargin{:});

v = p.Results.v;
lMax = p.Results.lMax;

%%
% Preallocate alm
alm = zeros(1,(lMax+1)^2);

nSide = sqrt(numel(v)/12);
for L = 0:lMax
  ndx = L^2+1;
  for M = -L:L
    alm(ndx+L+M) = sum(conj(ylm(nSide,L,M)).*v);
  end
end
alm = alm*4*pi/numel(v);

return