function cL = alm2spec(varargin)
% ALM2SPEC valuate angular power spectral density
%
%  cL = alm2spec(alm)
%  cL = alm2spec(alm, lMax)
%
%  alm     real valued coefficients
%  lMax    max order of harmonic to calculate (default 20)
%
%  cL      angular power spectrum
%
% Coefficients alm are arranged such that alm(k) is the coefficient
% corresponding to k = (L+1)^2 + M - L
%
% Example
%  % Power spectrum estimate of dummy data
%  pix = ylm(16,4,2) + ylm(16,7,2) + ylm(16,16,-4) + .1.*rand(1,12*16*16);
%  alm = pix2alm(16,pix,100);
%  cl = alm2spec(alm,100);
%  plot(cl);
%
% See also pix2alm, alm2pix
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: alm2spec.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.addRequired('alm',@qAlm);
p.addOptional('lMax',20,@qLMax);
p.parse(varargin{:});

alm = p.Results.alm;
lMax = p.Results.lMax;

%% Compute power spectrum
% Angular power spectrum for given L is sum over m of (|alm|^2/(2L+1))

% Preallocate
cL=zeros(1,lMax+1);

cL(1) = abs(alm(1))^2;
for L = 1:lMax
  % Coefficients for L run from L^2+1 to (L+1)^2
  cL(L+1) = sum(abs(alm(L^2+(1:(2*L+1)))).^2);
end
cL = cL./(2*(0:L)+1);

return