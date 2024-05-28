function v = alm2pix(nSide, alm, varargin)
% ALM2PIX Evaluate spherical harmonic function expansion on HEALPix pixels
%
% v = alm2pix( nSide, alm, 'Param1', Value1, 'Param2', Value2, ...)
%
% nSide   HEALPix resolution parameter (power of 2)
% alm     Spherical harmonic expansion coefficients
%
% Param   Value
% 'lmax'  max order of harmonic to calculate (default floor2*nSide/3)
% 'nest'  (true | {false})
%
% v       array of function values at pixels (numel(v) = 12*nSide^2)
%
% Coefficient of Y_{lm} is alm((l+1)*(l+1)+m-1)
% 
% Example
%  % Shows the gibbs-like error at the poles
%  hp3d(alm2pix(16,pix2alm(16,ones(1,12*16^2)))-1)
%
% Requires ylm
% See also pix2alm
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: alm2pix.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.addRequired('alm',@qALM);
p.addParamValue('lMax', 2*floor(nSide/3),@qLMax); %% SA corrected it on 11/03/2014 
%p.addParamValue('lMax', floor(2*nSide/3),@qLMax);
p.addParamValue('nest', false, @islogical)

p.parse(nSide, alm, varargin{:});

nSide = p.Results.nSide;
alm = p.Results.alm;
lMax = p.Results.lMax;
nest = p.Results.nest;

%% 
nPix = nSide2nPix(nSide);
v = zeros(1,nPix);
p = 0; % alm index
lMax = max(lMax,sqrt(numel(alm))-1); 
for L = 0:lMax
  for M = -L:L
    v = v + alm(1+p).*ylm(nSide,L,M,'nest',nest);
    p = p + 1;
  end
end

return

function tf = qLMax(lMax)
% QLMAX validate lMax
%
% lMax is a non-negative integer-valued scalar

tf = isnumeric(lMax) && isscalar(lMax) && floor(lMax) == lMax && lMax >= 0;

return
