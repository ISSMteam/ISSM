function v = ylm(nSide,L,M,varargin)
% YLM returns spherical harmonic basis Y_L^M on the HEALPix sphere
%
%  v = ylm(nSide,L,M,'Param1',Value1,'Param2',Value2,...);
%
%  nSide    HEALPix resolution parameter (power of 2)
%  L        spherical harmonic degree (0 <= L)
%  M        spherical harmonic order (-M <= L <= M)
%
%  Param    Value
%  'real'   return real values (true | {false})
%  'nest'   nested indexing flag (true | {false}) 
%
%  v        spherical harmonics evaluated at pixel centers
%
% Example:
%
% % Plot Ylm for (L,M) = (4,2) on nSide = 16 HEALPix sphere
% pix = ylm(16,4,2);
% hp3d(pix);
%
% Requires pix2ang
% Required by alm2pix, pix2alm
% See also hp3d
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: ylm.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.KeepUnmatched = false;
p.FunctionName = mfilename;

p.addRequired('nSide',@qNSide);
p.addRequired('L',@(x)(isnumeric(x) && round(x) == x && 0 <= x));
p.addRequired('M',@(x)(isnumeric(x) && round(x) == x));
p.addParamValue('real',false,@(x)(isscalar(x) && islogical(x)));
p.addParamValue('nest',false,@(x)(isscalar(x) && islogical(x)));

p.parse(nSide,L,M,varargin{:});
nSide = p.Results.nSide;
L = p.Results.L;
M = p.Results.M;
real = p.Results.real;
nest = p.Results.nest;

% Check m is in range
assert(abs(M)<=L,[mfilename ':M'],[mfilename ':M']);

% Setup HEALPix pixels to solve on
tp = pix2ang(nSide,'nest',nest);

% Evaluate for positive M, using symmetry to find value for negative M
signM = sign(M);
M = abs(M);

% Exploit HEALPix rings for legendre function calc
th = cellfun(@(tp)(tp(1)),tp);
[b, ~, n] = unique(th);
ringLmn = legendre(L,cos(b));
if L~=0
  ringLmn = squeeze(ringLmn(M+1,:,:));
end
Lmn = ringLmn(n);

% Calculate normalization constant
%a1 = (2*L+1);       % 4pi normalization 
a1 = (2*L+1)/(4*pi); % orthonormalization (default) 

% a2 = factorial(L-M)/factorial(L+M);
a2 = 1/prod((L-M+1):(L+M));

C = sqrt(a1*a2);

% Calculate complex Ymn
phi = cellfun(@(x)(x(2)),tp);
Ymn = C*Lmn.*exp(1i*M*phi);

% Use symmetry to find value for negative M
if 0 > signM
  Ymn = (-1)^M*conj(Ymn);
end

% Set v
v = Ymn;

% Convert to real form if true == real
if real && 0 ~= M
  if 0 < signM
    v = (Ymn+conj(Ymn))/sqrt(2);
  else % 0 > signM
    v = (-1)^M*(conj(Ymn)-Ymn)/sqrt(-2);
  end
end

return
