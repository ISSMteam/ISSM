function tp = pix2ang(nSide, varargin)
% PIX2ANG Find spherical coordinate location of HEALPix pixels
%
%  tp = pix2ang(nSide,pix,'Param1',Value1,...)
%
%  nSide    HEALPix resolution parameter (positive integer power of 2)
%  pix      (optional) array of pixel values to find coordinates of
%
%  Param    Values
%  'nest'   pixels are given in nested indexing (true | {false})
%
%  tp       cell array of [theta; phi] pixel locations (radians)
%
% Example
%   % Return theta/phi for all pix of nSide=4, nest scheme
%   tp = pix2ang(4,'nest')
%
% See also ang2pix, pix2vec, vec2pix
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: pix2ang.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


%% I/O Check
p = inputParser;
p.addRequired('nSide', @qNSide);
p.addOptional('pix',-1,@(n)(qPix(nSide,n)));
p.addParamValue('nest',false,@(n)(isscalar(n) && islogical(n)));

p.parse(nSide,varargin{:});

pix = p.Results.pix;
nest = p.Results.nest;

% Matlab indexes pixels from 1
if pix < 0
  pix = 1:nSide2nPix(nSide);
end

%% Work
% Save size
% Convert to HEALPix convention
% Get pixel location angles
% Restore size

sz = size(pix);
pix = reshape(pix - 1,1,[]);
if nest
  [tparray(1,:), tparray(2,:)]=f_pix2ang_nest(nSide,pix);
else
  [tparray(1,:), tparray(2,:)]=f_pix2ang_ring(nSide,pix);
end
tp=reshape(num2cell(tparray,1),sz);

return

function tf = qPix(nSide,n)
tf = ...
  isnumeric(n(:)) && ...
  all(n(:) == fix(abs(n(:)))) && ...
  all(n(:)<=nSide2nPix(nSide)); 
return