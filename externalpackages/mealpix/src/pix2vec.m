function varargout = pix2vec(nSide, varargin)
% PIX2VEC convert HEALPix pixel numbers to (x,y,z) unit vector
%
% [x,y,z] = pix2vec(nSide,pix,'Param1',Value1,...)
% xyz = pix2vec(nSide,pix,'Param1',Value1,...)
%
% nSide       HEALPix resolution parameter
% pix         (optional) numeric array of pixel numbers to find cartesian
%             coordinate vectors for 
%
% Param       Value
% 'nest'      indexing scheme is nested (true | {false}) 
%
%  x,y,z       size(pix) numeric array of unit vector coordinate
%              components corresponding to pixels pix
%  xyz         size(pix) cell array of cartesian unit vectors
%              corresponding to pixels pix
%
% Example
% % 
% nSide = 2^4;
% nPix = nSide2nPix(nSide);
% pix = randi(nPix,2,3);
% xyz = pix2vec(nSide,pix);
% xyzn = pix2vec(nSide,pix,'nest',true);
% 
% Requires nest2ring
% Required by angDist, hpInterp, inRing, queryDisc, ringNum
% See also ang2pix, pix2ang, vec2pix
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: pix2vec.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


%% IO Check
assert(1 == nargout || 3 == nargout,[mfilename ':nargout']);

%%% Setup parser
p = inputParser;
p.FunctionName = mfilename;
p.StructExpand = true;

% Argument declarations
p.addRequired('nSide',@qNSide);
p.addOptional('pix',-1,@(n)(isnumeric(n) && all(floor(n(:))==n(:))));
p.addParamValue('nest',false,@islogical);

%%% parse and validate input
p.parse(nSide,varargin{:});
nSide = p.Results.nSide;
pix = p.Results.pix;
nest = p.Results.nest;

%% Get pixel vector

if 1 == numel(pix) && pix < 0
  pix = 1:nSide2nPix(nSide);
end
assert(...
  all(0 < pix(:)) || all(pix(:) <= nSide2nPix(nSide)),...
  [mfilename ':pix']);

if nest
  pix = nest2ring(nSide,pix);
end
xyz = f_pix2vec_ring(nSide,pix-1);

%% Format return
switch nargout
  case 1,
    nSize = size(xyz);
    if 2 == numel(nSize)
      nSize = [nSize 1];
    end
    xyz = reshape(xyz,3,[]);
    xyz = num2cell(xyz,1);
    varargout{1} = reshape(xyz,nSize(2:end));
  case 3,
    nSize = size(xyz);
    if 2 == numel(nSize);
      nSize = [nSize 1];
    end
    sz = nSize(2:end);
    xyz = reshape(xyz,3,[]);
    varargout{1} = reshape(xyz(1,:),sz);
    varargout{2} = reshape(xyz(2,:),sz);
    varargout{3} = reshape(xyz(3,:),sz);
end

return