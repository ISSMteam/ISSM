function xyz = corners(nSide,varargin)
% CORNERS Find pixel corners
%
%  xyz = corners(nSide,nPix,'Param1',Value1,...);
%
%  nSide    HEALPix resolution parameter
%  nPix     (optional) pixel list. Default all pixels.
%  
%  Param    Value
%  'nest'   nested indexing (true | {false})
%
%  xyz    size(nPix) cell array of [3,1] cartesian vector pixel corner
%         locations
%
% Example:
% nSide = 2^2;
% % Corners for all 48 nSide=4 pixels in ring-indexed order
% xyzR = corners(nSide);
% size(xyzR)
% xyzR{1} 
% xyzR{end-1}
% 
% nPix = randi(nSide2nPix(nSide), 3, 4);
% % Corners for an assortment of pixels in nested index order
% xyzN = corners(nSide,nPix,'nest',true);
% % Convert to pixel list to ring indexing
% r = xyzR{nest2ring(nSide,nPix(2,3))}
% % Corners for same pixels requested via ring indexing
% n = xyzN{2,3}
% % Compare
% all(r == n)
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: corners.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

%% Validate and parse input
p = inputParser;
p.addRequired('nSide',@qNSide);
p.addOptional('nPix',[],@isnumeric);
p.addParamValue('nest',false,@islogical);

p.parse(nSide,varargin{:});
nPix = p.Results.nPix;
nest = p.Results.nest;

%% Get corners
if isempty(nPix)
  nPix = 1:nSide2nPix(nSide);
end
if nest
  % convert to ring indexing
  nPix = nest2ring(nSide,nPix);
end
nPix = nPix - 1;  % convert from MEALPix to HEALPix indexing
[~,xyzc]=f_pix2vec_ring(nSide,nPix(:));

xyz = reshape(xyzc,12,[]);
xyz = reshape(num2cell(xyz,1),size(nPix));
xyz = cellfun(@(x)(reshape(x,3,4)),xyz,'UniformOutput',false);

return