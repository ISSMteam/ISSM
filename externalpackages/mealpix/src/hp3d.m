function varargout = hp3d(vPix, varargin)
% HP3D Plots a HEALPix dataset on a 3D sphere
%
% h = hp3d(vPix,'Param1',Value1,'Param2',Value2,...);
%
% vPix   vector of values at HEALPix pixels
%
% Param      Value
% 'nest'     nest indexing (true | {false})
%
% h      (optional) array of patch handles for the sphere
%
% Example
%  % plot nested pixel numbers on a sphere
%  hp3d(1:768,'nest',true);
%
% Requires corners, nPix2nSide
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: hp3d.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.addRequired('vPix',@qVPix);
p.addParamValue('nest',false,@(x)(isscalar(x) && islogical(x)));

p.parse(vPix,varargin{:});
nest = p.Results.nest;

nSide = nPix2nSide(numel(vPix));
vPix = reshape(vPix,[],1);

%% Get pixel corners
c = corners(nSide,'nest',nest);

%% Plot patches
patchPlot = @(c,v)(...
  patch(c(1,:),c(2,:),c(3,:),v,'EdgeColor','none','FaceColor','flat'));
  % SA :: patch(c(1,:),c(2,:),c(3,:),v,'EdgeColor','black','FaceColor','white'));
h = cellfun(patchPlot,c,num2cell(reshape(vPix,1,[])));

% Cleanup the plot
%view(3);
view([90 90]);
axis equal off
colormap(jet(1024));
colorbar('southoutside');

% slower than opengl, but it doesnt crash X...
set(gcf, 'renderer', 'painters')

if 1 == nargout
  varargout{1} = h;
end

return

function tf = qVPix(vPix)
% QVPIX validate vPix
%
% vPix is a numeric vector of length 12*k^2 for k a non-negative integer
% power of 2

tf = isnumeric(vPix) && max(size(vPix)) == numel(vPix);

if tf
  nSide = floor(sqrt(numel(vPix)/12));
  ln2 = log2(nSide);
  tf = tf && ...
    numel(vPix) == 12*nSide^2 && ...
    floor(ln2) == ln2;
end

return
