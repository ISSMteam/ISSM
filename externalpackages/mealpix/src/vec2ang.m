function tp = vec2ang(xyz)
% VEC2ANG Convert cartesian direction vector to position angle on sphere
%
%  tp = vec2ang(xyz)
%
%  xyz    cell array of [3,1] cartesian coordinate direction vectors
%
% The direction vectors need not be normalized.
%
%  tp     cell array of [2,1] position angles ([theta; phi], radians) on
%         sphere 
%
% See also ang2vec
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: vec2ang.m 5831 2011-02-25 23:55:19Z lsf@GRAVITY.PSU.EDU $

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


if ~iscell(xyz)
  msgid = [mfilename ':xyz'];
  error(msgid,msgid);
end
tp = cellfun(@v2a,xyz,'UniformOutput',false);

return

function tp = v2a(xyz)
% dVec2Ang un-nomralized cartesian vector to position angle on sphere
%
% tp = v2a(xyz)
%
% xyz   [3,1]
%
% tp    [2,1] of [theta, phi]

assert(numel(xyz)==3,[mfilename ':v2a'],[mfilename 'v2a']);
xyz = xyz/norm(xyz);
tp = [acos(xyz(3)); atan2(xyz(2),xyz(1))];

return