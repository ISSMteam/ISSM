function xyz = ang2vec(tp)
% ANG2VEC Convert from spherical to cartesian coordinates
%
% xyz = ang2vec(tp)
%
%  tp        cell of [2,1] angular locations [theta; phi] in radians
%
%  xyz       size(tp) celll array of [3,1] vectors
%
% Example
% tp = [acos(2*rand(1,12)-1); 2*pi*rand(1,12)];
% tp = mat2cell(tp,2,ones(12,1));
% tp = reshape(tp,3,4);
% xyz = ang2vec(tp);
% xyz{2,3}
% tp{2,3}
%
% See also vec2ang
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: ang2vec.m 5831 2011-02-25 23:55:19Z lsf@GRAVITY.PSU.EDU $

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


%% IO check
error(nargchk(1,1,nargin,'struct'));
error(nargoutchk(1,1,nargout,'struct'));

xyz = cellfun(@dAng2Vec,tp,'UniformOutput',false);

return

function xyz = dAng2Vec(tp)
% dAng2Vec location on sphere to unit cartesian direction vector
%
% xyz = dAng2Vec(tp)
%
% tp   [theta, phi] in radians
%
% xyz   [3,1]

t = tp(1);
z = cos(t);
s = sin(t);
p = tp(2);
xyz = [cos(p)*s;sin(p)*s;z];

return