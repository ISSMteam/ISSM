function nPix = ring2nest(nSide,rPix)
% RING2NEST Convert MEALPix pixel numbers from ring to nest indexing
%
%  nPix = ring2nest(nSide,rPix)
%
%
%  nSide      HEALPix resolution parameter
%  rPix       ring indexed MEALPix pixel numbers
%
%  nPix       next indexed MEALPix pixel numbers
%
% See also nest2ring
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: ring2nest.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
error(nargchk(1,2,nargin,'struct'));
if ~qNSide(nSide)
  msgid = [mfilename ':nSide'];
  error(msgid,msgid);
end
if isempty(rPix)
  rPix = 1:(12*nSide^2);
elseif any(fix(abs(rPix(:))) ~= rPix(:)) || any(rPix(:)<1)
  msgid = [mfilename ':rPix'];
  error(msgid,msgid);
end

% Convert
sz = size(rPix);
rPix = rPix - 1;
nPix = f_ring2nest(nSide,reshape(rPix,[],1));
nPix = reshape(nPix+1,sz);

return