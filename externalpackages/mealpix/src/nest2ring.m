function rPix = nest2ring(nSide,nPix)
% NEST2RING Convert MEALPix pixel numbers from nest to ring indexing
%
%  rPix = nest2ring(nSide,nPix)
%
%  nSide     HEALPix resolution parameters
%  nPix      numeric array of MEALPix nest indexed pixel numbers
%
%  rPix      size(nPix) numeric array of ring indexed pixel numbers
%
% See also ring2nest
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: nest2ring.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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
p.addRequired('nPix',@(n)(qPix(nSide,n)));
p.parse(nSide,nPix);

%% Convert
% Save shape
% convert HEALPix indexing as vector
% convert to MEALPix indexing and restore shape 

sz = size(nPix);
rPix = f_nest2ring(nSide,nPix(:)-1);
rPix = reshape(rPix+1,sz);

return

function tf = qPix(nSide,nPix)
tf = isnumeric(nPix);
if tf
  nPix = nPix(:);
  tf = ...
    all(fix(abs(nPix)) == nPix) && ...
    all(0 < nPix) && ...
    all(nSide2nPix(nSide) >= nPix);
end
return