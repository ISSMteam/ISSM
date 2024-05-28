function res = nSide2res(nSide)
% NSIDE2RES Calculates angular resolution of a HEALPix map in arcsec 
%
% angres = nSide2res(nSide)
%
% nSide   HEALPix Resolution (or array of nSides)
%
% angres  Angular resolution of the given nSide in arcsec
%
% See also res2nSide
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: nSide2res.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


if ~(all(nSide > 0) && all(2.^fix(log2(nSide)) == nSide))
  msgid = [mfilename ':nSide'];
  error(msgid,msgid);
end
area = 4*pi*(180*3600/pi)^2;
nPix = nSide2nPix(nSide);
res = sqrt(area./nPix); 

return