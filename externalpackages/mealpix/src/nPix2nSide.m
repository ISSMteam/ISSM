function nSide = nPix2nSide(nPix)
% NPIX2NSIDE Convert number of pixels to number of facet side subdivisions
%
%  nSide = nPix2nSide(nPix);
%
%  nPix     number of pixels on sphere
%  nSide    sub-divisions of HEALPix facet side
%
% Example
% nPix = 12*[1:4].^2;
% nSide = nPix2nSide(nPix)
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: nPix2nSide.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

nSide = sqrt(nPix/12);
if any(nSide(:) ~= fix(nSide(:)))
  msgid = [mfilename ':nPix'];
  error(msgid,msgid);
end

return