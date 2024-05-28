function nPix = nSide2nPix(nSide)
% NSIDE2NPIX Convert resolution parameter to number of pixels
%
%  nPix = nSide2nPix(nSide)
%
%  nSide    integer power-2 sub-divisions of facet sides. May be an array.
%  nPix     size(nSide) number of pixels on sphere corresponding to nSide.
%
% HEALPix divides the sphere into 12 basic facets. Each facet is further
% divided into nSide sub-divisions for a total of 12*nSide^2 pixels. 
%
% Example
% nPix = nSide2nPix([1:4])
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: nSide2nPix.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

lg = log2(nSide(:));
if ~all(abs(lg) == fix(lg))
  msgid = [mfilename ':nSide'];
  error(msgid,msgid);
end
nPix = 12*nSide.^2;

return