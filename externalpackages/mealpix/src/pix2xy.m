function xy = pix2xy(nSide, nPix)
%% PIX2XY Find pixel cartesian face coordinates
%
%  xy = pix2xy(nSide,nPix)
%
%  nSide     HEALPix resolution parameter
%  nPix      nested indexing pixel numbers to find coordinates of
%
%  xy        size(nPix) cell array of pixel coordinates on face
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: pix2xy.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


if ~qNSide(nSide)
  msgid = [mfilename ':nSide'];
  error(msgid,msgid);
end  

if ~(all(nPix(:) > 0) && all(nPix(:) <= nSide2nPix(nSide)))
  msgid = [mfilename ':nPix'];
  error(msgid,msgid);
end 
sPix = size(nPix);
if 2 == numel(sPix)
  sPix = [sPix, 1];
end
ixy = zeros([2,sPix(:)']);

[pix2x,pix2y]=mk_pix2xy;

ipf = nPix(:)-1;

% content of the last 10 bits
ip_low = fix(mod(ipf,1024));

% truncation of the last 10 bits
ip_trunc = fix(ipf/1024);

% content of the next 10 bits
ip_med = fix(mod(ip_trunc,1024));

% content of the high weight 10 bits
ip_hi  = fix(ip_trunc/1024);
ixy(1,:) = fix(1024*pix2x(ip_hi+1) + 32*pix2x(ip_med+1) + pix2x(ip_low+1));
ixy(2,:) = fix(1024*pix2y(ip_hi+1) + 32*pix2y(ip_med+1) + pix2y(ip_low+1));

xy = reshape(num2cell(ixy,1),sPix);

return