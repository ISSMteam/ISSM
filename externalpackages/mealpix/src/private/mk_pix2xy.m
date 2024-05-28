function [pix2x,pix2y] = mk_pix2xy

%=======================================================================
% constructs the array giving x and y in the face from pixel number for the
% nested (quad-cube like) ordering of pixels 
%
% the bits corresponding to x and y are interleaved in the pixel number one
% breaks up the pixel number by even and odd bits 
%=======================================================================

% $Id: mk_pix2xy.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

pix2x=zeros(1024,1);
pix2y=zeros(1024,1);
for kpix=0:1023
  jpix = kpix;
  ix = 0;
  iy = 0;
  % bit position (in x and y)
  ip = 1;
  while (1);
    % go through all the bits
    if(jpix == 0)
      break
    end
    % bit value (in kpix), goes in ix
    id = mod(jpix,2);
    jpix = fix(jpix/2);
    ix = id.*ip+ix;
    % bit value (in kpix), goes in iy
    id = mod(jpix,2);
    jpix = fix(jpix/2);
    iy = id.*ip+iy;
    % next bit (in x and y)
    ip = 2*ip;
  end
  % in 0,31
  pix2x(kpix+1) = ix;
  % in 0,31
  pix2y(kpix+1) = iy;
end

return