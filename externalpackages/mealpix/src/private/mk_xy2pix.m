function [x2pix y2pix]=mk_xy2pix
%=======================================================================
%     sets the array giving the number of the pixel lying in (x,y)
%     x and y are in {1,128}
%     the pixel number is in {0,128**2-1}
%
%     if  i-1 = sum_p=0  b_p * 2^p
%     then ix = sum_p=0  b_p * 4^p
%          iy = 2*ix
%     ix + iy in {0, 128**2 -1}
%=======================================================================

% $Id: mk_xy2pix.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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


%for converting x,y into pixel numbers
x2pix=zeros(128,1);
y2pix=zeros(128,1);
for i = 1:128;
  j  = (i-1);
  k  = 0;
  ip = 1;
  while (1);
    if(j==0)
      x2pix(i) = k;
      y2pix(i) = 2*k;
      break;
    else
      id = mod(j,2);
      j  = fix(j/2);
      k  = ip*id+k;
      ip = 4*ip;
    end
  end
end
return