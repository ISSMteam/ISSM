function [iring, iphi] = correct_ring_phi(location, iring, iphi)

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

delta = zeros(size(iphi));
delta(iphi<0) = 1;
delta(iphi>4*iring) = -1;
nzMask = delta ~= 0;
if any(nzMask)
  if any(abs(location(nzMask)) ~= 1) 
    msgid = [mfilename,':location'];
    error(msgid,msgid);
  end
  iring(nzMask) = iring(nzMask) - location(nzMask).*delta(nzMask);
  iphi(nzMask)  = iphi(nzMask) + delta(nzMask).*(4*iring(nzMask));
end
return