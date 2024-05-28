function tf = qPix(pix)
% QPIX validate list of pixel values
%
% pix is a numeric vector with length 12*k^2 for integer k

% $Id: qPix.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

tf = isnumeric(pix) && numel(pix) == max(size(pix));
tf = tf && floor(sqrt(numel(pix)/12))^2*12 == numel(pix);

return