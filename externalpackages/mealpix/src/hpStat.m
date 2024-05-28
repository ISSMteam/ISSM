function stats = hpStat(v)
% HPSTAT Compute descriptive statistics of a pixel value list
%
%  stats = hpStat(v)
%
%  v         array of HEALPix pixel values
%
%  stats     struct
%            field     contents
%            ntot      total number of data points
%            nvalid    number n of valid data points
%            mind      minimum of valid data
%            maxd      maximum of valid data
%            average   average of valid points
%            absdev    absolute deviation
%            var       variance
%            stddev    standard deviation
%            skew      skewness factor
%            kurt      kurtosis factor
%
% Example
% % display stats of Y_3^1
% v = ylm(32,3,1,'real',true);
% stats = hpStat(abs(v))
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: hpStat.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

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

stats.ntot = numel(v);
nValid = sum(~isnan(v));
stats.nValid = nValid;
v = v(~isnan(v));
stats.mind = min(v);
stats.maxd = max(v);
stats.mean = mean(v);
stats.absdev = mean(abs(v-mean(v)));
stats.var = var(v);
stats.stddev = std(v);
stats.skew = skewness(v);
stats.kurt = kurtosis(v);

return