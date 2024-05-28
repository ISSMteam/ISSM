%% hpStat 
% Compute descriptive statistics of a pixel value list

%% Syntax
%  stats = hpStat(v)

%% Input Arguments
%  v         array of HEALPix pixel values

%% Output Arguments
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

%% Example
% display stats of Y_3^1
v = ylm(32,3,1,'real',true);
stats = hpStat(abs(v))

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.