%% MEALPix for the HEALPix User
%% MEALPix function names are derived from their HEALPix counterparts.
% MEALPix functions are typically camel-cased (e.g., <angDist_help.html
% |angDist|>, or <nSide2nPix_help.html |nSide2nPix|>). Underscores in HEALPix
% names are dropped in the conversion to MEALPix names. 

%% MEALPix functions accept array arguments
% MEALPix functions exploit Matlab(tm)'s vectorization facilities to
% allow for the processing of arrays of pixels, direction vectors and
% position angle vectors via a single function call. For example, the
% MEALPix implementation of <angDist_help.html |angDist|> accepts mixed arrays of
% cartesian direction vectors, angular position vectors, and ring or nest
% indexed pixels for either argument.  

%% MEALPix indexes pixels from 1
% HEALPix pixels are indexed from 0. In keeping with the Matlab standard
% MEALPix pixels are indexed from 1. MEALPix pixels are indexed from 1:
% i.e., MEALPix pixel number n corresponds to HEALPix pixel number n-1.  

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.