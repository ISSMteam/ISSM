%% Product Overview
% MEALPix is a native Matlab implementation of a substantial subset of the 
% <http://healpix.jpl.nasa.gov/ HEALPix> subroutine library. It supports
% the fast and accurate statistical and astrophysical analysis of massive
% full-sky data sets. HEALPix is supported by NASA and ESO, used
% extensively by the WMAP and PLANCK missions, and is part of the official 
% FITS World Coordinate System.
%
% MEALPix provides fully vectorized functions for 
%
% * Conversion between position angles, cartesian direction vectors, and
% HEALPix ring and nest pixel indexing schemes; 
%
% * Spherical harmonic resolution and analysis of functions on the sphere; 
%
% * Nearest pixel identification (pixels within a given disc, pixels in a
% ring, pixels in a HEALPix face). 
%
% MEALPix was created from the Fortran90 implementation
% <http://healpix.jpl.nasa.gov HEALPix> via a combination of hand coding
% and machine assisted conversion of the HEALPix F90 source code. Some
% routines were hand-coded based on the documented functionality of the F90
% subroutines on which they were based. Other routines were hand-coded via
% inspection of the HEALPix F90 implementation source. Still other routines
% were machine converted from F90 to Matlab(TM) using
% <http://www.mathworks.com/matlabcentral/fileexchange/5260 f2matlab>. In
% this latter case the resulting Matlab(TM) source was then hand-modified
% to vectorize the computations and take advantage of Matlab(TM) data types
% and other functions.

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.