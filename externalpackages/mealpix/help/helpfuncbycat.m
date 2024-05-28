%% Functions by Category
% MEALPix Toolbox
% Version 3.0
% 21 February 2011
 
%% Coordinate Conversions
% * <ang2vec_help.html |ang2vec|>    - Convert from spherical to cartesian coordinates
% * <vec2ang_help.html |vec2ang|>    - Convert cartesian direction vector to position angle on sphere

%% Coordinate-Pixel Conversions
% * <ang2pix_help.html |ang2pix|>    - Convert spherical coordinate location to HEALPix pixel number
% * <pix2ang_help.html |pix2ang|>    - Find spherical coordinate location of HEALPix pixels
% * <pix2vec_help.html |pix2vec|>    - Convert HEALPix pixel numbers to (x,y,z) unit vector
% * <ring2z_help.html |ring2z|>     - Find HEALPix ring number corresponding to a given z coordinate
% * <vec2pix_help.html |vec2pix|>    - Convert cartesian direction vectors to MEALPix pixel numbers
% * <z2ring_help.html |z2ring|>     - Find HEALPix ring numbers from direction vector z coordinate

%% Functions on the sphere
% * <alm2pix_help.html |alm2pix|>    - Evaluate spherical harmonic function expansion on HEALPix pixels
% * <angDist_help.html |angDist|>    - Computes angular distance on the unit sphere
% * <pix2alm_help.html |pix2alm|>    - Find spherical harmonic decomposition of function on sphere
% * <ylm_help.html |ylm|>        - returns spherical harmonic basis $Y_L^M$ on the HEALPix sphere

%% Map Statistics
% * <alm2spec_help.html |alm2spec|>   - valuate angular power spectral density
% * <hpStat_help.html |hpStat|>     - Compute descriptive statistics of a pixel value list

%% Pixel Indexing Conversions
% * <nest2ring_help.html |nest2ring|>  - Convert MEALPix pixel numbers from nest to ring indexing
% * <ring2nest_help.html |ring2nest|>  - Convert MEALPix pixel numbers from ring to nest indexing

%% Pixel Topology
% * <corners_help.html |corners|>    - Find pixel corners
% * <inRing_help.html |inRing|>     - return pixels in ring or ring section
% * <pix2xy_help.html |pix2xy|>     - Find pixel cartesian face coordinates
% * <queryDisc_help.html |queryDisc|>  - Find all pixels whose centers are within a specified disc
% * <xy2pix_help.html |xy2pix|>     - gives the pixel number ipix (NESTED) of ix, iy and face_num

%% Plotting
% * <hp3d_help.html |hp3d|>       - Plots a HEALPix dataset on a 3D sphere

%% Resolution Relations
% * <nPix2nSide_help.html |nPix2nSide|> - Convert number of pixels to number of facet side subdivisions
% * <nSide2nPix_help.html |nSide2nPix|> - Convert resolution parameter to number of pixels
% * <nSide2res_help.html |nSide2res|>  - Calculates angular resolution of a HEALPix map in arcsec 
% * <res2nSide_help.html |res2nSide|>  - Finds minimum base pixel sub-division for required resolution

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.
