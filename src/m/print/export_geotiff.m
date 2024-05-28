function export_geotiff(filename,ref)
%EXPORT_GEOTIF - export geotiff 
%
%   Usage:
%      export_geotif(filename,ref);
%
%      This function must be called after plotmodel
%      filname: no extension
%      ref:     UPS Greenland  EPSG:3413 (http://www.spatialreference.org/ref/epsg/3413/)
%               UPS Antarctica EPSG:3031 (http://www.spatialreference.org/ref/epsg/3031/)
%
%   Example:
%      export_geotif('Greenland','EPSG:3413');

%Get axis limits and convert to strings
XLIM = xlim(); x0 = num2str(XLIM(1)); x1 = num2str(XLIM(2));
YLIM = ylim(); y0 = num2str(YLIM(1)); y1 = num2str(YLIM(2));

%first export the figure
export_fig([filename '.jpg']);

%call gdal on this: 
system(['gdal_translate -a_srs ' ref ' -of GTiff -co "INTERLEAVE=PIXEL" -a_ullr ' x0 ' ' y1 ' ' x1 ' ' y0 ' ' filename '.jpg ' filename '.tif']);
system(['rm -rf ' filename '.jpg']);
