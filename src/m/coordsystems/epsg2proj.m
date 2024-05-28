function string = epsg2proj(epsg)
%EPSG2PROJ - uses gdalsrsinfo to provide PROJ.4 compatible string from 
%EPSG code
%
%   Usage:
%      proj4string = epsg2proj(4326);
%
%   Example:
%      proj4string = epsg2proj(4326);
%      return proj4string='+proj=longlat +datum=wgs84 +no_defs'
%

%Call PROJ library
[status, string]=system(['projinfo -o PROJ -q epsg:' num2str(epsg)]);

%Check status
if status~=0; error(string); end

%remove trailing  blanks
string = deblank(string);
