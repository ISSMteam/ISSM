function string = laea(lat,long)
%Lambert Azimuthal Equal Area projection at lat,long projection center. 
%
%   Usage:
%      string = laea(45,-90);
%
%   Example: 
%      string = laea(45,-90); 
%      return string='+proj=laea +lat_0=45 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'

	string=sprintf('+proj=laea +lat_0=%i +lon_0=%i +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs',lat,long);
