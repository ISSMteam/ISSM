function info=gdalinfo(imagename)
%GDALINFO - retrieve information from a geotiff or georeferenced image
%
%   Usage:
%      gdalinfo(imagename)
%
%

	%first, get pixel size: 
	[s,r]=system(sprintf('gdalinfo %s | command grep "Pixel Size"',imagename));
	if s, error('gdalinfo error message: could not run system command gdalinfo to find pixel size'); end
	d=sscanf(r,'Pixel Size = (%g,%g)');
	info.dx=abs(d(1)); info.dy=abs(d(2));

	%get upper left: 
	[s,r]=system(sprintf('gdalinfo %s | command grep "Upper Left"',imagename));
	if s, error('gdalinfo error message: could not run system command gdalinfo to find upper left'); end
	d=sscanf(r,'Upper Left ( %g,%g) %s');
	info.xmin=d(1);
	info.ymax=d(2);
	
	%get lower right: 
	[s,r]=system(sprintf('gdalinfo %s | command grep "Lower Right"',imagename));
	if s, error('gdalinfo error message: could not run system command gdalinfo to find lower right'); end
	d=sscanf(r,'Lower Right ( %g,%g) %s');
	info.xmax=d(1);
	info.ymin=d(2);

	%Size: 
	[s,r]=system(sprintf('gdalinfo %s | command grep "Size is"',imagename));
	if s, error('gdalinfo error message: could not run system command gdalinfo to find size'); end
	d=sscanf(r,'Size is %g, %g');
	info.nx=d(1);
	info.ny=d(2);

	%Name: 
	[pathstr,name,ext] = fileparts(imagename);
	info.rootname=name;
	info.ext=ext;
