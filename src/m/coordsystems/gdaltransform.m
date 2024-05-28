function [xout,yout] = gdaltransform(x,y,proj_in,proj_out)
%GDALTRANSFORM - switch from one projection system to another 
%
%	Usage:
%		[x,y] = gdaltransform(x1,y1,epsg_in, epsg_out);
%
%	Example:
%		[x,y] = gdaltransform(md.mesh.long,md.mesh.lat,'EPSG:4326','EPSG:3031')
%
%	For reference:
%		EPSG: 4326 (lat,long)
%		EPSG: 3411 (Greenland,  UPS 45W, 70N)
%		EPSG: 3031 (Antarctica, UPS 0E,  71S)
%
%		ll2xy default projection Antarctica:
%			+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
%		ll2xy default projection Greenland:
%			+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
%		Bamber's Greeland projection
%			+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
%
%	To get PROJ.4 string from EPSG, use gdalsrsinfo. Example:
%		gdalsrsinfo epsg:4326 | grep "PROJ.4" | sed "s/PROJ.4 : //"

	%give ourselves unique file names
	filename_in  = tempname();
	filename_out = tempname();

	fid=fopen(filename_in,'w');
	fprintf(fid,'%8g %8g\n',[x(:) y(:)]');
	fclose(fid);

	[s,r]=system(['gdaltransform -s_srs "',proj_in,'" -t_srs "',proj_out,'" < ' filename_in ' > ' filename_out]);
	if s~=0 | ~isempty(deblank(r)),
		error(r);
	end

	A=load(filename_out);

	%clean up
	delete(filename_in);
	delete(filename_out);

	xout=A(:,1);
	xout=reshape(xout,size(x));
	yout=A(:,2);
	yout=reshape(yout,size(y));
