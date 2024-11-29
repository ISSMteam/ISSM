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
	% Set the PROJ_LIB environment variable
	%setenv('PROJ_LIB', '/home/connord/issm/ISSM-Davis/externalpackages/proj/install/share/proj/');

	% Set the LD_LIBRARY_PATH environment variable
	%setenv('LD_LIBRARY_PATH', ['/issm/ISSM-Davis/externalpackages/proj/install/lib:', getenv('LD_LIBRARY_PATH')]);

	%give ourselves unique file names
	filename_in  = tempname();
	filename_out = tempname();

	fid=fopen(filename_in,'w');
	fprintf(fid,'%8g %8g\n',[x(:) y(:)]');
	fclose(fid);
	%proj_in = ['"', proj_in, '"'];
	%proj_out = ['"', proj_out, '"'];
	%[s,r] = system(['$ISSM_DIR/src/m/contrib/cdavis/gdaltransform_cdavis "', proj_in, '" "', proj_out, '" ', filename_in, ' ', filename_out]);
	% Call the shell script with the proper quoting
	%[s,r] = system(['$ISSM_DIR/src/m/contrib/cdavis/gdaltransform_cdavis ', proj_in, ' ', proj_out, ' ', filename_in, ' ', filename_out]);
	%[s,r]=system(['gdaltransform_cdavis ',proj_in,' ',proj_out,' ',filename_in,' ', filename_out]);
	
	%LAST WORKING PATH
    	%[s,r] =system(['$ISSM_DIR/src/m/contrib/cdavis/gdaltransform_cdavis ', proj_in, ' "', proj_out, '" ', filename_in, ' ', filename_out]);
	
	%[s,r]=system(['export $LD_LIBRARY_PATH = "/issm/ISSM-Davis/externalpackages/proj/install/lib:$LD_LIBRARY_PATH" && gdaltransform -s_srs "',proj_in,'" -t_srs "',proj_out,'" < ' filename_in ' > ' filename_out]);

	%[s,r]=system(['gdaltransform_cdavisd -s_srs "',proj_in,'" -t_srs "',proj_out,'" < ' filename_in ' > ' filename_out]);
	%[s,r]=system(['gdaltransform -s_srs "',proj_in,'" -t_srs "',proj_out,'" < ' filename_in ' > ' filename_out]);
	command = sprintf('$ISSM_DIR/src/m/contrib/cdavis/gdaltransform_cdavis "%s" "%s" %s %s', proj_in, proj_out, filename_in, filename_out);

	% Execute the system command
	[s, r] = system(command);
	if s~=0 | ~isempty(deblank(r)),
		error(r);
	end

	%Uncomment if ASCII error occurs
	A=load(filename_out,'-ASCII');
	%A=load(filename_out);

	%clean up
	delete(filename_in);
	delete(filename_out);

	xout=A(:,1);
	xout=reshape(xout,size(x));
	yout=A(:,2);
	yout=reshape(yout,size(y));
