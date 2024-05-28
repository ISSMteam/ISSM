function proj2shpprj(filename,proj)

	%filename: does it have an extension? remove it. 
	[path,name,ext] = fileparts(filename);

	%Figure out central meridian and latitude of origin: 
	[lat,long]=projlatlong(proj);

	fid=fopen([path '/' name '.prj'],'w');
	fprintf(fid,'PROJCS["Lambert_Azimuthal_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_unknown",SPHEROID["WGS84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_origin",%g],PARAMETER["central_meridian",%g],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]',lat,long);
	fclose(fid);
