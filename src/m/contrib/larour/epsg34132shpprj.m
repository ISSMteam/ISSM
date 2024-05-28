function epsg34132shpprj(filename)

	%filename: does it have an extension? remove it. 
	[path,name,ext] = fileparts(filename);

	fid=fopen([path '/' name '.prj'],'w');
	fprintf(fid,'PROJCS["WGS_84_NSIDC_Sea_Ice_Polar_Stereographic_North",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Stereographic_North_Pole"],PARAMETER["standard_parallel_1",70],PARAMETER["central_meridian",-45],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]');
	fclose(fid);
