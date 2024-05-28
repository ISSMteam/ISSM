function prjwrite(filename,proj)

	[path,name,ext] = fileparts(filename);
	projname=[path '/' name '.prj'];

	projdir='/Users/larour/issm-jpl/proj-group/qgis/Proj/';
	if strcmpi(proj,'EPSG:3031'),
		copyfile([projdir '/3031.prj'],projname,'f');
	elseif strcmpi(proj,'EPSG:3413'),
		copyfile([projdir '/3413.prj'],projname,'f');
	else 
		error('not supported yet');
	end
