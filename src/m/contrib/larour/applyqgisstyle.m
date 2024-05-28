function applyqgisstyle(filename,type); 

	%Getting extension out: 
	[path,name,ext]=fileparts(filename);

	%New name: 
	if isempty(path),
		newname=[name '.qml'];
	else
		newname=[path '/' name '.qml'];
	end

	if strcmpi(type,'mesh'),
		system(['cp /Users/larour/issm-jpl/usr/larour/Qgis/Mesh.qml ' newname]);
	elseif strcmpi(type,'line'),
		system(['cp /Users/larour/issm-jpl/usr/larour/Qgis/LineStyle.qml ' newname]);
	elseif strcmpi(type,'point'),
		system(['cp /Users/larour/issm-jpl/usr/larour/Qgis/SummitStyle.qml ' newname]);
	else
		error(['applyqgisstyle error message: ' type ' style not supported yet']);
	end
