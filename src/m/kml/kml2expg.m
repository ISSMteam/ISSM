function kml2expg(filename)

	[path,name,ext]=fileparts(filename);

	if strcmpi(ext,'.kmz'),
		eval(['!unzip ' filename]);
		eval(['!mv doc.kml ' name '.kml']);
		kml2exp([name '.kml'],[name '.exp']);
		expll2xy([name '.exp'],1);
	end

	kml2exp([name '.kml'],[name '.exp']);
	expll2xy([name '.exp'],1);
