function kml2exp(input,output)
%KML2EXP: transform kml file Argus exp file.
%
% Usage:    kmltoexp('temp.kml','temp2.exp')
%
%

%First, read polygon kml file.
structure=kml_shapefile(input);

%create exp file: 
domain=struct();
for i=1:length(structure),

	if isfield(structure,'name'),
		domain(end+1).name=structure(i).name;
	else
		domain(end+1).name='NaN';
	end

	domain(end).density=1;
	domain(end).x=structure(i).X;
	domain(end).y=structure(i).Y;
end
domain=domain(2:end);
expwrite(domain,output);
