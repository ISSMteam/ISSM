function dhdt=interpSmith2020(X,Y,string,varargin)
%INTERPSMITH2020 - interpolate Smith2020 data
%
%	Available data:
%		Filtered mass-change maps, for display only, units of m(ice-equivalent)/yr:
%		ais_floating_filt
%		ais_grounded_filt
%		gris_filt
%
%		Raw mass-change maps, suitable for generation of basin-by-basin mass-change estimates, units of m(ice-equivalent)/yr:
%		ais_floating
%		ais_grounded
%		gris
%
%   Usage:
%      [dataout] = interpSmith2020(X,Y,'ais_floating_filt')

options={'ais_floating_filt','ais_grounded_filt','gris_filt','ais_floating','ais_grounded','gris'};
tf=strcmp(string,options);

if ~any(tf)
	disp('String not available!');
	disp('The options are:');
	disp(options);
	error('String not available. See message above.');
end
	
switch oshostname(),
	case {'ronne'}
		if strcmp(string,'gris_filt') | strcmp(string,'gris')
			path='/home/ModelData/Greenland/DHDTSmith/';
		else
			path='/home/ModelData/Antarctica/DHDTSmith/';
		end
	case {'totten'}
		if strcmp(string,'gris_filt') | strcmp(string,'gris')
			path='/totten_1/ModelData/Greenland/DHDTSmith/';
		else
			path='/totten_1/ModelData/Antarctica/DHDTSmith/';
		end
	case {'recruta'}
		path='/home/santos/ModelData/ICESat1_ICESat2_mass_change/';
	otherwise
		error('machine not supported yet');
end

file=strcat(path,string,'.tif');

dhdt=interpFromGeotiff(file,X,Y);

end
