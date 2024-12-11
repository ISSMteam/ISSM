%{
% Explain
%  load netcdf format ISSM model (~ export_netCDF).
%
% Usage
%  md = loadmodel_netcdf('Antarctica.nc');
%
% Inputs
%
% Outputs
%}
function md = loadmodel_netcdf(fname)
	warning('WARNING: this script is under developement.');

	% initialize model
	md = model;

	% get information.
	info = ncinfo(fname);

	% get groups;
	GroupNames = {info.Groups.Name};

	for i = 1:length(GroupNames)
		GroupName = GroupNames{i};
		posGroup = strcmpi(GroupName, GroupNames); 

		disp(GroupName);
		for j = 1:length(info.Groups(posGroup).Variables)
			VarName = info.Groups(posGroup).Variables(j).Name;
			md.(GroupName).(VarName) = transpose(ncread(fname,[GroupName '/' VarName]));
		end
		break;
	end

	% figure out groups!

end
