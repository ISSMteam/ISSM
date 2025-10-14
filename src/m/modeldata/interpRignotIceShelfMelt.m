function output = interpRignotIceShelfMelt(X,Y,string)
%INTERPRIGNOTICESHELFMELT - interp melt rates from Rignot et al. 2013
%
%   Usage:
%      output = interpRignotIceShelfMelt(X,Y)

switch (oshostname())
	case {'ronne'}
		rignotmelt='/home/ModelData/Antarctica/RignotMeltingrate/Ant_MeltingRate.nc';
	case {'totten'}
		rignotmelt='/totten_1/ModelData/Antarctica/RignotMeltingrate/Ant_MeltingRate.nc';
	case {'amundsen.thayer.dartmouth.edu'}
		rignotmelt='/local/ModelData/AntarcticMeltRignot/Ant_MeltingRate.nc';
	case {'thwaites','larsen','murdo','astrid'}
		rignotmelt=['/u/astrid-r1b/ModelData/RignotAntarcticaMeltRates/Ant_MeltingRate.v2.nc'];
	otherwise
		error('hostname not supported yet');
end

if nargin==2,
	string = 'melt_actual';
end

disp(['   -- Rignot Ice Shelf Melt: loading ' string]);
xdata = double(ncread(rignotmelt,'xaxis'));
ydata = double(ncread(rignotmelt,'yaxis'));

disp(['   -- Rignot Ice Shelf Melt: loading' string]);
data  = double(ncread(rignotmelt,string))';

disp(['   -- Rignot Ice Shelf Melt: interpolating ' string]);
output = InterpFromGrid(xdata,ydata,data,X(:),Y(:));
output = reshape(output,size(X,1),size(X,2));
