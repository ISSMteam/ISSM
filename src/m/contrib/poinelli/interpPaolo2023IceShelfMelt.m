function output = interpPaolo2023IceShelfMelt(X,Y,varargin)
    %interpPaolo2023IceShelfMelt - interp melt rates from Paolo et al. 2023
    %
    %   Usage:
    %      output = interpPaolo2023IceShelfMelt(X,Y,initialDate,finalDate)
    %   
    %   Input:
    %       X               : md.mesh.x
    %       Y               : md.mesh.y
    %       initialDate     : Initial date in the form datetime(2020, 1, 1);
    %       finalDate       : Final date in the form datetime(2020, 1, 1);
    %
    % Citation: Paolo, F. S., Gardner, A. S., Greene, C. A., Nilsson, J., Schodlok, M. P., Schlegel, N.-J., and Fricker, H. A.: 
    % "Widespread slowdown in thinning rates of West Antarctic ice shelves", The Cryosphere, 17, 3409–3433, https://doi.org/10.5194/tc-17-3409-2023, 2023.

switch (oshostname())
    case {'thwaites','larsen','murdo','astrid'}
		meltrate=['/u/astrid-r1b/ModelData/RignotAntarcticaMeltRates/Ant_MeltingRate.v2.nc'];
	otherwise
		error('hostname not supported yet');
end

xdata = double(ncread(meltrate,'x'));
ydata = double(ncread(meltrate,'y'));

time = double(ncread(meltrate,'time'));
% Define the reference date
ref_date = datetime(1950, 1, 1);

% Convert the time vector to datetime format
time_datetime = ref_date + days(time);