function dataout = interpFromMEaSUREsGeotiff(X,Y,Tstart,Tend,varargin)
%interpFromMEaSUREsGeotiff: 
%	This function calls src/m/contrib/morlighem/modeldata/interpFromGeotiff.m for multiple times to load all avaliable 
%	tif data in  /totten_1/ModelData/Greenland/VelMEaSUREs/Jakobshavn_2008_2021/ within the given time period (in decimal years)
%	For some reason, each .tif file in this folder contains two sets of data, only the first dataset is useful
%
%   Usage:
%		 dataout = interpFromMEaSUREsGeotiff(X,Y,Tstart,Tend, varargin)
%
%	X, Y are the coordinates of the mesh 
%	Tstart and Tend decimal year of the start and end time
%
%   Example:
%			obsData = interpFromMEaSUREsGeotiff(md.mesh.x,md.mesh.y, tstart, tend);
%
%   Options:
%      - 'glacier':  which glacier to look for
options    = pairoptions(varargin{:});
glacier    = getfieldvalue(options,'glacier','Jakobshavn');

if strcmp(glacier, 'Jakobshavn')
	foldername = '/totten_1/ModelData/Greenland/VelMEaSUREs/Jakobshavn_2008_2021/';
elseif strcmp(glacier, 'Kangerlussuaq')
	foldername = '/totten_1/ModelData/Greenland/VelMEaSUREs/Kangerlussuaq_2006_2021/';
elseif strcmp(glacier, 'Store')
	foldername = '/totten_1/ModelData/Greenland/VelMEaSUREs/Store_2008_2021/';
elseif strcmp(glacier, 'Rink')
	foldername = '/totten_1/ModelData/Greenland/VelMEaSUREs/Rink_2008_2022/';
elseif strcmp(glacier, 'Upernavik')
	foldername = '/totten_1/ModelData/Greenland/VelMEaSUREs/Upernavik_2008_2022/';
elseif strcmp(glacier, 'Helheim')
	foldername = '/totten_1/ModelData/Greenland/VelMEaSUREs/Helheim_2008_2023/';
else
	error(['The velocity data for ', glacier, ' is not available, please download from NSIDC first.']);
end

% get the time info from file names
templist = dir([foldername,'*.meta']);
Ndata = length(templist);
dataTstart = zeros(Ndata,1);
dataTend = zeros(Ndata,1);

for i = 1:Ndata
	tempConv = split(templist(i).name, '_');
	% follow the naming convention
	dataPrefix(i) = join(tempConv(1:5), '_');
	dataTstart(i) = date2decyear(datenum(tempConv{3}));
	dataTend(i) = date2decyear(datenum(tempConv{4}));
end
disp(['  Found ', num2str(Ndata), ' records in ', foldername]);
disp(['    from ', datestr(decyear2date(min(dataTstart)),'yyyy-mm-dd'), ' to ', datestr(decyear2date(max(dataTend)),'yyyy-mm-dd') ]);


% find all the data files with Tstart<=t<=Tend
dataInd = (dataTend>=Tstart) & (dataTstart<=Tend);
disp([' For the selected period: ', datestr(decyear2date((Tstart)),'yyyy-mm-dd'), ' to ', datestr(decyear2date((Tend)),'yyyy-mm-dd'), ', there are ', num2str(sum(dataInd)), ' records' ]);

dataToLoad = dataPrefix(dataInd);
TstartToload = dataTstart(dataInd);
TendToload = dataTend(dataInd);

for i = 1:length(dataToLoad)
	dataout(i).vx = interpFromGeotiff([foldername, dataToLoad{i}, '_vx_v04.0.tif'], X, Y, 2e9);
	dataout(i).vy = interpFromGeotiff([foldername, dataToLoad{i}, '_vy_v04.0.tif'], X, Y, 2e9);
	dataout(i).vel = interpFromGeotiff([foldername, dataToLoad{i}, '_vv_v04.0.tif'], X, Y, -1);
	dataout(i).Tstart = TstartToload(i);
	dataout(i).Tend = TendToload(i);
end
