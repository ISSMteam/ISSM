function output = interpBedmap3(X,Y,string,method)
%INTERPBEDMAP3 - interpolate Bedmap3 data onto X and Y
%
%   available fields:
%      - surface_topography
%      - bed_uncertainty
%      - bed_topography
%      - mask (1 = grounded ice, 2 = transiently grounded ice shelf, 3 = floating ice shelf, 4 = rock)
%      - ice_thickness
%      - thickness_survey_count
%      - thickness_uncertainty

if nargin<3, string = 'bed_topography'; end
if nargin<4
	if strcmp(string,'mask') | strcmp(string,'source')
		method='nearest'; % default method
	else
		method='cubic'; % default method
	end
end

%Compatibility with BedMachine
if strcmp(string, 'bed');       string = 'bed_topography'; end
if strcmp(string, 'surface');   string = 'surface_topography'; end
if strcmp(string, 'thickness'); string = 'ice_thickness'; end

%List of common paths to try
paths = {...
	['/totten_1/ModelData/Antarctica/BedMap3/GRID/bedmap3.nc'],...
	['./bedmap3.nc'],...
	};
found = 0;
for i=1:numel(paths)
	if exist(paths{i},'file')
		ncfile = paths{i};
		found = 1;
		break;
	end
end
if ~found
	error(['Could not find Bedmap3.nc'])
end

xdata = double(ncread(ncfile,'x'));
ydata = double(ncread(ncfile,'y'));
offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
if isempty(posx), posx=numel(xdata); end
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
if isempty(posy), posy=numel(ydata); end
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

%Loading data
data  = double(ncread(ncfile,string,[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);

disp(['   -- Bedmap3: interpolating ' string]);
disp(['       -- Interpolation method: ' method]);
output = InterpFromGrid(xdata,ydata,data,double(X),double(Y),method); % now the interpolation method can be defined by the user

end
