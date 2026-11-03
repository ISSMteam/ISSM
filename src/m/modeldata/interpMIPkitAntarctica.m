function [dataout] = interpMIPkitAntarctica(X,Y,string,varargin),
%INTERPMIPKITANTARCTICA - interpolate ISMIP7 MIPkit Antarctica datasets
%
%   Available data:
%      - mouginot_basins        (2004, imbie.org/imbie-3/drainage-basins)
%      - icemask_greene         (time series, nsidc-0794)
%      - mask_bedmachine        (2007: from BedMachine Antarctica v4.1)
%      - surface_bedmachine     (2015: from BedMachine Antarctica v4.1)
%      - thickness_bedmachine   (2015: from BedMachine Antarctica v4.1)
%      - bed_bedmachine         (2015: from BedMachine Antarctica v4.1)
%      - mask_bedmap3           (2007, Bedmap3, doi:10.5285/2d0e4791-8e20-46a3-80e4-f5f6716025d2)
%      - surface_bedmap3        (2015, Bedmap3)
%      - thickness_bedmap3      (2015, Bedmap3)
%      - bed_bedmap3            (2015, Bedmap3)
%      - geoid                  (EIGEN-6C4, Forste et al 2014)
%      - vx_mosaic              (2010, nsidc-0484)
%      - vy_mosaic              (2010, nsidc-0484)
%      - vx_timeseries          (time series, nsidc-0720)
%      - vy_timeseries          (time series, nsidc-0720)
%      - dhdt_smith             (2019, Smith et al. 2020)
%      - dhdt_cpom              (time series, CPOM, cpom.ucl.ac.uk/csopr/icesheets3)
%      - geothermal_heat_flux1  (2021, Staal et al. 2021)
%      - geothermal_heat_flux2  (2018, Lösing & Ebbing 2021)
%
%   Usage:
%      [dataout] = interpMIPkitAntarctica(X, Y, string, [ncfile])
%
%   Examples:
%      md.geometry.bed = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'bed');
%      md.geometry.bed = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'bed', '../Data/AntarcticaObsISMIP7-v1.0.nc');

verbose=1;

if nargin<4
	%List of common paths to try
	filename = 'AntarcticaObsISMIP7-v1.2.nc';
	paths = {...
		['/totten_1/ModelData/ISMIP7/MIPkit/' filename ],...
		['./Data/' filename],...
		['./' filename],...
		};

	found = 0;
	for i=1:numel(paths)
		if exist(paths{i},'file')
			mipkitnc = paths{i}; found = 1; break;
		end
	end
	if ~found
		error(['Could not find ' filename '. You can add the path to the list or provide its path as 4th argument']);
	end
end

%Non timeseries, 1km resolution (x1km, y1km)
if verbose, disp(['   -- MIPkit: loading coordinates for ' string]); end
if ismember(string, {'dhdt_smith', 'geothermal_heat_flux1', 'geothermal_heat_flux2'})
	xdata = double(ncread(mipkitnc,'x1km'));
	ydata = double(ncread(mipkitnc,'y1km'));
	istimeseries = false;

%Non timeseries, high-resolution (x, y)
elseif ismember(string, {...
		'mouginot_basins', 'mask_bedmachine', 'mask_bedmap3', ...
		'surface_bedmachine', 'thickness_bedmachine', 'bed_bedmachine', ...
		'surface_bedmap3', 'thickness_bedmap3', 'bed_bedmap3', 'geoid',...
		'vx_mosaic', 'vy_mosaic'})
	xdata = double(ncread(mipkitnc,'x'));
	ydata = double(ncread(mipkitnc,'y'));
	istimeseries = false;

%Specific time series
elseif strcmp(string, 'icemask_greene')
	xdata = double(ncread(mipkitnc,'x'));
	ydata = double(ncread(mipkitnc,'y'));
	tdata = double(ncread(mipkitnc,'greene_mask_time'));
	istimeseries = true;
elseif strcmp(string, 'dhdt_cpom')
	xdata = double(ncread(mipkitnc,'x1km'));
	ydata = double(ncread(mipkitnc,'y1km'));
	tdata = double(ncread(mipkitnc,'cpom_dhdt_time'));
	istimeseries = true;
elseif ismember(string, {'vx_timeseries', 'vy_timeseries'})
	xdata = double(ncread(mipkitnc,'x'));
	ydata = double(ncread(mipkitnc,'y'));
	tdata = double(ncread(mipkitnc,'vel_time'));
	istimeseries = true;
else
	error(['data field ''' string ''' not supported']);
end

% Prepare subset based on provided X and Y
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

if verbose, disp(['   -- MIPkit: loading ' string]); end
if istimeseries
	id1t = 1;
	id2t = numel(tdata);
	data = permute(double(ncread(mipkitnc,string,[id1x id1y id1t],[id2x-id1x+1 id2y-id1y+1 id2t],[1 1 1])), [2 1 3]);
else
	data = double(ncread(mipkitnc,string,[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
end
xdata = xdata(id1x:id2x);
ydata = ydata(id1y:id2y);

%which interpolation type are we using?
if ismember(string, {'mouginot_basins','mask_bedmachine','mask_bedmap3','icemask_greene'})
	method = 'nearest';
else
	method = 'bilinear';
end

%Go ahead and interpolate
if verbose, disp(['   -- MIPkit: interpolating ' string ' (method: ' method ')']); end
if istimeseries
	if ~isvector(X) && ~isvector(Y)
		error('X and Y should be vectors');
	end
	tdata   = date2decyear(datenum(datetime('1900-01-01')+ days(tdata)));
	dataout = zeros(length(X)+1,numel(tdata));
	for i=1:numel(tdata)
		dataout(1:end-1,i) = InterpFromGrid(xdata, ydata, data(:,:,i), X, Y);
		dataout(end,    i) = tdata(i);
	end
else
	dataout = InterpFromGrid(xdata, ydata, data, X, Y);
end
