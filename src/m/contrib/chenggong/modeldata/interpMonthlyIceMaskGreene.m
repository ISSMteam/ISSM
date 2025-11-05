function icemask = interpMonthlyIceMaskGreene(X, Y, time, includedRockMask, ncdata)
%INTERPMONTHLYICEMASKGREENE - interpolate monthly reconstructed ice masks onto X and Y, within the given time period
%
%	 Usage:
%		distance = interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [md.timestepping.start_time, md.timestepping.final_time]);
%
%
%   - X, Y: coordinates of the mesh or grid
%	 - time: the starting and end point of the time series
%   - optional 4th input argument: a flag to cover rock mask by ice and merge it into the ice mask, 
%											so that there is no hole in the interior part of the domain, 
%											and the mask will only represent the levelset of the ice front 
%   - optional 5th input argument: path to the data file, by default it is the path on totten
%

% set icemask=-1 for the region with rocks
if nargin < 4
	includedRockMask = 1;
end
if nargin < 5
%	ncdata = '/totten_1/ModelData/Greenland/IceFrontsGreene/greenland_ice_masks_1972-2022_v1.nc';
	ncdata = '/totten_1/ModelData/Greenland/IceFrontsGreene/NSIDC-0793_19720915-20220215_V01.0.nc';
end

x = ncread(ncdata, 'x');
y = ncread(ncdata, 'y');
d = ncread(ncdata, 'time');
% convert t to decyear
t = date2decyear(datenum(datetime('1900-01-01')+days(d)));

offset=2;

% get x-index covers the domain
xmin=min(X(:)); xmax=max(X(:));
idx = sort(find((x>=xmin) & (x<=xmax)));
idx_min = max(idx(1)-offset, 1);
idx_max = min(idx(end)+offset, length(x));
x = x(idx_min:idx_max);

% get y-index covers the domain
ymin=min(Y(:)); ymax=max(Y(:));
idy = sort(find((y>=ymin) & (y<=ymax)));
idy_min = max(idy(1)-offset, 1);
idy_max = min(idy(end)+offset, length(y));
y = y(idy_min:idy_max);

% get time index
idt_min = max([find(t<=time(1), 1, 'last'), 1]);
idt_max = min([find(t>=time(end), 1, 'first'), length(t)]);
t = t(idt_min:idt_max);

% load icemask and rockmask from netCDF
ice = ncread(ncdata, 'ice_mask', [idx_min, idy_min, idt_min], [idx_max-idx_min+1, idy_max-idy_min+1, idt_max-idt_min+1], [1,1,1]);
rock = ncread(ncdata, 'rock_mask', [idx_min, idy_min], [idx_max-idx_min+1, idy_max-idy_min+1], [1,1]);

% merge ice and rock
if includedRockMask
	iceall = ice + rock;
else 
	iceall = ice;
end
% Convert to ice_levelset values
icemask = zeros(numel(X)+1, numel(t));
icemask(end,:) = t;
for i = 1:numel(t)
	icemask(1:end-1, i) = InterpFromGrid(x, y, double(1-2*iceall(:,:,i)'), X, Y,'nearest');
end
