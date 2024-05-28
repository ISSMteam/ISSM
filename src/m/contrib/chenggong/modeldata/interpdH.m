function dH = interpdH(X, Y, time, ncdata)
%INTERPDH - interpolate monthly reconstructed dH onto X and Y, within the given time period
%
%	 Usage:
%		dH = interpdH(md.mesh.x, md.mesh.y, [md.timestepping.start_time, md.timestepping.final_time]);
%
%   - X, Y: coordinates of the mesh or grid
%	 - time: the starting and end point of the time series
%   - optional 4th input argument: path to the data file, by default it is the path on totten
%

if nargin < 4
	ncdata = '/totten_1/ModelData/Greenland/DHKhan/dHdt_monthly_1km.nc';
end

x = ncread(ncdata, 'x');
y = ncread(ncdata, 'y');
t = ncread(ncdata, 'time');

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

% load dH
data = ncread(ncdata, 'dHdt', [idx_min, idy_min, idt_min], [idx_max-idx_min+1, idy_max-idy_min+1, idt_max-idt_min+1], [1,1,1]);

% Convert to ice_levelset values
dH = zeros(numel(X)+1, numel(t));
dH(end,:) = t;
for i = 1:numel(t)
	dH(1:end-1, i) = InterpFromGrid(x, y, data(:,:,i)', X, Y,'nearest');
end
