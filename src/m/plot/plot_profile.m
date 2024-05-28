function plot_profile(md,data,options,nlines,ncols,ii)
%PLOT_SECTION - plot a given field on a profile
%
%   Usage:
%      plot_profile(md,data,options,nlines,ncols,i)
%
%   See also: PLOTMODEL

%process model
[x_m y_m z_m elements_m is2d isplanet]=processmesh(md,[],options);
if is2d, error('only 3d model supported'); end

%Get number of curves and generate random colors
numcurves=size(data,2);
colorm=getfieldvalue(options,'colormap','lines');
color=eval([ colorm '(numcurves);']);
options=removefield(options,'colormap',0); %back to default colormap

%Get coordinates
location=getfieldvalue(options,'profile');
if ~isnumeric(location) | numel(location)~=2,
	error('location provided not supported (should be [x y])');
end
xprof=location(1);
yprof=location(2);

%Loop over number of curves
for i=1:numcurves,

	%Process data
	[datai datatype]=processdata(md,data(:,i),options);

	%resolution[z,data_interp]=ProfileValues(md,datai,xprof,yprof,resolution);
	if exist(options,'resolution'),
		resolution=getfieldvalue(options,'resolution');
	else %Default resolution
		resolution=[100];
		disp(['plot_profile warning: no resolution specified, using default: ' num2str(resolution) ]);
	end

	%Compute profile value
	[z,data_interp]=ProfileValues(md,datai,xprof,yprof,resolution);

	%plot profile
	subplot(nlines,ncols,ii)
	plot(data_interp,z,'color',color(i,:),'LineWidth',getfieldvalue(options,'linewidth',1),'LineStyle','-');
	hold on;
end

%apply options
options=addfielddefault(options,'title','Profile');
options=addfielddefault(options,'colorbar',0);
options=addfielddefault(options,'ylabel','z');
options=addfielddefault(options,'view',2);
applyoptions(md,[],options);
