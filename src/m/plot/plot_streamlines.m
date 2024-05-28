function plot_streamlines(md,options)
%PLOT_STREAMLINES - plot stream lines on a figure
%
%   Usage:
%      plot_streamlines(md,options)

%process data and model
[x y z index is2d isplanet]=processmesh(md,[],options);
[u datatype]=processdata(md,md.initialization.vx,options);
[v datatype]=processdata(md,md.initialization.vy,options);

%some checks
if ~is2d,
	disp('plot_streamlines error: streamlines option not supported for 3d plots. Project on a layer')
	return
end

%initialize flowpath
streamlines=getfieldvalue(options,'streamlines');
edgecolor=getfieldvalue(options,'streamlines_color','y');
if ischar(streamlines) & strcmpi(streamlines,'on');
	streamlines=60;
end
if iscell(streamlines)
	x0=[]; y0=[];
	for i=1:size(streamlines,2)
		coord=streamlines{i};
		x0=[x0;coord(1)]; y0=[y0;coord(2)];
	end
else
	x0=x(1:ceil(length(x)/streamlines):end);
	y0=y(1:ceil(length(x)/streamlines):end);
end

%Get flow lines
flowpath=flowlines(index,x,y,u,v,x0,y0);

%plot
hold on
for i=1:length(flowpath)
	patch('Xdata',[flowpath(i).x;NaN],'Ydata',[flowpath(i).y;NaN],'facecolor','none','edgecolor',edgecolor);
end
