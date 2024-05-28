function plot_tensor(md,options,width,i,type)
%PLOT_TENSOR - plot tensor components
%
%   Usage:
%      plot_tensor(md,options,width,i);
%
%   See also: PLOTMODEL

h=subplot(width,width,i); axis off; pos=get(h,'Position');

plot_options.offsetx=pos(1);
plot_options.offsety=pos(2);
plot_options.width=pos(3);
plot_options.height=pos(4);

%Figure out tensor type:
%FIXME does not work anymore
if strncmpi(type,'strain',6),
	tensor=md.results.strainrate;
elseif strncmpi(type,'stress',6),
	tensor=md.results.stress;
elseif strncmpi(type,'deviatoricstress',16),
	tensor=md.results.deviatoricstress;
else
	error('plot_tensor error message: unsupported type of tensor');
end

%Figure out type of plot being requested
if strncmpi(fliplr(type),fliplr('tensor'),6) | strcmpi(type,'strainrate') | strcmpi(type,'deviatoricstress') | strcmpi(type,'stress'),
	plot_tensor_components(md,options,width,i,tensor,type,plot_options);
	return;
elseif strncmpi(fliplr(type),fliplr('principal'),9),
	plot_tensor_principal(md,options,width,i,tensor,type,plot_options);
	return;
elseif strncmpi(fliplr(type(1:end-1)),fliplr('principalaxis'),13),
	plot_tensor_principalaxis(md,options,width,i,tensor,type,plot_options);
	return;
else
	error('plot_tensor error message: unsurported type of plot');
end
