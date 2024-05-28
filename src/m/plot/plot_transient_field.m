function plot_transient_field(md,options,width,i,data)
%PLOT_TRANSIENT_FIELD - plot transient results
%
%   Usage:
%      plot_transient_field(md,options,width,i,data);
%
%   See also: PLOTMODEL

%Check that they are transient results
if (~isfield(md.results,'TransientSolution')),
	error('plot_transient_field error message: no transient results in the model');
end

%Figure out the iterations to plot and check if it is possible
transient=md.results.TransientSolution;
maxiteration=size(transient,2);
steps=getfieldvalue(options,'steps',1:1:maxiteration);

if max(steps)>maxiteration | min(steps)<1,
	error(['plot_transient_field error message: problem with the steps requested, must be an interger between 0 and ' num2str(maxiteration)]);
end
subplotwidth=ceil(sqrt(length(steps)));

%Figure out the field to plot

field=getfieldvalue(options,'field','Vel');

%process mes only once
[x y z elements is2d isplanet]=processmesh(md,[],options);

%plot data for all steps
for i=1:length(steps),

	%process data and change title if needed
	[data datatype]=processdata(md,transient(steps(i)).(field),options);
	options=changefieldvalue(options,'title',[field ' at time ' num2str(transient(steps(i)).time/md.constants.yts) ' a']);

	%create plot of step i
	subplot(subplotwidth,subplotwidth,i);
	plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options)
	applyoptions(md,data,options);

end
