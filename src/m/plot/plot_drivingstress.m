function plot_drivingstress(md,options,width,i)
%PLOT_DRIVINGSTRESS - plot driving stress
%
%   Usage:
%      plot_drivingstress(md,options,width,i);
%
%   See also: PLOTMODEL, PLOT_UNIT, PLOT_MANAGER

%get driving stress
[sx sy s]=drivingstress(md);

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[dstress datatype]=processdata(md,s,options);
dstress=dstress/1000;

%plot mesh quivervel
subplot(width,width,i); 
plot_unit(x,y,z,elements,dstress,is2d,isplanet,datatype,options)

%apply options
options=addfielddefault(options,'title','Driving stress [kPa]');
options=addfielddefault(options,'view',2);
applyoptions(md,dstress,options);
