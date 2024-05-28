function plot_thermaltransient_results(md,options,width,i)
%PLOT_THERMALTRANSIENT_RESULTS - plot  results of a thermal transient solution
%
%   Usage:
%      plot_thermaltransient_results(md,options,width,i);
%
%   See also: PLOTMODEL

string='plot(md';
for i=1:length(md.thermaltransient_results),
	string=[string ',''data'',''thermaltransient_results(' num2str(i) ').temperature'',''view'',3,''title'',''Temperature at time ' num2str(md.thermaltransient_results(i).time) ' a'''];
end
string=[string ',''figure'',1,''colorbar#all'',''on'',''view'',3,''fontsize'',' num2str(options.fontsize) ',''fontweight'',' options.fontweight ');'];
eval(string);
clear string;

string='plot(md';
for i=2:length(md.thermaltransient_results),
	string=[string ',''data'',md.thermaltransient_results(' num2str(i) ').temperature-md.thermaltransient_results(' num2str(i-1) ').temperature,''view'',3,''title'',''Delta temperature at time ' num2str(md.thermaltransient_results(i).time) ' a'''];
end
string=[string ',''figure'',2,''colorbar#all'',''on'',''fontsize'',' num2str(options.fontsize) ',''fontweight'',' options.fontweight ');'];
eval(string);
clear string;
