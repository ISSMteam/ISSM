function plot_transient_results(md,options,width,i)
%PLOT_TRANSIENT_RESULTS - plot transient results
%
%   Usage:
%      plot_transient_results(md,options,width,i);
%
%   See also: PLOTMODEL

fontsize=getfieldvalue(options,'fontsize',14);
fontweight=getfieldvalue(options,'fontweight','n');

%Prepare window distribution
%Get screen geometry
mp = get(0, 'MonitorPositions');
%Build window sizes
if size(mp,1)>=2        %several monitors, use the second one
	bdwidth=mp(2,1)+5; topbdwidth=mp(2,2)+20; W=mp(2,3)/3; H=mp(2,4)/2;
else                    %only one monitor
	bdwidth=5;         topbdwidth=20;         W=mp(1,3)/3; H=mp(1,4)/2;
end
pos1=[bdwidth  H+bdwidth  W-2*bdwidth  H-bdwidth-topbdwidth];
pos2=pos1+[W 0 0 0]; pos3=pos1+[2*W 0 0 0]; pos4=pos1+[0 -H 0 0]; pos5=pos1+[W -H 0 0]; pos6=pos1+[2*W -H 0 0];
%Create windows
figure(1);close;
figure('Position',pos1); figure('Position',pos2);figure('Position',pos3);figure('Position',pos4);figure('Position',pos5);figure('Position',pos6);

string='plotmodel(md';
for i=1:length(md.results.transient),
	string=[string ',''data'',md.results.transient(' num2str(i) ').thickness,''title'',''Thickness at time ' num2str(md.results.transient(i).time) ' a'''];
end
string=[string ',''figure'',1,''colorbar#all'',''on'',''fontsize'',' num2str(fontsize) ',''fontweight'',''' fontweight ''');'];
eval(string);
clear string;

string='plotmodel(md';
for i=1:length(md.results.transient),
	string=[string ',''data'',md.results.transient(' num2str(i) ').vel,''view'',3,''title'',''Velocity at time ' num2str(md.results.transient(i).time) ' a'''];
end
string=[string ',''figure'',2,''colorbar#all'',''on'',''fontsize'',' num2str(fontsize) ',''fontweight'',''' fontweight  ''');'];
eval(string);
clear string;

if dimension(md.mesh)==3,
	string='plotmodel(md';
	for i=1:length(md.results.transient),
		string=[string ',''data'',md.results.transient(' num2str(i) ').temperature,''view'',3,''title'',''Temperature at time ' num2str(md.results.transient(i).time) ' a'''];
	end
	string=[string ',''figure'',3,''colorbar#all'',''on'',''view'',3,''fontsize'',' num2str(fontsize) ',''fontweight'',''' fontweight  ''');'];
	eval(string);
	clear string;
end

string='plotmodel(md';
for i=2:length(md.results.transient),
	string=[string ',''data'',md.results.transient(' num2str(i) ').thickness-md.results.transient(' num2str(i-1) ').thickness,''title'',''Delta thickness at time ' num2str(md.results.transient(i).time) ' a'''];
end
string=[string ',''figure'',4,''colorbar#all'',''on'',''fontsize'',' num2str(fontsize) ',''fontweight'',''' fontweight  ''');'];
eval(string);
clear string;

string='plotmodel(md';
for i=2:length(md.results.transient),
	string=[string ',''data'',md.results.transient(' num2str(i) ').vel-md.results.transient(' num2str(i-1) ').vel,''view'',3,''title'',''Delta velocity at time ' num2str(md.results.transient(i).time) ' a'''];
end
string=[string ',''figure'',5,''colorbar#all'',''on'',''fontsize'',' num2str(fontsize) ',''fontweight'',''' fontweight  ''');'];
eval(string);
clear string;

if dimension(md.mesh)==3,
	string='plotmodel(md';
	for i=2:length(md.results.transient),
		string=[string ',''data'',md.results.transient(' num2str(i) ').temperature-md.results.transient(' num2str(i-1) ').temperature,''view'',3,''title'',''Delta temperature at time ' num2str(md.results.transient(i).time) ' a'''];
	end
	string=[string ',''figure'',6,''colorbar#all'',''on'',''fontsize'',' num2str(fontsize) ',''fontweight'',''' fontweight  ''');'];
	eval(string);
	clear string;
end
