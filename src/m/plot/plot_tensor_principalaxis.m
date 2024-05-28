function plot_tensor_principalaxis(md,options,width,i,tensor,type,plot_options)
%PLOT_TENSOR_PRINCIPALAXIS - plot ytensor principal axis
%
%   Usage:
%      plot_tensor_principalaxis(md,options,width,i);
%
%   See also: PLOTMODEL

%prepare subplot
subplot(width,width,i); 

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

if dimension(md.mesh)==2,
	eval(['Vx=tensor.principalaxis' type(end) '(:,1); Vy=tensor.principalaxis' type(end) '(:,2);'])
	eval(['value=tensor.principalvalue' type(end) ';']);
	[Vx datatype]=processdata(md,Vx,options);
	[Vy datatype]=processdata(md,Vy,options);
	[value datatype]=processdata(md,value,options);
else
	eval(['Vx=tensor.principalaxis' type(end) '(:,1); Vy=tensor.principalaxis' type(end) '(:,2); Vz=tensor.principalaxis' type(end) '(:,3);'])
	[Vx datatype]=processdata(md,Vx,options);
	[Vy datatype]=processdata(md,Vy,options);
	[Vz datatype]=processdata(md,Vz,options);
	[value datatype]=processdata(md,value,options);
end

%take the center of each element if ~isonnode
if datatype==1,
	x=mean(x(elements'))'; y=mean(y(elements'))'; z=mean(z(elements'))';
end

%plot quivers
if dimension(md.mesh)==2,

	%density
	if exist(options,'density')
		density=getfieldvalue(options,'density');
		x=x(1:density:end);
		y=y(1:density:end);
		Vx=Vx(1:density:end);
		Vy=Vy(1:density:end);
		value=value(1:density:end);
	end

	%scaling:
	delta=((min(x)-max(x))^2+(min(y)-max(y))^2)/numel(x);
	scale=0.5/max(sqrt((Vx.^2+Vy.^2)/delta));
	Vx=scale*Vx; Vy=scale*Vy;

	pos=find(value>=0);
	q1=quiver(x(pos),y(pos),Vx(pos),Vy(pos),'Color','r','ShowArrowHead','off','AutoScale','off');
	hold on
	pos=find(value<0);
	q2=quiver(x(pos),y(pos),Vx(pos),Vy(pos),'Color','b','ShowArrowHead','off','AutoScale','off');

else
	%density
	if exist(options,'density')
		density=getfieldvalue(options,'density');
		x=x(1:density:end);
		y=y(1:density:end);
		z=z(1:density:end);
		Vx=Vx(1:density:end);
		Vy=Vy(1:density:end);
		Vz=Vz(1:density:end);
		value=value(1:density:end);
	end

	%scaling:
	delta=((min(x)-max(x))^2+(min(y)-max(y))^2)/numel(x);
	scale=0.5/max(sqrt((Vx.^2+Vy.^2)/delta));
	Vx=scale*Vx; Vy=scale*Vy; Vz=scale*Vz;

	pos=find(value>=0);
	q1=quiver3(x(pos),y(pos),z(pos),Vx(pos),Vy(pos),Vz(pos),'Color','r','ShowArrowHead','off','AutoScale','off');
	hold on
	pos=find(value<0);
	q2=quiver3(x(pos),y(pos),z(pos),Vx(pos),Vy(pos),Vz(pos),'Color','b','ShowArrowHead','off','AutoScale','off');
end

%legend
if strcmpi(type(1:6),'strain')
	legend([q1 q2],'extension','compression')
elseif strcmpi(type(1:6),'stress')
	legend([q1 q2],'compression','traction')
end

%apply options
strings=strsplit_strict(type,'_');
string=strings{1};
options=addfielddefault(options,'title',[upper(string(1)) string(2:end) ' principal axis ' type(end)]);
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
