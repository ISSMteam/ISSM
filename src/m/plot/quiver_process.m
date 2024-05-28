function [quivers,palette]=quiver_process(x,y,u,v,options)
%QUIVER_PROCESS - process data for color quiver plot
%
%   Usage:
%      [quivers,palette]=quiver_process(x,y,u,v,options)

%keep only non NaN elements
pos=find(~isnan(x) & ~isnan(y) & ~isnan(u) & ~isnan(v));
x=x(pos); y=y(pos);
u=u(pos); v=v(pos);

%get Norm Min and Max
Norm=sqrt(u.^2+v.^2);
Min=min(Norm);
Max=max(Norm);

%Scale data
scalingfactor=getfieldvalue(options,'scaling',0.40);
if strcmpi(getfieldvalue(options,'autoscale','on'),'off'),
	delta=((min(x)-max(x))^2+(min(y)-max(y))^2)/numel(x);
	u=scalingfactor*sqrt(delta)*u./Norm;
	v=scalingfactor*sqrt(delta)*v./Norm;
else
	delta=((min(x)-max(x))^2+(min(y)-max(y))^2)/numel(x);
	u=scalingfactor*sqrt(delta)*u./max(Norm);
	v=scalingfactor*sqrt(delta)*v./max(Norm);
end

%number of colors?
colorlevels=getfieldvalue(options,'colorlevels',30);
if isnumeric(colorlevels),
	if isnan(colorlevels),
		numcolors=30;
	else
		numcolors=colorlevels;
	end
	levels=round_ice(linspace(Min,Max,numcolors+1),2);
else
	levels=zeros(1,length(colorlevels)+2);
	levels(1)=Min;
	for i=1:length(colorlevels)
		levels(i+1)=colorlevels{i};
	end
	levels(end)=Max;
	levels=sort(unique(levels));
	numcolors=length(levels)-1;
end

%create colorind for colors
colorind=ones(length(u),1);
for i=1:numcolors
	pos=find((Norm>=levels(i)) & (Norm<=levels(i+1)) );
	colorind(pos)=i;
end
colorind(find(Norm>levels(end)))=numcolors;

%build output
quivers=struct('x',x,'y',y,'u',u,'v',v,'levels',levels,'colorind',colorind,'numcolors',numcolors);

%set the colormap 
if numcolors==2;
	%blue and red
	palette=colormap([0 0 1;1 0 0]);
elseif numcolors==3,
	%blue yellow and red
	palette=colormap([0 0 1;1 1 0;1 0 0]);
else
	%use turbo
	palette=colormap(turbo(numcolors));
end
