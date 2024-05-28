function plot_quiver3(x,y,z,u,v,w,options),
%PLOT_QUIVER3 - 3d quiver plot with colors
%
%   to be perfected tomorrow
%
%   Usage:
%      plot_quiver3(x,y,z,u,v,w,options)
%
%   Example:
%      plot_quiver(md.mesh.x,md.mesh.y,md.mesh.z,md.initialization.vx,md.initialization.vy,md.initialization.vz,options);

%keep only non NaN elements
pos=find(~isnan(x) & ~isnan(y) & ~isnan(z) & ~isnan(u) & ~isnan(v) & ~isnan(w));
x=x(pos); y=y(pos); z=z(pos);
u=u(pos); v=v(pos); w=w(pos);

%get norm Min and Max
Norm=sqrt(u.^2+v.^2+w.^2);
Min=min(Norm);
Max=max(Norm);

%process options: scaling factor?
scalingfactor=getfieldvalue(options,'scaling',0.40);

%number of colors?
colorlevels=getfieldvalue(options,'colorlevels',NaN);
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

%set the colormap 
if numcolors==2;
	%blue and red
	c=[0 0 1;1 0 0];
elseif numcolors==3,
	%blue yellow and red
	c=[0 0 1;1 1 0;1 0 0];
else
	%use turbo
	c=colormap(turbo(numcolors));
end

%Scale data
if strcmpi(getfieldvalue(options,'autoscale','on'),'off'),
	delta=((min(x)-max(x))^2+(min(y)-max(y))^2)/numel(x);
	u=scalingfactor*sqrt(delta)*u./Norm;
	v=scalingfactor*sqrt(delta)*v./Norm;
else
	delta=((min(x)-max(x))^2+(min(y)-max(y))^2)/numel(x);
	u=scalingfactor*sqrt(delta)*u./max(Norm);
	v=scalingfactor*sqrt(delta)*v./max(Norm);
end

%loop over the number of colors
hold on
h=[];
for i=1:numcolors
	pos=find( (Norm>=levels(i)) & (Norm<=levels(i+1)) );
	hprime=quiver3(x(pos),y(pos),z(pos),u(pos),v(pos),w(pos),'Color',c(i,:),'ShowArrowHead','on','AutoScale','off');
	h=[h;hprime];
end

%take care of colorbar
if  ~strcmpi(getfieldvalue(options,'colorbar','on'),'off'),

	%build ticks
	hcb=colorbar('peer',gca,'location','EastOutside');
	ticklabel=cell(1,length(levels));
	for i=1:length(levels),
		ticklabel{i}=num2str(round_ice(levels(i),3));
	end
	tickpos=1:numcolors+1;

	%remove ticks if to many have been created
	proportion=round(length(levels)/10);
	if proportion>1,
		ticklabel=ticklabel(1:proportion:end);
		tickpos=tickpos(1:proportion:end);
	end

	%draw colorbar
	set(hcb,'YTickLabel',ticklabel,'YTick',tickpos);
	%position
	if exist(options,'colorbarpos'),
		set(hcb,'Position',getfieldvalue(options,'colorbarpos'));
	end
	%fontsize
	fontsize=getfieldvalue(options,'fontsize',14);
	set(hcb,'FontSize',fontsize);
end
