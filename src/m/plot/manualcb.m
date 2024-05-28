function h=manualcb(zmin,zmax,cmap,varargin)
%MANUALCB - custom colorbar
%
%   Usage:
%      manualcb(min,max,colormap,options)
%
%   Available options:
%      - 'fontsize'    : default is 12
%      - 'fontcolor'   : default is 'k'
%      - 'smallbars'   : bars next to each tick (default is false)
%      - 'position'    : colorbar position in normalized units
%      - 'orientation' : 'vertical' (default) or 'horizontal'
%      - 'title'       : colorbar title
%      - 'xlabel'      : colorbar x-label
%      - 'ylabel'      : colorbar y-label
%      - 'tick'        : specified values of tick labels
%      - 'ticksep'     : spacing between ticks
%      - 'inverttickposition' : put ticks on the left hand side for vertical cb

%check inputs
if nargin<3,
	help manualcb
	error('bad usage');
end
if zmin>zmax,
	error('zmin should be smaller than zmax');
end

%Get plot axes
mainaxes = gca;

%process options
options = pairoptions(varargin{:});
if exist(options,'tick') & exist(options,'ticksep'),
	error('only one of tick or ticksep can be specified');
end
fontsize  = getfieldvalue(options,'fontsize',12);
fontcolor = getfieldvalue(options,'fontcolor','k');
smallbars = getfieldvalue(options,'smallbars',false);

%Colorbar position
if ~exist(options,'position'),
	position = plotboxpos;
	xstart   = position(1)+position(3)+0.01;
	ystart   = position(2);
	width    = .02;
	height   = position(4);
else
	position = getfieldvalue(options,'position');
	xstart = position(1);
	ystart = position(2);
	width  = position(3);
	height = position(4);
end
axes('Units','normalized','Position',[xstart ystart width height],'XTickLabel','','YTickLabel','','Visible','on');
xlim([0 1]);
ylim([0 1]);

%Prepare ticks
if ~exist(options,'log'),
	deltaz = getfieldvalue(options,'ticksep',dtick(zmax-zmin));
	ztick  = getfieldvalue(options,'tick',(deltaz*ceil(zmin/deltaz)):deltaz:zmax);
	if (any(ztick>zmax) | any(ztick<zmin)),
		error('one or more specified tick values falls outside of [zmin,zmax]');
	end
	ytick  = (ztick-zmin)/(zmax-zmin);
else
	%old method
	ztick = getfieldvalue(options,'tick',round( logspace(log10(zmin),log10(zmax),8) ));
	ytick = linspace(0,1,numel(ztick));

	%New method
	test=logspace(-10,10,21);
	if zmax<0 %Negative log (RARE!!)
		pos=find(test>=-zmax & test<=-zmin);
		ztick= -test(pos);
		ytick= (log(ztick) - log(zmin))/(log(zmax) - log(zmin));
		ztick
		ytick
	else
		pos=find(test>=zmin & test<=zmax);
		ztick= test(pos);
		ytick= (log(ztick) - log(zmin))/(log(zmax) - log(zmin));
	end
end

%Display colorbar
hold on
numcolors=size(cmap,1);
if 0,
	%disappears somtimes
	if strcmpi(getfieldvalue(options,'orientation','vertical'),'vertical'),
		image_rgb = ind2rgb(repmat((1:numcolors)',1,10),cmap);
	else
		image_rgb = ind2rgb(repmat((1:numcolors),10,1),cmap);
	end

	imagesc([0 1],[0 1],image_rgb);
else
	%Creates triangles when exported as pdf
	if strcmpi(getfieldvalue(options,'orientation','vertical'),'vertical'),
		for i=1:numcolors,
			patch([0,0,1,1],[(i-1)/numcolors,i/numcolors,i/numcolors,(i-1)/numcolors],0,'FaceColor',cmap(i,:),'Clipping','off','EdgeColor','none')
		end
	else
		for i=1:numcolors,
			patch([(i-1)/numcolors,i/numcolors,i/numcolors,(i-1)/numcolors],[0,0,1,1],0,'FaceColor',cmap(i,:),'Clipping','off','EdgeColor','none')
		end
	end
end
patch([0,0,1,1],[0,1,1,0],fontcolor,'FaceColor','none','Clipping','off','Edgecolor',fontcolor)

%Add ticks
if strcmpi(getfieldvalue(options,'orientation','vertical'),'vertical'),
	%Use FOR LOOP otherwise numbers are not correcly centered
	if getfieldvalue(options,'inverttickposition',0)==1,
		for i=1:length(ytick), text(-0.5,ytick(i),num2str(ztick(i)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',fontsize,'Color',fontcolor); end
	else
		for i=1:length(ytick), text(1.5,ytick(i),num2str(ztick(i)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',fontsize,'Color',fontcolor); end
	end
	if smallbars,
		for i=1:numel(ztick)
			patch([0.8 1.0],[ytick(i) ytick(i)],fontcolor,'Edgecolor',fontcolor)
			patch([0.0 0.2],[ytick(i) ytick(i)],fontcolor,'Edgecolor',fontcolor)
		end
	end
else
	%Use FOR LOOP otherwise numbers are not correcly centered
	for i=1:length(ytick), text(ytick(i),-0.5,num2str(ztick(i)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',fontsize,'Color',fontcolor); end
	if smallbars,
		for i=1:numel(ztick)
			patch([ytick(i) ytick(i)],[0.8 1.0],[ytick(i) ytick(i)],fontcolor,'Edgecolor',fontcolor)
			patch([ytick(i) ytick(i)],[0.0 0.2],[ytick(i) ytick(i)],fontcolor,'Edgecolor',fontcolor)
		end
	end
end

if exist(options,'title'),
	title(getfieldvalue(options,'title'),'FontSize',getfieldvalue(options,'titlefontsize',fontsize),'Color',fontcolor);
end
if exist(options,'ylabel'),
	if strcmpi(getfieldvalue(options,'orientation','vertical'),'horizontal'),
		th=title(getfieldvalue(options,'title'),'FontSize',fontsize,'Color',fontcolor);
		set(th,'Position',[ytick(end)+0.075,-0.3]);
	elseif getfieldvalue(options,'inverttickposition',0)==1,
		text(1.9,.7,getfieldvalue(options,'ylabel'),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',fontsize,'Color',fontcolor,'rotation',90);
	else
		ylabel(getfieldvalue(options,'ylabel'),'FontSize',fontsize,'Color',fontcolor);
	end
end
	
%Back to original axes
h=gca;
if getfieldvalue(options,'showregion',0)==0,
	%Do it this way in order to preserve the figure visibility
	set(gcf,'CurrentAxes',mainaxes);
end

function delta = dtick(range)
%Tick intervals
m = 10^floor(log10(range));
p = ceil(range/m);
if p <= 1,     delta = .1*m;
elseif p == 2, delta = .2*m;
elseif p <= 5, delta = .5*m;
else           delta = m;
end
