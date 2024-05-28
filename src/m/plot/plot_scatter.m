function plot_scatter(x,y,level,varargin),
%PLOT_SCATTER - scatter plot
%
%   Usage:
%      plot_scatter(x,y,level,options);
%
%   Available options:
%      'caxis'      : default is full range
%      'MarkerSize' : default is 3
%      'subset'     : only plot the indices provided
%      'Line'       : use line instead of circles
%      'Cutoff'     : cut the line if the distance between 2 points is
%                     greater than Cutoff (default is 1000)

if nargin == 4,
	options = varargin{1};
else
	options=pairoptions(varargin{:}); 
end

%check input
if numel(x)~=numel(y) | numel(x)~=numel(level),
	error('x, y and data should have the same size');
end

if exist(options,'subset'),
	pos=getfieldvalue(options,'subset');
	x=x(pos);
	y=y(pos);
	level=level(pos);
end

%Some processing
Min=min(level);
Max=max(level);
if exist(options,'caxis');
	range=getfieldvalue(options,'caxis');
	Min=min(range);
	Max=max(range);
end
Siz=length(level);
nlab=10;
%Min=0;
%Max=1300;

%OK, should we create a new colorbar for the occasion?
if isempty(findobj(gcf,'tag','TMW_COLORBAR')) && isempty(findobj(gcf,'Type','Colorbar')),
	alreadyplot=false;
else
	alreadyplot=true;
end

%generate levels
if (alreadyplot),
	phch = get(findall(gcf,'type','image','tag','TMW_COLORBAR'),{'parent'});
	if ~isempty(phch),
		h    = phch{1};
		ylim=get(h,'YLim');
	else
		%R2014b +
		h = findobj(gcf,'Type','Colorbar');
		ylim = h.Limits;
	end
	palette=turbo();%colormap();
	numcolors=size(palette,1);
	levels=round_ice(linspace(ylim(1),ylim(2),numcolors+1),2);
else
	palette=getcolormap(options);
	colormap(palette);
	numcolors=size(palette,1);
	levels=round_ice(linspace(Min,Max,numcolors+1),3);
end

colorind=ones(Siz,1);
for i=1:numcolors
	pos=find((level>=levels(i)) & (level<=levels(i+1)) );
	colorind(pos)=i;
end
colorind(find(level>levels(end)))=numcolors;

%loop over the number of colors
hold on
hp=[];
if ~exist(options,'line'),
	for i=1:numcolors
		pos=find(colorind==i);
	%	hprime=plot3(x(pos),y(pos),ones(size(x(pos))),...
		hprime=plot(x(pos),y(pos),...
			'o','MarkerSize',getfieldvalue(options,'MarkerSize',3),'MarkerEdgeColor',palette(i,:),...
			'MarkerFaceColor',palette(i,:));
		hp=[hp;hprime];
	end
else
	distances = sqrt( (x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2);
	pos=find(distances>getfieldvalue(options,'Cutoff',1000));
	x(pos,:)=NaN;
	y(pos,:)=NaN;
	for j=1:numcolors;
		pos=find(colorind==j);
		if(~isempty(pos) & pos(1)==1), pos(1)=[]; end
		if ~isempty(pos),
			tempx = [x(pos-1) x(pos) NaN(size(pos))]';
			tempy = [y(pos-1) y(pos) NaN(size(pos))]';
			line(tempx(1:end-1),tempy(1:end-1),'color',palette(j,:),'linewidth',getfieldvalue(options,'LineWidth',2));
		end
	end
end

%Stop MATLAB's default interactivity
disableDefaultInteractivity(gca);

if ~alreadyplot,
	% format the colorbar
	h    = colorbar;
	caxis([min(levels) max(levels)]);
end
