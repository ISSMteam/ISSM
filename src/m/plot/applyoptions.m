function applyoptions(md,data,options)
%APPLYOPTIONS - apply the options to current plot
%
%   Usage:
%      applyoptions(md,data,options)
%
%   See also: PLOTMODEL, PARSE_OPTIONS

%fontsize
fontsize=getfieldvalue(options,'fontsize',14);

%fontweight
fontweight=getfieldvalue(options,'fontweight','normal');

%title
if exist(options,'title')
	titlevalue=getfieldvalue(options,'title');
	if iscell(titlevalue),
		title(titlevalue,'FontSize',fontsize,'FontWeight',fontweight);
	else
		if ~isnan(titlevalue),
			title(titlevalue,'FontSize',fontsize,'FontWeight',fontweight);
		end
	end
end

%xlabel, ylabel and zlabel
if exist(options,'xlabel');
	xlabel(getfieldvalue(options,'xlabel'),'FontSize',fontsize,'FontWeight',fontweight);
end
if exist(options,'ylabel');
	ylabel(getfieldvalue(options,'ylabel'),'FontSize',fontsize,'FontWeight',fontweight);
end
if exist(options,'zlabel');
	zlabel(getfieldvalue(options,'zlabel'),'FontSize',fontsize,'FontWeight',fontweight);
end

%xticks, yticks and zticks
if exist(options,'xtick'), set(gca,'XTick',getfieldvalue(options,'xtick')); end
if exist(options,'ytick'), set(gca,'YTick',getfieldvalue(options,'ytick')); end
if exist(options,'ztick'), set(gca,'ZTick',getfieldvalue(options,'ztick')); end

%view 
if dimension(md.mesh)==3 && ~(exist(options,'layer') || exist(options,'depthaverage'))
	view(getfieldvalue(options,'view',3));
else
	view(getfieldvalue(options,'view',2));
end

%axis
set(gca,'FontSize',getfieldvalue(options,'axisfontsize',fontsize));;
if exist(options,'axis')
	axisopts = getfieldvalue(options,'axis');
	if ischar(axisopts) & any(axisopts==' ');
		axisopts = strsplit(axisopts);
		axis(axisopts{:});
	else
		axis(axisopts);
	end
else
	if strcmp(domaintype(md.mesh),'3D'),
		if ~(exist(options,'layer') || exist(options,'depthaverage'))
			axis auto tight
		else
			axis tight equal
		end
	elseif strcmp(domaintype(md.mesh),'2Dvertical'),
		axis auto tight
	elseif strcmp(domaintype(md.mesh),'3Dsurface'),
		axis auto tight

	elseif strcmp(domaintype(md.mesh),'2Dhorizontal'),
		axis tight equal;
	else
		error('type of domain not supported');
	end
end

%box
if exist(options,'box')
	box(getfieldvalue(options,'box'));
end

%xlim, ylim and zlim
if exist(options,'xlim');
	xlim(getfieldvalue(options,'xlim'));
end
if exist(options,'ylim');
	ylim(getfieldvalue(options,'ylim'));
end
if exist(options,'zlim');
	zlim(getfieldvalue(options,'zlim'));
end

%latlon
%Must be done here (before xlim and ylim??) so that it uses the same xlim and ylim as plot_overlay
%these are changed by axis that follows
if ~strcmpi(getfieldvalue(options,'latlon','off'),'off')
	latlonoverlay(md,options);
end

%Basinzoom
if exist(options,'basin');
	basinzoom(options);
end

%Zoom
if exist(options,'zoom');
	zoom(getfieldvalue(options,'zoom',2));
end

%ShowBasins
if strcmpi(getfieldvalue(options,'showbasins','off'),'on')
	showbasins(options);
end

%Caxis
if exist(options,'caxis'),
	caxis(getfieldvalue(options,'caxis'));
end

%shading
if exist(options,'shading'),
	shading(getfieldvalue(options,'shading'));
end

%grid
if exist(options,'grid'),
	if strcmpi(getfieldvalue(options,'grid'),'on'),
		grid on;
	end
end

%colormap
c = getcolormap(options);
h = colormap(gca,c);

%wrapping
if exist(options,'wrapping'),
	if ~exist(options,'colormap'),
		h=turbo();
	end
	colormap(repmat(h,getfieldvalue(options,'wrapping',1),1));
end

%colorbar
if getfieldvalue(options,'colorbar',1)==1,
	if exist(options,'colorbarcornerposition'),
		c=colorbar(getfieldvalue(options,'colorbarcornerposition'),'peer',gca);
	elseif exist(options,'colorbarpos') & ischar(getfieldvalue(options,'colorbarpos')),
		c=colorbar(getfieldvalue(options,'colorbarpos'));
	else 
		c=colorbar('peer',gca);
	end
	set(c,'FontSize',getfieldvalue(options,'colorbarfontsize',fontsize),'YColor',getfieldvalue(options,'FontColor','k'));
	if exist(options,'wrapping')
		lim=get(c,'Ylim');
		lim=[lim(1) lim(1)+(lim(2)-lim(1))/getfieldvalue(options,'wrapping')];
		set(c,'Ylim',lim);
	end
	if exist(options,'colorbarpos') & isnumeric(getfieldvalue(options,'colorbarpos')),
		set(c,'Position',getfieldvalue(options,'colorbarpos'));
	end
	if exist(options,'log'),
		nlab=length(get(c,'YTick'));
		logvalue=getfieldvalue(options,'log');

		scaleminmax=caxis;
		Min=min(scaleminmax);
		Max=max(scaleminmax);
		set(c,'YLim',[Min Max]); % set colorbar limits
		set(c,'YTick',linspace(Min,Max,nlab));     % set tick mark locations

		labels = cell(1,nlab);
		tick_vals = linspace(Min,Max,nlab);
		tick_vals = exp(log(logvalue)*tick_vals);
		warning off MATLAB:log:logOfZero;
		for i = 1:nlab
			labels{i} = sprintf('%-3.4g',round_ice(tick_vals(i),2));
			%labels{i} = sprintf('%-.4g',round_ice(tick_vals(i),2));
		end
		warning on MATLAB:log:logOfZero;
		set(c,'YTickLabel',labels);
	end 
 	if exist(options,'cbYLim'); 
		set(c,'YLim',getfieldvalue(options,'cbYLim'));
	end
	if exist(options,'colorbartitle'),
		set(get(c,'title'),'FontSize',getfieldvalue(options,'colorbarfontsize',fontsize),'String',getfieldvalue(options,'colorbartitle'),...
			'Color',getfieldvalue(options,'FontColor','k'));
	end
	if exist(options,'colorbarYLabel'),
		set(get(c,'Ylabel'),'FontSize',getfieldvalue(options,'colorbarfontsize',fontsize),'String',getfieldvalue(options,'colorbarYLabel'),...
			'Color',getfieldvalue(options,'FontColor','k'));
	end
	if exist(options,'colorbarwidth'),
		posaxes=get(gca,'Position');
		alpha=getfieldvalue(options,'colorbarwidth',1);
		position=get(c,'Position');
		dx=position(3);
		newdx=dx*alpha;
		position(1)=position(1)+(dx-newdx)/2;
		position(3)=newdx;
		set(c,'Position',position);
		set(gca,'Position',posaxes);
	end
	if exist(options,'colorbarheight'),
		posaxes=get(gca,'Position');
		alpha=getfieldvalue(options,'colorbarheight',1);
		position=get(c,'Position');
		dy=position(4);
		newdy=dy*alpha;
		position(2)=position(2)+(dy-newdy)/2;
		position(4)=newdy;
		set(c,'Position',position);
		set(gca,'Position',posaxes);
	end
	if exist(options,'cbYTickLabel');
		tick_vals=getfieldvalue(options,'cbYTickLabel');
		if ~isnumeric(tick_vals) & strcmp(tick_vals,'on')
			tick_vals=get(c,'YTick')';
			if exist(options,'log')
				logval= getfieldvalue(options,'log');
				for i= 1:numel(tick_vals)
					tick_vals(i)= logval^(tick_vals(i));
				end
			elseif numel(tick_vals) == 3
				tick_vals=[tick_vals(1); mean(tick_vals(1:2)); tick_vals(2); ...
					mean(tick_vals(2:3)); tick_vals(3)];
				set(c,'YTick',tick_vals);
			end
		else
			if exist(options,'log')
				logvalue=getfieldvalue(options,'log');
				set(c,'YTick',log(tick_vals)./log(logvalue));
			else
				set(c,'YTick',tick_vals);
			end
		end
		labels = cell(1,numel(tick_vals));
		for i = 1:numel(tick_vals)
			labels{i} = num2str(tick_vals(i));
		end
		set(c,'YTickLabel',labels);
	end

elseif getfieldvalue(options,'colorbar',1)==0,
	colorbar('off');
else
	%do nothing
end

%area
if exist(options,'area'),
	antzoom(getfieldvalue(options,'area'));
end

%expdisp
if exist(options,'expdisp'),
	filename=(getfieldvalue(options,'expdisp'));
	style=(getfieldvalue(options,'expstyle'));
	linewidth=(getfieldvalue(options,'linewidth',1));
	for i=1:length(getfieldvalue(options,'expdisp')),
		filenamei=filename{i};
		stylei=style{i};
		if length(linewidth)==1,
			linewidthi=linewidth;
		else
			linewidthi=linewidth{i};
		end
		expdisp(filenamei,'linestyle',stylei,'linewidth',linewidthi,'multiplier',getfieldvalue(options,'unit',1));
	end
end

%shpdisp
if exist(options,'shpdisp'),
	filename=(getfieldvalue(options,'shpdisp'));
	style=(getfieldvalue(options,'shpstyle',{'r.-'}));
	linewidth=(getfieldvalue(options,'linewidth',1));
	for i=1:length(getfieldvalue(options,'shpdisp')),
		filenamei=filename{i};
		stylei=style{i};
		if length(linewidth)==1,
			linewidthi=linewidth;
		else
			linewidthi=linewidth{i};
		end
		%shpdisp(filenamei,'linestyle',stylei,'linewidth',linewidthi,'multiplier',getfieldvalue(options,'unit',1));
		shpdisp(filenamei,1,stylei,linewidthi,getfieldvalue(options,'unit',1));
	end
end
if exist(options,'contours'),

	hold on;
	contours=getfieldvalue(options,'contours');
	style=getfieldvalue(options,'contourstyle',{'-'});
	linewidth=getfieldvalue(options,'linewidth',{1});
	color=getfieldvalue(options,'contourcolor',{'r'});
	contourheight=getfieldvalue(options,'contourheight',1); 

	radius=md.solidearth.planetradius;
	ratio=1+(contourheight*1000/radius);

	
	if ~isa(contours,'cell'),
		contours={contours};
	end
	nc=length(contours);
	if ~isa(style,'cell'), error('contour style should be a cell array'); end
	if ~isa(linewidth,'cell'), error('contour line width should be a cell array'); end
	if ~isa(color,'cell'), error('contour color should be a cell array'); end

	for i=1:length(contours),
		ci=contours{i};
		if length(style)==1, sti=style{1}; else sti=style{i}; end
		if length(color)==1, coli=color{1}; else coli=color{i}; end
		if length(linewidth)==1, li=linewidth{1}; else li=linewidth{i}; end

		for j=1:length(ci),
			cijx=ci(j).x*ratio;
			cijy=ci(j).y*ratio;
			cijz=ci(j).z*ratio;

			plot3(cijx,cijy,cijz,'LineWidth',li,'LineStyle',sti,'Color',coli);
		end
	end

end

%shpdisp3d
if exist(options,'shpdisp3d'),
	filename=(getfieldvalue(options,'shpdisp3d'));
	style=(getfieldvalue(options,'shpstyle',{'r.-'}));
	linewidth=(getfieldvalue(options,'linewidth',1));
	for i=1:length(getfieldvalue(options,'shpdisp3d')),
		filenamei=filename{i};
		stylei=style{i};
		if length(linewidth)==1,
			linewidthi=linewidth;
		else
			linewidthi=linewidth{i};
		end
		shpdisp3d(filenamei,'figure',1,'style',stylei,'linewidth',linewidthi);
	end
end

%text (default value is empty, not NaN...)
if exist(options,'text');
	textstring=getfieldvalue(options,'text');
	textweight=getfieldvalue(options,'textweight','b');
	textsize=getfieldvalue(options,'textsize');
	textcolor=getfieldvalue(options,'textcolor');
	textposition=getfieldvalue(options,'textposition');
	textrotation=getfieldvalue(options,'textrotation');
	text3d=getfieldvalue(options,'text3d',0);
	for i=1:length(getfieldvalue(options,'text'));
		textstringi=textstring{i};
		textweighti=textweight{i};
		textsizei=textsize{i};
		textcolori=textcolor{i};
		textpositioni=textposition{i};
		textrotationi=textrotation{i};
		if ~text3d,
			h=text(textpositioni(1),textpositioni(2),textstringi,'FontSize',textsizei,'FontWeight',textweighti,'Color',textcolori,'Rotation',textrotationi);
		else
			h=text(textpositioni(1),textpositioni(2),textpositioni(3),textstringi,'FontSize',textsizei,'FontWeight',textweighti,'Color',textcolori,'Rotation',textrotationi);
		end
		if strcmpi(getfieldvalue(options,'textclip','on'),'on'),
			set(h,'Clipping','on'); %prevent text from appearing outside of the box
		end
	end
end

%north arrow
if exist(options,'northarrow'),
	northarrow(getfieldvalue(options,'northarrow'));
end

%curved arrow
if exist(options,'curvedarrow'),
	curvedoptions=getfieldvalue(options,'curvedarrow');
	curvedarrow(curvedoptions{:});
end

%Scale ruler
if exist(options,'scaleruler'),
	scaleruler(options);
end

%streamlines
if exist(options,'streamlines'),
	plot_streamlines(md,options);
end

%contours
if exist(options,'contourlevels'),
	plot_contour(md,data,options);
end

%coastlines
if (strcmpi(getfieldvalue(options,'coastlines','off'),'on') | ...
	strcmpi(getfieldvalue(options,'coastlines','off'),'on'))
	plot_coastlines(md.mesh,options);
end

%YTickLabel
if exist(options,'yticklabel'),
	set(gca,'YTickLabel',getfieldvalue(options,'YTickLabel'));
end

%XTickLabel
if exist(options,'xticklabel'),
	set(gca,'XTickLabel',getfieldvalue(options,'XTickLabel'));
end

%xtick
if exist(options,'xtick'),
	set(gca,'xtick',getfieldvalue(options,'xtick'));
end

%ytick
if exist(options,'ytick'),
	set(gca,'ytick',getfieldvalue(options,'ytick'));
end

%Axis positions
if exist(options,'offsetaxispos'),
	offset=getfieldvalue(options,'offsetaxispos');
	P=get(gca,'pos');
	P(1)=P(1)+offset(1);
	P(2)=P(2)+offset(2);
	P(3)=P(3)+offset(3);
	P(3)=P(4)+offset(4);
	set(gca,'pos',P);
end
if exist(options,'axispos'),
	Axis=getfieldvalue(options,'axispos');
	hold on
	set(gca,'pos',Axis);
end

%showregion
if strcmpi(getfieldvalue(options,'showregion','off'),'on'),
	%Keep pointer of main axis
	maingca=gca;
	%get inset relative position (x,y,width,height)
	insetpos=getfieldvalue(options,'insetpos',[0.02 0.70 0.18 0.18]);
	%get current plos position
	cplotpos=get(maingca,'pos');
	%compute inset position
	PosInset=[cplotpos(1)+insetpos(1)*cplotpos(3),cplotpos(2)+insetpos(2)*cplotpos(4), insetpos(3)*cplotpos(3), insetpos(4)*cplotpos(4)];
	axes('pos',PosInset);
	axis equal off
	%box off
	if md.mesh.epsg==3413,
		A=expread('/u/astrid-r1b/ModelData/Exp/GreenlandBoxFront.exp');
		[A.x A.y]=ll2xy(A.x,A.y,+1,45,70);
		A.x = A.x(1:30:end);
		A.y = A.y(1:30:end);
	elseif md.mesh.epsg==3031,
		A=expread('/u/astrid-r1b/ModelData/Exp/Antarctica.exp');
	else
		error('applyoptions error message: md.mesh.epsg not defined');
	end
	offset=3*10^4;
	Ax=[min(A.x)-offset max(A.x)+offset];
	Ay=[min(A.y)-offset max(A.y)+offset];
	%if we are zooming on a basin, don't take the mesh for the boundaries!
	if exist(options,'basin'),
		[mdx mdy]=basinzoom(options);
	elseif exist(options,'xlim') | exist(options,'ylim'),
		mdx=getfieldvalue(options,'xlim');
		mdy=getfieldvalue(options,'ylim');
	else
		mdx=[min(md.mesh.x)-offset max(md.mesh.x)+offset];
		mdy=[min(md.mesh.y)-offset max(md.mesh.y)+offset];
	end
	line(A.x,A.y,ones(size(A.x)),'color','b');
	patch([Ax(1)  Ax(2)  Ax(2)  Ax(1) Ax(1)],[Ay(1)  Ay(1)  Ay(2)  Ay(2) Ay(1)],[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1,'FaceLighting','none')
	patch([mdx(1) mdx(2) mdx(2) mdx(1)],[mdy(1) mdy(1) mdy(2) mdy(2)],ones(4,1),'EdgeColor',[0 0 0],'FaceColor','r','FaceAlpha',0.5)
	colorbar('off');
	%back to main gca
	set(gcf,'CurrentAxes',maingca)
end

%flag edges of a partition
%
% TODO:
% - Figure out how to expand string representing an object (e.g. 
%	'md.qmu.variables.thickness.partition') to get the actual value *without* 
%	using `eval` and like functions (md.('token').('token') works, but how do 
%	we access md given the string 'md').
if exist(options,'partition')
	partition=getfieldvalue(options,'partition','');
	if isvector(partition) & isnumeric(partition) % partition option is vector
		% do nothing: we already have a partition vector
	% elseif ischar(partition) & ~isempty(partition) % partition option is string
	% 	% expand string
	% 	partition=eval(partition);
	% 	class(partition);
	% 	if ~(isvector(partition) & isnumeric(partition))
	% 		error('String passed to ''partition'' option does not represent a partition vector');
	% 	end
	else
		error('If ''partition'' option is supplied, it should be a partition vector object');
	end
	mdp=md;
	[xsegments ysegments]=flagedges(mdp.mesh.elements,mdp.mesh.x,mdp.mesh.y,partition);
	xsegments=xsegments*getfieldvalue(options,'unit',1);
	ysegments=ysegments*getfieldvalue(options,'unit',1);
	color=getfieldvalue(options,'partitionedgescolor','r-');
	linewidth=getfieldvalue(options,'linewidth',1);
	hold on;
	for i=1:length(xsegments),
		if (isnumeric(color))
			h=plot(xsegments(i,:),ysegments(i,:),'Color',color,'LineWidth',linewidth);
		else
			plot(xsegments(i,:),ysegments(i,:),color,'LineWidth',linewidth);
		end
	end
end

%Scatter
if exist(options,'scatter')
	data=getfieldvalue(options,'scatter');
	hold on
	plot_scatter(data(:,1),data(:,2),data(:,3),options);
end

%backgroundcolor
set(gca,'color',getfieldvalue(options,'backgroundcolor','none'));

%lighting
if strcmpi(getfieldvalue(options,'light','off'),'on'),
	set(gca,'FaceLighting','gouraud','FaceColor','interp','AmbientStrength',0.5);
	light('Position',[0 0.1 0.1],'Style','infinite');
end

%cloud of points: 
if exist(options,'cloud'),
	field=getfieldvalue(options,'cloud');
	x=field(:,1);
	y=field(:,2);
	%unit multiplier:
	if exist(options,'unit'),
		unit=getfieldvalue(options,'unit');
		x=x*unit;
		y=y*unit;
	end
	hold on,p=plot(x,y,'k.');
	markersize=getfieldvalue(options,'markersize',12);
	color=getfieldvalue(options,'cloudcolor','k');
	set(p,'Color',color);
	set(p,'MarkerSize',markersize);
end

%========================%
%OK VERY LAST STEP: INSET|
%========================%
if exist(options,'inset'),

	%Keep pointer of main axis
	maingca=gca;
	%get inset relative position (x,y,width,height)
	insetpos=getfieldvalue(options,'insetpos',[0.56 0.55 0.35 0.35]);
	%get current plot position
	cplotpos=get(gca,'pos');

	X1=getfieldvalue(options,'insetx',xlim);
	Y1=getfieldvalue(options,'insety',ylim);

	for i=1:length(getfieldvalue(options,'insetx')),
		if length(insetpos)==4,
			insetposi=insetpos;
		else
			insetposi=insetpos{i};
		end
		PosInseti=[cplotpos(1)+insetposi(1)*cplotpos(3),cplotpos(2)+insetposi(2)*cplotpos(4), insetposi(3)*cplotpos(3), insetposi(4)*cplotpos(4)];
		%show pos
		if iscell(X1),
			X1i=X1{i};
		else
			X1i=X1;
		end
		if iscell(Y1),
			Y1i=Y1{i};
		else
			Y1i=Y1;
		end
		if strcmpi(getfieldvalue(options,'showinset','off'),'on')
			line(X1i([1 2 2 1 1]),Y1i([1 1 2 2 1]),zeros(1,5),'Color','k','LineWidth',2);
		end

		%Get current figure
		ax1=gca;

		%plot inset
		axes('pos',PosInseti);
		copyobj(get(ax1,'children'),gca);
		patch('Faces',[1 2 3 4 1],'Vertices',[X1i([1 2 2 1])' Y1i([1 1 2 2])'],'FaceColor','None','EdgeColor','k','LineWidth',2);

		%apply options
		options=removefield(options,'text',0);
		options=removefield(options,'title',0);
		options=removefield(options,'xlabel',0);
		options=removefield(options,'ylabel',0);
		options=removefield(options,'inset',0);
		options=removefield(options,'offsetaxispos',0);
		options=removefield(options,'showregion',0);
		options=changefieldvalue(options,'colorbar',0);
		options=changefieldvalue(options,'latlon','off');
		options=changefieldvalue(options,'axis','equal off');
		options=changefieldvalue(options,'xlim',X1i);
		options=changefieldvalue(options,'ylim',Y1i);
		applyoptions(md,data,options);

		%back to main gca
		set(gcf,'CurrentAxes',maingca)
	end
end

%clipping off always: 
set(gca,'Clipping',getfieldvalue(options,'clipping','on'));
