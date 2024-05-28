function latlonoverlay(md,options)
%LATLONOVERLAY - overlay latitude and longitude lines on current figure
%
%   latstep,lonstep, in latitude and longitude degreees, between two latitudinal, longitudinal profiles.
%   color: [1 1 1] for example
%   resolution: profile resolution ( in lat,lon degrees) 
%   gap: gap (in meters) to plug lat,lon degree numbers;
%
%   Usage:
%      latlonoverlay(options)

%get options
latlon=getfieldvalue(options,'latlon');
numbering=getfieldvalue(options,'latlonnumbering','off');
latlonclick=getfieldvalue(options,'latlonclick',0);
fontsize=getfieldvalue(options,'fontsize',16);

%recover arguments (set default parameters if needed)
%1: latlon
if ~iscell(latlon),
	if ischar(latlon) & strcmpi(latlon,'on'),
		%defaults
		latstep=3; lonstep=3;
		resolution=0.1;
		color=[1 0 1];
	else return; end
else
	if length(latlon)<2
		error('latlonoverlay error message: at least 2 arguments are required, or use ''on'' option.');
	end
	if length(latlon)>3, color=latlon{4};      else color=[1 1 1]; end
	if length(latlon)>2, resolution=latlon{3}; else resolution=0.1;end
	latstep=latlon{1};
	lonstep=latlon{2};
end

%2: numbering
if ~iscell(numbering) & isnan(numbering),
	numbering=false;
else
	if ~iscell(numbering),
		if strcmpi(char(numbering),'on'),
			%defaults
			latgap=2; longap=2;
			colornumber=color;
			latangle=0; lonangle=0;
			numbering=true;
		else
			numbering=false;
		end
	else
		latgap=numbering{1}; longap=numbering{2};
		colornumber=numbering{3};
		latangle=numbering{4}; lonangle=numbering{5};
		numbering=true;
	end
end

%what are the x and y limits
xlimits=getfieldvalue(options,'xlim',xlim);
ylimits=getfieldvalue(options,'ylim',ylim);

%lat
for lat=-90:latstep:90
	longitudes=0:resolution:360;
	latitudes =lat*ones(size(longitudes));

	if md.mesh.epsg==3413,
		if lat<0, continue; end
		[x,y]=ll2xy(latitudes,longitudes,+1,45,70);
	elseif md.mesh.epsg==3031,
		if lat>0, continue; end
		[x,y]=ll2xy(latitudes,longitudes,-1, 0,71);
	else error('field md.mesh.epsg not supported yet'); end

	pos=find(x<=xlimits(2) & x>=xlimits(1) & y<=ylimits(2) & y>=ylimits(1));
	if length(pos)<=1, continue; end
	x=x(pos);y=y(pos);
	l=line(x,y,'Color',color);

	if numbering
		ind=length(x)-2*latgap;
		if (ind<=0), continue; end
		xcorner=x(ind);            ycorner=y(ind);
		xcorner2=x(max(ind-10,1)); ycorner2=y(max(ind-10,1));

		if (xcorner>xlimits(1) & xcorner<xlimits(2) & ycorner>ylimits(1) & ycorner<ylimits(2)),
			angle=mod((180)/pi*atan2((ycorner2-ycorner),(xcorner2-xcorner))+latangle,360);
			if lat<0, label=[num2str(abs(lat)) '^{\circ} S'];
			else      label=[num2str(abs(lat)) '^{\circ} N']; end
			th=text(xcorner,ycorner,label);
			set(th,'Color',colornumber,'Rotation',angle,'FontSize',fontsize,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on');

			%erase line and redraw it in two parts, to leave space for latitude number
			delete(l);
			line(x(1:ind-latgap),y(1:ind-latgap),'Color',color);hold on;
			line(x(ind+latgap:end),y(ind+latgap:end),'Color',color);
			set(gcf,'InvertHardcopy','off');
		end

	end
end

%lon
for lon=-180:lonstep:180

	if md.mesh.epsg==3413,
		latitudes =0:resolution:90;
		longitudes=lon*ones(size(latitudes));
		[x,y]=ll2xy(latitudes,longitudes,+1,45,70);
	elseif md.mesh.epsg==3031,
		latitudes =-90:resolution:0;
		longitudes=lon*ones(size(latitudes));
		[x,y]=ll2xy(latitudes,longitudes,-1, 0,71);
	else
		error('field md.mesh.epsg not supported yet'); 
	end

	pos=find(x<=xlimits(2) & x>=xlimits(1) & y<=ylimits(2) & y>=ylimits(1));
	if length(pos)<=1, continue; end
	x=x(pos);y=y(pos);
	l=line(x,y,'Color',color);

	if numbering,
		ind=length(x)-2*longap;
		if (ind<=0), continue; end
		xcorner=x(ind);            ycorner=y(ind);
		xcorner2=x(max(ind-10,1)); ycorner2=y(max(ind-10,1));

		if (xcorner>xlimits(1) & xcorner<xlimits(2) & ycorner>ylimits(1) & ycorner<ylimits(2)),
			angle=mod((180)/pi*atan2((ycorner2-ycorner),(xcorner2-xcorner))+lonangle,360);
			if lon<0, label=[num2str(abs(lon)) '^{\circ} W'];
			else      label=[num2str(abs(lon)) '^{\circ} E']; end
			th=text(xcorner,ycorner,label);
			set(th,'Color',colornumber,'Rotation',angle,'FontSize',fontsize,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on');

			%erase line and redraw it in two parts, to leave space for latitude number
			delete(l);
			line(x(1:ind-longap),y(1:ind-longap),'Color',color);hold on;
			line(x(ind+longap:end),y(ind+longap:end),'Color',color);
		end

	end
end

%Back to original limits
xlim(xlimits);
ylim(ylimits);
