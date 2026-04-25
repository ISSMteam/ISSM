function latlongrid(varargin)
%LATLONGRID - make lat lon grid on figure
%
%   Usage:
%      latlongrid(varargin)
%
%   Supported options:
%      - 'linecolor': line color (default is [0 0 0])
%      - 'latlonres': lat/lon resolution (default is 0.01 deg)
%		 - 'sgn': 1 : north latitude (default is mer=45 lat=70)
%					-1 : south latitude (default is mer=0  lat=71)
%      - 'lat': series of latitude to display
%      - 'lon': series of longitude to display

%process options
options = pairoptions(varargin{:});
linecolor = getfieldvalue(options,'linecolor',[0 0 0]);
latlonres = getfieldvalue(options,'latlonres',0.01);
sgn = getfieldvalue(options,'sgn', 1);  

%Some parameters
numpts = 1000;

%Get plot extrema
XLIM=xlim(); YLIM=ylim();
xmin = XLIM(1); xmax = XLIM(2);
ymin = YLIM(1); ymax = YLIM(2);

%Get lat/lon "extremes" (warning: they may be INSIDE the plot, which is why we use meshgrid)
[allX   allY]   = meshgrid(linspace(xmin,xmax,numpts),linspace(ymin,ymax,numpts));
if sgn == 1
	disp('== Northern hemisphere ==')
	delta = 45;
	slat = 70;
else
	delta = 0;
	slat = 71;
	disp('== Southern hemisphere ==')
end
[alllat alllon] = xy2ll(allX(:),allY(:),sgn,delta,slat);

%Get lat list (easy)
latmin = min(alllat); latmax = max(alllat);
latrange = latmax - latmin;
dlat = dtick(latrange);
latlist = (dlat*ceil(latmin/dlat)):dlat:latmax;
latlist_long = latmin:latlonres:latmax;

%Get lon list (hard because it jumps from -180 to 180)
sortedlons = sort(alllon);
pos = find(diff(sortedlons)>10);
if isempty(pos)
	lonmin = min(alllon); lonmax = max(alllon);
	lonrange = lonmax-lonmin;
	dlon = dtick(lonrange);
	lonlist = [(dlon*ceil(lonmin/dlon)):dlon:lonmax];
	lonlist_long = lonmin:latlonres:lonmax;
else
	assert(numel(pos)==1);
	lonmin1 = sortedlons(1);     lonmax1 = sortedlons(pos);
	lonmin2 = sortedlons(pos+1); lonmax2 = 180;%sortedlons(end);
	lonrange = lonmax1-lonmin1 + lonmax2-lonmin2;
	dlon = dtick(lonrange);
	lonlist = [(dlon*ceil(lonmin2/dlon)):dlon:lonmax2  (dlon*ceil(lonmin1/dlon)):dlon:lonmax1];
	lonlist_long = [lonmin2:latlonres:lonmax2 lonmin1:latlonres:lonmax1];
end

%Overwrite if need be
if exist(options,'lat'); latlist = getfieldvalue(options,'lat'); end
if exist(options,'lon'); lonlist = getfieldvalue(options,'lon'); end

%draw lats and lon
hold on
tic
for thislat=latlist
	lon=lonlist_long;
	lat=thislat*ones(size(lon));
	[xll,yll]=ll2xy(lat,lon,sgn,delta,slat);
	pos = find(xll>xmin & xll<xmax & yll>ymin & yll<ymax);
	if isempty(pos); continue; end
	pos = find(xll<xmin | xll>xmax | yll<ymin | yll>ymax);
	xll(pos) = NaN;
	yll(pos) = NaN;
	p=plot(xll,yll);set(p,'Color',linecolor,'LineWidth',1,'LineStyle','-')
%	uistack(p, 'bottom')
end
for thislon=lonlist
	lat=latlist_long;
	lon=thislon*ones(size(lat));
	[xll,yll]=ll2xy(lat,lon,sgn,delta,slat);
	pos = find(xll>xmin & xll<xmax & yll>ymin & yll<ymax);
	if isempty(pos); continue; end
	pos = find(xll<xmin | xll>xmax | yll<ymin | yll>ymax);
	xll(pos) = NaN;
	yll(pos) = NaN;
	p=plot(xll,yll);set(p,'Color',linecolor,'LineWidth',1,'LineStyle','-')
%	uistack(p, 'bottom')
end

%Now change axis ticks Let's assume we have lat on x axis, but that may change...
ax = gca;

%x-axis compute how many ticks we would have it we display longitude
[xlat xlon] = xy2ll(linspace(xmin,xmax,numpts),linspace(ymin,ymin,numpts),sgn,delta,slat);
sortedlons = sort(xlon);
pos = find(diff(sortedlons)>10);
if isempty(pos)
	lonmin = min(xlon); lonmax = max(xlon);
	lontick = lonlist(find(lonlist>lonmin & lonlist<lonmax));
else
	assert(numel(pos)==1);
	lonmin1 = sortedlons(1);     lonmax1 = sortedlons(pos);
	lonmin2 = sortedlons(pos+1); lonmax2 = 180;%sortedlons(end);
	lontick = [lonlist(fliplr(find(lonlist>lonmin1 & lonlist<lonmax1))) ...
		180 ...
		lonlist(fliplr(find(lonlist>lonmin2 & lonlist<lonmax2)))];
end

%x-axis compute how many ticks we would have it we display latitudes
latmin = min(xlat); latmax = max(xlat);
lattick = latlist(find(latlist>latmin & latlist<latmax));

if numel(lontick) > numel(lattick)
	%x-axis -> lonngitude
	xtick   = interp1(xlon,linspace(xmin,xmax,numpts),lontick,'linear','extrap');
	labels  = compose(['%g' char(176) 'E'],lontick);
	pos = find(lontick<0);
	labels(pos) = compose(['%g' char(176) 'W'],abs(lontick(pos)));
	if numel(xtick)>1 & xtick(2)-xtick(1)<0
		xtick   = fliplr(xtick);
		labels  = fliplr(labels);
	end
	ax.XTick      = xtick;
	ax.XTickLabel = labels;

	%y-axis -> latitude
	[ylat ylon] = xy2ll(linspace(xmin,xmin,numpts),linspace(ymin,ymax,numpts),sgn,delta,slat);
	latmin = min(ylat); latmax = max(ylat);
	lonmin = min(ylon); lonmax = max(ylon);
	lattick = latlist(find(latlist>latmin & latlist<latmax));
	ytick   = interp1(ylat,linspace(ymin,ymax,numpts),lattick);
	labels  = compose(['%g' char(176) 'N'],lattick);
	pos = find(lattick<0);
	labels(pos) = compose(['%g' char(176) 'S'],abs(lattick(pos)));
	if ylat(2)-ylat(1)<0
		ytick   = fliplr(ytick);
		labels  = fliplr(labels);
	end
	ax.YTick      = ytick;
	ax.YTickLabel = labels;
else
	%x-axis -> latitude
	lattick = latlist(find(latlist>latmin & latlist<latmax));
	xtick   = interp1(xlat,linspace(xmin,xmax,numpts),lattick);
	labels  = compose(['%g' char(176) 'N'],lattick);
	pos = find(lattick<0);
	labels(pos) = compose(['%g' char(176) 'S'],abs(lattick(pos)));
	if xlat(2)-xlat(1)<0
		xtick   = fliplr(xtick);
		labels  = fliplr(labels);
	end
	ax.XTick      = xtick;
	ax.XTickLabel = labels;

	%y-axis -> longitude
	[ylat ylon] = xy2ll(linspace(xmin,xmin,numpts),linspace(ymin,ymax,numpts),sgn,delta,slat);
	latmin = min(ylat); latmax = max(ylat);
	lonmin = min(ylon); lonmax = max(ylon);
	lontick = lonlist(find(lonlist>lonmin & lonlist<lonmax));
	ytick   = interp1(ylon,linspace(ymin,ymax,numpts),lontick);
	labels  = compose(['%g' char(176) 'E'],lontick);
	pos = find(lontick<0);
	labels(pos) = compose(['%g' char(176) 'W'],abs(lontick(pos)));
	if ylon(2)-ylon(1)<0
		ytick   = fliplr(ytick);
		labels  = fliplr(labels);
	end
	ax.YTick      = ytick;
	ax.YTickLabel = labels;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta = dtick(range)
	%Tick intervals
	m = 10^floor(log10(range));
	p = ceil(range/m);
	if p <= 1,     delta = .1*m;
	elseif p == 2, delta = .2*m;
	elseif p <= 5, delta = .5*m;
	else           delta = m;
	end
