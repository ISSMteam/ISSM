function [vxout vyout errxout erryout stdxout stdyout]= interpMouginotAntTimeSeries1973to2018(X,Y,T)
%INTERPMOUGINOTANTTIMESERIES1973TO2018 - interpolate observed (time series) velocities 
%
%   Inputs
%      X,Y: spatial (scatter) coordinates
%      T: time (indexed by YEAR1 (YEAR2 is optional); see below) 
%
%   Outputs
%      vxout,vyout: interpolated velocities at X,Y, for each time requested in T
%
%   Available time series:
%
%          YEAR1  YEAR2
%    1     1973   1975
%    2     1973   1984
%    3     1973   1988
%    4     1984   1988
%    5     1986   1988
%    6     1988   1990
%    7     1991   1992
%    8     1995   1996
%    9     2000   2001
%   10     2002   2003
%   11     2003   2004
%   12     2005   2006
%   13     2006   2007
%   14     2007   2008
%   15     2008   2009
%   16     2009   2010
%   17     2010   2011
%   18     2011   2012
%   19     2012   2013
%   20     2013   2014
%   21     2014   2015
%   22     2015   2016
%   23     2016   2017
%   24     2017   2018
%
%   Usage:
%      T refers to YEAR1, but the user can also use YEAR2 (e.g., the "1973" case in YEAR1).
%  
%      Then, these codes generate the same results:
%
%      [vxout vyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1986; 1991; 1995; 2000]);
%      [vxout vyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1986 1988; 1991 1992; 1995 1996; 2000 2001]);
%
%      Another examples:
%      [vxout vyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1973 1975; 1973 1988; 1991 1992; 2011 2012]);
%      [vel]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1986; 1991; 1995; 2000]);
%      [vxout vyout errxout erryout stdxout stdyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1986; 1991; 1995; 2000]);

%read data
switch (oshostname()),
	case {'ronne'}
		nc = '/home/ModelData/Antarctica/MouginotVel/ASE_TimeSeries_1973-2018.nc';
	case {'totten'}
		nc = '/totten_1/ModelData/Antarctica/MouginotVel/ASE_TimeSeries_1973-2018.nc';
	otherwise
		error('hostname not supported yet');
end

xdata = double(ncread(nc,'x'));
ydata = double(ncread(nc,'y'));
year1 = ncread(nc,'YEAR1');
year2 = ncread(nc,'YEAR2');

% get the positions related to T
if nargin==3
	% initial checks %{{{
	if ~isvector(T)
		error('Size of input T not supported!');
	end
	T=T(:);
	if size(T,2)==1 & any(T(:,1)==1973),
		disp(' ');
		disp('   Found year=1973 in T (array). Please, specify the data series using a second index.');
		disp('   Data available for 1973:');
		disp('      1973   1975');
		disp('      1973   1984');
		disp('      1973   1988');
		disp(' ');
		disp('   Usage:');
		disp('      [vxout vyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1973 1975; 1973 1988; 1991 1992; 2011 2012])');
		disp(' ');
		error('   Change input T before continuing.');
	end %}}}
	pos = [];
	for i=1:size(T,1),
		flag = (T(i,1)==year1);
		if size(T,2)==2, % ok, check both indexes (year1 and year2)
			flag = (T(i,1)==year1).*(T(i,2)==year2);
		end
		pos = [pos; find(flag)];
	end
	% check again {{{
	if length(pos)~=size(T,1) | length(unique(pos))~=length(pos),
		disp(' ');
		disp('   Time resquested does not exist in data set or is repeated!');
		disp('   Data resquested:');
		for i=1:length(T(:,1)),
			str = ['      ' int2str(T(i,1)) '   '];
			if size(T,2)==2, % ok, check both indexes (year1 and year2)
				str = [str int2str(T(i,2))];
			end
			disp(str);
		end
		disp(' ');
		disp('   Data available (24 series):');
		for i=1:length(year1),
			str = ['      ' int2str(year1(i)) '   ' int2str(year2(i))];
			disp(str);
		end
		disp(' ');
		disp('   Usage:');
		disp('      [vxout vyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1986; 1991; 1995; 2000])');
		disp('      [vxout vyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1986 1988; 1991 1992; 1995 1996; 2000 2001])');
		disp('      [vxout vyout]= interpMouginotAntTimeSeries1973to2018(md.mesh.x,md.mesh.y,[1973 1975; 1973 1988; 1991 1992; 2011 2012])');
		disp(' ');
		error('   Change input T before continuing.');
	end%}}}
elseif nargin<3,
	pos = 1:24; % all available data		
else
	error('nargin not supported yet!');
end
if nargout~=1 & nargout~=2 & nargout~=6
	error('nargout not supported!');
end


% get the spatial positions
offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

disp(['   -- Mouginot Time Series 1973 to 2018: loading velocities']);
vxdata = [];
vydata = [];
if nargout==6 % it includes ERRX, ERRY, STDX and STDY
	errxdata = [];
	errydata = [];
	stdxdata = [];
	stdydata = [];
end
for i=1:length(pos), 
	disp(['      step = ' int2str(i) '/' int2str(length(pos)) ', position = ' int2str(pos(i)) ', year = '  int2str(year1(pos(i))) ' - ' int2str(year2(pos(i)))]);
	vx = double(ncread(nc,'VX',[id1x id1y pos(i)],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]));
	vy = double(ncread(nc,'VY',[id1x id1y pos(i)],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]));
	vxdata(:,:,i) = permute(vx,[2 1 3]);
	vydata(:,:,i) = permute(vy,[2 1 3]);
	if nargout==6 % it includes ERRX, ERRY, STDX and STDY
		errx = double(ncread(nc,'ERRX',[id1x id1y pos(i)],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]));
		erry = double(ncread(nc,'ERRY',[id1x id1y pos(i)],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]));	
		stdx = double(ncread(nc,'STDX',[id1x id1y pos(i)],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]));
		stdy = double(ncread(nc,'STDY',[id1x id1y pos(i)],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]));	
		errxdata(:,:,i) = permute(errx,[2 1 3]);
		errydata(:,:,i) = permute(erry,[2 1 3]);
		stdxdata(:,:,i) = permute(stdx,[2 1 3]);
		stdydata(:,:,i) = permute(stdy,[2 1 3]);
	end
end
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);

disp(['   -- Mouginot Time Series 1973 to 2018: interpolating']);
vxout = [];
vyout = [];
if nargout==6 % it includes ERRX, ERRY, STDX and STDY
	errxout = [];
	erryout = [];
	stdxout = [];
	stdyout = [];
end
for i=1:length(pos),
	disp(['      step = ' int2str(i) '/' int2str(length(pos)) ', position = ' int2str(pos(i)) ', year = '  int2str(year1(pos(i))) ' - ' int2str(year2(pos(i)))]);
	vxout = [vxout InterpFromGrid(xdata,ydata,vxdata(:,:,i),double(X),double(Y))];
	vyout = [vyout InterpFromGrid(xdata,ydata,vydata(:,:,i),double(X),double(Y))];
	if nargout==6 % it includes ERRX, ERRY, STDX and STDY
		errxout = [errxout InterpFromGrid(xdata,ydata,errxdata(:,:,i),double(X),double(Y))];
		erryout = [erryout InterpFromGrid(xdata,ydata,errydata(:,:,i),double(X),double(Y))];
		stdxout = [stdxout InterpFromGrid(xdata,ydata,stdxdata(:,:,i),double(X),double(Y))];
		stdyout = [stdyout InterpFromGrid(xdata,ydata,stdydata(:,:,i),double(X),double(Y))];
	end
end

%return vel if only one output is requested
if nargout==1,
	vxout = sqrt(vxout.^2+vyout.^2);
end
