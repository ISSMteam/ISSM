function [dh_raw_out dh_fil_out T_out] = interpPaolo2015(X,Y,T,method)
%INTERPPAOLO2015 - interpolate observed (time series) height change [m]
%
%   Time series are average height changes [m] with respect to 1994 every three months (72 time steps)
%
%   Inputs
%      X,Y: spatial (scatter) coordinates
%      T: time (see below the available years) 
%      ATTENTION: it is assumed that X and Y come in Polar Stereographic Projection (Std Latitude: 71S Meridian: 0E)
%
%   Outputs
%      dh_raw_out: interpolated raw height change at X,Y, for each time requested in T
%      dh_fil_out: interpolated filtered height change at X,Y, for each time requested in T
%      T_out: time related to dh_raw_out and dh_fil_out (see below)
%
%   Available time series:
%
% 		 1		 1994.038
% 		 2		 1994.285
% 		 3		 1994.534
% 		 4		 1994.786
% 		 5		 1995.038
% 		 6		 1995.285
% 		 7		 1995.534
% 		 8		 1995.786
% 		 9		 1996.038
% 		10		 1996.287
% 		11		 1996.536
% 		12		 1996.787
% 		13		 1997.038
% 		14		 1997.285
% 		15		 1997.534
% 		16		 1997.786
% 		17		 1998.038
% 		18		 1998.285
% 		19		 1998.534
% 		20		 1998.786
% 		21		 1999.038
% 		22		 1999.285
% 		23		 1999.534
% 		24		 1999.786
% 		25		 2000.038
% 		26		 2000.287
% 		27		 2000.536
% 		28		 2000.787
% 		29		 2001.038
% 		30		 2001.285
% 		31		 2001.534
% 		32		 2001.786
% 		33		 2002.038
% 		34		 2002.285
% 		35		 2002.534
% 		36		 2002.786
% 		37		 2003.038
% 		38		 2003.285
% 		39		 2003.534
% 		40		 2003.786
% 		41		 2004.038
% 		42		 2004.287
% 		43		 2004.536
% 		44		 2004.787
% 		45		 2005.038
% 		46		 2005.285
% 		47		 2005.534
% 		48		 2005.786
% 		49		 2006.038
% 		50		 2006.285
% 		51		 2006.534
% 		52		 2006.786
% 		53		 2007.038
% 		54		 2007.285
% 		55		 2007.534
% 		56		 2007.786
% 		57		 2008.038
% 		58		 2008.287
% 		59		 2008.536
% 		60		 2008.787
% 		61		 2009.038
% 		62		 2009.285
% 		63		 2009.534
% 		64		 2009.786
% 		65		 2010.038
% 		66		 2010.285
% 		67		 2010.534
% 		68		 2010.786
% 		69		 2011.038
% 		70		 2011.285
% 		71		 2011.534
% 		72		 2011.786
%
%
%   Usage:
%      % Get data at specific time:
%      % In this example, T_out = [2006.038; 2007.038; 2008.038].
%      [dh_raw_out dh_fil_out T_out] = interpPaolo2015(md.mesh.x, md.mesh.y, [2006.038; 2007.038; 2008.038]);
%
% 		 % Get all data in the provided years:
%      % In this example, T_out = [2006.038; 2006.285; 2006.534; 2006.786; 2007.038; 2007.285; 2007.534; 2007.786]. 
%      [dh_raw dh_fil T_out] = interpPaolo2015(md.mesh.x, md.mesh.y, [2006; 2007]);
%
% 		 % Get all data set:
%      % In this example, T_out = [1994.038; ... ; 2011.786]. (all available time)
%      [dh_raw dh_fil T_out] = interpPaolo2015(md.mesh.x, md.mesh.y);
%
%
%   Info from ice_shelf_dh_v1.h5:	
%      The dataset is a rectangular grid (480 points in x, 80 points in y) with x- and y-axes being longitude and latitude, respectively.
%      Longitude/latitude coordinates refer to the center of the grid cells.
%      The grid has a resolution of lon x lat: 0.75 x 0.25 deg (~27 km at latitude -71).
%
%
%   Data are (grids in HDF5, ice_shelf_dh_v1.h5):
%      time         : time coordinate [year; 72 values at 3-month time step]
%      lon          : x-coordinate [degrees east; range 0/360]
%      lat          : y-coordinate [degrees north; range -82/-62]
%      height_raw   : Raw time series of height change [m]     
%      height_filt  : Filtered time series of height change [m]
%      height_err   : 2-standard-error time series [m]
%

if nargin>4 | nargin<2,
	error('nargin not supported yet!');
end

% read data
switch (oshostname()),
	case {'ronne'}
		h5 = '/home/ModelData/Antarctica/Paolo2015/ice_shelf_dh_v1.h5';
	otherwise
		error('hostname not supported yet');
end

disp(['   -- Paolo''s Time Series 1994 to 2012: loading data set']);
t_data = h5read(h5,'/time');
lat_data = h5read(h5,'/lat');
lon_data = h5read(h5,'/lon');
dh_raw_data = h5read(h5,'/height_raw');
dh_fil_data = h5read(h5,'/height_filt');

% set interpolation method
if nargin<4,
	method = 'linear'; % default method
end

% get the positions related to T
if nargin<3,
	pos = 1:length(t_data); % all available data		
else
	% initial check %{{{
	if size(T,2)>1 | size(T,1)<1 | size(T,2)<1,
		error('Size of input T not supported!');
	end 
	if size(X,1)>1 & size(X,2)>1
		error('Size of input X not supported! X and Y should be vectors');
	end
	%}}}
	% Loop over T
	pos = [];
	epsilon = 5e-4;
	for i=1:length(T),
		% find specific time
		flag = (T(i)-epsilon<t_data & T(i)+epsilon>t_data);
		if ~any(flag), 
			% ok, find the time related to the requested year
			flag = (T(i)==floor(t_data));
		end
		if ~any(flag)
			error(['requested time (' num2str(T(i)) ') not found in data set'])
		end
		pos = [pos; find(flag)];
	end
	% Check if there is repeated positions
	posunique = unique(pos);
	if length(posunique)~=length(pos),
		disp('   WARNING: found repeated positions in requested time');
	end
end

% convert x/y to lat/lon:
[LAT, LON] = xy2ll(X,Y,-1); % attention: it is assumed that X and Y comes in Polar Stereographic Projection (Std Latitude: 71S Meridian: 0E)
posLON = find(LON<0);
LON(posLON) =360+LON(posLON);

disp(['   -- Paolo''s Time Series 1994 to 2012: interpolating in Lat/Long grid']);
dh_raw_out = [];
dh_fil_out = [];
for i=1:length(pos),
	disp(['      step = ' int2str(i) '/' int2str(length(pos)) ', position = ' int2str(pos(i)) ', year = '  num2str(t_data(pos(i)))]);
	dh_raw_out = [dh_raw_out InterpFromGrid(lat_data(1,:),lon_data(:,1),dh_raw_data(:,:,pos(i)),LAT,LON,method)];
	dh_fil_out = [dh_fil_out InterpFromGrid(lat_data(1,:),lon_data(:,1),dh_fil_data(:,:,pos(i)),LAT,LON,method)];
end

T_out = t_data(pos);
