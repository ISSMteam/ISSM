function dhdt=interpShepherd2019(X,Y,string,varargin)
%INTERPSHEPHERD2019 - interpolate Shepherd2019 data
%
%   Available data:
%      1.  dhdt_1992_1996
%      2.  dhdt_1997_2001
%      3.  dhdt_2002_2006
%      4.  dhdt_2007_2011
%      5.  dhdt_2012_2016
%      6.  dhdt_1992_2017
%      7.  uncert_1992_1996
%      8.  uncert_1997_2001
%      9.  uncert_2002_2006
%      10.  uncert_2007_2011
%      11.  uncert_2012_2016
%      12.  uncert_1992_2017
%
%   Usage:
%      [dataout] = interpShepherd2019(X,Y,'dhdt_1992_2017')

options={'dhdt_1992_1996','dhdt_1997_2001','dhdt_2002_2006','dhdt_2007_2011','dhdt_2012_2016','dhdt_1992_2017',...
			'uncert_1992_1996','uncert_1997_2001','uncert_2002_2006','uncert_2007_2011','uncert_2012_2016','uncert_1992_2017'};
tf=strcmp(string,options);

if ~any(tf)
	disp('String not available!');
   disp('The options are:');
   disp(options);
   error('String not available. See message above.');
end

switch oshostname(),
   case {'ronne'}
		nc='/home/ModelData/Antarctica/DHDTShepherd/antarctic_dhdt_5km_grid_1992_2017.nc';
	case {'totten'}
		nc='/totten_1/ModelData/Antarctica/CPOM_dhdt/antarctic_dhdt_5km_grid_1992_2017.nc';
	case {'recruta'}
		nc='/home/santos/ModelData/CPOM_dhdt_shepherd_2019/antarctic_dhdt_5km_grid_1992_2017.nc';
	otherwise
      error('machine not supported yet');
end

if nargin==3,
   method='linear';% default
else
   method=varargin{1};
end

xdata=double(ncread(nc,'x'));
ydata=double(ncread(nc,'y'));
data=double(ncread(nc,string))';

dhdt=InterpFromGrid(xdata,ydata,data,X,Y,method);

end
