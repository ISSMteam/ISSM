function output = interpMAR(X,Y,string,loadyear)
%%%%%% %Interpolate from MAR (Fettweis et al 2017)

%Available outputs
	% -LAT
	% -LON
	% -ME					(Melt)
	% -MSK_BAM01		(Mask Bamber et al 2001 5x5 km)
	% -MSK_BAM13		(Mask Bamber et al 2013 1x1 km)
	% -MSK_MAR			(MAR 10x10 km Ice Mask)
	% -RF					(Rainfall)
	% -RU					(Runoff)
	% -SF					(Snowfall)
	% -SMB				(Surface Mass Balance (without corrections)
	% -SMBCORR			(Surface Mass Balance (with corrections)
	% -SRF_BAM01		(Bamber et al 2001 5x5km Surface height)
	% -SRF_BAM13		(Bamber et al 2013 1x1 km Surface height)
	% -SRF_MAR			(MAR 10x10 km Surface height)
	% -ST					(Surface temperature)
	% -ST					(Surface temperature (with corrections)
	% -SU					(Sublimation/Evaporation)
	% -TIME				(Time)
	% -X					(x)
	% -Y					(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==2,
	string = 'SMBCORR';
	loadyear = 2014;
end

isverbose=0;

%Define data source
ncpath='/Users/henning/issm/trunk-jpl/projects/ModelData/MAR3/MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc';

% Build grid from Bamber2001
xdata = -800000:5000:700000;
ydata = -3400000:5000:-600000;

%Convert to MAR projections
if isverbose, disp('   -- MAR: converting coordinates'); end
[LAT,LON] = xy2ll(double(X(:)),double(Y(:)),+1,45,70); %convert model mesh xy to Bamber's projection
[xMAR,yMAR] = ll2xy(LAT,LON,+1,39,71); %convert from lat/long Bamber's to MAR projection xy

if isverbose, disp('   -- MAR: loading data'); end
data  = double(ncread(ncpath,string));

%Define what year to load data for 1979-2014
startyear=1979;
year=loadyear-startyear+1;

%Get data from given year
data2=squeeze(data(:,:,year));

if isverbose, disp('   -- MAR: interpolating data to grid'); end
output = InterpFromGrid(xdata,ydata,data2',xMAR,yMAR);
output = reshape(output,size(X,1),size(X,2));

