function output = interpBamber2013(X,Y,string),
%INTERPBAMBER2013 - interpolate Bamber 2013 data
%
%   Available data:
%      BedrockElevation
%      SurfaceElevation
%      IceThickness
%      SurfaceRMSE
%      BedrockError
%      LandMask (Land mask, 0=ocean, 1=land, 2=ice sheet, 3=non-Greenlandic land, 4=ice shelf)
%      NumberAirbornePoints
%      Geoid
%      BedrockChangeMask
%      IceShelfSourceMask
%      BedrockElevation_unprocessed
%      IceThickness_unprocessed
%      BathymetryDataMask

switch oshostname(),
	case {'murdo','thwaites','astrid'}
		bamber2013nc='/u/astrid-r1b/morlighe/issmjpl/proj-morlighem/DatasetGreenland/Data/Bamber2013/Greenland_bedrock_topography_V3.nc';
	case {'ronne'}
		bamber2013nc='/home/ModelData/Greenland/Bamber2013/Greenland_bedrock_topography_V3.nc';
	case {'totten'}
		bamber2013nc='/totten_1/ModelData/Greenland/Bamber2013/Greenland_bedrock_topography_V3.nc';
	otherwise
		error('machine not supported yet');
end
verbose = 0;

if nargin==2,
	string = 'BedrockElevation';
end

%Convert to Bamber's projections
if verbose, disp('   -- Bamber2013: converting coordinates'); end
[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
[x3971,y3971] = ll2xy(LAT,LON  ,+1,39,71);

if verbose, disp('   -- Bamber2013: loading coordinates'); end
xdata = double(ncread(bamber2013nc,'projection_x_coordinate'));%*1000;
ydata = double(ncread(bamber2013nc,'projection_y_coordinate'));%*1000;

if verbose, disp(['   -- Bamber2013: loading ' string]); end
data  = double(ncread(bamber2013nc,string))';
if verbose, disp(['   -- Bamber2013: interpolating ' string]); end
if strcmpi(string,'LandMask');
	output = InterpFromGrid(xdata,ydata,data,x3971,y3971,'nearest');
else
	output = InterpFromGrid(xdata,ydata,data,x3971,y3971);
end
output = reshape(output,size(X,1),size(X,2));
