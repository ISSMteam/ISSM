function output = interpGriggs2013(X,Y,string),

disp('============================================');
disp(' ');
disp('WARNING: interpBamber2013 should now be used');
disp(' ');
disp('============================================');
error('interpBamber2013 should now be used');
griggs2013nc='/u/astrid-r1b/morlighe/issmjpl/proj-morlighem/DatasetGreenland/Data/Griggs2012/Greenland_bedrock_topography_and_geometry_062012_JGriggs.nc';
verbose = 0;

if nargout==2,
	string = 'BedrockElevation';
end

%Convert to Bamber's projections
if verbose, disp('   -- Griggs2013: converting coordinates'); end
[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
[x3971,y3971] = ll2xy(LAT,LON  ,+1,39,71);

if verbose, disp('   -- Griggs2013: loading coordinates'); end
xdata = double(ncread(griggs2013nc,'projection_x_coordinate'))*1000;
ydata = double(ncread(griggs2013nc,'projection_y_coordinate'))*1000;

if verbose, disp(['   -- Griggs2013: loading ' string]); end
data  = double(ncread(griggs2013nc,string))';
if verbose, disp(['   -- Griggs2013: interpolating ' string]); end
output = InterpFromGrid(xdata,ydata,data,x3971,y3971);
output = reshape(output,size(X,1),size(X,2));
