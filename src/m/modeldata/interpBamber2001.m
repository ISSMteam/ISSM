function [bedout thicknessout] = interpBamber2001(X,Y),

switch oshostname(),
	case {'murdo','thwaites','astrid'}
		bamber2001bedpath ='/u/astrid-r1b/ModelData/BamberDEMGreenland5km/bedrock.mat';
		bamber2001thxpath ='/u/astrid-r1b/ModelData/BamberDEMGreenland5km/thickness.mat';
	case {'ronne'}
		bamber2001bedpath ='/home/ModelData/Greenland/Bamber2001/bedrock.mat';
		bamber2001thxpath ='/home/ModelData/Greenland/Bamber2001/thickness.mat';
	case {'totten'}
		bamber2001bedpath ='/totten_1/ModelData/Greenland/Bamber2001/bedrock.mat';
		bamber2001thxpath ='/totten_1/ModelData/Greenland/Bamber2001/thickness.mat';
	otherwise
		error('machine not supported yet');
end

verbose = 0;

%Convert to Bamber's projections
if verbose, disp('   -- Bamber2001: converting coordinates'); end
[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
[x3971,y3971] = ll2xy(LAT,LON  ,+1,39,71);

if verbose, disp('   -- Bamber2001: loading bed'); end
load(bamber2001bedpath);
if verbose, disp('   -- Bamber2001: interpolating bed'); end
bedout = InterpFromGrid((x_m(1:end-1)+x_m(2:end))/2,(y_m(1:end-1)+y_m(2:end))/2,bedrock,x3971,y3971);
bedout = reshape(bedout,size(X,1),size(X,2));

if nargout>1
	if verbose, disp('   -- Bamber2001: loading thickness'); end
	load(bamber2001thxpath);
	if verbose, disp('   -- Bamber2001: interpolating thickness'); end
	thicknessout = InterpFromGrid((x_m(1:end-1)+x_m(2:end))/2,(y_m(1:end-1)+y_m(2:end))/2,thickness,x3971,y3971);
	thicknessout = reshape(thicknessout,size(X,1),size(X,2));
end
