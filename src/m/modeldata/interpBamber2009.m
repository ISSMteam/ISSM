function demout = interpBamber2009(X, Y)
%INTERPBAMBER2009 - interpolate surface dem of Bamber 2009
%
%   Surface dem Nominal year 2004 (WGS84, no firn correction)
%
%   Usage:
%      demout = interpBamber2009(X, Y)

switch oshostname(),
	case {'totten'}
		bamber2009path ='/totten_1/ModelData/Antarctica/Bamber2009DEM/krigged_dem_nsidc.mat';
	otherwise
		error('machine not supported yet');
end

%Convert to Bamber's projections
%disp('   -- Bamber2009: loading dem'); 
load(bamber2009path);

disp('   -- Bamber2009: interpolating dem (WGS84)');
demout = InterpFromGrid(x, y, surfacedem, X, Y);
