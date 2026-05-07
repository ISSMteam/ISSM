function H = HeightAboveFloatation(md,step)
%HEIGHTABOVEFLOATATION - returns height above hydrostatic floatation per element within the ice levelset.
% H<0 will be returned for elements which are floating. 
%
% USAGE:
%      H = HeightAboveFloatation(md)          % calculate HAF from md.geometry
%      H = HeightAboveFloatation(md,[])       % calculate HAF from md.geometry
%      H = HeightAboveFloatation(md,10)       % calculate HAF from md.results.TransientSolution(10) -- step 10 of transient solution
%
% INPUT:
%      md		% ISSM model containing the geometry
%		 step		% index of md.results.TransientSolution(step) from which to pull the geometry
%
% OUTPUT:
%      H       % height above floatation over each element
%
% SEE ALSO: VolumeAboveFloatation

% process inputs
if (nargin<2 | isempty(step))
	istransientstep = 0;
else
	istransientstep = 1;
end

% define element index
switch class(md.mesh)
	case 'mesh2d'
		index = md.mesh.elements;
	case 'mesh3dprisms'
		index = md.mesh.elements2d;
	otherwise 
		error('not supported yet');
end

% define density parameters
rho_ice   = md.materials.rho_ice; % density of glacial ice (kg m^-3)
rho_water = md.materials.rho_water; % density of ocean water (kg m^-3)

% compute geometry on the elements and define ice levelset
if ~istransientstep
	base           = mean(md.geometry.base(index),2);
	surface        = mean(md.geometry.surface(index),2);
	bed     = mean(md.geometry.bed(index),2);
	ice_levelset   = md.mask.ice_levelset;
else
	if isprop(md.results.TransientSolution(step),'MaskIceLevelset')
		ice_levelset   = md.results.TransientSolution(step).MaskIceLevelset;
	else
		ice_levelset   = md.mask.ice_levelset;
	end
	base           = mean(md.results.TransientSolution(step).Base(index),2);
	surface        = mean(md.results.TransientSolution(step).Surface(index),2);
	if isprop(md.results.TransientSolution(step),'Bed')
		bed  = mean(md.results.TransientSolution(step).Bed(index),2);
	else
		bed  = mean(md.geometry.bed(index),2);
	end
end

% compute height above floatation
H = surface - base + min(rho_water/rho_ice*bed,0);

% mask out elements which are outside the ice levelset
pos = find(min(ice_levelset(index),[],2)>0);
H(pos) = 0;
