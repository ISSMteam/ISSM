function V = VolumeAboveFloatation(md,step,flags)
%VOLUMEABOVEFLOATATION - returns volume above floatation
%
%   Usage:
%      V = VolumeAboveFloatation(md)          % uses model fiels alone
%      V = VolumeAboveFloatation(md,10)       % Will look at step 10 of transient solution
%      V = VolumeAboveFloatation(md,10,flags) % Will look at step 10 of transient solution, only flaged elements

%Special case if 3d
if isa(md.mesh,'mesh3dprisms')
	index = md.mesh.elements2d;
	x = md.mesh.x2d;
	y = md.mesh.y2d;
elseif isa(md.mesh,'mesh2d'),
	index = md.mesh.elements;
	x = md.mesh.x;
	y = md.mesh.y;
else
	error('not supported yet');
end

%1. get some parameters
rho_ice   = md.materials.rho_ice;
rho_water = md.materials.rho_water;

%2. compute averages
if nargin==1
	base           = mean(md.geometry.base(index),2);
	surface        = mean(md.geometry.surface(index),2);
	bathymetry     = mean(md.geometry.bed(index),2);
	ice_levelset   = md.mask.ice_levelset;
	ocean_levelset = md.mask.ocean_levelset;
else
	if isprop(md.results.TransientSolution(step),'MaskIceLevelset')
		ice_levelset   = md.results.TransientSolution(step).MaskIceLevelset;
	else
		ice_levelset   = md.mask.ice_levelset;
	end
   ocean_levelset = md.results.TransientSolution(step).MaskOceanLevelset;
   base           = mean(md.results.TransientSolution(step).Base(index),2);
   surface        = mean(md.results.TransientSolution(step).Surface(index),2);
	if isprop(md.results.TransientSolution(step),'Bed')
		bathymetry  = mean(md.results.TransientSolution(step).Bed(index),2);
	else
		 bathymetry  = mean(md.geometry.bed(index),2);
	 end
end

%3. get areas of all triangles
areas = GetAreas(index,x,y);

%4. Compute volume above floatation
V = areas.*(surface-base+min(rho_water/rho_ice*bathymetry,0.));

%5. take out the ones that are outside of levelset or floating
pos = find(min(ice_levelset(index),[],2)>0 | min(ocean_levelset(index),[],2)<0);
V(pos) = 0;

%In case we are only looking at one portion of the domain...
if nargin==3
	V(find(~flags)) = 0;
end

%sum individual contributions
V = sum(V);
