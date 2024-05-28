function md = transientrestart(md,step)
%TRANSIENTRESTART - reinitialize model from last transient step
%
%   Usage:
%      md = transientrestart(md)
%      md = transientrestart(md,step)
%
%   By default, transientrestart will use the last step provided in md.results.TransientSolution

%Get result and save it again
if nargin==1,
	step = numel(md.results.TransientSolution);
end
if step<1,
	error('step needs to be >0');
elseif step>numel(md.results.TransientSolution)
	error(['md.results.TransientSolution has only ' num2str(numel(md.results.TransientSolution)) ' steps']);
end
results = md.results.TransientSolution(step);

newname = ['TransientSolution' num2str(numel(fields(md.results))+1)];
if isfield(md.results,newname)
	error(['Cannot save ' newname ' in md.results']);
else
	disp(['Moving results to ' newname]);
	md.results.(newname) = md.results.TransientSolution;
	md.results.TransientSolution  = struct();
end

%Change time
md.timestepping.start_time = results.time;

%Change initialization fields
if isfield(results,'Vx'),          md.initialization.vx=results.Vx; end
if isfield(results,'Vy'),          md.initialization.vy=results.Vy; end
if isfield(results,'Vz'),          md.initialization.vz=results.Vz; end
if isfield(results,'Vel'),         md.initialization.vel=results.Vel; end
if isfield(results,'Temperature'), md.initialization.temperature=results.Temperature; end
if isfield(results,'Pressure'),    md.initialization.pressure=results.Pressure; end
if isfield(results,'Waterfraction'),md.initialization.waterfraction=results.Waterfraction; end
if isfield(results,'Watercolumn'), md.initialization.watercolumn=results.Watercolumn; end
if isfield(results,'Enthalpy'),    md.initialization.enthalpy=results.Enthalpy; end
if isfield(results,'DebrisThickness'),md.initialization.debris=results.DebrisThickness; end

%Deal with new geometry
if isfield(results,'Base') & isfield(results,'Thickness'),
	base=results.Base;
	thickness=results.Thickness;
	if isa(md.mesh,'mesh3dprisms')
		md.mesh.z=base+thickness./md.geometry.thickness.*(md.mesh.z-md.geometry.base);
	elseif isa(md.mesh,'mesh2dvertical')
		md.mesh.y=base+thickness./md.geometry.thickness.*(md.mesh.y-md.geometry.base);
	end
	md.geometry.base=base;
	md.geometry.thickness=thickness;
	md.geometry.surface=md.geometry.base+md.geometry.thickness;
end

%Update mask
if isfield(results,'MaskOceanLevelset'),
	md.mask.ocean_levelset = results.MaskOceanLevelset;
end
if isfield(results,'MaskIceLevelset'),
	md.mask.ice_levelset = results.MaskIceLevelset;
end
%Update mask in case of old model version, before name change groundedice -> ocean levelset
if isfield(results,'MaskGroundediceLevelset'),
	md.mask.ocean_levelset = results.MaskGroundediceLevelset;
end
