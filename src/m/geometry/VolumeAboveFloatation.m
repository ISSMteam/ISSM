function V = VolumeAboveFloatation(md,step,flags)
%VOLUMEABOVEFLOATATION - returns volume above floatation
% Sums per-element volume above floatation over grounded ice, considering only positive height above floatation.
%
% USAGE:
%      V = VolumeAboveFloatation(md)          % calculate VAF from md.geometry
%      V = VolumeAboveFloatation(md,10)       % calculate VAF from md.results.TransientSolution(10) -- step 10 of transient solution 
%      V = VolumeAboveFloatation(md,10,flags) % calculate VAF from md.results.TransientSolution(10), only flagged elements
%      V = VolumeAboveFloatation(md,[],flags) % calculate VAF from md.geometry, only flagged elements
%
% INPUT:
%      md      % ISSM model containing the geometry
%      step    % index of md.results.TransientSolution(step) from which to pull the geometry
%      flags   % boolean vector of length md.mesh.numberofelements flagging which elements to integrate over
%
% OUTPUT:
%      V       % volume above floatation over each element
%
% SEE ALSO: HeightAboveFloatation, GetAreas

% process inputs
if nargin<2
	step = []; % calculate VAF from md.geometry
end
if nargin<3
	flags = []; % integrate over entire domain
end
if (nargin<2 | isempty(step))
	istransientstep = 0;
else
	istransientstep = 1;
end

% define element index
switch class(md.mesh)
	case 'mesh2d'
		index = md.mesh.elements;
		x = md.mesh.x;
		y = md.mesh.y;
	case 'mesh3dprisms'
		index = md.mesh.elements2d;
		x = md.mesh.x2d;
		y = md.mesh.y2d;
	otherwise
		error('not supported yet');
end

% retrieve ice and ocean levelsets
if ~istransientstep
	ice_levelset   = md.mask.ice_levelset;
	ocean_levelset = md.mask.ocean_levelset;
else
	if isprop(md.results.TransientSolution(step),'MaskIceLevelset')
		ice_levelset   = md.results.TransientSolution(step).MaskIceLevelset;
	else
		ice_levelset   = md.mask.ice_levelset;
	end
   ocean_levelset = md.results.TransientSolution(step).MaskOceanLevelset;
end

% compute height above floatation per element
H = HeightAboveFloatation(md,step); % (m)

% mask out elements which are outside the levelset or are floating
pos = find(min(ice_levelset(index),[],2)>0 | min(ocean_levelset(index),[],2)<0); % index of elements outside ice or inside ocean mask
H(pos) = 0; 
H = max(H,0); % Ensure that negative height above floatation is not integrated 

% compute volume above floatation per element
areas = GetAreas(index,x,y); % (m^2)
V = areas.*H; % (m^3)

% apply flags
V(find(~flags)) = 0;

% sum over all elements
V = sum(V);
