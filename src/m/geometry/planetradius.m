function radius=planetradius(planet)
%PLANETRADIUS - return planet radius according to planetary body name.
%
%   Usage:
%      radius=planetradius(planet);
%
%   Examples:
%      earthradius=planetradius('earth');

if strcmpi(planet,'earth'),
	radius=6.371012*10^6;
elseif strcmpi(planet,'europa'),
	radius=1.5008*10^6;
else 
	error(['planet type ' planet ' not supported yet!']);
end
