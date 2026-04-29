function N = effectivepressure(md)
% EFFECTIVEPRESSURE - Calculate the effective basal pressure N from md.geometry and effective pressure coupling rule in md.friction
% USEAGE:
%   N = effectivepressure(md); 
% INPUT:
%   md      ISSM model from which to calculate pressure
% OUTPUT:
%   N			effective pressure at base (Pa)
%
%   See also: BASALSTRESS

sealevel = 0; % sea level reference elevation(m)

p_ice   = md.constants.g * md.materials.rho_ice * md.geometry.thickness; % ice pressure (Pa)
p_water = md.constants.g * md.materials.rho_water * (sealevel-md.geometry.base); % water pressure (Pa)

% calculate effective pressure using coupling definition in md.friction (Pa)
switch(md.friction.coupling)
   case 0 % uniform sheet (negative water pressure ok, default)
      N = p_ice-p_water;
   case 1 % effective pressure is equal to the overburden pressure
      N = p_ice;
   case 2 % uniform sheet (water pressure >= 0)
      p_water=max(p_water,0.0);
      N = p_ice-p_water;
   case 3 % Use effective pressure prescrived in md.friction.effective_pressure
      N = max(md.friction.effective_pressure, 0);
   case 4 % Use effective pressure dynamically calculated by the hydrology model (i.e., fully coupled)
      error('md.friction.coupling=4 is not supported yet.');
   otherwise
      error('not supported yet');
end
