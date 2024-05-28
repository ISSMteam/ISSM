function tprofile=robintemperature(heatflux,accumrate,thickness,surftemp,z)
%ROBINTEMPERATURE - compute vertical temperature profile of an ice sheet (Robin, 1955)
%
%	This routine computes the vertical temperature profile of an 
%	ice sheet according to the solution of Robin (1955), neglecting 
%	friction and horizontal advection.  The solution is thus most
%  appropriate at an ice divide.
%
%	The coordinate system for the solution runs from z=0 at the base 
%	to z=H at the surface of the ice.
%
%	Parameters (SI units):
%		-heatflux	Geothermal heat flux (W m^-2)
%		-accumrate	Surface accumulation rate (m s^-1 ice equivalent)
%		-thickness	Ice thickness (m)
%		-surftemp	Surface temperature (K)
%		-z				Vertical position at which to calculate temperature
%						(z can be a scalar or a vector)
%
%	Returns a vector the same length as z containing the temperature in K
%
%	Usage:
%		tprofile=robintemperature(heatflux,accumrate,thickness,surftemp,z)

%checks
if nargin~=5
	help robintemperature
	error('bad usage - wrong number of arguments.')
end

%some constants (from Holland and Jenkins, 1999)
alphaT=1.14e-6; % thermal diffusivity (m^2 s^-1)
c=2009; % specific heat capacity (J kg^-1 K^-1)
rho=917; % ice density (kg m^-3)

%create vertical coordinate variable
zstar=sqrt(2*alphaT*thickness./accumrate);

tprofile=surftemp+sqrt(2*thickness*pi./accumrate/alphaT).*(-heatflux)/2/rho/c.*(erf(z./zstar)-erf(thickness./zstar));

% difference between surface and base temperature for check (Cuffey2010 p412):
% tprofile-surftemp
