function smb=PattynSMB(md,Tf)
%PATTYNSMB- Compute SMB over Antarctica (from Pattyn 2006, pg. 18, "GRANTISM: An ExcelTM model for Greenland and 
%           Antarctic ice-sheet response to climate changes")
%
%   Usage:
%      md=PattynSMB(md,Tf)
%
%      where Tf is a background forcing temperature ("an anomalous temperature relative to the present conditions)
%
%
%   See also: SETICESHELFBC, SETMARINEICESHEETBC

 
% Tma    : Mean annual surface temperature in [deg C]
% Tms    : Mean summer temperature in [deg C]
% h      : Surface/bedrock elevation (I assume in meters but paper does not specify)
% phi    : Latitude in degrees SOUTH
% lambda : Longitude in degrees WEST
% Tf     : Background forcing temperature ("an anomalous temperature relative to the present conditions)
% ACCdot : Accumulation rate in units of [m/a] ice equivalent
% ABLdot : Surface ablation rate in [m/a] ice equivalent

	%Double check lat and long exist:
	if isnan(md.mesh.lat),
		error('PattynSMB error message: md.mesh.lat field required');
	end

	%Delta Temperature forcing 
	if (nargin==1),
		Tf=0; 
	end

	% Calculate mean annual surface temperature, Eqn (11)
	% Here, -0.012 is the atmospheric Lapse rate from sea level in deg/m.
	% It is multiplied by surface elevation from sea level
	Tma = -15.15 - 0.012*md.geometry.surface;

	%% Calculate summer temperature, Eqn (12)
	%% No melting at PIG in mean conditions - need about 6 degress Tf to start having a negative yearly SMB
	Tms = 16.81 - 0.00692*md.geometry.surface - 0.27937*abs(md.mesh.lat) + Tf;

	%% Calculate Accumulation perturbation with Tf forcing, Eqn (9)
	ACCdot = 2.5*2.^((Tma+Tf)/10.) - 2.5*2.^(Tma/10.);

	%% Calculate Ablation, Eqn (10) (use for both Antarctica & Greenland), max melt is 10m/a
	ABLdot=zeros(size(ACCdot));
	pos=Tms>=0;
	ABLdot(pos)=min(1.4*Tms(Tms>=0),10);

	smb=ACCdot-ABLdot;
end
