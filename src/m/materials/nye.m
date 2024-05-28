function rigidity = nye(temperature,ice_type)
%NYE - Nye viscosity coefficient
%
%   Compute rigidity of ice (either CO2 or H2O) for a given temperature
%   rigidity (in s^(1/n)Pa) is the flow law parameter in the flow law
%   sigma=B*e(1/n) (Nye, p2000).  temperature is in Kelvin degrees
%
%   Usage:
%      rigidity=nye(temperature,ice_type) % ice_type = 1: CO2 ice // ice_type = 2: H2O ice

	% Beyond-melting-point cases
	warning OFF BACKTRACE
	if (ice_type==1)
		if (any(temperature>200&temperature<220))
			warning('nye.m: CO2 ICE - POSSIBLE MELTING. Some temperature values are between 200K and 220K.\nLook at indexes: %s', mat2str(find(temperature>200 & temperature<220))');
		end
		if (any(temperature>=220))
			warning('nye.m: CO2 ICE - GUARANTEED MELTING. Some temperature values are beyond 220K.\nLook at indexes: %s', mat2str(find(temperature>=220))');
		end
	elseif ((ice_type==2)&&(any(temperature>273.15)))
		warning('nye.m: H2O ICE - GUARANTEED MELTING. Some temperature values are beyond 273.15K.\nLook at indexes: %s', mat2str(find(temperature>273.15))');
	end

	% Coefficients
	Rg=8.3144598;       % J mol^-1 K^-1

	if(ice_type==1)     % CO2 ice
		A_const     = 10^(13.0);    % s^-1 MPa
		Q           = 66900;        % J mol^-1
		n           = 8;            % Glen's exponent
	elseif(ice_type==2) % H2O ice
		A_const     = 9e4;          % s^-1 MPa
		Q           = 60000;        % J mol^-1
		n           = 3;            % Glen's exponent
	else
		error('Ice type not supported');
	end

	% Arrhenius Law
	A=A_const*exp(-Q./(temperature*Rg));  % s^-1 MPa
	rigidity=A.^(-1/n)*1e6;               % s^(1/n) Pa

end
