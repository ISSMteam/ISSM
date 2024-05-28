%OLD materials class definition

classdef materials
	properties (SetAccess=public) 
		rho_ice                    = 0.;
		rho_water                  = 0.;
		rho_freshwater             = 0.;
		mu_water                   = 0.;
		heatcapacity               = 0.;
		latentheat                 = 0.;
		thermalconductivity        = 0.;
		meltingpoint               = 0.;
		beta                       = 0.;
		mixed_layer_capacity       = 0.;
		thermal_exchange_velocity  = 0.;
		rheology_B   = NaN;
		rheology_n   = NaN;
		rheology_Z   = NaN;
		rheology_law = '';
	end
end
