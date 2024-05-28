%DIAGNOSTIC class definition
%
%   Usage:
%      diagnostic=diagnostic();

classdef diagnostic
	properties (SetAccess=public) 
		spcvx                    = NaN;
		spcvy                    = NaN;
		spcvz                    = NaN;
		restol                   = 0;
		reltol                   = 0;
		abstol                   = 0;
		isnewton                 = 0;
		FSreconditioning     = 0;
		viscosity_overshoot      = 0;
		icefront                 = NaN;
		maxiter                  = 0;
		shelf_dampening          = 0;
		vertex_pairing           = NaN;
		penalty_factor           = NaN;
		rift_penalty_lock        = NaN;
		rift_penalty_threshold   = 0;
		referential              = NaN;
		loadingforce             = NaN;
		requested_outputs        = NaN;
	end
end
