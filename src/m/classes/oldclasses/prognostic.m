%MASSTRANSPORT class definition
%
%   Usage:
%      prognostic=prognostic();

classdef prognostic
	properties (SetAccess=public) 
		 spcthickness           = NaN;
		 min_thickness          = 0;
		 hydrostatic_adjustment = 0;
		 stabilization          = 0;
		 vertex_pairing         = NaN;
		 penalty_factor         = 0;
		 requested_outputs      = NaN;
	end
end
