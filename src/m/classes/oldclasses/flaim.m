%FLAIM class definition
%
%   Usage:
%      flaim=flaim();

classdef flaim
	properties (SetAccess=public) 
		targets            = ''
		tracks             = ''
		flightreqs         = struct()
		criterion          = NaN
		gridsatequator     = 200000
		usevalueordering   = true
		split_antimeridian = true
		solution           = ''
		quality            = 0
		path_optimize      = false
		opt_ndir           = 1
		opt_dist           = 25
		opt_niter          = 30000
	end
end
