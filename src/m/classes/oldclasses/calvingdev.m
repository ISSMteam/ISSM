%CALVINGDEV class definition
%
%   Usage:
%      calvingdev=calvingdev();

classdef calvingdev
	properties (SetAccess=public) 
		stress_threshold_groundedice = 0.;
		stress_threshold_floatingice = 0.;
		meltingrate   = NaN;
	end
end
