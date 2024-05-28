%HYDROLOGY class definition
%
%   Usage:
%      hydrology=hydrology();

classdef hydrology
	properties (SetAccess=public) 
		spcwatercolumn = NaN;
		n              = 0;
		CR             = 0;
		p              = 0;
		q              = 0;
		kn             = 0;
		stabilization  = 0;
	end
end
