%SURFACEFORCINGS Class definition
%
%   Usage:
%      smb=surfaceforcings();

classdef surfaceforcings
	properties (SetAccess=public) 
		precipitation             = NaN;
		mass_balance              = NaN;
		ispdd                     = 0;
		issmbgradients            = 0;
		isdelta18o                = 0;
		href                      = NaN;
		smbref                    = NaN;
		b_pos                     = NaN;
		b_neg                     = NaN;
		monthlytemperatures       = NaN;
		delta18o                  = NaN;
		delta18o_surface          = NaN;
		temperatures_presentday   = NaN;
		temperatures_lgm          = NaN;
		precipitations_presentday = NaN;
		desfac                    = 0.5;
		s0p                       = 0;
	end
end
