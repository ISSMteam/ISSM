function plot_none(md,options,nlines,ncols,i)
%PLOT_NONE - plot nothing, just apply options
%
%   Usage:
%      plot_mesh(md,options,nlines,ncols,i);
%
%   See also: PLOTMODEL
return;

options=addfielddefault(options,'colorbar','none');
options=addfielddefault(options,'map','none');
options=addfielddefault(options,'axis','equal');

if exist(options,'overlay'),
	plot_overlay(md,'none',options,nlines,ncols,i);
	return;
end

%apply options
applyoptions(md,[],options);
