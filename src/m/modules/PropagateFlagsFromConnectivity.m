function [pool] = PropagateFlagsFromConnectivity(connectivity,pool,index,flags);
%PROPAGATEFLAGSFROMCONNECTIVITY - Propagate flags onto mesh, element by element, using connectivity
%
%   Usage:
%      [pool] = PropagateFlagsFromConnectivity(connectivity,pool,index,flags);

% Check usage
if nargin~=4
	help PropagateFlagsFromConnectivity
	error('Wrong usage (see above)');
end

% Call mex module
[pool] = PropagateFlagsFromConnectivity_matlab(connectivity,pool,index,flags);
