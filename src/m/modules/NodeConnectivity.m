function connectivity = NodeConnectivity(elements,numnodes);
%NODECONNECTIVITY - Build node connectivity from elements
%
%   Usage:
%      connectivity = NodeConnectivity(elements,numnodes);

% Check usage
if nargin~=2
	help NodeConnectivity
	error('Wrong usage (see above)');
end

% Call mex module
connectivity = NodeConnectivity_matlab(elements,numnodes);
