function elementconnectivity = ElementConnectivity(elements,nodeconnectivity);
%ELEMENTCONNECTIVITY - Build element connectivity using node connectivity and elements
%
%   Usage:
%      elementconnectivity = ElementConnectivity(elements,nodeconnectivity);

% Check usage
if nargin~=2
	help ElementConnectivity
	error('Wrong usage (see above)');
end

% Call mex module
elementconnectivity = ElementConnectivity_matlab(elements,nodeconnectivity);
