function value = IssmConfig(string);
%ISSMCONFIG - return the value of an ISSM compile-time configuration option
%
%   Usage:
%      value = IssmConfig('string');

% Check usage
if nargin~=1
	help IssmConfig
	error('Wrong usage (see above)');
end

% Call mex module
value = IssmConfig_matlab(string);
