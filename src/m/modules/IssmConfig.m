function value = IssmConfig(string);
%ISSMCONFIG
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
