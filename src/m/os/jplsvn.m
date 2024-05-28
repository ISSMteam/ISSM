function JPL_SVN=jplsvn()
%ISSMDIR - Get JPL_SVN environment variable
%
%   Usage:
%      JPL_SVN=jplsvn()

if ~ispc(),
	JPL_SVN =getenv('JPL_SVN');
else
	JPL_SVN =getenv('JPL_SVN_WIN');
end

if (isempty(JPL_SVN)),
	error('jplsvn error message: ''JPL_SVN'' environment variable is empty! You should define JPL_SVN in your .cshrc or .bashrc');
end
