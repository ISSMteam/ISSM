function bool=isexp(filename)
%ISEXP - is a file an expfile? 
%
%   Usage:
%      isexp(filename);
%
%   See also EXPREAD, EXPDOC

[path,name,ext]=fileparts(filename); 
if strcmpi(ext,'.exp'),
	bool=1;
else 
	bool=0;
end
