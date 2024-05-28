function bool=isshp(filename)
%ISSHP - is a file an shpfile? 
%
%   Usage:
%      isshp(filename);
%
%   See also SHPREAD, SHPDOC

[path,name,ext]=fileparts(filename); 
if strcmpi(ext,'.shp'),
	bool=1;
else 
	bool=0;
end
