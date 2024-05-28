function ISSM_DIR=issmdir()
%ISSMDIR - Get ISSM_DIR environment variable
%
%   Usage:
%      ISSM_DIR=issmdir()

%Initialize output ISSM_DIR
ISSM_DIR='';

%Get ISSM_DIR from function path (we do not want to force users to edit their bashrc)
path = mfilename('fullpath');

%issmdir might be in bin,
slash = filesep();
pos   = strfind(path,['bin' slash 'issmdir']);
if ~isempty(pos),
	ISSM_DIR=path(1:pos-1);
else
	pos=strfind(path,['src' slash 'm' slash 'os' slash 'issmdir']);
	if ~isempty(pos),
		ISSM_DIR=path(1:pos-1);
	end
end

if isempty(ISSM_DIR),
	error('Could not determine the location of ISSM...');
end
