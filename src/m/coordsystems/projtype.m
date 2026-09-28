function type=projtype(proj)
%PROJTYPE - extract the projection type (+proj value) from a PROJ.4 string
%
%   Usage:
%      type=projtype(proj);

	ind=findstr(proj,'+proj');
	proj=proj(ind+1:end);
	ind=findstr(proj,'+');
	proj=proj(1:ind-2);
	ind=findstr(proj,'=');
	type=proj(ind+1:end);
