function name=sessionname()
%SESSIONNAME - Get screen session name
%
%   Usage:
%      name=sessionname()


	[status,styname]=system('export | grep STY'); 

	sessionname=styname(17:end-2);
	if isempty(sessionname),
		name='none';
		return;
	end
	index=findstr(sessionname,'.');
	name=sessionname(index+1:end);
