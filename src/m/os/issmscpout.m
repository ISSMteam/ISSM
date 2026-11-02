function issmscpout(host, path, login, port, packages, transfermethod)
%ISSMSCPOUT - Send files to host
%
%   Usage:
%      issmscpout(host,path,login,port,packages,transfermethod)
%
%      transfermethod: 'auto'         - use rsync, fall back to scp if the remote host
%                                       does not have rsync installed (default)
%                      'rsync'        - use rsync only
%                      'scp'          - use scp only, with escaped brackets \{\}
%                      'scp-noescape' - use scp only, with unescaped brackets {}
%                                       (for hosts whose shell must not see backslashes)

%get hostname
hostname=oshostname();

%how do we move the files over?
if nargin<6 || isempty(transfermethod)
	transfermethod = 'auto';
end

%if hostname and host are the same, do a simple symlinks
if strcmpi(host,hostname)

	%Process both paths and add \ if there are any white spaces
	here = replace(pwd(), ' ', '\ ');

	for i=1:numel(packages)
		system(['rm -rf ' path '/' packages{i} ]);
		system(['ln -s ' here '/' packages{i} ' ' path]);
	end
	return;
end

%General case: this is not a local machine
usescp = strcmpi(transfermethod,'scp') | strcmpi(transfermethod,'scp-noescape');
if ~usescp & ~strcmpi(transfermethod,'rsync') & ~strcmpi(transfermethod,'auto')
	error(['issmscpout error message: unknown transfermethod ''' transfermethod '''']);
end

%no point in trying rsync if this machine does not have it
if ~usescp & ~haversync()
	if strcmpi(transfermethod,'rsync')
		error('issmscpout error message: rsync requested but not found on this machine');
	end
	usescp = 1;
end

if ~usescp
	[status,cmdout]=rsyncout(host,login,port,path,packages);

	%rsync returns 127 when the remote host does not have rsync installed.
	%Any other non-zero status is a genuine transfer error and must not
	%trigger a second connection (it would prompt for credentials twice).
	if status==127 & strcmpi(transfermethod,'auto')
		warning(['issmscpout: rsync not found on ' host ', falling back to scp']);
		usescp = 1;
	end
end

if usescp
	[status,cmdout]=scpout(host,login,port,path,packages,strcmpi(transfermethod,'scp-noescape'));
end

%check that the transfer worked
if status~=0
	error(['issmscpout error message: ' cmdout])
end

end

function tf=haversync() % {{{
%HAVERSYNC - is rsync available on this machine? (probed once per session)

	persistent found;

	if isempty(found)
		if ispc
			[status,~]=system('where rsync');
		else
			[status,~]=system('command -v rsync');
		end
		found = (status==0);
	end
	tf = found;
end % }}}
function [status,cmdout]=rsyncout(host,login,port,path,packages) % {{{
%RSYNCOUT - push packages with rsync, passing every package as a source

	if port
		remote = [login '@localhost:'];
		sshopt = [' -e ''ssh -p ' num2str(port) ''''];
	else
		remote = [login '@' host ':'];
		sshopt = '';
	end

	sources = '';
	for i=1:numel(packages)
		sources = [sources ' ''' packages{i} ''''];
	end

	cmd = ['rsync --progress -a' sshopt sources ' ' remote '''' path ''''];
	[status,cmdout]=system(cmd,'-echo');
end % }}}
function [status,cmdout]=scpout(host,login,port,path,packages,noescape) % {{{
%SCPOUT - push packages with scp, relying on the local shell to expand \{\}

	if numel(packages)==1
		fileliststr=packages{1};
	else
		fileliststr='\{';
		for i=1:numel(packages)-1
			fileliststr=[fileliststr packages{i} ','];
		end
		fileliststr=[fileliststr packages{end} '\}'];

		%remove backslashes if noescape is set
		if noescape
			fileliststr = [fileliststr(2:end-2) fileliststr(end)];
		end
	end

	cmdout='';
	if port
		[status]=system(['scp -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path],'-echo');
		if status~=0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path],'-echo');
		end
	else
		[status]=system(['scp ' fileliststr ' ' login '@' host ':' path],'-echo');
		if status~=0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT ' fileliststr ' ' login '@' host ':' path],'-echo');
		end
	end
end % }}}
