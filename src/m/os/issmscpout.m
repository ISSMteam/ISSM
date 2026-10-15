function issmscpout(host, path, login, port, packages, bracketstyle, transfermethod)
%ISSMSCPOUT - Send files to host
%
%   Usage:
%      issmscpout(host,path,login,port,packages,bracketstyle,transfermethod)
%
%      bracketstyle:    1 - \{\}    (escaped; default)
%                       2 - {}      (not escaped)
%                       Only used by the scp method: rsync sends the file list
%                       over the protocol and never relies on the remote shell.
%
%      transfermethod: 'auto'  - use rsync, fall back to scp if the remote host
%                                does not have rsync installed (default)
%                      'rsync' - use rsync only
%                      'scp'   - use scp only

%get hostname
hostname=oshostname();

%does machine require escaped brackets?
if nargin<6 || isempty(bracketstyle)
	bracketstyle = 1;
end

%how do we move the files over?
if nargin<7 || isempty(transfermethod)
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
usescp = strcmpi(transfermethod,'scp');
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
	[status,cmdout]=scpout(host,login,port,path,packages,bracketstyle);
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
%RSYNCOUT - push packages with rsync, handing the list over the protocol

	%--files-from takes the file names from a local file, so the remote shell
	%never sees (and never has to expand) a list. Names are relative to './'.
	listfile = [tempname '.lst'];
	fid = fopen(listfile,'w');
	if fid==-1
		error(['issmscpout error message: could not create temporary file ' listfile]);
	end
	for i=1:numel(packages)
		fprintf(fid,'%s\n',packages{i});
	end
	fclose(fid);

	if port
		cmd = ['rsync -a --files-from=''' listfile ''' -e ''ssh -p ' num2str(port) ''' ./ ' login '@localhost:''' path ''''];
	else
		cmd = ['rsync -a --files-from=''' listfile ''' ./ ' login '@' host ':''' path ''''];
	end
	[status,cmdout]=system(cmd);

	delete(listfile);
end % }}}
function [status,cmdout]=scpout(host,login,port,path,packages,bracketstyle) % {{{
%SCPOUT - push packages with scp, relying on the local shell to expand \{\}

	if numel(packages)==1
		fileliststr=packages{1};
	else
		fileliststr='\{';
		for i=1:numel(packages)-1
			fileliststr=[fileliststr packages{i} ','];
		end
		fileliststr=[fileliststr packages{end} '\}'];

		%remove backslashes if bracketstyle is 2
		if bracketstyle==2
			fileliststr = [fileliststr(2:end-2) fileliststr(end)];
		end
	end

	cmdout='';
	if port
		[status]=system(['scp -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path]);
		if status~=0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path]);
		end
	else
		[status]=system(['scp ' fileliststr ' ' login '@' host ':' path]);
		if status~=0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT ' fileliststr ' ' login '@' host ':' path]);
		end
	end
end % }}}
