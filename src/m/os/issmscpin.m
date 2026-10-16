function issmscpin(host, login, port, path, filelist, transfermethod)
%ISSMSCPIN - Get filelist from host
%
%   Usage:
%      issmscpin(host,login,port,path,filelist,transfermethod)
%
%      transfermethod: 'auto'         - use rsync, fall back to scp if the remote host
%                                       does not have rsync installed (default)
%                      'rsync'        - use rsync only
%                      'scp'          - use scp only, with escaped brackets \{\}
%                      'scp-noescape' - use scp only, with unescaped brackets {}
%                                       (for hosts whose shell must not see backslashes)

%how do we move the files over?
if nargin<6 || isempty(transfermethod)
	transfermethod = 'auto';
end

%get hostname
hostname=oshostname();

%if hostname and host are the same, do a simple copy
if strcmpi(hostname,host)
	for i=1:numel(filelist)
		status = copyfile([path '/' filelist{i} ], '.');
		if(~status) warning(['Could not copy ' path '/' filelist{i} ]); end
	end
else
	usescp = strcmpi(transfermethod,'scp') | strcmpi(transfermethod,'scp-noescape');
	if ~usescp & ~strcmpi(transfermethod,'rsync') & ~strcmpi(transfermethod,'auto')
		error(['unknown transfermethod ''' transfermethod '''']);
	end

	%no point in trying rsync if this machine does not have it
	if ~usescp & ~haversync()
		if strcmpi(transfermethod,'rsync')
			error('rsync requested but not found on this machine');
		end
		usescp = 1;
	end

	if ~usescp
		[status,cmdout]=rsyncin(host,login,port,path,filelist);

		%rsync returns 127 when the remote host does not have rsync installed.
		%Any other non-zero status is a genuine transfer error and must not
		%trigger a second connection (it would prompt for credentials twice).
		if status==127 & strcmpi(transfermethod,'auto')
			warning(['rsync not found on ' host ', falling back to scp']);
			usescp = 1;
		end
	end

	if usescp
		[status,cmdout]=scpin(host,login,port,path,filelist,strcmpi(transfermethod,'scp-noescape'));
	end

	%check that the transfer worked
	if status ~= 0
		error(['issmscpin error message: ' cmdout])
	end
	for i=1:numel(filelist)
		if ~exist(['./' filelist{i}])
			warning(['could not transfer ' filelist{i}]);
		end
	end
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
function [status,cmdout]=rsyncin(host,login,port,path,filelist) % {{{
%RSYNCIN - pull filelist with rsync, handing the list over the protocol

	%--files-from takes the file names from a local file, so the remote shell
	%never sees (and never has to expand) a list. Names are relative to path.
	listfile = [tempname '.lst'];
	fid = fopen(listfile,'w');
	if fid==-1
		error(['could not create temporary file ' listfile]);
	end
	for i=1:numel(filelist)
		fprintf(fid,'%s\n',filelist{i});
	end
	fclose(fid);

	if port
		cmd = ['rsync --progress -a --files-from=''' listfile ''' -e ''ssh -p ' num2str(port) ''' ' login '@localhost:''' path '/'' ./'];
	else
		cmd = ['rsync --progress -a --files-from=''' listfile ''' ' login '@' host ':''' path '/'' ./'];
	end
	[status,cmdout]=system(cmd);

	delete(listfile);
end % }}}
function [status,cmdout]=scpin(host,login,port,path,filelist,noescape) % {{{
%SCPIN - pull filelist with scp, relying on the remote shell to expand \{\}

	if numel(filelist)==1
		fileliststr=filelist{1};
	else
		fileliststr='\{';
		for i=1:numel(filelist)-1
			fileliststr=[fileliststr filelist{i} ','];
		end
		fileliststr=[fileliststr filelist{end} '\}'];

		%remove \ if noescape is set
		if noescape
			fileliststr = [fileliststr(2:end-2) '}'];
		end
	end

	cmdout='';
	if port
		[status]=system(['scp -P ' num2str(port) ' ' login '@localhost:' path '/' fileliststr ' ./']);
		if status ~= 0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT -P ' num2str(port) ' ' login '@localhost:' path '/' fileliststr ' ./']);
		end
	else
		[status]=system(['scp ' login '@' host ':' path '/' fileliststr ' ./']);
		if status ~= 0
			%List expansion is a bashism. Try again with '-OT'.
			disp('issmscp info: standard scp command did not work, trying again with -OT flag');
			[status,cmdout]=system(['scp -OT ' login '@' host ':' path '/' fileliststr ' ./']);
		end
	end
end % }}}
