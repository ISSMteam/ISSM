function hostname=oshostname()
%OSHOSTNAME - Determine hostname, irrespective of os type
%
%   Usage:
%      hostname=oshostname();

%Initialize output
hostname = '';

%First try using java (most stable way)
if usejava('jvm'),
	hostname = char(getHostName(java.net.InetAddress.getLocalHost));
end

%Method 2: use system command (MATLAB bug includes what's in the clipboard)
if isempty(hostname),
	%See http://www.mathworks.com/help/matlab/ref/system.html "tips" section
	%We need to add < /dev/null otherwise what is in the clipboard is added
	[status,hostname]=system('hostname < /dev/null');
end

%Method 3, last chance
if isempty(hostname),
	if ispc % If OS is MinGW, $COMPUTERNAME and $HOSTNAME are identical
		hostname = getenv('COMPUTERNAME');
	else
		hostname = getenv('HOSTNAME');
	end
end

% Take out minus signs
hostname = strrep(hostname,'-','');

% Trim and lower case
hostname = strtrim(lower(hostname));

% Check that machine name is not empty
if isempty(hostname),
	error('Cannot determine machine name');
end
