function v = sshversion(inputs)
%SSHVERSION get ssh version, using "ssh" on unix.
%
% Usage
%  v = sshversion;
%}

%Get ssh version
[status, stdout] = system('ssh -V');

if status ~= 0
	error(sprintf('ERROR: check error message: %s', stdout));
end

% search date.
pattern = '(?<month>\d{1,2}) (?<day>[A-Za-z]{3}) (?<year>\d{4})';
match = regexp(stdout, pattern, 'match');

% return ssh version information.
v = struct('version',stdout,...
	'date',datetime(match,'InputFormat','dd MMM yyyy'));

end
