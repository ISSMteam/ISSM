function p = recursivepath(d)
%RECURSIVEPATH - generate paths in a directory
%
%   this routine is equivalent to Matlab's genpath except that it skips CVS and
%   .svn directories
%
%   Usage:
%      p = recursivepath(d)

%initialize path to be returned
p = '';
sep=pathsep;  %directory separator

% Generate path based on given root directory
files=dir(d);
if isempty(files)
	return
end

% Add d to the path even if it is empty.
p = [p d sep];

% set logical vector for subdirectory entries in d
isdir = logical(cat(1,files.isdir));

% Recursively goes through the subdirectories of d
dirs=files(isdir); % select only directory entries from the current listing
for i=1:length(dirs)
	dirname=dirs(i).name;
	if ~strcmp(dirname,'.')    & ...
		~strcmp(dirname,'..')   & ...
		~strcmp(dirname,'.svn') & ...
		~strcmp(dirname,'CVS')  & ...
		~strncmp(dirname,'@',1) & ... %Method directories not allowed in MATLAB path
		~strcmp(dirname,'private')    %private directories not allowed in MATLAB path

		p = [p recursivepath(fullfile(d,dirname))];
	end
end
