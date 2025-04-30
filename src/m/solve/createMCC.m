function createMCC(filename)
%CREATEMCC - create mcc file from MATLAB file
%
%   takes an existing matlab script and compiles it to an executable file.
%   A number of files are produced, including 'run_MCCexecutable.sh' and 'MCCexecutable,'
%   and are saved into a directory './mccfiles' which is cleared if it already exists and is
%   created if it does not.
%
%   USAGE:
%      $ matlab 
%      >> createMCC(filename);
%      >> exit
%      $ ./mccfiles/run_MCCexecutable.sh /nasa/matlab/2022b/
%
%   INPUT
%      filename .m file to be turned into an executable
%
%   OUTPUT
%      no direct output is produced, however CREATEMCC creates a number of files in ./mccfiles

% check if mccfiles directory exists in current directory
if exist('./mccfiles','dir')
	system('rm ./mccfiles/* ');
else
	mkdir('./mccfiles');
end

%Get dependencies
files = matlab.codetools.requiredFilesAndProducts(filename);

%Creaste long string for command, ignore Matlab's statistical toolbox license
deps = [];
for i=1:numel(files)
   if contains(files{i},'normfit_issm.m')
      continue
   elseif contains(files{i},'dakota_moments.m')
      continue
   elseif contains(files{i},'dakota_out_parse.m')
      continue
   else
      deps = [deps ' ' files{i}];
   end
end

%Create command
command = ['mcc -m ' deps ' -o MCCexecutable'];

%Create executable
cd ./mccfiles
system(command);
cd ..
