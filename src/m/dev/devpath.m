% clear the last warning to focus on the warnings of the ISSM path
lastwarn(''); 

%Recover ISSM_DIR , or if on a Windows machine, ISSM_DIR_WIN
if ispc,
	ISSM_DIR=getenv('ISSM_DIR_WIN');
else
	ISSM_DIR=getenv('ISSM_DIR');
end

if (isempty(ISSM_DIR)),
	error('''ISSM_DIR'' environment variable is empty! You should define ISSM_DIR in your .cshrc or .bashrc!');
end

%Now add all issm code paths necessary to run issm smoothly. 
%We capture the error output, so that we can warn the user to update 
%the variable ISSM_DIR in this file, in case it is not correctly setup. 

%ISSM path
addpath([ISSM_DIR '/src/m/os/']);       %load recursivepath
addpath([ISSM_DIR '/lib']);             %load MEX files
addpath(recursivepath([ISSM_DIR '/src/m']));
addpath(recursivepath([ISSM_DIR '/externalpackages/scotch']));
addpath(recursivepath([ISSM_DIR '/externalpackages/canos']));
addpath(recursivepath([ISSM_DIR '/externalpackages/kml']));
addpath(recursivepath([ISSM_DIR '/externalpackages/export_fig']));
addpath(recursivepath([ISSM_DIR '/externalpackages/googleearthtoolbox']));
addpath(recursivepath([ISSM_DIR '/externalpackages/howatmask']));
addpath(recursivepath([ISSM_DIR '/externalpackages/dem']));
addpath(recursivepath([ISSM_DIR '/externalpackages/mealpix']));
addpath(recursivepath([ISSM_DIR '/externalpackages/pcatool']));

%Check on any warning messages that might indicate that the paths were not correct. 
if ~isempty(lastwarn),
	fprintf('\n  Error trying to setup ''ISSM'' code paths. Try and update the ISSM_DIR variable in your .cshrc or .bashrc!\n');
	fprintf('  ''ISSM'' will not  work at all until this is resolved\n\n');
else
	fprintf('\n  ISSM development path correctly loaded\n\n');
end

warning ('off','all');
addpath([ISSM_DIR '/lib-precompiled']); %load MEX files (precompiled; remove after MEX file compilation is supported on Silicon-based Macs)
warning ('on','all');
clear ISSM_DIR;

%disable matlab bell!
beep off;
