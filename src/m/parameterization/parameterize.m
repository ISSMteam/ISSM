function md=parameterize(md,parametername)
%PARAMETERIZE - parameterize a model
%
%   from a parameter MATLAB file, start filling in all the @model fields that were not 
%   filled in by the mesh.m and mask.m @model methods.
%   Warning: the parameter file must be able to be run in MATLAB
%
%   Usage:
%      md=parameterize(md,parametername)
%
%   Example:
%      md=parameterize(md,'Square.par');

%some checks
if ~exist(parametername),
	error(['parameterize error message: file ' parametername ' not found!']);
end

%Try and run parameter file.
temporaryname=['TemporaryParameterFile' num2str(feature('GetPid')) ];
copyfile(parametername,[temporaryname '.m']);

%WARNING: this is a bug of matlab: the TemporaryParameterFile must be cleared
%otherwise matlab keeps the previous version of this file which is not what
%we want!!!
eval(['clear ' temporaryname]);

try,
	eval(temporaryname);
	delete([temporaryname '.m']);
catch me,
	delete([temporaryname '.m']);

	%copy error message
	me2=struct('message',me.message,'stack',me.stack);

	%rename parameter file
	for i=1:length(me2.stack)-1,
		if strcmp(me2.stack(i).name,temporaryname)
			me2.stack(i).file = strrep(me2.stack(i).file,temporaryname,parametername);
			me2.stack(i).name = parametername;
		end
		if strcmp(me2.stack(i).name,'parameterize'),
			%remove error (eval(temporaryname);) misleading
			me2.stack(i)=[];
		end
	end

	%throw error message
	rethrow(me2);
end

%Name and notes
if isempty(md.miscellaneous.name), 
	[path,root,ext]=fileparts(parametername);
	md.miscellaneous.name=root; 
end
if isempty(md.miscellaneous.notes), 
	md.miscellaneous.notes=['Model created by using parameter file: ' parametername ' on: ' datestr(now)];
end
