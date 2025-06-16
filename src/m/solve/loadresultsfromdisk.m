function md=loadresultsfromdisk(md,filename)
%LOADRESULTSFROMDISK - load results of solution sequence from disk file "filename"            
%
%   Usage:
%      md=loadresultsfromdisk(md,filename);

%check number of inputs/outputs
if ((nargin~=2) | (nargout~=1)),
	help loadresultsfromdisk;
	error('loadresultsfromdisk: error message.');
end

if ~md.qmu.isdakota
	%Check that file exists
	if ~exist(filename,'file')
		error(sprintf(['\n'...
			'=========================================================================\n'...
			'   Binary file ' filename ' not found                                    \n'... 
			'                                                                         \n'...
			'   This typically results from an error encountered during the simulation\n'...
			'   Please check for error messages above or in the outlog                \n'...
			'=========================================================================\n'...
			]));
		return;
	end

	%initialize md.results if not a structure yet
	if ~isstruct(md.results)
		md.results=struct();
	end

	%load results onto model
	structure=parseresultsfromdisk(md,filename,~md.settings.io_gather);
	if isempty(fieldnames(structure))
		error(['No result found in binary file ' filename '. Check for solution crash.']);
	end
	if isempty(structure(1).SolutionType),
		if ~isempty(structure(end).SolutionType)
			structure(1).SolutionType=structure(end).SolutionType;
		else
			warning(['Cannot find a solution type in the results! Ascribing one: ''NoneSolution''.']);
			structure(1).SolutionType='NoneSolution';
		end
	end
	md.results.(structure(1).SolutionType)=structure;

	%recover solution_type from results
	md.private.solution=structure(1).SolutionType;

	%read log files onto fields (only keep the first 1000 lines!)
	if exist([md.miscellaneous.name '.errlog'],'file')
		md.results.(structure(1).SolutionType)(1).errlog=char(textread([md.miscellaneous.name '.errlog'],'%s',1000,'delimiter','\n'));
	else
		md.results.(structure(1).SolutionType)(1).errlog='';
	end

	if exist([md.miscellaneous.name '.outlog'],'file')
		md.results.(structure(1).SolutionType)(1).outlog=char(textread([md.miscellaneous.name '.outlog'],'%c',4000,'delimiter','\n'));
	else
		md.results.(structure(1).SolutionType)(1).outlog='';
	end

	if ~isempty(md.results.(structure(1).SolutionType)(1).errlog)
		disp(['loadresultsfromdisk info message: error during solution. Check your errlog and outlog model fields']);
	end

%postprocess qmu results if necessary
else
	md=postqmu(md);
end
