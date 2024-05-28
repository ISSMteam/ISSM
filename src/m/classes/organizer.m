%ORGANIZER class definition
%
%   Supported options:
%      repository: directory where all models will be saved
%      prefix:     prefix for saved model names
%      steps:      requested steps
%      color:      color of step title (default is '41;37')
%
%   Usage:
%      org = organizer(varargin)
%
%   Examples:
%      org = organizer('repository','Models/','prefix','AGU2015','steps',0);  %build an empty organizer object with a given repository

classdef organizer < handle
	properties (SetAccess=private) 
		% {{{
		currentstep   =0;
	end
    properties (SetAccess=public) 
		repository    ='';
		prefix        ='';
		color         ='';
		steps         =[];
		skipio        = false;
		requestedsteps=[0];
		%}}}
	end
	methods
		function org=organizer(varargin) % {{{

			%process options
			options=pairoptions(varargin{:});

			%Get prefix
			prefix=getfieldvalue(options,'prefix','model_');
			if ~ischar(prefix),                            error('prefix is not a string'); end
			if ~strcmp(regexprep(prefix,'\s+',''),prefix), error('prefix should not have any white space'); end
			org.prefix=prefix;

			%Get repository
			repository=getfieldvalue(options,'repository','./');
			if ~ischar(repository),        error('repository is not a string'); end
			if exist(repository,'dir')~=7, error(['Directory ' repository ' not found']), end
			org.repository=repository;

			%Color
			org.color=getfieldvalue(options,'color','41;37');

			%Get steps
			org.requestedsteps=getfieldvalue(options,'steps',0);

			%Skip io?
			org.skipio=getfieldvalue(options,'skipio',0);
		end
		%}}}
		function reset(org,varargin) % {{{
			
			%process options
			options=pairoptions(varargin{:});
		
			%Zero out some fields: 
			org.currentstep   =0;
			org.steps=[];
	
			%get requested step: 
			org.requestedsteps=getfieldvalue(options,'steps',0);
			
		end
		%}}}
		function disp(org) % {{{
			disp(sprintf('   Repository: ''%s''',org.repository));
			disp(sprintf('   Prefix:     ''%s''\n',org.prefix));
			disp(sprintf('   Color:      ''%s''\n',org.color));
			disp(sprintf('   skipio:     %i\n',org.skipio));
			if isempty(org.steps)
				disp('   no step');
			else
				for i=1:length(org.steps),
					disp(sprintf('   step #%2i: ''%s''',org.steps(i).id,org.steps(i).string));
				end
			end
		end
		%}}}
		function md=load(org,string),% {{{

			%Get model path
			if ~ischar(string), error('argument provided is not a string'); end
			path=[org.repository '/' org.prefix string];

			%Skip if requested
			if org.skipio,
				disp(['WARNING: Skipping loading ' path]);
				md = evalin('base', 'md');
				return;
			end

			%figure out if the model is there
			if exist(path,'file'),
				path=path;
			elseif exist([path '.mat'],'file'),
				path=[path '.mat'];
			else
				error(['Could not find ' path ]);
			end

			struc=load(path,'-mat');
			name=char(fieldnames(struc));
			md=struc.(name);
			if nargout,
				varargout{1}=md;
			end
		end%}}}
		function md=loadmodel(org,string),% {{{

			%Get model path
			if ~ischar(string), error('argument provided is not a string'); end
			path=[org.repository '/' org.prefix string];

			%Skip if requested
			if org.skipio,
				disp(['WARNING: Skipping loading ' path]);
				md = evalin('base', 'md');
				return;
			end

			%figure out if the model is there, otherwise, we have to use the default path supplied by user.
			if exist(path,'file') | exist([path '.mat'],'file'),
				md=loadmodel(path);
				return;
			end

			%If we are here, the data has not been found. 
			error(['Could not find ' path ]);
		end%}}}
		function loaddata(org,string),% {{{

			%Get model path
			if ~ischar(string), error('argument provided is not a string'); end
			path=[org.repository '/' org.prefix string];

			%figure out if the data is there, otherwise, we have to use the default path supplied by user.
			if exist(path,'file'),
				path=path;
			elseif exist([path '.mat'],'file'),
				path=[path '.mat'];
			else
				error(['Could not find ' path ]);
			end
			if exist(path,'file')
				evalin('caller',['load -mat ' path]);
				return;
			end

			%If we are here, the data has not been found. 
			error(['Could not find ' path ]);
		end%}}}
		function loaddatanoprefix(org,string),% {{{

			%Get model path
			if ~ischar(string), error('argument provided is not a string'); end
			path=[org.repository '/' string];

			%figure out if the data is there, otherwise, we have to use the default path supplied by user.
			if exist(path,'file'),
				path=path;
			elseif exist([path '.mat'],'file'),
				path=[path '.mat'];
			else
				error(['Could not find ' path ]);
			end
			if exist(path,'file')
				evalin('caller',['load -mat ' path]);
				return;
			end

			%If we are here, the data has not been found. 
			error(['Could not find ' path ]);
		end%}}}
		function bool=perform(org,varargin) % {{{

			bool=false;
			
			%group,string are the variable arguments length: 
			if nargin==2,
				string=varargin{1};
			elseif nargin==3,
				string=sprintf('%s.%s',varargin{1},varargin{2});
			end

			%Some checks
			if ~ischar(string),                            error('Step provided should be a string'); end
			if ~strcmp(regexprep(string,'\s+',''),string), error('Step provided should not have any white space'); end
			if (org.currentstep>0 & ismember({string},{org.steps.string})) 
				error(['Step ' string ' already present. Change name']); 
			end

			%Add step
			org.steps(end+1).id=length(org.steps)+1;
			org.steps(end).string=string;
			org.currentstep=org.currentstep+1;

			%if requestedsteps = 0, print all steps in org 
			if any(org.requestedsteps==0),
				if org.currentstep==1,
					disp(sprintf('   prefix: %s',org.prefix));
				end
				disp(sprintf('   step #%2i : %s',org.steps(org.currentstep).id,org.steps(org.currentstep).string));
			end

			%Ok, now if currentstep is a member of steps, return true
			if ismember(org.currentstep,org.requestedsteps),
				if usejava('desktop'),
					disp(sprintf('\n   step #%i : %s\n',org.steps(org.currentstep).id,org.steps(org.currentstep).string));
				else
					%Print on a red background
					fprintf(['\n\033[' org.color 'm   step #' num2str(org.steps(org.currentstep).id) ': ' org.steps(org.currentstep).string '   \033[0m\n\n']);
				end
				bool=true;

				%last check: is this step locked? 
				string=org.steps(org.currentstep).string; 
				if length(string)>7,
					if strcmpi(string(end-5:end),'Locked'),
						error('organizer: you are  trying to run a locked step! Unlock it first!');
					end
				end
			end
		end%}}}
		function savemodel(org,md) % {{{

			%check
			if (org.currentstep==0), error('Cannot save model because organizer (org) is empty! Make sure you did not skip any perform call'); end
			if (org.currentstep>length(org.steps)), error('Cannot save model because organizer (org) is not up to date!'); end

			name=[org.repository '/' org.prefix org.steps(org.currentstep).string ];
			disp(['saving model as: ' name]);

			%Skip if requested
			if org.skipio,
				disp(['WARNING: Skipping saving ' name]);
				return;
			end

			%check that md is a model
			if ~isa(md,'model') & ~isa(md,'sealevelmodel'), warning('second argument is not a model'); end
			if (org.currentstep>length(org.steps)), error(['organizer error message: element with id ' num2str(org.currentstep) ' not found']); end

			%save model
			save(name,'md','-v7.3');
		end%}}}
		function savedata(org,varargin) % {{{

			%check
			if (org.currentstep==0), error('Cannot save data because organizer (org) is empty! Make sure you did not skip any perform call'); end
			if (org.currentstep>length(org.steps)), error('Cannot save data because organizer (org) is not up to date!'); end

			name=[org.repository '/' org.prefix org.steps(org.currentstep).string ];
			disp(['saving data in: ' name]);

			%Skip if requested
			if org.skipio,
				disp(['WARNING: Skipping saving ' name]);
				return;
			end

			%check that md is a model
			if (org.currentstep>length(org.steps)), error(['organizer error message: element with id ' num2str(org.currentstep) ' not found']); end

			%list of variable names: 
			variables='';
			for i=2:nargin, 
				variables=[variables ',' '''' inputname(i) ''''];
				eval([inputname(i) '= varargin{' num2str(i-1) '};']);
			end
			eval(['save(''' name '''' variables ',''-v7.3'');']);
		end%}}}
		function savedatanoprefix(org,varargin) % {{{

			%check
			if (org.currentstep==0), error('Cannot save data because organizer (org) is empty! Make sure you did not skip any perform call'); end
			if (org.currentstep>length(org.steps)), error('Cannot save data because organizer (org) is not up to date!'); end

			name=[org.repository '/' org.steps(org.currentstep).string ];
			disp(['saving data in: ' name]);

			%Skip if requested
			if org.skipio,
				disp(['WARNING: Skipping saving ' name]);
				return;
			end

			%check that md is a model
			if (org.currentstep>length(org.steps)), error(['organizer error message: element with id ' num2str(org.currentstep) ' not found']); end

			%list of variable names: 
			variables='';
			for i=2:nargin, 
				variables=[variables ',' '''' inputname(i) ''''];
				eval([inputname(i) '= varargin{' num2str(i-1) '};']);
			end
			eval(['save(''' name '''' variables ',''-v7.3'');']);
		end%}}}
	end
end
