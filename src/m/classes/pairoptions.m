%PAIROPTIONS class definition
%
%   Usage:
%      pairoptions=pairoptions();
%      pairoptions=pairoptions('module',true,'solver',false);

classdef pairoptions < matlab.mixin.Copyable
	properties (SetAccess = private,GetAccess = private) 
		list = cell(0,3);
	end
	methods
		function self = pairoptions(varargin) % {{{

			%initialize list
			if nargin==0,
				%Do nothing,
			else
				self=buildlist(self,varargin{:});
			end
		end % }}}
		function self = buildlist(self,varargin) % {{{
		%BUILDLIST - build list of obj from input

			%check length of input
			if mod((nargin-1),2)
				error('Invalid parameter/value pair arguments') 
			end
			numoptions = (nargin-1)/2;

			%Allocate memory
			self.list=cell(numoptions,3);

			%go through varargin and build list of obj
			for i=1:numoptions
				if ischar(varargin{2*i-1})
					self.list{i,1}=varargin{2*i-1};
					self.list{i,2}=varargin{2*i};
					self.list{i,3}=false; %used?
				else
					%option is not a string, ignore it
					disp(['WARNING: option number ' num2str(i) ' is not a string, it will be ignored']);
					self.list(i,:)=[];
					continue
				end
			end
		end % }}}
		function self = addfield(self,field,value) % {{{
			if ischar(field),
				self.list{end+1,1} = field;
				self.list{end,2}   = value;
				self.list{end,3}   = false;
			end
		end % }}}
		function self = addfielddefault(self,field,value) % {{{
		%ADDFIELDDEFAULT - add a field to an options list if it does not exist
			if ischar(field),
				if ~exist(self,field),
					self.list{end+1,1} = field;
					self.list{end,2}   = value;
					self.list{end,3}   = true;  %It is a default so user will not be notified if not used
				end
			end
		end % }}}
		function obj2 = AssignObjectFields(options,obj2) % {{{
		%ASSIGNOBJECTFIELDS - assign object fields from options
			listproperties=properties(obj2);
			for i=1:size(options.list,1),
				fieldname=options.list{i,1};
				fieldvalue=options.list{i,2};
				if ismember(fieldname,listproperties),
					obj2.(fieldname)=fieldvalue;
				else
					disp(['WARNING: ''' fieldname ''' is not a property of ''' class(obj2) '''']);
				end
			end
		end % }}}
		function self = changefieldvalue(self,field,newvalue) % {{{
		%CHANGEOPTIONVALUE - change the value of an option in an option list

			%track occurrence of field
			lines=find(strcmpi(self.list(:,1),field));

			%replace value
			if isempty(lines),
				%add new field if not found
				self=addfield(self,field,newvalue);
				self.list{end,3}=true; % do not notify user if unused
			else
				for i=1:length(lines),
					self.list{lines(i),2}=newvalue;
				end
			end
		end % }}}
		function self = deleteduplicates(self,warn) % {{{
		%DELETEDUPLICATES - delete duplicates in an option list

			%track the first occurrence of each option
			[dummy lines]=unique(self.list(:,1),'first');
			clear dummy

			%warn user if requested
			if warn,
				numoptions=size(self.list,1);
				for i=1:numoptions,
					if ~ismember(i,lines),
						disp(['WARNING: option ' self.list{i,1} ' appeared more than once. Only its first occurrence will be kept'])
					end
				end
			end

			%remove duplicates from the options list
			self.list=self.list(lines,:);
		end % }}}
		function displayunused(self) % {{{
			%DISPLAYUNUSED - display unused options

			numoptions=size(self.list,1);
			for i=1:numoptions,
				if ~self.list{i,3},
					disp(['WARNING: option ' self.list{i,1} ' was not used'])
				end
			end
		end % }}}
		function disp(self) % {{{
			if ~isempty(self.list),
				disp(sprintf('   list: (%ix%i)\n',size(self.list,1),size(self.list,2)));
				for i=1:size(self.list,1),
					if ischar(self.list{i,2}),
						disp(sprintf('     field: %-10s value: ''%s''',self.list{i,1},self.list{i,2}));
					elseif isnumeric(self.list{i,2}) & length(self.list{i,2})==1,
						disp(sprintf('     field: %-10s value: %g',self.list{i,1},self.list{i,2}));
					elseif isnumeric(self.list{i,2}) & length(self.list{i,2})==2,
						disp(sprintf('     field: %-10s value: [%g %g]',self.list{i,1},self.list{i,2}));
					else
						disp(sprintf('     field: %-10s value: (%ix%i)',self.list{i,1},size(self.list{i,2},1),size(self.list{i,2},2)));
					end
				end
			else
				disp(sprintf('   list: empty'));
			end
		end % }}}
		function bool = exist(self,field) % {{{
		%EXIST - check if the option exists

			%some argument checking: 
			if ((nargin~=2) | (nargout~=1)),
				error('exist error message: bad usage');
			end
			if ~ischar(field),
				error('exist error message: field should be a string');
			end

			%Recover option
			pos=find(strcmpi(field,self.list(:,1)));
			if ~isempty(pos),
				bool=true;
				self.list{pos,3}   = true;  %It is a default so user will not be notified if not used
			else
				bool=false;
			end
		end % }}}
		function num = fieldoccurrences(self,field), % {{{
		%FIELDOCCURRENCES - get number of occurrence of a field

			%check input 
			if ~ischar(field),
				error('fieldoccurrences error message: field should be a string');
			end

			%get number of occurrence
			num=sum(strcmpi(field,self.list(:,1)));
		end % }}}
		function value = getfieldvalue(self,field,varargin), % {{{
		%GETFIELDVALUE - get the value of an option
		%
		%   Usage:
		%      value=getfieldvalue(self,field,varargin)
		%
		%   Find an option value from a field. A default option
		%   can be given in input if the field does not exist
		%
		%   Examples:
		%      value=getfieldvalue(options,'caxis');
		%      value=getfieldvalue(options,'caxis',[0 2]);

			%some argument checking: 
			if nargin~=2 && nargin~=3,
				help getfieldvalue
				error('getfieldvalue error message: bad usage');
			end

			if ~ischar(field),
				error('getfieldvalue error message: field should be a string');
			end

			%Recover option
			pos=find(strcmpi(self.list(:,1),field));
			if ~isempty(pos),
				value=self.list{pos(1),2}; % ignore extra entry
				self.list{pos(1),3}=true;  % option used
				return;
			end

			%The option has not been found, output default if provided
			if nargin==3,
				value=varargin{1};
			else
				error(['error message: field ' field ' has not been provided by user (and no default value has been specified)'])
			end
		end % }}}
		function values = getfieldvalues(self,field,varargin), % {{{
		%GETFIELDVALUES - get the value of an option (if the option is repeated, return multiple values)
		%
		%   Usage:
		%      values=getfieldvalues(self,field,varargin)
		%
		%   Find all option values from a field. Default options
		%   can be given in input if the field does not exist
		%
		%   Examples:
		%      values=getfieldvalues(options,'caxis');
		%      values=getfieldvalues(options,'caxis',{[0 2],[3 4]});

			%some argument checking: 
			if nargin~=2 && nargin~=3,
				help getfieldvalues
				error('getfieldvalues error message: bad usage');
			end

			if ~ischar(field),
				error('getfieldvalues error message: field should be a string');
			end

			%Recover options
			pos=find(strcmpi(self.list(:,1),field));
			if ~isempty(pos),
				values={};
				for i=1:length(pos),
					values{i}=self.list{pos(i),2};
					self.list{pos(i),3}=true; % option used
				end
				return;
			end

			%The option has not been found, output default if provided
			if nargin==3,
				values=varargin{1};
			else
				error(['error message: field ' field ' has not been provided by user (and no default value has been specified)'])
			end
		end % }}}
		function self = removefield(self,field,warn)% {{{
		%REMOVEFIELD - delete a field in an option list
		%
		%   Usage:
		%      self=removefield(self,field,warn)
		%
		%   if warn==1 display an info message to warn user that
		%   some of their options have been removed.

			%check is field exist
			if exist(self,field),

				%find where the field is located
				lines=find(~strcmpi(self.list(:,1),field));

				%remove duplicates from the options list
				self.list=self.list(lines,:);

				%warn user if requested
				if warn
					disp(['removefield info: option ' field ' has been removed from the list of options.'])
				end
			end
		end % }}}
		function marshall(self,fid)% {{{

			for i=1:size(self.list,1),
				name  = self.list{i,1};
				value = self.list{i,2};

				%Write option value
				if (isnumeric(value) & numel(value)==1),
					WriteData(fid,'','name',['md.' name],'data',value,'format','Double');
				elseif ischar(value),
					WriteData(fid,'','name',['md.' name],'data',value,'format','String');
				else
					error(['Cannot marshall option ' name ': format not supported yet']);
				end
			end
		end % }}}
	end
end
