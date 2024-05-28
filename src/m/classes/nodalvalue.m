%NODALVALUE class definition
%
%   Usage:
%      nodalvalue=nodalvalue();
%      nodalvalue=nodalvalue('name','SealevelchangeSNodalValue',...
%                    'definitionstring','Outputdefinition1', 
%                    'model_string','SealevelchangeS',
%                    'node',1);

classdef nodalvalue
	properties (SetAccess=public)
		%nodalvalue
		name              = '';
		definitionstring   = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-10]'
		model_string      = ''; %string for field that is being retrieved
		node             = NaN; %for which node are we retrieving the value?
	end
	
	methods
		function self = nodalvalue(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.model_string=getfieldvalue(options,'model_string');
				self.node=getfieldvalue(options,'node',NaN);

			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('nodalvalue error message: ''name'' field should be a string!');
			end
			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);
			md = checkfield(md,'fieldname','self.node','field',self.node,'values',[1:md.mesh.numberofvertices]);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   Nodalvalue:\n'));

			fielddisplay(self,'name','identifier for this nodalvalue response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'model_string','string for field that is being retrieved');
			fielddisplay(self,'node','vertex index at which we retrieve the value');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'data',self.name,'name','md.nodalvalue.name','format','String');
			WriteData(fid,prefix,'data',self.definitionstring,'name','md.nodalvalue.definitionenum','format','String');
			WriteData(fid,prefix,'data',self.model_string,'name','md.nodalvalue.model_enum','format','String');
			WriteData(fid,prefix,'data',self.node,'name','md.nodalvalue.node','format','Integer');

		end % }}}
	end
end
