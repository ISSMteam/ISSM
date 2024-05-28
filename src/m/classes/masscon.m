%MASSCON class definition
%
%   Usage:
%      masscon=masscon();
%      masscon=masscon('name','MassCon58','definitionstring','Outputdefinition1',... %name of a North-East Greenland JPL MassCon
%                    'levelset',level);
% 
%   where level is a levelset vectorial field.
%
%   See also: MISFIT, MASSCONAXPBY, REGIONALOUTPUT

classdef masscon
	properties (SetAccess=public)
		%masscon
		name              = '';
		definitionstring  = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-10]'
		levelset          = NaN; %levelset vectorial field which identifies the boundaries of the masscon
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.levelset)
				self.levelset=project3d(md,'vector',self.levelset,'type','node');
			end
		end % }}}
		function self = masscon(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.levelset=getfieldvalue(options,'levelset');
				self.definitionstring=getfieldvalue(options,'definitionstring');

			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('masscon error message: ''name'' field should be a string!');
			end

			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end

			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);
			md = checkfield(md,'fieldname','self.levelset','field',self.levelset,'timeseries',1,'NaN',1,'Inf',1);

		end % }}}
		function md = disp(self) % {{{

			disp(sprintf('   Masscon:\n'));

			fielddisplay(self,'name','identifier for this masscon response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'levelset','levelset vectorial field which identifies the boundaries of the masscon');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'data',self.name,'name','md.masscon.name','format','String');
			WriteData(fid,prefix,'data',self.definitionstring,'name','md.masscon.definitionstring','format','String');
			WriteData(fid,prefix,'data',self.levelset,'name','md.masscon.levelset','format','DoubleMat','mattype',1);

		end % }}}
	end
end
