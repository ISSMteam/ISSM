%RADAR class definition
%
%   Usage:
%      radar=radar();
%      radar=radar('name','Radar1',...
%						 'definitionname','Outputdefinition1',...
%						 'ice_period', ones(md.mesh.numberofvertices,1));

classdef radar
	properties (SetAccess=private)  
		%radarattenuation
		name						 = '';
		definitionstring		 = '';
	end

	methods
		function self = extrude(self,md) % {{{
			return;
		end % }}}
		function self = radar(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			
			if ~ischar(self.name),
				error('radar error message: ''name'' field should be a string!');
			end	
			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'field',self.definitionstring,'values',OutputdefinitionStringArray);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Radar:\n'));

			fielddisplay(self,'name','identifier for this radar response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-100]''');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
		WriteData(fid,prefix,'data',self.name,'name','md.radar.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring, 'name','md.radar.definitionstring','format','String');
	end % }}}
	end
end
