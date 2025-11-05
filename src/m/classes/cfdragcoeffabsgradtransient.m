%CFDRAGCOEFFABSGRADTRANSIENT class definition
%
%   Usage:
%      cfdragcoeffabsgradtransient=cfdragcoeffabsgradtransient();
%      cfdragcoeffabsgradtransient=cfdragcoeffabsgradtransient('name','SurfaceAltimetry',...
%                    'definitionstring','Outputdefinition1',... 
%                    'weights',ones(md.mesh.numberofvertices+1,1));
%
%

classdef cfdragcoeffabsgradtransient
	properties (SetAccess=public)
		%cfdragcoeffabsgrad
		name               = '';
		definitionstring   = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-100]'
		weights            = NaN; %weight coefficients for every vertex
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.weights)
				self.weights=project3d(md,'vector',self.weights,'type','node');
			end
		end % }}}
		function self = cfdragcoeffabsgradtransient(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.weights=getfieldvalue(options,'weights',NaN);
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('cfdragcoeffabsgrad error message: ''name'' field should be a string!');
			end
			OutputdefinitionStringArray={};
			for i=1:2000
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);
			md = checkfield(md,'fieldname','self.weights','field',self.weights,'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   cfdragcoeffabsgradtransient:\n'));

			fielddisplay(self,'name','identifier for this cfdragcoeffabsgradtransient response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'weights','weights (at vertices) to apply to the cfdragcoeffabsgradtransient');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

		WriteData(fid,prefix,'data',self.name,'name','md.cfdragcoeffabsgradtransient.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring,'name','md.cfdragcoeffabsgradtransient.definitionstring','format','String');
		WriteData(fid,prefix,'data',self.weights,'name','md.cfdragcoeffabsgradtransient.weights','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		end % }}}
	end
end
