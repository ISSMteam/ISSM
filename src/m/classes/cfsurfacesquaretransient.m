%MISFIT class definition
%
%   Usage:
%      cfsurfacesquaretransient=cfsurfacesquaretransient();
%      cfsurfacesquaretransient=cfsurfacesquaretransient('name','SurfaceAltimetry',...
%                    'definitionstring','Outputdefinition1',... 
%							'model_string','Surface',...
%                    'observations',[md.geometry.surface;0],...
%                    'weights',ones(md.mesh.numberofvertices+1,1));
%
%

classdef cfsurfacesquaretransient
	properties (SetAccess=public)
		%cfsurfacesquaretransient
		name                = '';
		definitionstring    = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-100]'
		model_string        = ''; %string for field that is modeled
		observations        = NaN;%observed field that we compare the model against
		weights             = NaN;%weight coefficients for every vertex
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.weights)
				self.weights=project3d(md,'vector',self.weights,'type','node');
			end
			if ~isnan(self.observations)
				self.observations=project3d(md,'vector',self.observations,'type','node');
			end
		end % }}}
		function self = cfsurfacesquaretransient(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name                = getfieldvalue(options,'name','');
				self.definitionstring    = getfieldvalue(options,'definitionstring');
				self.model_string        = getfieldvalue(options,'model_string');
				self.observations        = getfieldvalue(options,'observations',NaN);
				self.weights             = getfieldvalue(options,'weights',NaN);
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('cfsurfacesquaretransient error message: ''name'' field should be a string!');
			end
			OutputdefinitionStringArray={};
			for i=1:2000
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);
			md = checkfield(md,'fieldname','self.observations','field',self.observations,'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','self.weights','field',self.weights,'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   cfsurfacesquaretransient:\n'));

			fielddisplay(self,'name','identifier for this cfsurfacesquaretransient response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'model_string','string for field that is modeled');
			fielddisplay(self,'observations','observed field that we compare the model against');
			fielddisplay(self,'weights','weights (at vertices) to apply to the cfsurfacesquaretransient');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

		WriteData(fid,prefix,'data',self.name,'name','md.cfsurfacesquaretransient.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring,'name','md.cfsurfacesquaretransient.definitionstring','format','String');
		WriteData(fid,prefix,'data',self.model_string,'name','md.cfsurfacesquaretransient.model_string','format','String');
		WriteData(fid,prefix,'data',self.observations,'name','md.cfsurfacesquaretransient.observations','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'data',self.weights,'name','md.cfsurfacesquaretransient.weights','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		end % }}}
	end
end
