%MISFIT class definition
%
%   Usage:
%      misfit=misfit();
%      misfit=misfit('name','SurfaceAltimetry',...
%                    'definitionstring','Outputdefinition1',... 
%							'model_string','Surface',...
%                    'observation_string','SurfaceObservations',...
%                    'observation',md.geometry.surface,...
%                    'timeinterpolation','nearestneighbor',...
%                    'local',1,...
%                    'weights',ones(md.mesh.numberofvertices,1),...
%                    'weights_string','WeightsSurfaceObservations');
%
%

classdef misfit
	properties (SetAccess=public)
		%misfit
		name               = '';
		definitionstring   = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-100]'
		model_string       = ''; %string for field that is modeled
		observation        = NaN; %observed field that we compare the model against
		observation_string = ''; %string for observed field.
		timeinterpolation  = '';
		local              = 1;
		weights            = NaN; %weight coefficients for every vertex
		weights_string     = ''; %string to identify this particular set of weights
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.weights)
				self.weights=project3d(md,'vector',self.weights,'type','node');
			end
			if ~isnan(self.observation)
				self.observation=project3d(md,'vector',self.observation,'type','node');
			end
		end % }}}
		function self = misfit(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.model_string=getfieldvalue(options,'model_string');
				self.observation=getfieldvalue(options,'observation',NaN);
				self.observation_string=getfieldvalue(options,'observation_string');
				self.local=getfieldvalue(options,'local',1);
				self.timeinterpolation=getfieldvalue(options,'timeinterpolation','nearestneighbor');
				self.weights=getfieldvalue(options,'weights',NaN);
				self.weights_string=getfieldvalue(options,'weights_string','');

			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.local=1;
			self.timeinterpolation='nearestneighbor';
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('misfit error message: ''name'' field should be a string!');
			end
			OutputdefinitionStringArray={};
			for i=1:2000
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);

			if ~ischar(self.timeinterpolation),
				error('misfit error message: ''timeinterpolation'' field should be a string!');
			end
			md = checkfield(md,'fieldname','self.observation','field',self.observation,'timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','self.timeinterpolation','field',self.timeinterpolation,'values',{'nearestneighbor'});
			md = checkfield(md,'fieldname','self.weights','field',self.weights,'timeseries',1,'NaN',1,'Inf',1);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   Misfit:\n'));

			fielddisplay(self,'name','identifier for this misfit response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'model_string','string for field that is modeled');
			fielddisplay(self,'observation','observed field that we compare the model against');
			fielddisplay(self,'observation_string','observation string');
			fielddisplay(self,'local','is the response local to the elements, or global? (default is 1)''');
			fielddisplay(self,'timeinterpolation','interpolation routine used to interpolate misfit between two time steps (default is ''nearestneighbor''');
			fielddisplay(self,'weights','weights (at vertices) to apply to the misfit');
			fielddisplay(self,'weights_string','string for weights for identification purposes');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

		WriteData(fid,prefix,'data',self.name,'name','md.misfit.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring,'name','md.misfit.definitionstring','format','String');
		WriteData(fid,prefix,'data',self.model_string,'name','md.misfit.model_string','format','String');
		WriteData(fid,prefix,'data',self.observation,'name','md.misfit.observation','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'data',self.observation_string,'name','md.misfit.observation_string','format','String');
		WriteData(fid,prefix,'data',self.local,'name','md.misfit.local','format','Integer');
		WriteData(fid,prefix,'data',self.timeinterpolation,'name','md.misfit.timeinterpolation','format','String');
		WriteData(fid,prefix,'data',self.weights,'name','md.misfit.weights','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'data',self.weights_string,'name','md.misfit.weights_string','format','String');

		end % }}}
	end
end
