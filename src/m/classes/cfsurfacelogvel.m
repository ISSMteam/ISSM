%MISFIT class definition
%
%   Usage:
%      cfsurfacelogvel=cfsurfacelogvel();
%      cfsurfacelogvel=cfsurfacelogvel('name','SurfaceAltimetry',...
%                    'definitionstring','Outputdefinition1',... 
%                    'observation_string','SurfaceObservations',...
%                    'observation',md.geometry.surface,...
%                    'weights',ones(md.mesh.numberofvertices,1),...
%                    'weights_string','WeightsSurfaceObservations',...
%							'datatime',time);
%
%

classdef cfsurfacelogvel
	properties (SetAccess=public)
		%cfsurfacelogvel
		name               = '';
		definitionstring   = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-100]'
		vxobs					 = NaN; %observed field that we compare the model against
		vxobs_string		 = ''; %string for observed field.
		vyobs			       = NaN; %observed field that we compare the model against
		vyobs_string		 = ''; %string for observed field.
		weights            = NaN; %weight coefficients for every vertex
		weights_string     = ''; %string to identify this particular set of weights
		datatime				 = 0; %time in years from start that the data is from 
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.weights)
				self.weights=project3d(md,'vector',self.weights,'type','node');
			end
			if ~isnan(self.vxobs)
				self.vxobs=project3d(md,'vector',self.vxobs,'type','node');
			end
		end % }}}
		function self = cfsurfacelogvel(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.vxobs=getfieldvalue(options,'vxobs',NaN);
				self.vyobs=getfieldvalue(options,'vyobs',NaN);
				self.vxobs_string=getfieldvalue(options,'vxobs_string');			
				self.vyobs_string=getfieldvalue(options,'vyobs_string');
				self.weights=getfieldvalue(options,'weights',NaN);
				self.weights_string=getfieldvalue(options,'weights_string','');
				self.datatime = getfieldvalue(options, 'datatime');

			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.datatime = 0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('cfsurfacelogvel error message: ''name'' field should be a string!');
			end
			OutputdefinitionStringArray={};
			for i=1:2000
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);

			md = checkfield(md,'fieldname','self.vxobs','field',self.vxobs,'timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','self.weights','field',self.weights,'timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','self.datatime','field',self.datatime,'<=',md.timestepping.final_time);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   cfsurfacelogvel:\n'));

			fielddisplay(self,'name','identifier for this cfsurfacelogvel response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-10]''');
			fielddisplay(self,'vxobs','observed field that we compare the model against');
			fielddisplay(self,'vxobs_string','observation string');
			fielddisplay(self,'weights','weights (at vertices) to apply to the cfsurfacelogvel');
			fielddisplay(self,'weights_string','string for weights for identification purposes');
			fielddisplay(self,'datatime','time to compare data to model for misfit');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

		WriteData(fid,prefix,'data',self.name,'name','md.cfsurfacelogvel.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring,'name','md.cfsurfacelogvel.definitionstring','format','String');
		WriteData(fid,prefix,'data',self.vxobs,'name','md.cfsurfacelogvel.vxobs','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
		WriteData(fid,prefix,'data',self.vxobs_string,'name','md.cfsurfacelogvel.vxobs_string','format','String');
		WriteData(fid,prefix,'data',self.vyobs,'name','md.cfsurfacelogvel.vyobs','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
		WriteData(fid,prefix,'data',self.vyobs_string,'name','md.cfsurfacelogvel.vyobs_string','format','String');
		WriteData(fid,prefix,'data',self.weights,'name','md.cfsurfacelogvel.weights','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'data',self.weights_string,'name','md.cfsurfacelogvel.weights_string','format','String');
		WriteData(fid,prefix,'data',round(self.datatime*md.constants.yts),'name','md.cfsurfacelogvel.datatime','format','Double');
		end % }}}
	end
end
