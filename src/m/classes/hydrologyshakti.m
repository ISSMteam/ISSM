%HYDROLOGYSHAKTI class definition
%
%   Usage:
%      hydrologyshakti=hydrologyshakti();

classdef hydrologyshakti
	properties (SetAccess=public) 
		head            = NaN;
		gap_height      = NaN;
		gap_height_min  = 1e-3;
		gap_height_max  = 1.;
		bump_spacing    = NaN;
		bump_height     = NaN;
		englacial_input = NaN;
		moulin_input    = NaN;
		reynolds        = NaN;
		spchead         = NaN;
		neumannflux     = NaN;
		relaxation      = 0;
		storage         = NaN;
		requested_outputs = {};
	end
	methods
		function self = extrude(self,md) % {{{
			self.head = project3d(md, 'vector', self.head , 'type', 'node');
			self.gap_height = project3d(md, 'vector', self.gap_height, 'type', 'element');
			self.bump_spacing = project3d(md, 'vector', self.bump_spacing, 'type', 'element');
			self.bump_height = project3d(md, 'vector', self.bump_height, 'type', 'element');
			self.englacial_input = project3d(md, 'vector', self.englacial_input, 'type', 'node');
			self.moulin_input = project3d(md, 'vector', self.moulin_input, 'type', 'node');
			self.reynolds = project3d(md, 'vector', self.reynolds, 'type', 'element');
			self.neumannflux = project3d(md, 'vector', self.neumannflux, 'type', 'element');
			self.spchead = project3d(md, 'vector', self.spchead, 'type', 'node');
		end % }}}
		function self = hydrologyshakti(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(self,varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'HydrologyHead','HydrologyGapHeight','EffectivePressure','HydrologyBasalFlux','DegreeOfChannelization'};
		end % }}}    

		function self = setdefaultparameters(self) % {{{
			% Set under-relaxation parameter to be 1 (no under-relaxation of nonlinear iteration)	
			self.gap_height_min  = 1e-3;
			self.gap_height_max  = 1.;
			self.relaxation=1;
			self.storage=0;
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologyShaktiAnalysis',analyses)
				return;
			end

			md = checkfield(md,'fieldname','hydrology.head','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.gap_height','>=',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.gap_height_min','>=',0,'numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.gap_height_max','>=',0,'numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.bump_spacing','>',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.bump_height','>=',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.englacial_input','NaN',1,'Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.moulin_input','NaN',1,'Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.reynolds','>',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.neumannflux','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.spchead','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.relaxation','>=',0);	
			md = checkfield(md,'fieldname','hydrology.storage','>=',0,'universal',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   hydrologyshakti solution parameters:'));
			fielddisplay(self,'head','subglacial hydrology water head (m)');
			fielddisplay(self,'gap_height','height of gap separating ice to bed (m)');
			fielddisplay(self,'gap_height_min','minimum allowed gap height (m)');
			fielddisplay(self,'gap_height_max','maximum allowed gap height (m)');
			fielddisplay(self,'bump_spacing','characteristic bedrock bump spacing (m)');
			fielddisplay(self,'bump_height','characteristic bedrock bump height (m)');
			fielddisplay(self,'englacial_input','liquid water input from englacial to subglacial system (m/yr)');
			fielddisplay(self,'moulin_input','liquid water input from moulins (at the vertices) to subglacial system (m^3/s)');
			fielddisplay(self,'reynolds','Reynolds'' number');
			fielddisplay(self,'neumannflux','water flux applied along the model boundary (m^2/s)');
			fielddisplay(self,'spchead','water head constraints (NaN means no constraint) (m)');
			fielddisplay(self,'relaxation','under-relaxation coefficient for nonlinear iteration');
			fielddisplay(self,'storage','englacial storage coefficient (void ratio)');
			fielddisplay(self,'requested_outputs','additional outputs requested');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.hydrology.model','data',3,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','head','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','gap_height','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','gap_height_min','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','gap_height_max','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_spacing','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_height','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','englacial_input','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','moulin_input','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','reynolds','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','neumannflux','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','spchead','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','relaxation','format','Double');
			if(size(self.storage,1)==md.mesh.numberofvertices | size(self.storage,1)==md.mesh.numberofvertices+1 | (size(self.storage,1)==md.mesh.numberofelements && size(self.storage,2)>1))
				mattype=1; tsl = md.mesh.numberofvertices;
			else
				mattype=2; tsl = md.mesh.numberofelements;
			end
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','storage','format','DoubleMat','mattype',mattype,'timeserieslength',tsl+1,'yts',md.constants.yts);

			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];  %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}
	end
end

