%MISMIP BASAL FORCINGS class definition
%
%   Usage:
%      mismipbasalforcings=mismipbasalforcings();

classdef mismipbasalforcings
	properties (SetAccess=public) 
		groundedice_melting_rate  = NaN;
		meltrate_factor           = NaN;
		threshold_thickness       = 0.;
		upperdepth_melt           = 0.;
		geothermalflux            = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1); 
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','node','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = mismipbasalforcings(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(mismipbasalforcings(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%default values for melting parameterization
			self.meltrate_factor        = 0.2;
			self.threshold_thickness    = 75;
			self.upperdepth_melt        = -100;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.ismasstransport==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'size','universal','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'size','universal','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',1);
			end
			if ismember('ThermalAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.isthermal==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'size','universal','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			end
			if isnan(md.geometry.bed),
				md = checkmessage(md,['requesting mismip basal melting parameterization, but bathymetry is absent!']);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   MISMIP+ basal melt parameterization:'));

			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]');
			fielddisplay(self,'meltrate_factor','Melt-rate rate factor [1/yr] (sign is opposite to MISMIP+ benchmark to remain consistent with ISSM convention of positive values for melting)');
			fielddisplay(self,'threshold_thickness','threshold thickness for saturation of basal melting [m]');
			fielddisplay(self,'upperdepth_melt','depth above which the melt rate is zero [m]');
			fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;
			if yts~=365.2422*24.*3600.
				disp('WARNING: value of yts for MISMIP+ runs different from ISSM default!');
			end

			WriteData(fid,prefix,'name','md.basalforcings.model','data',3,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','name','md.basalforcings.groundedice_melting_rate','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','name','md.basalforcings.geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','meltrate_factor','format','DoubleMat','mattype',1,'scale',1./yts)
			WriteData(fid,prefix,'object',self,'fieldname','threshold_thickness','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','upperdepth_melt','format','Double')
		end % }}}
	end
end
