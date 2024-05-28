%LINEAR BASAL FORCINGS class definition
%
%   Usage:
%      linearbasalforcings=linearbasalforcings();

classdef linearbasalforcings
	properties (SetAccess=public)
		deepwater_melting_rate    = 0.;
		upperwater_melting_rate   = 0.;
		deepwater_elevation       = 0.;
		upperwater_elevation      = 0.;
		groundedice_melting_rate  = NaN;
		perturbation_melting_rate = NaN;
		geothermalflux            = NaN;
	end
	methods
		function self = linearbasalforcings(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(linearbasalforcings(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.perturbation_melting_rate=project3d(md,'vector',self.perturbation_melting_rate,'type','node','layer',1);
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1);
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','node','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end
			% TODO: Should we be setting self.geothermalflux here (we do so in linearbasalforcings.py)?

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%default values for melting parameterization
			self.deepwater_melting_rate = 50.;
			self.upperwater_melting_rate = 0.;
			self.deepwater_elevation    = -800.;
			self.upperwater_elevation   = -400.;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if numel(md.basalforcings.perturbation_melting_rate)>1
				md = checkfield(md,'fieldname','basalforcings.perturbation_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			end

			if ismember('MasstransportAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.ismasstransport==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0,'singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.upperwater_melting_rate','>=',0,'singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.deepwater_elevation','<','basalforcings.upperwater_elevation','singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.upperwater_elevation','<=',0,'singletimeseries',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0,'singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.upperwater_melting_rate','>=',0,'singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.deepwater_elevation','<','basalforcings.upperwater_elevation','singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.upperwater_elevation','<=',0,'singletimeseries',1);
			end
			if ismember('ThermalAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.isthermal==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0,'singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.upperwater_melting_rate','>=',0,'singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.deepwater_elevation','<','basalforcings.upperwater_elevation','singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.upperwater_elevation','<=',0,'singletimeseries',1);
				md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   basal forcings parameters:'));

			fielddisplay(self,'deepwater_melting_rate','basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]');
			fielddisplay(self,'deepwater_elevation','elevation of ocean deepwater [m]');
			fielddisplay(self,'upperwater_melting_rate','basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]');
			fielddisplay(self,'upperwater_elevation','elevation of ocean upperwater [m]');
			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]');
			fielddisplay(self,'perturbation_melting_rate','basal melting rate perturbation added to computed melting rate (positive if melting) [m/yr]');
			fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','perturbation_melting_rate','format','DoubleMat','name','md.basalforcings.perturbation_melting_rate','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','name','md.basalforcings.groundedice_melting_rate','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','name','md.basalforcings.geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','deepwater_melting_rate','format','DoubleMat','mattype',3,'timeserieslength',2,'name','md.basalforcings.deepwater_melting_rate','scale',1./yts,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','deepwater_elevation','format','DoubleMat','mattype',3,'name','md.basalforcings.deepwater_elevation','yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','upperwater_melting_rate','format','DoubleMat','mattype',3,'timeserieslength',2,'name','md.basalforcings.upperwater_melting_rate','scale',1./yts,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','upperwater_elevation','format','DoubleMat','mattype',3,'name','md.basalforcings.upperwater_elevation','yts',md.constants.yts)
		end % }}}
	end
end
