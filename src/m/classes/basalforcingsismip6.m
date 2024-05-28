%ISMIP6 BASAL FORCINGS class definition
%
%   Usage:
%      basalforcingsismip6=basalforcingsismip6();

classdef basalforcingsismip6
	properties (SetAccess=public) 
		num_basins                = 0;
		basin_id                  = NaN;
		gamma_0                   = 0.;
		tf                        = NaN;
		tf_depths                 = NaN;
		delta_t                   = NaN;
		islocal                   = 0;
		geothermalflux            = NaN;
		groundedice_melting_rate  = NaN;
		melt_anomaly              = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.basin_id=project3d(md,'vector',self.basin_id,'type','element','layer',1);
			%self.tf=project3d(md,'vector',self.tf,'type','element','layer',1);
			%self.delta_t=project3d(md,'vector',self.delta_t,'type','element','layer',1);
			self.tf=project3d(md,'vector',self.tf,'type','node');
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','element','layer',1); %bedrock only gets geothermal flux
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1);
			self.melt_anomaly=project3d(md,'vector',self.melt_anomaly,'type','element','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = basalforcingsismip6(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=setdefaultparameters(self);
					self=structtoobj(self,varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = initialize(self,md) % {{{
			if self.gamma_0 == 0,
				self.gamma_0 = 14477;
				disp('      no basalforcings.gamma_0 specified: value set to 14477 m/yr');
			end
			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.gamma_0 = 14477; %m/yr
			self.islocal = false;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','basalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','basalforcings.basin_id','Inf',1,'>=',0,'<=',md.basalforcings.num_basins,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','basalforcings.gamma_0','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','basalforcings.tf_depths','NaN',1,'Inf',1,'size',[1,NaN],'<=',0);
			md = checkfield(md,'fieldname','basalforcings.delta_t','NaN',1,'Inf',1,'numel',md.basalforcings.num_basins,'size',[1,md.basalforcings.num_basins]);
			md = checkfield(md,'fieldname','basalforcings.islocal','values',[0 1]);
			md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'>=',0,'timeseries',1);
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			if length(md.basalforcings.melt_anomaly)>1,
				md = checkfield(md,'fieldname','basalforcings.melt_anomaly','NaN',1,'Inf',1,'timeseries',1);
			end

			md = checkfield(md,'fieldname','basalforcings.tf','size',[1,1,numel(md.basalforcings.tf_depths)]);
			for i=1:numel(md.basalforcings.tf_depths)
				md = checkfield(md,'fieldname',['basalforcings.tf{' num2str(i) '}'],'field',md.basalforcings.tf{i},'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'>=',0,'timeseries',1);
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   ISMIP6 basal melt rate parameterization:'));
			fielddisplay(self,'num_basins','number of basins the model domain is partitioned into [unitless]');
			fielddisplay(self,'basin_id','basin number assigned to each node (unitless)');
			fielddisplay(self,'gamma_0','melt rate coefficient (m/yr)');
			fielddisplay(self,'tf_depths','elevation of vertical layers in ocean thermal forcing dataset');
			fielddisplay(self,'tf','thermal forcing (ocean temperature minus freezing point) (degrees C)');
			fielddisplay(self,'delta_t','Ocean temperature correction per basin (degrees C)');
			fielddisplay(self,'islocal','boolean to use the local version of the ISMIP6 melt rate parameterization (default false)');
			fielddisplay(self,'geothermalflux','geothermal heat flux (W/m^2)');
			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) (m/yr)');
			fielddisplay(self,'melt_anomaly','floating ice basal melt anomaly (m/yr)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',7,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.basalforcings.basin_id','format','IntMat','mattype',2);   %0-indexed
			WriteData(fid,prefix,'object',self,'fieldname','gamma_0','format','Double','scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','tf_depths','format','DoubleMat','name','md.basalforcings.tf_depths');
			WriteData(fid,prefix,'object',self,'fieldname','tf','format','MatArray','name','md.basalforcings.tf','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','delta_t','format','DoubleMat','name','md.basalforcings.delta_t','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','islocal','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','name','md.basalforcings.geothermalflux','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','melt_anomaly','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);

		end % }}}
	end
end
