%PICO BASAL FORCINGS class definition
%
%   Usage:
%      basalforcingspico=basalforcingspico();

classdef basalforcingspico
	properties (SetAccess=public) 
		num_basins                = 0;
		basin_id                  = NaN;
		maxboxcount               = 0;
		overturning_coeff         = NaN;
		gamma_T                   = 0.;
		farocean_temperature      = NaN;
		farocean_salinity         = NaN;
		isplume                   = 0;
		geothermalflux            = NaN;
		groundedice_melting_rate  = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.basin_id=project3d(md,'vector',self.basin_id,'type','element','layer',1);
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','element','layer',1); %bedrock only gets geothermal flux
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1);
		end % }}}
		function self = basalforcingspico(varargin) % {{{
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
			if isnan(self.maxboxcount),
	            self.maxboxcount = 5;
		   			disp('      no maximum number of boxes set, setting value to 5');
		   end
			if isnan(self.overturning_coeff)
				self.overturning_coeff = 1e6*ones(md.mesh.numberofvertices,1); %m^3/s
				disp('      no overturning strength set, setting value to 1e6');
			end
			if isnan(self.gamma_T)
				self.gamma_T = 2e-5; %m/s
				disp('      no turbulent temperature exchange velocity set, setting value to 2e-5');
			end
			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.maxboxcount       = 5;
			self.overturning_coeff = 1e6; %m^3/s
			self.gamma_T           = 2e-5; %m/s
			self.isplume           = false;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

				md = checkfield(md,'fieldname','basalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
				md = checkfield(md,'fieldname','basalforcings.basin_id','Inf',1,'>=',0,'<=',md.basalforcings.num_basins,'size',[md.mesh.numberofelements 1]);
				md = checkfield(md,'fieldname','basalforcings.maxboxcount','numel',1,'NaN',1,'Inf',1,'>',0);
				if numel(self.overturning_coeff)==1
					md = checkfield(md,'fieldname','basalforcings.overturning_coeff','numel',1,'NaN',1,'Inf',1,'>',0);
				else
					md = checkfield(md,'fieldname','basalforcings.overturning_coeff','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1,'>',0);
				end
				md = checkfield(md,'fieldname','basalforcings.gamma_T','numel',1,'NaN',1,'Inf',1,'>',0);
				md = checkfield(md,'fieldname','basalforcings.farocean_temperature','NaN',1,'Inf',1,'size',[md.basalforcings.num_basins+1 NaN]);
				md = checkfield(md,'fieldname','basalforcings.farocean_salinity','NaN',1,'Inf',1,'>',0,'size',[md.basalforcings.num_basins+1 NaN]);
				md = checkfield(md,'fieldname','basalforcings.isplume','values',[0 1]);
				md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'>=',0,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   PICO basal melt rate parameterization:'));
			fielddisplay(self,'num_basins','number of basins the model domain is partitioned into [unitless]');
			fielddisplay(self,'basin_id','basin number assigned to each node [unitless]');
			fielddisplay(self,'maxboxcount','maximum number of boxes initialized under all ice shelves');
			fielddisplay(self,'overturning_coeff','overturning strength [m^3/s]');
			fielddisplay(self,'gamma_T','turbulent temperature exchange velocity [m/s]');
			fielddisplay(self,'farocean_temperature','depth averaged ocean temperature in front of the ice shelf for basin i [K]');
			fielddisplay(self,'farocean_salinity','depth averaged ocean salinity in front of the ice shelf for basin i [psu]');
			fielddisplay(self,'isplume','boolean to use buoyant plume melt rate parameterization from Lazeroms et al., 2018 (default false)');
			fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]');
			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',5,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','maxboxcount','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','overturning_coeff','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','gamma_T','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','farocean_temperature','format','DoubleMat','name','md.basalforcings.farocean_temperature','timeserieslength',md.basalforcings.num_basins+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','farocean_salinity','format','DoubleMat','name','md.basalforcings.farocean_salinity','timeserieslength',md.basalforcings.num_basins+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.basalforcings.basin_id','format','IntMat','mattype',2);   %Change to 0-indexing
			WriteData(fid,prefix,'object',self,'fieldname','isplume','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','name','md.basalforcings.geothermalflux','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);

		end % }}}
	end
end
