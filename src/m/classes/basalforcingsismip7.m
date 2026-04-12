%ISMIP7 BASAL FORCINGS class definition
%
%   Usage:
%      basalforcingsismip7=basalforcingsismip7();

classdef basalforcingsismip7
	properties (SetAccess=public) 
		num_basins                = 0;
		basin_id                  = 0;
		gamma                     = 0;
		coriolis_f                = NaN;

		salinity                  = NaN;
		tf                        = NaN;
		tf_depths                 = NaN;
		
		geothermalflux            = NaN;
		groundedice_melting_rate  = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			%self.tf=project3d(md,'vector',self.tf,'type','element','layer',1);
			%self.delta_t=project3d(md,'vector',self.delta_t,'type','element','layer',1);
			self.tf=project3d(md,'vector',self.tf,'type','node');
			self.salinity=project3d(md,'vector',self.salinity,'type','node');
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','element','layer',1); %bedrock only gets geothermal flux
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1);
		end % }}}
		function self = basalforcingsismip7(varargin) % {{{
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
			%Update fixed-coriolis parameter
			if isnan(md.mesh.lat) | isempty(md.mesh.lat),
				disp('      no md.mesh.lat specified.');
				if md.mesh.epsg == 3031 % For Antarctica
					[lat, lon] = xy2ll(md.mesh.x,md.mesh.y,-1);
				elseif md.mesh.epsg == 3413 % For Greenland
					[lat, lon] = xy2ll(md.mesh.x,md.mesh.y,1);
				else
					error('      md.mesh.lat not specified and cannot be calculated from md.mesh.epsg');
				end
			else
				lat = md.mesh.lat;
			end
			omega=7.2921e-5; %angular velocity of the Earth (rad/s)
			self.coriolis_f=2*omega*sin(lat/180*pi);

			if self.gamma == 0,
				self.gamma = 14477;
				disp('      no basalforcings.gamma specified: value set to 14477 m/yr');
			end
			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.gamma = 0.0; % ?
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','basalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','basalforcings.basin_id','Inf',1,'>=',0,'<=',md.basalforcings.num_basins,'size',[md.mesh.numberofelements 1]);

			md = checkfield(md,'fieldname','basalforcings.gamma','numel',1,'NaN',1,'Inf',1,'>',0);

			md = checkfield(md,'fieldname','basalforcings.coriolis_f','size',[md.mesh.numberofvertices, 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','basalforcings.tf_depths','NaN',1,'Inf',1,'size',[1,NaN],'<=',0);
			md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'>=',0,'timeseries',1);
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);

			md = checkfield(md,'fieldname','basalforcings.tf','size',[1,1,numel(md.basalforcings.tf_depths)]);
			md = checkfield(md,'fieldname','basalforcings.salinity','size',[1,1,numel(md.basalforcings.tf_depths)]);
			for i=1:numel(md.basalforcings.tf_depths)
				md = checkfield(md,'fieldname',['basalforcings.tf{' num2str(i) '}'],'field',md.basalforcings.tf{i},'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'>=',0,'timeseries',1);
				md = checkfield(md,'fieldname',['basalforcings.salinity{' num2str(i) '}'],'field',md.basalforcings.salinity{i},'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'>=',0,'timeseries',1);
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   ISMIP7 basal melt rate parameterization:'));
			fielddisplay(self,'num_basins','[TODO] number of basins the model domain is partitioned into [unitless]');
			fielddisplay(self,'basin_id','[TODO] basin number assigned to each node (unitless)');
			fielddisplay(self,'gamma','melt rate coefficient (m/yr)');
			fielddisplay(self,'tf_depths','elevation of vertical layers in ocean thermal forcing dataset');
			fielddisplay(self,'tf','thermal forcing (ocean temperature minus freezing point) (degrees C)');
			fielddisplay(self,'salinity','salinity (psu)');
			fielddisplay(self,'coriolis_f','Coriolis parameter (s^-1)');
			fielddisplay(self,'geothermalflux','geothermal heat flux (W/m^2)');
			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) (m/yr)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',10,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.basalforcings.basin_id','format','IntMat','mattype',2);   %0-indexed
			WriteData(fid,prefix,'object',self,'fieldname','gamma','format','Double','scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','coriolis_f','format','DoubleMat','name','md.basalforcings.coriolis_f','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','tf_depths','format','DoubleMat','name','md.basalforcings.tf_depths');
			WriteData(fid,prefix,'object',self,'fieldname','tf','format','MatArray','name','md.basalforcings.tf','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','salinity','format','MatArray','name','md.basalforcings.salinity','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','name','md.basalforcings.geothermalflux','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			
		end % }}}
	end
end
