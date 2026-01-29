%BECKMANNGOOSSE BASAL FORCINGS class definition
%
%   Usage:
%      basalforcingsbeckmanngoosse=basalforcingsbeckmanngoosse();

classdef basalforcingsbeckmanngoosse
	properties (SetAccess=public) 
		groundedice_melting_rate  = NaN;
		geothermalflux            = NaN;
		meltrate_factor           = 0.;
		ocean_temp                = 0.;
		ocean_salinity            = NaN;
		ocean_thermalforcing      = NaN;
		isthermalforcing          = 0;
	end
	methods
		function createxml(self,fid) % {{{
			fprintf(fid, '\n\n');
			fprintf(fid, '%s\n', '<!-- basalforcings -->');
		        fprintf(fid,'%s%s%s%s%s\n%s\n%s\n','<parameter key ="geothermalflux" type="',class(self.geothermalflux),'" default="',num2str(self.geothermalflux),'">', '     <section name="basalforcings" />','     <help> geothermal heat flux [W/m^2] </help>','</parameter>');
			fprintf(fid,'%s%s%s%s%s\n%s\n%s\n%s\n','<parameter key ="melting_rate" type="',class(self.melting_rate),'" default="',num2str(self.melting_rate),'">','     <section name="basalforcings" />','     <help> basal melting rate (positive if melting) [m/yr] </help>','</parameter>');
			fprintf(fid,'%s%s%s%s%s\n%s\n%s\n%s\n','<parameter key ="ocean_temp" type="',class(self.ocean_temp),'" default="',num2str(self.ocean_temp),'">','     <section name="basalforcings" />','     <help> ocean_temp [degC] </help>','</parameter>');
		        fprintf(fid,'%s%s%s%s%s\n%s\n%s\n%s\n','<parameter key ="ocean_salinity" type="',class(self.ocean_salinity),'" default="',num2str(self.ocean_salinity),'">','     <section name="basalforcings" />','     <help> ocean_salinity [psu] </help>','</parameter>');
		        fprintf(fid,'%s%s%s%s%s\n%s\n%s\n%s\n','<parameter key ="ocean_thermalforcing" type="',class(self.ocean_thermalforcing),'" default="',num2str(self.ocean_thermalforcing),'">','     <section name="basalforcings" />','     <help> ocean_thermalforcing [K] </help>','</parameter>');
        	end % }}}
		function self = extrude(self,md) % {{{
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1); 
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','node','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = basalforcingsbeckmanngoosse(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(basalforcingsbeckmanngoosse(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end
			if isnan(self.ocean_temp),
				self.ocean_temp=-1.7*ones(md.mesh.numberofvertices,1);
				disp('      no basalforcings.ocean_temp specified: values set as -1.7degC');
			end
			if isnan(self.ocean_salinity),
				self.ocean_salinity=35.0*ones(md.mesh.numberofvertices,1);
				disp('      no basalforcings.ocean_salinity specified: values set as 35 psu');
			end


		end % }}}
		function self = setdefaultparameters(self) % {{{

			%default values for melting parameterization
			self.meltrate_factor        = 0.5;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if(self.isthermalforcing==0)
				if(any(~isnan(self.ocean_thermalforcing)))
					error('basalforcingsbeckmanngoosse has isthermalforcing=0 but has specified entries for ocean_thermalforcing (should be unused)');
				end
			elseif(self.isthermalforcing==1)
				if(any(self.ocean_temp~=0) || any(~isnan(self.ocean_salinity)))
					error('basalforcingsbeckmanngoosse has isthermalforcing=1 but has specified entries for ocean_temp or ocean_salinity (should be unused)');
				end
				if(strcmp(class(md.frontalforcings),'frontalforcingsrignot') && ~isequal(md.frontalforcings.thermalforcing,self.ocean_thermalforcing))
					error('basalforcing.ocean_thermalforcing is not consistent with frontalforcings.thermalforcing')
				end
				if(strcmp(class(md.frontalforcings),'frontalforcingsrignotarma') && any(~isnan(self.ocean_thermalforcing)))
					error('basalforcings has specified entries for ocean_thermalforcing, but frontalforcingsrignotarma already calculates thermalforcing internally')
				end
			end


			if ismember('MasstransportAnalysis',analyses) & ~(solution=='TransientSolution' & md.transient.ismasstransport==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'size','universal','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','basalforcings.isthermalforcing','values',[0 1]);
				if(self.isthermalforcing==0)
					md = checkfield(md,'fieldname','basalforcings.ocean_temp','NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','basalforcings.ocean_salinity','NaN',1,'Inf',1,'timeseries',1);	
				elseif(self.isthermalforcing==1 && ~strcmp(class(md.frontalforcings),'frontalforcingsrignotarma'))
					md = checkfield(md,'fieldname','basalforcings.ocean_thermalforcing','NaN',1,'Inf',1,'timeseries',1);	
				end
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'size','universal','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','basalforcings.isthermalforcing','values',[0 1]);
				if(self.isthermalforcing==0)
					md = checkfield(md,'fieldname','basalforcings.ocean_temp','NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','basalforcings.ocean_salinity','NaN',1,'Inf',1,'timeseries',1);	
				elseif(self.isthermalforcing==1 && ~strcmp(class(md.frontalforcings),'frontalforcingsrignotarma'))
					md = checkfield(md,'fieldname','basalforcings.ocean_thermalforcing','NaN',1,'Inf',1,'timeseries',1);	
				end
			end
			if ismember('ThermalAnalysis',analyses) & ~(solution=='TransientSolution' & md.transient.isthermal==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0);
				md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'size','universal','NaN',1,'Inf',1);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Beckmann & Goosse (2003) basal melt parameterization:'));

			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) (m/yr)');
			fielddisplay(self,'geothermalflux','geothermal heat flux (W/m^2)');
			fielddisplay(self,'meltrate_factor','Melt-rate rate factor');
			fielddisplay(self,'isthermalforcing','whether to use ocean_temp and ocean_salinity (0) or ocean_thermalforcing (1) (default:0)');
			fielddisplay(self,'ocean_temp','ocean temperature (degC)');
			fielddisplay(self,'ocean_salinity','ocean salinity (psu)');
			fielddisplay(self,'ocean_thermalforcing','ocean thermalforcing (K)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=365.2422*24.0*3600.0;

			floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
			T_f=(0.0939 - 0.057 * md.basalforcings.ocean_salinity + 7.64e-4 * md.geometry.base);
			ocean_heat_flux = 1.5 * md.materials.rho_water * md.materials.mixed_layer_capacity * md.materials.thermal_exchange_velocity  * (md.basalforcings.ocean_temp - T_f);
			floatingice_melting_rate=ocean_heat_flux/(md.materials.latentheat*md.materials.rho_ice);


			WriteData(fid,prefix,'name','md.basalforcings.model','data',8,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','name','md.basalforcings.groundedice_melting_rate','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','name','md.basalforcings.geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1);
			WriteData(fid,prefix,'object',self,'fieldname','meltrate_factor','format','DoubleMat','mattype',1,'name','md.basalforcings.meltrate_factor');
			WriteData(fid,prefix,'object',self,'fieldname','isthermalforcing','format','Boolean');
			if(self.isthermalforcing==0)
				WriteData(fid,prefix,'object',self,'fieldname','ocean_temp','format','DoubleMat','name','md.basalforcings.ocean_temp','mattype',1,'timeserieslength',md.mesh.numberofvertices+1);
				WriteData(fid,prefix,'object',self,'fieldname','ocean_salinity','format','DoubleMat','name','md.basalforcings.ocean_salinity','mattype',1,'timeserieslength',md.mesh.numberofvertices+1);
			elseif(self.isthermalforcing==1)
				WriteData(fid,prefix,'object',self,'fieldname','ocean_thermalforcing','format','DoubleMat','name','md.basalforcings.ocean_thermalforcing','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end
		end % }}}
	end
end
