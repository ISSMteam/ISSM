%LINEAR BASAL FORCINGS ARMA class definition
%
%   Usage:
%      linearbasalforcingsarma=linearbasalforcingsarma();

classdef linearbasalforcingsarma
	
	properties (SetAccess=public) 
		num_basins                = 0;
		num_params                = 0;
		num_breaks                = 0;
		polynomialparams          = NaN;
      datebreaks                = NaN;
      ar_order                  = 0;
      ma_order                  = 0;
      arma_timestep             = 0;
      arlag_coefs               = NaN;
      malag_coefs               = NaN;
      basin_id                  = NaN;
		groundedice_melting_rate  = NaN;
		deepwater_elevation       = NaN;
		upperwater_melting_rate   = NaN;
		upperwater_elevation      = NaN;
		geothermalflux            = NaN;
	end
	methods
		function self = linearbasalforcingsarma(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1); 
			self.deepwater_elevation=project3d(md,'vector',self.deepwater_elevation,'type','node','layer',1); 
			self.upperwater_melting_rate=project3d(md,'vector',self.upperwater_melting_rate,'type','node','layer',1); 
			self.upperwater_elevation=project3d(md,'vector',self.upperwater_elevation,'type','node','layer',1); 
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','node','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end
			if isnan(self.trend)
            self.trend = zeros(1,self.num_basins); %no trend in SMB
            disp('      basalforcings.trend (trend) not specified: value set to 0');
         end
         if (self.ar_order==0)
            self.ar_order = 1; %dummy 1 value for autoregression
            self.arlag_coefs      = zeros(self.num_basins,self.ar_order); %autoregression coefficients all set to 0
            disp('      basalforcings.ar_order (order of autoregressive model) not specified: order of autoregressive model set to 0');
         end
         if (self.ma_order==0)
            self.ma_order = 1; %dummy 1 value for moving-average
            self.arlag_coefs      = zeros(self.num_basins,self.ma_order); %moving-average coefficients all set to 0
            disp('      basalforcings.ma_order (order of moving-average model) not specified: order of moving-average model set to 0');
         end
         if (self.arma_timestep==0)
            self.arma_timestep = md.timestepping.time_step; %ARMA model has no prescribed time step
            disp('      basalforcings.arma_timestep (timestep of ARMA model) not specified: set to md.timestepping.time_step');
         end
         if isnan(self.arlag_coefs)
            self.arlag_coefs = zeros(self.num_basins,self.ar_order); %autoregression model of order 0
            disp('      basalforcings.arlag_coefs (AR lag coefficients) not specified: order of autoregressive model set to 0');
         end

		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.ar_order    = 0.0; %autoregression model of order 0
         self.ma_order    = 0.0; %moving-average model of order 0
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.ismasstransport==0),
				nbas = md.basalforcings.num_basins;
				nprm = md.basalforcings.num_params;
				nbrk = md.basalforcings.num_breaks;
				md = checkfield(md,'fieldname','basalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
				md = checkfield(md,'fieldname','basalforcings.num_params','numel',1,'NaN',1,'Inf',1,'>',0);
				md = checkfield(md,'fieldname','basalforcings.num_breaks','numel',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.deepwater_elevation','NaN',1,'Inf',1,'size',[1,md.basalforcings.num_basins],'numel',md.basalforcings.num_basins);
				md = checkfield(md,'fieldname','basalforcings.upperwater_melting_rate','NaN',1,'Inf',1,'>=',0,'size',[1,md.basalforcings.num_basins],'numel',md.basalforcings.num_basins);
				md = checkfield(md,'fieldname','basalforcings.upperwater_elevation','NaN',1,'Inf',1,'<=',0,'size',[1,md.basalforcings.num_basins],'numel',md.basalforcings.num_basins);
            md = checkfield(md,'fieldname','basalforcings.basin_id','Inf',1,'>=',0,'<=',md.basalforcings.num_basins,'size',[md.mesh.numberofelements,1]);
            if(nbas>1 && nbrk>=1 && nprm>1)
					md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1,nprm],'numel',nbas*(nbrk+1)*nprm);
				elseif(nbas==1)
					md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
				elseif(nbrk==0)
					md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nprm],'numel',nbas*(nbrk+1)*nprm);
				elseif(nprm==1)
					md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
				end
				md = checkfield(md,'fieldname','basalforcings.ar_order','numel',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','basalforcings.ma_order','numel',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md,'fieldname','basalforcings.arma_timestep','numel',1,'NaN',1,'Inf',1,'>=',md.timestepping.time_step); %moving-average time step cannot be finer than ISSM timestep
            md = checkfield(md,'fieldname','basalforcings.arlag_coefs','NaN',1,'Inf',1,'size',[md.basalforcings.num_basins,md.basalforcings.ar_order]);
            md = checkfield(md,'fieldname','basalforcings.malag_coefs','NaN',1,'Inf',1,'size',[md.basalforcings.num_basins,md.basalforcings.ma_order]);
				if(nbrk>0)
					md = checkfield(md,'fieldname','basalforcings.datebreaks','NaN',1,'Inf',1,'size',[nbas,nbrk]);
				elseif(numel(md.basalforcings.datebreaks)==0 || all(isnan(md.basalforcings.datebreaks)))
					;
				else
					error('md.basalforcings.num_breaks is 0 but md.basalforcings.datebreaks is not empty');
         end
			end
			if ismember('BalancethicknessAnalysis',analyses),
				error('not implemented yet!');
			end
			if ismember('ThermalAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.isthermal==0),
				error('not implemented yet!');
			end
			if numel(md.basalforcings.geothermalflux)>1
            md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0);
         end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   autoregression linear basal forcings parameters:'));
			disp(sprintf('   autoregressive model is applied for deepwater_melting_rate'));

			fielddisplay(self,'num_basins','number of different basins [unitless]');
         fielddisplay(self,'basin_id','basin number assigned to each element [unitless]');
         fielddisplay(self,'num_breaks','number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)');
         fielddisplay(self,'num_params','number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)');
         fielddisplay(self,'polynomialparams','coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders');
         disp(sprintf('%51s  ex: polyparams=cat(3,intercepts,trendlinearcoefs,trendquadraticcoefs)',' '));
         fielddisplay(self,'datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]');
			fielddisplay(self,'ar_order','order of the autoregressive model [unitless]');
         fielddisplay(self,'ar_order','order of the moving-average model [unitless]');
         fielddisplay(self,'arma_timestep','time resolution of the ARMA model [yr]');
         fielddisplay(self,'arlag_coefs','basin-specific vectors of AR lag coefficients [unitless]');
         fielddisplay(self,'malag_coefs','basin-specific vectors of MA lag coefficients [unitless]');
			fielddisplay(self,'deepwater_elevation','basin-specific elevation of ocean deepwater [m]');
			fielddisplay(self,'upperwater_melting_rate','basin-specific basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]');
			fielddisplay(self,'upperwater_elevation','basin-specific elevation of ocean upperwater [m]');
			fielddisplay(self,'groundedice_melting_rate','node-specific basal melting rate (positive if melting) [m/yr]');
			fielddisplay(self,'geothermalflux','node-specific geothermal heat flux [W/m^2]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;
			nbas = md.basalforcings.num_basins;
         nprm = md.basalforcings.num_params;
         nper = md.basalforcings.num_breaks+1;
         % Scale the parameters %
         polyparamsScaled   = md.basalforcings.polynomialparams;
         polyparams2dScaled = zeros(nbas,nper*nprm);
         if(nprm>1)
            % Case 3D %
            if(nbas>1 && nper>1)
               for(ii=[1:nprm])
                  polyparamsScaled(:,:,ii) = polyparamsScaled(:,:,ii)*((1/yts)^(ii));
               end
               % Fit in 2D array %
               for(ii=[1:nprm])
                  jj = 1+(ii-1)*nper;
                  polyparams2dScaled(:,jj:jj+nper-1) = polyparamsScaled(:,:,ii);
               end
            % Case 2D and higher-order params at increasing row index %
            elseif(nbas==1)
               for(ii=[1:nprm])
                  polyparamsScaled(ii,:) = polyparamsScaled(ii,:)*((1/yts)^(ii));
               end
               % Fit in row array %
               for(ii=[1:nprm])
                  jj = 1+(ii-1)*nper;
                  polyparams2dScaled(1,jj:jj+nper-1) = polyparamsScaled(ii,:);
               end
            % Case 2D and higher-order params at incrasing column index %
            elseif(nper==1)
               for(ii=[1:nprm])
                  polyparamsScaled(:,ii) = polyparamsScaled(:,ii)*((1/yts)^(ii));
               end
               % 2D array is already in correct format %
               polyparams2dScaled = polyparamsScaled;
            end
         else
				polyparamsScaled   = polyparamsScaled*(1/yts);
            % 2D array is already in correct format %
            polyparams2dScaled = polyparamsScaled;
         end
         if(nper==1) %a single period (no break date)
            dbreaks = zeros(nbas,1); %dummy
         else
            dbreaks = md.basalforcings.datebreaks;
         end

			WriteData(fid,prefix,'name','md.basalforcings.model','data',9,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','name','md.basalforcings.groundedice_melting_rate','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','name','md.basalforcings.geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','num_params','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','num_breaks','format','Integer');
         WriteData(fid,prefix,'object',self,'fieldname','ar_order','format','Integer');
         WriteData(fid,prefix,'object',self,'fieldname','ma_order','format','Integer');
         WriteData(fid,prefix,'object',self,'fieldname','arma_timestep','format','Double','scale',yts);
         WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.basalforcings.basin_id','format','IntMat','mattype',2); %0-indexed
         WriteData(fid,prefix,'object',self,'fieldname','arlag_coefs','format','DoubleMat','name','md.basalforcings.arlag_coefs','yts',yts);	
         WriteData(fid,prefix,'object',self,'fieldname','malag_coefs','format','DoubleMat','name','md.basalforcings.malag_coefs','yts',yts);	
			WriteData(fid,prefix,'data',polyparams2dScaled,'name','md.basalforcings.polynomialparams','format','DoubleMat');
			WriteData(fid,prefix,'data',dbreaks,'name','md.basalforcings.datebreaks','format','DoubleMat','scale',yts);
			WriteData(fid,prefix,'object',self,'fieldname','deepwater_elevation','format','DoubleMat','name','md.basalforcings.deepwater_elevation');
			WriteData(fid,prefix,'object',self,'fieldname','upperwater_melting_rate','format','DoubleMat','name','md.basalforcings.upperwater_melting_rate','scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','upperwater_elevation','format','DoubleMat','name','md.basalforcings.upperwater_elevation');
		end % }}}
	end
end
