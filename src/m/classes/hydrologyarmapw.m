%HYDROLOGYARMAPW class definition
%
%   Usage:
%      hydrologyarmapw=hydrologyarmapw();

classdef hydrologyarmapw
	properties (SetAccess=public) 
		num_basins               = 0;
      num_params               = 0;
      num_breaks               = 0;
		basin_id                 = NaN;
      monthlyfactors           = NaN;
		polynomialparams         = NaN;
		ar_order                 = 0;
      ma_order                 = 0;
      arma_timestep            = 0;
		arlag_coefs              = NaN;
      malag_coefs              = NaN;
		datebreaks               = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.basin_id    = project3d(md,'vector',self.basin_id,'type','element');
		end % }}}
		function self = hydrologyarmapw(varargin) % {{{
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
			list = {'FrictionWaterPressure'};
		end % }}}    
		function self = setdefaultparameters(self) % {{{
			%No default parameters
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologyArmapwAnalysis',analyses)
				return;
			end

			nbas  = md.hydrology.num_basins;
         nprm  = md.hydrology.num_params;
         nbrk  = md.hydrology.num_breaks;

			md = checkfield(md,'fieldname','hydrology.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','hydrology.num_breaks','numel',1,'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','hydrology.num_params','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','hydrology.basin_id','Inf',1,'>=',0,'<=',md.hydrology.num_basins,'size',[md.mesh.numberofelements,1]);
			
			% Check if monthly factors are provided %
			if(numel(md.hydrology.monthlyfactors)>1 || ~isnan(md.hydrology.monthlyfactors))
				md = checkfield(md,'fieldname','hydrology.monthlyfactors','NaN',1,'Inf',1,'size',[md.hydrology.num_basins,12]);
				isseasonality = false;
				for(rr=[1:md.hydrology.num_basins])
					for(cc=[1:12])
						if(md.hydrology.monthlyfactors(rr,cc)~=1)
							isseasonality = true;
						end
					end
				end
				if(isseasonality && md.timestepping.time_step>=1)
					error('md.timestepping.time_step is too large to use hydrologyarmapw() with monthlyfactors');
				end
			end

			if(nbas>1 && nbrk>=1 && nprm>1)
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1,nprm],'numel',nbas*(nbrk+1)*nprm);
         elseif(nbas==1)
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
         elseif(nbrk==0)
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nprm],'numel',nbas*(nbrk+1)*nprm);
         elseif(nprm==1)
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
         end

			md = checkfield(md,'fieldname','hydrology.ar_order','numel',1,'NaN',1,'Inf',1,'>=',0);
         md = checkfield(md,'fieldname','hydrology.ma_order','numel',1,'NaN',1,'Inf',1,'>=',0);
         md = checkfield(md,'fieldname','hydrology.arma_timestep','numel',1,'NaN',1,'Inf',1,'>=',md.timestepping.time_step); %ARMA time step cannot be finer than ISSM timestep
         md = checkfield(md,'fieldname','hydrology.arlag_coefs','NaN',1,'Inf',1,'size',[md.hydrology.num_basins,md.hydrology.ar_order]);
         md = checkfield(md,'fieldname','hydrology.malag_coefs','NaN',1,'Inf',1,'size',[md.hydrology.num_basins,md.hydrology.ma_order]);

			if(nbrk>0)
            md = checkfield(md,'fieldname','hydrology.datebreaks','NaN',1,'Inf',1,'size',[nbas,nbrk]);
         elseif(numel(md.hydrology.datebreaks)==0 || all(isnan(md.hydrology.datebreaks)))
            ;
         else
            error('md.hydrology.num_breaks is 0 but md.hydrology.datebreaks is not empty');
         end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   hydrologyarmapw'));
			disp(sprintf('   subglacial water pressure is calculated as Pw=monthlyfactor[month]*(rho_water*g*bed+Pw_arma) where Pw_arma is the perturbation calculated as an ARMA process'));
			disp(sprintf('   polynomialparams includes the constant, linear trend, quadratic trend, etc. of the ARMA process'));
			disp(sprintf('   arlag_coefs and malag_coefs include the coefficients of the ARMA process'));
			fielddisplay(self,'num_basins','number of different basins');
			fielddisplay(self,'basin_id','basin number assigned to each element');
			fielddisplay(self,'num_breaks','number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)');
         fielddisplay(self,'num_params','number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)');
			fielddisplay(self,'monthlyfactors','monthly multiplicative factor on the subglacial water pressure, specified per basin (size:[num_basins,12])');
			fielddisplay(self,'polynomialparams','coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders');
         disp(sprintf('%51s  ex: polyparams=cat(3,intercepts,trendlinearcoefs,trendquadraticcoefs)',' '));
         fielddisplay(self,'datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]');
         fielddisplay(self,'ar_order','order of the autoregressive model [unitless]');
         fielddisplay(self,'ma_order','order of the moving-average model [unitless]');
         fielddisplay(self,'arma_timestep','time resolution of the autoregressive model [yr]');
         fielddisplay(self,'arlag_coefs','basin-specific vectors of AR lag coefficients [unitless]');
         fielddisplay(self,'malag_coefs','basin-specific vectors of MA lag coefficients [unitless]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;
         %%% Deal with polynomial %%%
         nbas  = md.hydrology.num_basins;
         nprm  = md.hydrology.num_params;
         nper  = md.hydrology.num_breaks+1;
         % Scale the parameters %
         polyparamsScaled   = md.hydrology.polynomialparams;
         polyparams2dScaled = zeros(nbas,nper*nprm);
         if(nprm>1)
            % Case 3D %
            if(nbas>1 && nper>1)
               for(ii=[1:nprm])
                  polyparamsScaled(:,:,ii) = polyparamsScaled(:,:,ii)*((1/yts)^(ii-1));
               end
               % Fit in 2D array %
               for(ii=[1:nprm])
                  jj = 1+(ii-1)*nper;
                  polyparams2dScaled(:,jj:jj+nper-1) = polyparamsScaled(:,:,ii);
               end
            % Case 2D and higher-order params at increasing row index %
            elseif(nbas==1)
               for(ii=[1:nprm])
                  polyparamsScaled(ii,:) = polyparamsScaled(ii,:)*((1/yts)^(ii-1));
               end
               % Fit in row array %
               for(ii=[1:nprm])
                  jj = 1+(ii-1)*nper;
                  polyparams2dScaled(1,jj:jj+nper-1) = polyparamsScaled(ii,:);
               end
            % Case 2D and higher-order params at incrasing column index %
            elseif(nper==1)
               for(ii=[1:nprm])
                  polyparamsScaled(:,ii) = polyparamsScaled(:,ii)*((1/yts)^(ii-1));
               end
               % 2D array is already in correct format %
               polyparams2dScaled = polyparamsScaled;
            end
         else
            % 2D array is already in correct format and no need for scaling %
            polyparams2dScaled = polyparamsScaled;
         end
         if(nper==1) %a single period (no break date)
            dbreaks = zeros(nbas,1); %dummy
         else
            dbreaks = md.hydrology.datebreaks;
         end

			% If no monthlyfactors provided: set them all to 1 %
			if(numel(md.hydrology.monthlyfactors)==1)
				tempmonthlyfactors = ones(nbas,12);
			else
				tempmonthlyfactors = md.hydrology.monthlyfactors;
			end

			WriteData(fid,prefix,'name','md.hydrology.model','data',7,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','num_breaks','format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','num_params','format','Integer');
         WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','basin_id','data',self.basin_id-1,'name','md.hydrology.basin_id','format','IntMat','mattype',2); %0-indexed
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ar_order','format','Integer');
         WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','ma_order','format','Integer');
         WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','arma_timestep','format','Double','scale',yts);
         WriteData(fid,prefix,'data',polyparams2dScaled,'name','md.hydrology.polynomialparams','format','DoubleMat');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','arlag_coefs','format','DoubleMat','name','md.hydrology.arlag_coefs','yts',yts);
         WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','malag_coefs','format','DoubleMat','name','md.hydrology.malag_coefs','yts',yts);
			WriteData(fid,prefix,'data',dbreaks,'name','md.hydrology.datebreaks','format','DoubleMat','scale',yts);
			WriteData(fid,prefix,'data',tempmonthlyfactors,'name','md.hydrology.monthlyfactors','format','DoubleMat');
			WriteData(fid,prefix,'data',{'FrictionWaterPressure'},'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}
	end
end

