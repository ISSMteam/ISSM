%FRONTAL FORCINGS rignot arma class definition
%
%   Usage:
%      frontalforcingsrignotarma=frontalforcingsrignotarma();

classdef frontalforcingsrignotarma
	properties (SetAccess=public) 
		num_basins               = 0;
		num_params               = 0;
		num_breaks               = 0;
		polynomialparams         = NaN;
		datebreaks               = NaN;
		ar_order                 = 0;
		ma_order                 = 0;
		arma_timestep            = 0;
		arlag_coefs              = NaN;
		malag_coefs              = NaN;
		monthlyvals_intercepts   = NaN;
		monthlyvals_trends       = NaN;
		monthlyvals_numbreaks    = 0;
		monthlyvals_datebreaks   = NaN;
		basin_id                 = NaN;
		subglacial_discharge     = NaN;
		isdischargearma          = 0;
		sd_ar_order              = 0;
		sd_ma_order              = 0;
		sd_arma_timestep         = 0;
		sd_arlag_coefs           = NaN;
		sd_malag_coefs           = NaN;
		sd_monthlyfrac           = NaN;
		sd_num_breaks            = 0;
		sd_num_params            = 0;
		sd_polynomialparams      = NaN;
		sd_datebreaks            = NaN;
	end
	methods
		function self = frontalforcingsrignotarma(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('frontalforcingsrignotarma');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
		    % nothing for now
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.basin_id             = NaN;
			self.num_basins           = 0;
			self.subglacial_discharge = NaN;
			self.ar_order             = 0.0; %autoregression model of order 0
			self.ma_order             = 0.0; %moving-average model of order 0

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
         %Early return
         if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

         nbas  = md.frontalforcings.num_basins;
         nprm  = md.frontalforcings.num_params;
         nbrk  = md.frontalforcings.num_breaks;
         nMbrk = md.frontalforcings.monthlyvals_numbreaks;
         md = checkfield(md,'fieldname','frontalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
         md = checkfield(md,'fieldname','frontalforcings.num_params','numel',1,'NaN',1,'Inf',1,'>',0);
         md = checkfield(md,'fieldname','frontalforcings.num_breaks','numel',1,'NaN',1,'Inf',1,'>=',0);
         md = checkfield(md,'fieldname','frontalforcings.basin_id','Inf',1,'>=',0,'<=',md.frontalforcings.num_basins,'size',[md.mesh.numberofelements 1]);
         if(nbas>1 && nbrk>=1 && nprm>1)
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1,nprm],'numel',nbas*(nbrk+1)*nprm);
         elseif(nbas==1)
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
         elseif(nbrk==0)
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nprm],'numel',nbas*(nbrk+1)*nprm);
         elseif(nprm==1)
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
         end
         md = checkfield(md,'fieldname','frontalforcings.ar_order','numel',1,'NaN',1,'Inf',1,'>=',0);
         md = checkfield(md,'fieldname','frontalforcings.ma_order','numel',1,'NaN',1,'Inf',1,'>=',0);
         md = checkfield(md,'fieldname','frontalforcings.arma_timestep','numel',1,'NaN',1,'Inf',1,'>=',md.timestepping.time_step); %ARMA time step cannot be finer than ISSM timestep
         md = checkfield(md,'fieldname','frontalforcings.arlag_coefs','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,md.frontalforcings.ar_order]);
         md = checkfield(md,'fieldname','frontalforcings.malag_coefs','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,md.frontalforcings.ma_order]);
         if(nbrk>0)
            md = checkfield(md,'fieldname','frontalforcings.datebreaks','NaN',1,'Inf',1,'size',[nbas,nbrk]);
         elseif(numel(md.frontalforcings.datebreaks)==0 || all(isnan(md.frontalforcings.datebreaks)))
            ;
         else
            error('md.frontalforcings.num_breaks is 0 but md.frontalforcings.datebreaks is not empty');
         end

			%%% Check if some monthly forcings are provided %%%
			if((numel(md.frontalforcings.monthlyvals_intercepts)>1 || ~isnan(md.frontalforcings.monthlyvals_intercepts)) || (numel(md.frontalforcings.monthlyvals_trends)>1 || ~isnan(md.frontalforcings.monthlyvals_trends)) || (numel(md.frontalforcings.monthlyvals_datebreaks)>1 || ~isnan(md.frontalforcings.monthlyvals_datebreaks)))
				isMonthly = true;
			else
				isMonthly = false;
			end
			if(numel(md.frontalforcings.monthlyvals_trends)>1 || ~isnan(md.frontalforcings.monthlyvals_trends))
				isMonthlyTrend = true;
			else
				isMonthlyTrend = false;
			end
			if(isMonthly)
            md = checkfield(md,'fieldname','frontalforcings.monthlyvals_numbreaks','numel',1,'NaN',1,'Inf',1,'>=',0);
			   if(nbas>1 && nMbrk>=1)
               md = checkfield(md,'fieldname','frontalforcings.monthlyvals_intercepts','NaN',1,'Inf',1,'size',[nbas,12,nMbrk+1],'numel',nbas*(nMbrk+1)*12);
			   	if(isMonthlyTrend)
			   		md = checkfield(md,'fieldname','frontalforcings.monthlyvals_trends','NaN',1,'Inf',1,'size',[nbas,12,nMbrk+1],'numel',nbas*(nMbrk+1)*12);
			   	end
			   elseif(nbas==1)
               md = checkfield(md,'fieldname','frontalforcings.monthlyvals_intercepts','NaN',1,'Inf',1,'size',[nMbrk+1,12],'numel',nbas*(nMbrk+1)*12);
			   	if(isMonthlyTrend)
			   		md = checkfield(md,'fieldname','frontalforcings.monthlyvals_trends','NaN',1,'Inf',1,'size',[nMbrk+1,12],'numel',nbas*(nMbrk+1)*12);
               end
			   elseif(nMbrk==0)
               md = checkfield(md,'fieldname','frontalforcings.monthlyvals_intercepts','NaN',1,'Inf',1,'size',[nbas,12],'numel',nbas*(nMbrk+1)*12);
			   	if(isMonthlyTrend)
			   		md = checkfield(md,'fieldname','frontalforcings.monthlyvals_trends','NaN',1,'Inf',1,'size',[nbas,12],'numel',nbas*(nMbrk+1)*12);
               end
			   end
			   if(nMbrk>0)
               md = checkfield(md,'fieldname','frontalforcings.monthlyvals_datebreaks','NaN',1,'Inf',1,'size',[nbas,nMbrk]);
            elseif(numel(md.frontalforcings.monthlyvals_datebreaks)==0 || all(isnan(md.frontalforcings.monthlyvals_datebreaks)))
               ;
            else
               error('md.frontalforcings.monthlyvals_numbreaks is 0 but md.frontalforcings.monthlyvals_datebreaks is not empty');
            end
		   end

			%%% Checking subglacial discharge %%%
			md = checkfield(md,'fieldname','frontalforcings.isdischargearma','values',[0 1]);
			if(~self.isdischargearma)
				md = checkfield(md,'fieldname','frontalforcings.subglacial_discharge','>=',0,'NaN',1,'Inf',1,'timeseries',1);
			else
				sdnbrk  = md.frontalforcings.sd_num_breaks; 
				sdnprm  = md.frontalforcings.sd_num_params;
				md = checkfield(md,'fieldname','frontalforcings.sd_ar_order','numel',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','frontalforcings.sd_ma_order','numel',1,'NaN',1,'Inf',1,'>=',0);
         	md = checkfield(md,'fieldname','frontalforcings.sd_arma_timestep','numel',1,'NaN',1,'Inf',1,'>=',max(1,md.timestepping.time_step)); %ARMA time step cannot be finer than ISSM timestep and annual timestep
         	md = checkfield(md,'fieldname','frontalforcings.sd_arlag_coefs','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,md.frontalforcings.sd_ar_order]);
         	md = checkfield(md,'fieldname','frontalforcings.sd_malag_coefs','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,md.frontalforcings.sd_ma_order]);
         	md = checkfield(md,'fieldname','frontalforcings.sd_monthlyfrac','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,12]);
				if(any(abs(sum(self.sd_monthlyfrac,2)-1)>1e-3))
					error('the 12 entries for each basin of md.frontalforcings.sd_monthlyfrac should add up to 1');
				end
	         md = checkfield(md,'fieldname','frontalforcings.sd_num_params','numel',1,'NaN',1,'Inf',1,'>',0);
		      md = checkfield(md,'fieldname','frontalforcings.sd_num_breaks','numel',1,'NaN',1,'Inf',1,'>=',0);
	         if(nbas>1 && sdnbrk>=1 && sdnprm>1)
	            md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[nbas,sdnbrk+1,sdnprm],'numel',nbas*(sdnbrk+1)*sdnprm);
	         elseif(nbas==1)
	            md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[sdnprm,sdnbrk+1],'numel',nbas*(sdnbrk+1)*sdnprm);
	         elseif(sdnbrk==0)
	            md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[nbas,sdnprm],'numel',nbas*(sdnbrk+1)*sdnprm);
	         elseif(sdnprm==1)
	            md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[nbas,sdnbrk+1],'numel',nbas*(sdnbrk+1)*sdnprm);
	         end
	         if(sdnbrk>0)
	            md = checkfield(md,'fieldname','frontalforcings.sd_datebreaks','NaN',1,'Inf',1,'size',[nbas,sdnbrk]);
	         elseif(numel(md.frontalforcings.sd_datebreaks)==0 || all(isnan(md.frontalforcings.sd_datebreaks)))
	            ;
	         else
	            error('md.frontalforcings.sd_num_breaks is 0 but md.frontalforcings.sd_datebreaks is not empty');
	         end

  			end

      end % }}}
		function disp(self) % {{{
			disp(sprintf('   Frontalforcings parameters:'));
			fielddisplay(self,'num_basins','number of different basins');
         fielddisplay(self,'basin_id','basin number assigned to each element [unitless]');
			fielddisplay(self,'num_breaks','number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)');
			fielddisplay(self,'num_params','number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)');
         fielddisplay(self,'polynomialparams','coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders');
			disp(sprintf('%51s  ex: polyparams=cat(3,intercepts,trendlinearcoefs,trendquadraticcoefs)',' '));
         fielddisplay(self,'datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]');
         fielddisplay(self,'ar_order','order of the autoregressive model [unitless]');
         fielddisplay(self,'ma_order','order of the moving-average model [unitless]');
         fielddisplay(self,'arma_timestep','time resolution of the autoregressive model [yr]');
         fielddisplay(self,'arlag_coefs','basin-specific vectors of AR lag coefficients [unitless]');
         fielddisplay(self,'malag_coefs','basin-specific vectors of MA lag coefficients [unitless]');
         fielddisplay(self,'monthlyvals_intercepts','intercept of basin-specific piecewise-linear monthly values of TF added at corresponding month (default: all 0) [°C]');
         fielddisplay(self,'monthlyvals_trends','trends in basin-specific piecewise-linear monthly values of TF added at corresponding month (default: all 0) [°C/yr]');
         fielddisplay(self,'monthlyvals_numbreaks','number of breakpoints in the piecewise-linear functions of monthly values');
         fielddisplay(self,'monthlyvals_datebreaks','dates at which the breakpoints in the piecewise-linear monthly values occur (1 row per basin)');
         fielddisplay(self,'isdischargearma','whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used)');
         fielddisplay(self,'subglacial_discharge','sum of subglacial discharge for each basin [m^3/d]');
			disp(sprintf('%51s  if isdischargearma is 1, sd_variables are used (sd arma model variable: sum of subglacial discharge for each basin [m^3/d])',' '));
         fielddisplay(self,'sd_ar_order','order of the subglacial discharge autoregressive model [unitless]');
         fielddisplay(self,'sd_ma_order','order of the subglacial discharge moving-average model [unitless]');
         fielddisplay(self,'sd_arma_timestep','time resolution of the subglacial discharge autoregressive model [yr]');
         fielddisplay(self,'sd_arlag_coefs','basin-specific vectors of AR lag coefficients for subglacial discharge [unitless]');
         fielddisplay(self,'sd_malag_coefs','basin-specific vectors of MA lag coefficients for subglacial discharge [unitless]');
         fielddisplay(self,'sd_monthlyfrac','basin-specific vectors of 12 values with fraction of the annual discharge occuring every month [unitless]');
			fielddisplay(self,'sd_num_params','number of different parameters in the subglacial discharge piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)');
			fielddisplay(self,'sd_num_breaks','number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods)');
         fielddisplay(self,'sd_datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]');
         fielddisplay(self,'sd_polynomialparams','coefficients for the sd_polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders');
 		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			%%% Deal with polynomial %%%
			nbas  = md.frontalforcings.num_basins;
			nprm  = md.frontalforcings.num_params;
			nper  = md.frontalforcings.num_breaks+1;
			% Scale the parameters %
			polyparamsScaled   = md.frontalforcings.polynomialparams;
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
				dbreaks = md.frontalforcings.datebreaks;
			end

			%%% Deal with monthly effects %%%
			nMper = md.frontalforcings.monthlyvals_numbreaks+1;
			if((numel(md.frontalforcings.monthlyvals_intercepts)<=1 && isnan(md.frontalforcings.monthlyvals_intercepts)))
				interceptsM = zeros(nbas,12); %monthly intercepts not provided, set to 0
				trendsM     = zeros(nbas,12); %set monthly trends also to 0
			else
				interceptsM3d = md.frontalforcings.monthlyvals_intercepts;
				if((numel(md.frontalforcings.monthlyvals_trends)<=1 && isnan(md.frontalforcings.monthlyvals_trends)))
					trendsM3d = 0*interceptsM3d; %monthly trends not provided, set to 0
				else
					trendsM3d = md.frontalforcings.monthlyvals_trends;
				end
			end
			% Create 2D arrays from 3D arrays if needed %
			if(nMper>1 && (numel(md.frontalforcings.monthlyvals_intercepts)>1 || ~isnan(md.frontalforcings.monthlyvals_intercepts)))
				interceptsM = zeros(nbas,12*nMper);
				trendsM     = zeros(nbas,12*nMper);
				for(ii=[1:nMper])
					jj = 1+(ii-1)*12;
					interceptsM(:,jj:jj+12-1) = interceptsM3d(:,:,ii);
					trendsM(:,jj:jj+12-1)     = trendsM3d(:,:,ii);
				end
			elseif(nMper==1 && (numel(md.frontalforcings.monthlyvals_intercepts)>1 || ~isnan(md.frontalforcings.monthlyvals_intercepts))) 
				interceptsM = interceptsM3d;
				trendsM     = trendsM3d;
			end
			if(nMper==1) %a single period (no break date)
				dMbreaks = zeros(nbas,1); %dummy
			else
				dMbreaks = md.frontalforcings.monthlyvals_datebreaks;
			end

			%%% Deal with the subglacial discharge polynomial %%%
			if(self.isdischargearma)
				sdnprm  = md.frontalforcings.sd_num_params;
				sdnper  = md.frontalforcings.sd_num_breaks+1;
				sdpolyparamsScaled   = md.frontalforcings.sd_polynomialparams;
	         sdpolyparams2dScaled = zeros(nbas,sdnper*sdnprm);
	         if(sdnprm>1)
	            % Case 3D %
	            if(nbas>1 && sdnper>1)
	               for(ii=[1:sdnprm])
	                  sdpolyparamsScaled(:,:,ii) = sdpolyparamsScaled(:,:,ii)*((1/yts)^(ii-1));
	               end
	               % Fit in 2D array %
	               for(ii=[1:sdnprm])
	                  jj = 1+(ii-1)*sdnper;
	                  sdpolyparams2dScaled(:,jj:jj+sdnper-1) = sdpolyparamsScaled(:,:,ii);
	               end
	            % Case 2D and higher-order params at increasing row index %
	            elseif(nbas==1)
	               for(ii=[1:sdnprm])
	                  sdpolyparamsScaled(ii,:) = sdpolyparamsScaled(ii,:)*((1/yts)^(ii-1));
	               end
	               % Fit in row array %
	               for(ii=[1:sdnprm])
	                  jj = 1+(ii-1)*sdnper;
	                  sdpolyparams2dScaled(1,jj:jj+sdnper-1) = sdpolyparamsScaled(ii,:);
	               end
	            % Case 2D and higher-order params at incrasing column index %
	            elseif(sdnper==1)
	               for(ii=[1:sdnprm])
	                  sdpolyparamsScaled(:,ii) = sdpolyparamsScaled(:,ii)*((1/yts)^(ii-1));
	               end
	               % 2D array is already in correct format %
						sdpolyparams2dScaled = sdpolyparamsScaled;
	            end
	         else
					% 2D array is already in correct format and no need for scaling %
               sdpolyparams2dScaled = sdpolyparamsScaled;
	         end
	         if(sdnper==1) %a single period (no break date)
	            sd_dbreaks = zeros(nbas,1); %dummy
	         else
	            sd_dbreaks = md.frontalforcings.sd_datebreaks;
	         end
			end

			WriteData(fid,prefix,'name','md.frontalforcings.parameterization','data',3,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','num_breaks','format','Integer');
			WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','num_params','format','Integer');
         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','ar_order','format','Integer');
         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','ma_order','format','Integer');
         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','arma_timestep','format','Double','scale',yts);
         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','basin_id','data',self.basin_id-1,'name','md.frontalforcings.basin_id','format','IntMat','mattype',2); %0-indexed
         WriteData(fid,prefix,'data',polyparams2dScaled,'name','md.frontalforcings.polynomialparams','format','DoubleMat');
         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','arlag_coefs','format','DoubleMat','name','md.frontalforcings.arlag_coefs','yts',yts);
         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','malag_coefs','format','DoubleMat','name','md.frontalforcings.malag_coefs','yts',yts);
         WriteData(fid,prefix,'data',dbreaks,'name','md.frontalforcings.datebreaks','format','DoubleMat','scale',yts);
			WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','monthlyvals_numbreaks','format','Integer');
         WriteData(fid,prefix,'data',dMbreaks,'name','md.frontalforcings.monthlyvals_datebreaks','format','DoubleMat','scale',yts);
         WriteData(fid,prefix,'data',interceptsM,'name','md.frontalforcings.monthlyvals_intercepts','format','DoubleMat');
         WriteData(fid,prefix,'data',trendsM,'name','md.frontalforcings.monthlyvals_trends','format','DoubleMat','scale',1/yts);
			WriteData(fid,prefix,'object',self,'fieldname','isdischargearma','format','Boolean');
			if(self.isdischargearma==0)
				WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','subglacial_discharge','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			else
				WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_num_breaks','format','Integer');
				WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_num_params','format','Integer');
	         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_ar_order','format','Integer');
	         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_ma_order','format','Integer');
	         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_arma_timestep','format','Double','scale',yts);
				WriteData(fid,prefix,'data',sdpolyparams2dScaled,'name','md.frontalforcings.sd_polynomialparams','format','DoubleMat');
	         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_arlag_coefs','format','DoubleMat','name','md.frontalforcings.sd_arlag_coefs','yts',yts);
	         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_malag_coefs','format','DoubleMat','name','md.frontalforcings.sd_malag_coefs','yts',yts);
	         WriteData(fid,prefix,'data',sd_dbreaks,'name','md.frontalforcings.sd_datebreaks','format','DoubleMat','scale',yts);
	         WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_monthlyfrac','format','DoubleMat','name','md.frontalforcings.sd_monthlyfrac','yts',yts);
			end
		end % }}}
	end
end
