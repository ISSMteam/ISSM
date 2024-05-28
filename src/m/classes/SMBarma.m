%SMBarma Class definition
%
%   Usage:
%      SMBarma=SMBarma();

classdef SMBarma
	properties (SetAccess=public)
		num_basins        = 0;
		num_breaks        = 0;
		num_params        = 0;
		arma_timestep     = 0;
		ar_order			 = 0;
		arlag_coefs       = NaN;
		ma_order			 = 0;
		malag_coefs       = NaN;
		polynomialparams  = NaN;
		datebreaks        = NaN;
		basin_id			 = NaN;
		lapserates        = NaN;
		elevationbins     = NaN;
		refelevation      = NaN;
		steps_per_step    = 1;
		averaging			= 0;
		requested_outputs = {};
	end
	methods
		function self = SMBarma(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			%Nothing for now
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{
			if (self.ar_order==0)
				self.ar_order = 1; %dummy 1 value for autoregression
				self.arlag_coefs      = zeros(self.num_basins,self.ar_order); %autoregression coefficients all set to 0 
				disp('      smb.ar_order (order of autoregressive model) not specified: order of autoregressive model set to 0');
			end
			if (self.ma_order==0)
				self.ma_order = 1; %dummy 1 value for moving-average
				self.malag_coefs      = zeros(self.num_basins,self.ma_order); %moving-average coefficients all set to 0 
				disp('      smb.ma_order (order of moving-average model) not specified: order of moving-average model set to 0');
			end
			if (self.arma_timestep==0)
				self.arma_timestep = md.timestepping.time_step; %ARMA model has no prescribed time step
				disp('      smb.arma_timestep (timestep of ARMA model) not specified: set to md.timestepping.time_step');
			end
			if isnan(self.arlag_coefs)
				self.arlag_coefs = zeros(self.num_basins,self.ar_order); %autoregression model of order 0 
				disp('      smb.arlag_coefs (AR lag coefficients) not specified: order of autoregressive model set to 0');
			end
			if isnan(self.malag_coefs)
				self.malag_coefs = zeros(self.num_basins,self.ma_order); %autoregression model of order 0 
				disp('      smb.malag_coefs (MA lag coefficients) not specified: order of moving-average model set to 0');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.ar_order    = 0.0; %autoregression model of order 0
			self.ma_order    = 0.0; %moving-average model of order 0
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				nbas = md.smb.num_basins;
				nprm = md.smb.num_params;
				nbrk = md.smb.num_breaks;
				md = checkfield(md,'fieldname','smb.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
				md = checkfield(md,'fieldname','smb.num_params','numel',1,'NaN',1,'Inf',1,'>',0);
				md = checkfield(md,'fieldname','smb.num_breaks','numel',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.basin_id','Inf',1,'>=',0,'<=',nbas,'size',[md.mesh.numberofelements,1]);
				if(nbas>1 && nbrk>=1 && nprm>1)
					md = checkfield(md,'fieldname','smb.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1,nprm],'numel',nbas*(nbrk+1)*nprm); 
				elseif(nbas==1)
					md = checkfield(md,'fieldname','smb.polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
				elseif(nbrk==0)
					md = checkfield(md,'fieldname','smb.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nprm],'numel',nbas*(nbrk+1)*nprm);
				elseif(nprm==1)
					md = checkfield(md,'fieldname','smb.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1],'numel',nbas*(nbrk+1)*nprm);
				end
				md = checkfield(md,'fieldname','smb.ar_order','numel',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.ma_order','numel',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.arma_timestep','numel',1,'NaN',1,'Inf',1,'>=',md.timestepping.time_step); %arma time step cannot be finer than ISSM timestep
				md = checkfield(md,'fieldname','smb.arlag_coefs','NaN',1,'Inf',1,'size',[nbas,md.smb.ar_order]);
				md = checkfield(md,'fieldname','smb.malag_coefs','NaN',1,'Inf',1,'size',[nbas,md.smb.ma_order]);
				
				if(nbrk>0)
					md = checkfield(md,'fieldname','smb.datebreaks','NaN',1,'Inf',1,'size',[nbas,nbrk]);
				elseif(numel(md.smb.datebreaks)==0 || all(isnan(md.smb.datebreaks)))
					;
				else
					error('md.smb.num_breaks is 0 but md.smb.datebreaks is not empty');
				end
				if (any(isnan(md.smb.refelevation)==0) || numel(md.smb.refelevation)>1)
					md = checkfield(md,'fieldname','smb.refelevation','NaN',1,'Inf',1,'>=',0,'size',[1,nbas],'numel',nbas);
				end
				nbas     = size(md.smb.lapserates,1);
				nbins    = size(md.smb.lapserates,2);
				ntmlapse = size(md.smb.lapserates,3);
				if(ntmlapse>1 && ntmlapse~=12)
					error('3rd dimension of md.smb.lapserates must be of size 1 or 12 (for monthly lapse rates)');
				end
				if (any(isnan(reshape(md.smb.lapserates,[1,nbas*nbins*ntmlapse]))==0) || numel(md.smb.lapserates)>1)
					md = checkfield(md,'fieldname','smb.lapserates','NaN',1,'Inf',1,'size',[nbas,nbins,ntmlapse],'numel',nbas*nbins*ntmlapse);
					md = checkfield(md,'fieldname','smb.elevationbins','NaN',1,'Inf',1,'size',[nbas,max(1,nbins-1),ntmlapse],'numel',nbas*max(1,nbins-1)*ntmlapse);
					if(issorted(md.smb.elevationbins,2)==0)
						error('md.smb.elevationbins should have rows in order of increasing elevation');
					end
				elseif (isnan(md.smb.elevationbins(1,1,1))==0 || numel(md.smb.elevationbins)>1)
					%elevationbins specified but not lapserates: this will inevitably lead to inconsistencies
					nbas     = size(md.smb.elevationbins,1);
					nbins    = size(md.smb.elevationbins,2);
					nbins    = nbins+1;
					ntmlapse = size(md.smb.elevationbins,3);
					md = checkfield(md,'fieldname','smb.lapserates','NaN',1,'Inf',1,'size',[nbas,max(1,nbins-1),ntmlapse],'numel',nbas*nbins*ntmlapse);
					md = checkfield(md,'fieldname','smb.elevationbins','NaN',1,'Inf',1,'size',[nbas,max(1,nbins-1),ntmlapse],'numel',nbas*max(1,nbins-1)*ntmlapse);
				end
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));
			fielddisplay(self,'num_basins','number of different basins [unitless]');
			fielddisplay(self,'basin_id','basin number assigned to each element [unitless]');
			fielddisplay(self,'num_breaks','number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)');
			fielddisplay(self,'num_params','number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)');
			fielddisplay(self,'polynomialparams','coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders');
			disp(sprintf('%51s  ex: polyparams=cat(3,intercepts,trendlinearcoefs,trendquadraticcoefs)',' '));
			fielddisplay(self,'datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]');
			fielddisplay(self,'ar_order','order of the autoregressive model [unitless]');
			fielddisplay(self,'ma_order','order of the moving-average model [unitless]');
			fielddisplay(self,'arma_timestep','time resolution of the ARMA model [yr]');
			fielddisplay(self,'arlag_coefs','basin-specific vectors of AR lag coefficients [unitless]');
			fielddisplay(self,'malag_coefs','basin-specific vectors of MA lag coefficients [unitless]');
			fielddisplay(self,'lapserates','basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin, 1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying values [m ice eq yr^-1 m^-1] (default: no lapse rate)');
			fielddisplay(self,'elevationbins','basin-specific separations between elevation bins, 1 row per basin, 1 column per limit between bins, dimension 3 can be of size 12 to prescribe monthly varying values [m] (default: no basin separation)');
			fielddisplay(self,'refelevation','basin-specific reference elevations at which SMB is calculated, and from which SMB is downscaled using lapserates (default: basin mean elevation) [m]');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;
			nbas = md.smb.num_basins;
			nprm = md.smb.num_params;
			nper = md.smb.num_breaks+1;

			templapserates    = md.smb.lapserates;
			tempelevationbins = md.smb.elevationbins;
			temprefelevation  = md.smb.refelevation;
			nbas     = size(md.smb.lapserates,1);
			nbins    = size(md.smb.lapserates,2);
			ntmlapse = size(md.smb.lapserates,3);
			if(any(isnan(reshape(md.smb.lapserates,[1,nbas*nbins*ntmlapse]))))
				templapserates = zeros(md.smb.num_basins,2,12);
				disp('      smb.lapserates not specified: set to 0');
			   tempelevationbins = zeros(md.smb.num_basins,1,12); %dummy elevation bins
			elseif(ntmlapse==1)
				templapserates    = repmat(templapserates,1,1,12); %same values each month
				tempelevationbins = repmat(tempelevationbins,1,1,12); %same values each month
			end
			nbas     = size(templapserates,1);
			nbins    = size(templapserates,2);
			ntmlapse = size(templapserates,3);
			if(any(isnan(md.smb.refelevation)))
				temprefelevation = zeros(1,md.smb.num_basins);
				areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);
				for ii=1:md.smb.num_basins
					indices = find(md.smb.basin_id==ii);
					elemsh  = zeros(numel(indices),1);
					for jj=1:numel(indices)
						elemsh(jj) = mean(md.geometry.surface(md.mesh.elements(indices(jj),:)));
					end
					temprefelevation(ii) = sum(areas(indices).*elemsh)/sum(areas(indices));
				end
				if(any(reshape(templapserates,[1,nbas*nbins*12])~=0))
					disp('      smb.refelevation not specified: Reference elevations set to mean surface elevation of basins');
				end
			end
			temp2dlapserates    = zeros(nbas,nbins*12);
			temp2delevationbins = zeros(nbas,max(1,nbins-1)*12);
			for(ii=[1:12])
				jj = 1+(ii-1)*nbins;
				temp2dlapserates(:,jj:jj+nbins-1)    = templapserates(:,:,ii);
				kk = 1+(ii-1)*(nbins-1);
				temp2delevationbins(:,kk:kk+nbins-2) = tempelevationbins(:,:,ii);
			end

			% Scale the parameters %
			polyparamsScaled   = md.smb.polynomialparams;
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
				dbreaks = md.smb.datebreaks;
			end

			WriteData(fid,prefix,'name','md.smb.model','data',13,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','num_params','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','num_breaks','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ar_order','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ma_order','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','arma_timestep','format','Double','scale',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','basin_id','data',self.basin_id-1,'name','md.smb.basin_id','format','IntMat','mattype',2); %0-indexed
			WriteData(fid,prefix,'data',polyparams2dScaled,'name','md.smb.polynomialparams','format','DoubleMat');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','arlag_coefs','format','DoubleMat','name','md.smb.arlag_coefs','yts',yts); 
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','malag_coefs','format','DoubleMat','name','md.smb.malag_coefs','yts',yts);
			WriteData(fid,prefix,'data',dbreaks,'name','md.smb.datebreaks','format','DoubleMat','scale',yts);
			WriteData(fid,prefix,'data',temp2dlapserates,'format','DoubleMat','name','md.smb.lapserates','scale',1./yts,'yts',yts);
			WriteData(fid,prefix,'data',temp2delevationbins,'format','DoubleMat','name','md.smb.elevationbins');
			WriteData(fid,prefix,'data',temprefelevation,'format','DoubleMat','name','md.smb.refelevation');
			WriteData(fid,prefix,'data',nbins,'format','Integer','name','md.smb.num_bins');
			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];											%remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)];	%add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
	end
end
