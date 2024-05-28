%STOCHASTICFORCING class definition
%
%   Usage:
%      stochasticforcing=stochasticforcing();

classdef stochasticforcing
	properties (SetAccess=public)
		isstochasticforcing	= 0;
		fields					= NaN;
		defaultdimension		= 0;
		default_id				= NaN;
		covariance				= NaN;
		timecovariance			= NaN;
		stochastictimestep   = 0;
		randomflag				= 1;
	end
	methods
		function self = stochasticforcing(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.default_id = project3d(md,'vector',self.default_id,'type','element');
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.isstochasticforcing	= 0; %stochasticforcing is turned off by default
			self.randomflag				= 1; %true randomness is implemented by default
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~self.isstochasticforcing, return; end

			num_fields = numel(self.fields);
			if(self.stochastictimestep==0)
				md.stochasticforcing.stochastictimestep = md.timestepping.time_step; %by default: stochastictimestep set to ISSM time step
			end


			if(numel(size(self.covariance)==3))
				numtcovmat = numel(self.covariance(1,1,:)); %number of covariance matrices in time
				lsCovmats = {};
				for ii=[1:numtcovmat] %loop over 3rd dimension
					lsCovmats{ii} = self.covariance(:,:,ii);
      			%Check that covariance matrix is positive definite
      			try
      				chol(self.covariance(:,:,ii));
      			catch
      				error('an entry in md.stochasticforcing.covariance is not positive definite');
					end
				end
			elseif(numel(size(self.covariance)==2))
				numtcovmat = 1; %number of covariance matrices in time
				lsCovmats = {self.covariance};
   			%Check that covariance matrix is positive definite
   			try
   				chol(self.covariance);
   			catch
   				error('md.stochasticforcing.covariance is not positive definite');
   			end
			end

			%Check that all fields agree with the corresponding md class and if any field needs the default params
			checkdefaults	= false; %need to check defaults only if one of the fields does not have its own dimensionality
			structstoch		= structstochforcing();
			% Check if hydrolgyarmapw is used
			if(strcmp(class(md.hydrology),'hydrologyarmapw') && md.transient.ishydrology==1)
				ispwHydroarma = 1;
			else
				ispwHydroarma = 0;
			end
			for field=self.fields
				%Checking agreement of classes
				if(contains(field,'SMBarma'))
					mdname = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					if~(isequal(class(md.smb),char(mdname)))
						error('md.smb does not agree with stochasticforcing field %s', char(field));
					end
				end
				if(contains(field,'SMBforcing'))
					mdname = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					if~(isequal(class(md.smb),char(mdname)))
						error('md.smb does not agree with stochasticforcing field %s', char(field));
					end
				end
				if(contains(field,'FrontalForcings'))
					mdname = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					if~(isequal(class(md.frontalforcings),char(mdname)))
						error('md.frontalforcings does not agree with stochasticforcing field %s', char(field));
					end
				end
				if(contains(field,'Calving'))
					mdname = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					if~(isequal(class(md.calving),char(mdname)))
						error('md.calving does not agree with stochasticforcing field %s', char(field));
					end
				end
				if(contains(field,'BasalforcingsFloatingice'))
					mdname = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					if~(isequal(class(md.basalforcings),char(mdname)))
						error('md.basalforcings does not agree with stochasticforcing field %s', char(field));
					end
				end
				if(contains(field,'BasalforcingsSpatialDeepwaterMeltingRate'))
					mdname = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					if~(isequal(class(md.basalforcings),char(mdname)))
						error('md.basalforcings does not agree with stochasticforcing field %s', char(field));
					end
				end
				if(contains(field,'BasalforcingsDeepwaterMeltingRatearma'))
					mdname = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					if~(isequal(class(md.basalforcings),char(mdname)))
						error('md.basalforcings does not agree with stochasticforcing field %s', char(field));
					end
				end
				if(contains(field,'WaterPressure'))

					mdnames = structstoch.mdnames(find(strcmp(structstoch.fields,char(field))));
					found   = 0;
					for(ii=[1:numel(mdnames)])
						if(isequal(class(md.friction),char(mdnames{ii}))) found=1; end
					end
					if(found==0)
						error('md.friction does not agree with stochasticforcing field %s', char(field));
					end
					if(strcmp(class(md.friction),'friction') || strcmp(class(md.friction),'frictionschoof') || strcmp(class(md.friction),'frictioncoulomb'))
						if(md.friction.coupling~=0 && md.friction.coupling~=1 && md.friction.coupling~=2)
							error('stochasticforcing field %s is only implemented for cases md.friction.coupling 0 or 1 or 2', char(field));
						end
					end
					if(strcmp(class(md.friction),'friction'))
						if(any(md.friction.q==0))
							error('stochasticforcing field %s requires non-zero q exponent',char(field));
						end
					end
				end
				%Checking for specific dimensions
				if ~(strcmp(field,'SMBarma') || strcmp(field,'FrontalForcingsRignotarma') || strcmp(field,'BasalforcingsDeepwaterMeltingRatearma') || strcmp(field,'FrontalForcingsSubglacialDischargearma') || ((strcmp(field,'FrictionWaterPressure') && ispwHydroarma)))
					checkdefaults = true; %field with non-specific dimensionality
				end
			end
			%Retrieve all the field dimensionalities
			dimensions = self.defaultdimension*ones(1,num_fields);
			indSMBarma = -1; %about to check for index of SMBarma
			indTFarma  = -1; %about to check for index of FrontalForcingsRignotarma
			indSdarma  = -1; %about to check for index of SubglacialDischargearma
			indBDWarma = -1; %about to check for index of BasalforcingsDeepwaterMeltingRatearma
			indPwarma  = -1; %about to check for index of hydrologyarmapw
			

			if any(contains(self.fields,'SMBarma'))
				indSMBarma = find(contains(self.fields,'SMBarma')); %index of SMBarma, now check for consistency with other arma timesteps 
				dimensions(indSMBarma) = md.smb.num_basins;
				if(md.smb.arma_timestep<self.stochastictimestep)
					error('SMBarma cannot have a timestep shorter than stochastictimestep');
				end
			end
			if any(contains(self.fields,'FrontalForcingsRignotarma'))
				indTFarma	= find(contains(self.fields,'FrontalForcingsRignotarma')); %index of TFarma, now check for consistency with other arma timesteps 
				dimensions(indTFarma) = md.frontalforcings.num_basins;
				if(md.frontalforcings.arma_timestep<self.stochastictimestep)
					error('FrontalForcingsRignotarma cannot have a timestep shorter than stochastictimestep');
				end
			end
			if any(contains(self.fields,'FrontalForcingsSubglacialDischargearma'))
				indSdarma	= find(contains(self.fields,'FrontalForcingsSubglacialDischargearma')); %index of Sdarma, now check for consistency with other arma timesteps 
				dimensions(indSdarma) = md.frontalforcings.num_basins;
				if(md.frontalforcings.sd_arma_timestep<self.stochastictimestep)
					error('FrontalForcingsSubglacialDischargearma cannot have a timestep shorter than stochastictimestep');
				end
			end
			if any(contains(self.fields,'BasalforcingsDeepwaterMeltingRatearma'))
				indBDWarma	= find(contains(self.fields,'BasalforcingsDeepwaterMeltingRatearma')); %index of BDWarma, now check for consistency with other arma timesteps 
				dimensions(indBDWarma) = md.basalforcings.num_basins;
				if(md.basalforcings.arma_timestep<self.stochastictimestep)
					error('BasalforcingsDeepwaterMeltingRatearma cannot have a timestep shorter than stochastictimestep');
				end
			end
			if (any(contains(self.fields,'FrictionWaterPressure')) && ispwHydroarma)
				indPwarma	= find(contains(self.fields,'FrictionWaterPressure')); %index of Pwarma, now check for consistency with other arma timesteps 
				dimensions(indPwarma) = md.hydrology.num_basins;
				if(md.hydrology.arma_timestep<self.stochastictimestep)
					error('hydrologyarmapw cannot have a timestep shorter than stochastictimestep');
				end
			end
			size_tot = sum(dimensions);

			%%% Check consistency between ARMA models %%%
			if(indBDWarma~=-1)
				if(indPwarma~=-1)
					if(md.basalforcings.arma_timestep~=md.hydrology.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indBDWarma-1)):sum(dimensions(1:indBDWarma)),1+sum(dimensions(1:indPwarma-1)):sum(dimensions(1:indPwarma))),1,[]);
							if any(crossentries~=0)
								error('BasalforcingsDeepwaterMeltingRatearma and hydrologyarmapw have different arma_timestep and non-zero covariance');
							end
						end
					end
				elseif(indSdarma~=-1)
					if(md.frontalforcings.sd_arma_timestep~=md.basalforcings.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indSdarma-1)):sum(dimensions(1:indSdarma)),1+sum(dimensions(1:indBDWarma-1)):sum(dimensions(1:indBDWarma))),1,[]);
							if any(crossentries~=0)
								error('FrontalForcingsSubglacialDischargearma and BasalforcingsDeepwaterMeltingRatearma have different arma_timestep and non-zero covariance');
							end
						end
					end
				elseif(indSMBarma~=-1)
					if(md.smb.arma_timestep~=md.basalforcings.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indSMBarma-1)):sum(dimensions(1:indSMBarma)),1+sum(dimensions(1:indBDWarma-1)):sum(dimensions(1:indBDWarma))),1,[]);
							if any(crossentries~=0)
								error('SMBarma and BasalforcingsDeepwaterMeltingRatearma have different arma_timestep and non-zero covariance');
							end
						end
					end
				elseif(indTFarma~=-1)
					if(md.frontalforcings.arma_timestep~=md.basalforcings.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indTFarma-1)):sum(dimensions(1:indTFarma)),1+sum(dimensions(1:indBDWarma-1)):sum(dimensions(1:indBDWarma))),1,[]);
							if any(crossentries~=0)
								error('FrontalForcingsRignotarma and BasalforcingsDeepwaterMeltingRatearma have different arma_timestep and non-zero covariance');
							end
						end
					end
				end
			elseif(indPwarma~=-1)
				if(indSdarma~=-1)
					if(md.frontalforcings.sd_arma_timestep~=md.hydrology.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indSdarma-1)):sum(dimensions(1:indSdarma)),1+sum(dimensions(1:indPwarma-1)):sum(dimensions(1:indPwarma))),1,[]);
							if any(crossentries~=0)
								error('FrontalForcingsSubglacialDischargearma and hydrologyarmapw have different arma_timestep and non-zero covariance');
							end
						end
					end
				elseif(indSMBarma~=-1)
					if(md.smb.arma_timestep~=md.hydrology.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indSMBarma-1)):sum(dimensions(1:indSMBarma)),1+sum(dimensions(1:indPwarma-1)):sum(dimensions(1:indPwarma))),1,[]);
							if any(crossentries~=0)
								error('SMBarma and hydrologyarmapw have different arma_timestep and non-zero covariance');
							end
						end
					end
				elseif(indTFarma~=-1)
					if(md.frontalforcings.arma_timestep~=md.hydrology.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indTFarma-1)):sum(dimensions(1:indTFarma)),1+sum(dimensions(1:indPwarma-1)):sum(dimensions(1:indPwarma))),1,[]);
							if any(crossentries~=0)
								error('FrontalForcingsRignotarma and hydrologyarmapw have different arma_timestep and non-zero covariance');
							end
						end
					end
				end
			elseif(indSdarma~=-1)
				if(indSMBarma~=-1)
					if(md.smb.arma_timestep~=md.frontalforcings.sd_arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indSMBarma-1)):sum(dimensions(1:indSMBarma)),1+sum(dimensions(1:indSdarma-1)):sum(dimensions(1:indSdarma))),1,[]);
							if any(crossentries~=0)
								error('SMBarma and FrontalForcingsSubglacialDischargearma have different arma_timestep and non-zero covariance');
							end
						end
					end
				elseif(indTFarma~=-1)
					if(md.frontalforcings.sd_arma_timestep~=md.frontalforcings.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indSdarma-1)):sum(dimensions(1:indSdarma)),1+sum(dimensions(1:indTFarma-1)):sum(dimensions(1:indTFarma))),1,[]);
							if any(crossentries~=0)
								error('FrontalForcingsRignotarma and FrontalForcingsSubglacialDischargearma have different arma_timestep and non-zero covariance');
							end
						end
					end
				end
			elseif(indSMBarma~=-1)
				if(indTFarma~=-1)
					if(md.smb.arma_timestep~=md.frontalforcings.arma_timestep)
						for(ii=[1:numel(lsCovmats)])
							covm = lsCovmats{ii};
							crossentries = reshape(covm(1+sum(dimensions(1:indSMBarma-1)):sum(dimensions(1:indSMBarma)),1+sum(dimensions(1:indTFarma-1)):sum(dimensions(1:indTFarma))),1,[]);
							if any(crossentries~=0)
								error('SMBarma and FrontalForcingsRignotarma have different arma_timestep and non-zero covariance');
							end
						end
					end
				end
			end
			%%% End of consistency checks between ARMA models %%%

			md = checkfield(md,'fieldname','stochasticforcing.isstochasticforcing','values',[0 1]);
			md = checkfield(md,'fieldname','stochasticforcing.fields','numel',num_fields,'cell',1,'values',supportedstochforcings());
			md = checkfield(md,'fieldname','stochasticforcing.covariance','NaN',1,'Inf',1,'size',[size_tot,size_tot,numtcovmat]); %global covariance matrix
			md = checkfield(md,'fieldname','stochasticforcing.stochastictimestep','NaN',1,'Inf',1,'>=',md.timestepping.time_step);
			md = checkfield(md,'fieldname','stochasticforcing.randomflag','numel',[1],'values',[0 1]);
			if(numtcovmat>1) %check the time steps at which each covariance matrix starts to be applied
				md = checkfield(md,'fieldname','stochasticforcing.timecovariance','NaN',1,'Inf',1,'>=',md.timestepping.start_time,'<=',md.timestepping.final_time,'size',[1,numtcovmat]);
			end
			if(checkdefaults) %need to check the defaults
				md = checkfield(md,'fieldname','stochasticforcing.defaultdimension','numel',1,'NaN',1,'Inf',1,'>',0);
				md = checkfield(md,'fieldname','stochasticforcing.default_id','Inf',1,'>=',0,'<=',self.defaultdimension,'size',[md.mesh.numberofelements,1]);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   stochasticforcing parameters:'));
			fielddisplay(self,'isstochasticforcing','is stochasticity activated?');
			fielddisplay(self,'fields','fields with stochasticity applied, ex: [{''SMBarma''}], or [{''SMBforcing''},{''DefaultCalving''}]');
			fielddisplay(self,'defaultdimension','dimensionality of the noise terms (does not apply to fields with their specific dimension)');
			fielddisplay(self,'default_id','id of each element for partitioning of the noise terms (does not apply to fields with their specific partition)');
			fielddisplay(self,'covariance',{'covariance matrix for within- and between-fields covariance (units must be squared field units)','multiple matrices can be concatenated along 3rd dimension to apply different covariances in time'}); 
			fielddisplay(self,'timecovariance','starting dates at which covariances apply (only applicabe if multiple covariance matrices are prescribed)'); 
			fielddisplay(self,'stochastictimestep','timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step)');
			fielddisplay(self,'randomflag','whether to apply real randomness (true) or pseudo-randomness with fixed seed (false)');
			disp('Available fields:');
			for field=supportedstochforcings()
				fprintf('   %s\n',string(field));
			end
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;
			num_fields = numel(self.fields);

			WriteData(fid,prefix,'object',self,'fieldname','isstochasticforcing','format','Boolean');
			if ~self.isstochasticforcing 
				return
			else

				if(self.stochastictimestep==0)
					disp('      stochasticforcing.stocahstictimestep not specified: set to md.timestepping.time_step');
					self.stochastictimestep = md.timestepping.time_step; %by default: stochastictimestep set to ISSM time step
				end

				% Check if hydrolgyarmapw is used
				if(strcmp(class(md.hydrology),'hydrologyarmapw') && md.transient.ishydrology==1)
         	   ispwHydroarma = 1;
         	else
         	   ispwHydroarma = 0;
         	end
				%Retrieve dimensionality of each field
				dimensions = self.defaultdimension*ones(1,num_fields);
				ind = 1;
				for field=self.fields
					%Checking for specific dimensions
					if(strcmp(field,'SMBarma'))
						dimensions(ind) = md.smb.num_basins;
					end
					if(strcmp(field,'FrontalForcingsRignotarma'))
						dimensions(ind) = md.frontalforcings.num_basins;
					end
					if(strcmp(field,'FrontalForcingsSubglacialDischargearma'))
						dimensions(ind) = md.frontalforcings.num_basins;
					end
					if(strcmp(field,'BasalforcingsDeepwaterMeltingRatearma'))
						dimensions(ind) = md.basalforcings.num_basins;
					end
					if(strcmp(field,'BasalforcingsDeepwaterMeltingRatearma'))
						dimensions(ind) = md.basalforcings.num_basins;
					end
					if(strcmp(field,'FrictionWaterPressure') && ispwHydroarma)
						dimensions(ind) = md.hydrology.num_basins;
				   end
					ind = ind+1;
				end

   			if(numel(size(self.covariance))==3)
   				[nrow,ncol,numtcovmat] = size(self.covariance);
   				lsCovmats = {};
   				for ii=[1:numtcovmat] %loop over 3rd dimension
   					lsCovmats{ii} = self.covariance(:,:,ii);
   				end
					if(md.timestepping.interp_forcing==1)
					   disp('WARNING: md.timestepping.interp_forcing is 1, but be aware that there is no interpolation between covariance matrices');
					   disp('         the changes between covariance matrices occur at the time steps specified in md.stochasticforcing.timecovariance');
					end
				elseif(numel(size(self.covariance)==2))
   				[nrow,ncol] = size(self.covariance);
   				numtcovmat = 1; %number of covariance matrices in time
   				lsCovmats = {self.covariance};
   			end
   
				%Scaling covariance matrix (scale column-by-column and row-by-row)
				scaledfields = {'BasalforcingsDeepwaterMeltingRatearma','BasalforcingsSpatialDeepwaterMeltingRate','DefaultCalving','FloatingMeltRate','SMBarma','SMBforcing'}; %list of fields that need scaling *1/yts
				tempcovariance2d = zeros(numtcovmat,sum(nrow*ncol)); %covariance matrices in 2d array
				% Loop over covariance matrices %
				for kk=[1:numtcovmat]
					kkcov = self.covariance(:,:,kk); %extract covariance at index kk
					% Loop over the fields %
					for i=1:num_fields
						if any(strcmp(scaledfields,self.fields(i)))
							inds = [1+sum(dimensions(1:i-1)):1:sum(dimensions(1:i))];
							for row=inds %scale rows corresponding to scaled field
								kkcov(row,:) = 1./yts*kkcov(row,:);
							end
							for col=inds %scale columns corresponding to scaled field
								kkcov(:,col) = 1./yts*kkcov(:,col);
							end
						end
					end
					% Save scaled covariance %
					for rr=[1:nrow]
						ind0 = 1+(rr-1)*ncol;
						tempcovariance2d(kk,ind0:ind0+ncol-1) = kkcov(rr,:);
					end
				end
				%Set dummy default_id vector if defaults not used
				if isnan(self.default_id)
					self.default_id = zeros(md.mesh.numberofelements,1);
				end
				%Set dummy timecovariance vector if a single covariance matrix is used
				if(numtcovmat==1)
					self.timecovariance = [md.timestepping.start_time];
				end

				WriteData(fid,prefix,'data',num_fields,'name','md.stochasticforcing.num_fields','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','fields','format','StringArray');
				WriteData(fid,prefix,'data',dimensions,'name','md.stochasticforcing.dimensions','format','IntMat');
				WriteData(fid,prefix,'object',self,'fieldname','default_id','data',self.default_id-1,'format','IntMat','mattype',2); %0-indexed
				WriteData(fid,prefix,'object',self,'fieldname','defaultdimension','format','Integer');
				WriteData(fid,prefix,'data',numtcovmat,'name','md.stochasticforcing.num_timescovariance','format','Integer');
				WriteData(fid,prefix,'data',tempcovariance2d,'name','md.stochasticforcing.covariance','format','DoubleMat');
				WriteData(fid,prefix,'object',self,'fieldname','timecovariance','format','DoubleMat','scale',yts);
				WriteData(fid,prefix,'object',self,'fieldname','stochastictimestep','format','Double','scale',yts);
				WriteData(fid,prefix,'object',self,'fieldname','randomflag','format','Boolean');
			end
		end % }}}
	end
end
function list = supportedstochforcings() % {{{
	% Defines list of fields supported
	% by the class md.stochasticforcing

	list = structstochforcing();
	list = list.fields;
end % }}}
function structure = structstochforcing() % {{{
	% Defines structure with list of fields
	% supported and corresponding md names
	structure.fields = {...
		'BasalforcingsDeepwaterMeltingRatearma',...
		'BasalforcingsSpatialDeepwaterMeltingRate',...
		'DefaultCalving',...
		'FloatingMeltRate',...
		'FrictionWaterPressure',...
		'FrictionWaterPressure',...
		'FrictionWaterPressure',...
		'FrontalForcingsRignotarma',...
		'FrontalForcingsSubglacialDischargearma',...
		'SMBarma',...
		'SMBforcing'
		};
	structure.mdnames = {...
		'linearbasalforcingsarma',...
		'spatiallinearbasalforcings',...
		'calving',...
		'basalforcings',...
		'friction',...
		'frictioncoulomb',...
		'frictionschoof',...
		'frontalforcingsrignotarma',...
		'frontalforcingsrignotarma',...
		'SMBarma',...
		'SMBforcing'
	};
end % }}}
