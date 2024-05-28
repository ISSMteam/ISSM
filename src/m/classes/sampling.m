%SAMPLING class definition
%
%   Usage:
%      sampling=sampling();

classdef sampling
	properties (SetAccess=public)
		kappa             = NaN;
		tau               = 0;
		beta              = NaN;
		phi               = NaN;
		alpha             = 0;
		robin             = 0;
		seed              = 0;
		requested_outputs = {};
	end
	methods
		function self = sampling(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{

			disp(sprintf('   Sampling parameters:'));

			disp(sprintf('\n      %s','Parameters of PDE operator (kappa^2 I-Laplacian)^(alpha/2)(tau):'));
			fielddisplay(self,'kappa','coefficient of the identity operator');
			fielddisplay(self,'tau','scaling coefficient of the solution');
			fielddisplay(self,'alpha','exponent in PDE operator, (default: 2.0, BiLaplacian covariance operator)');

			disp(sprintf('\n      %s','Parameters of Robin boundary conditions nabla () \cdot normvec + beta ():'));
			fielddisplay(self,'robin','Apply Robin boundary conditions (1 if applied and 0 for homogenous Neumann boundary conditions) (default: 0)');
			fielddisplay(self,'beta','Coefficient in Robin boundary conditions (to be defined for robin = 1)');

			disp(sprintf('\n      %s','Parameters for first-order autoregressive process (X_t = phi X_{t-1} + noise) (if transient):'));
			fielddisplay(self,'phi','Temporal correlation factor (|phi|<1 for stationary process, phi = 1 for random walk process) (default 0)');

			disp(sprintf('\n      %s','Other parameters of stochastic sampler:'));
			fielddisplay(self,'seed','Seed for pseudorandom number generator (given seed if >=0 and random seed if <0) (default: -1)');
			fielddisplay(self,'requested_outputs','additional outputs requested (not implemented yet)');

		end % }}}','
		function self = setdefaultparameters(self) % {{{

			%Apply Robin boundary conditions
			self.robin=0;

			%Exponent in fraction SPDE (default: 2, biLaplacian covariance
			%operator)
			self.alpha=2; % Default 

			%Seed for pseudorandom number generator (default: -1, for random seed)
			self.seed=-1;

			%default output
			self.requested_outputs={'default'};

		end % }}}
		function list = defaultoutputs(self,md) % {{{

			list = {};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SamplingAnalysis',analyses), return; end

			md = checkfield(md,'fieldname','sampling.kappa','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1],'>',0);
			md = checkfield(md,'fieldname','sampling.tau','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1],'>',0);
			md = checkfield(md,'fieldname','sampling.robin','numel',1,'values',[0 1]);
			if(md.sampling.robin)
				md = checkfield(md,'fieldname','sampling.beta','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1],'>',0);
			end
			md = checkfield(md,'fieldname','sampling.alpha','NaN',1,'Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','sampling.seed','NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','sampling.requested_outputs','stringrow',1);

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'object',self,'fieldname','kappa','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','tau','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','beta','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','phi','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','alpha','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','robin','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','seed','format','Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.sampling.requested_outputs','format','StringArray');
		end % }}}
		function md = setparameters(self,md,lc,sigma) % {{{

			nu = self.alpha-1;
			KAPPA = sqrt(8*nu)./lc;
			TAU = sqrt(gamma(nu)./(gamma(self.alpha)*(4*pi)*KAPPA.^(2*nu).*sigma.^2));
			md.sampling.kappa = KAPPA.*ones(md.mesh.numberofvertices,1);
			md.sampling.tau = TAU.*ones(md.mesh.numberofvertices,1);

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejsdouble(fid,[modelname '.sampling.kappa'],self.kappa);
			writejsdouble(fid,[modelname '.sampling.tau'],self.tau);
			writejsdouble(fid,[modelname '.sampling.beta'],self.beta);
			writejsdouble(fid,[modelname '.sampling.phi'],self.phi);
			writejsdouble(fid,[modelname '.sampling.alpha'],self.alpha);
			writejsdouble(fid,[modelname '.sampling.robin'],self.robin);
			writejsdouble(fid,[modelname '.sampling.seed'],self.seed);
			writejscellstring(fid,[modelname '.sampling.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
