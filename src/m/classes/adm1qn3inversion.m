%ADM1QN3INVERSION class definition
%
%   Usage:
%      adm1qn3inversion=adm1qn3inversion();

classdef adm1qn3inversion
	properties (SetAccess=public) 
		iscontrol                   = 0
		maxsteps                    = 0
		maxiter                     = 0
		dxmin                       = 0
		dfmin_frac                  = 0
		gttol                       = 0

	end
	methods
		function self = adm1qn3inversion(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(adm1qn3inversion(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
		function self = setdefaultparameters(self) % {{{


			%parameter to be inferred by control methods (only
			%drag and B are supported yet)
			%self.control_parameters={'FrictionCoefficient'};

			%number of iterations
			self.maxsteps=20;
			self.maxiter=40;

			%m1qn3 parameters
			self.dxmin      = 0.1;
			self.dfmin_frac = 1.;
			self.gttol      = 1e-4;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~self.iscontrol, return; end

			if ~IssmConfig('_HAVE_M1QN3_'),
				md = checkmessage(md,['M1QN3 has not been installed, ISSM needs to be reconfigured and recompiled with M1QN3']);
			end

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.maxsteps','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.maxiter','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.dxmin','numel',1,'>',0);
         md = checkfield(md,'fieldname','inversion.dfmin_frac','numel',1,'>=',0., '<=', 1.);
			md = checkfield(md,'fieldname','inversion.gttol','numel',1,'>',0);

	end % }}}
		function disp(self) % {{{
			disp(sprintf('   adm1qn3inversion parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'maxsteps','maximum number of iterations (gradient computation)');
			fielddisplay(self,'maxiter','maximum number of Function evaluation (forward run)');
			fielddisplay(self,'dxmin','convergence criterion: two points less than dxmin from eachother (sup-norm) are considered identical');
         fielddisplay(self,'dfmin_frac','expected reduction of J during the first step (e.g., 0.3=30% reduction in cost function)');
			fielddisplay(self,'gttol','convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','iscontrol','format','Boolean');
			WriteData(fid,prefix,'name','md.inversion.type','data',4,'format','Integer');
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxsteps','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','dxmin','format','Double');
         WriteData(fid,prefix,'object',self,'class','inversion','fieldname','dfmin_frac','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','gttol','format','Double');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.inversion.iscontrol'],self.iscontrol);
			writejsdouble(fid,[modelname '.inversion.maxsteps'],self.maxsteps);
			writejsdouble(fid,[modelname '.inversion.maxiter'],self.maxiter);
			writejsdouble(fid,[modelname '.inversion.dxmin'],self.dxmin);
         writejsdouble(fid,[modelname '.inversion.dfmin_frac'],self.dfmin_frac);
			writejsdouble(fid,[modelname '.inversion.gttol'],self.gttol);

		end % }}}
	end
end
