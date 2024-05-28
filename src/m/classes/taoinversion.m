%TAOINVERSION class definition
%
%   Usage:
%      taoinversion=taoinversion();

classdef taoinversion
	properties (SetAccess=public) 
		iscontrol                   = 0
		incomplete_adjoint          = 0
		control_parameters          = NaN
		maxsteps                    = 0
		maxiter                     = 0
		fatol                       = 0
		frtol                       = 0
		gatol                       = 0
		grtol                       = 0
		gttol                       = 0
		algorithm                   = ''
		cost_functions              = NaN
		cost_functions_coefficients = NaN
		min_parameters              = NaN
		max_parameters              = NaN
		vx_obs                      = NaN
		vy_obs                      = NaN
		vz_obs                      = NaN
		vel_obs                     = NaN
		thickness_obs               = NaN
		surface_obs               = NaN
	end
	methods
		function self = extrude(self,md) % {{{
			self.vx_obs=project3d(md,'vector',self.vx_obs,'type','node');
			self.vy_obs=project3d(md,'vector',self.vy_obs,'type','node');
			self.vel_obs=project3d(md,'vector',self.vel_obs,'type','node');
			self.thickness_obs=project3d(md,'vector',self.thickness_obs,'type','node');
			if numel(self.cost_functions_coefficients)>1,self.cost_functions_coefficients=project3d(md,'vector',self.cost_functions_coefficients,'type','node');end;
			if numel(self.min_parameters)>1,self.min_parameters=project3d(md,'vector',self.min_parameters,'type','node');end;
			if numel(self.max_parameters)>1,self.max_parameters=project3d(md,'vector',self.max_parameters,'type','node');end;
		end % }}}
		function self = taoinversion(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(taoinversion(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%default is incomplete adjoint for now
			self.incomplete_adjoint=1;

			%parameter to be inferred by control methods (only
			%drag and B are supported yet)
			self.control_parameters={'FrictionCoefficient'};

			%number of iterations and steps
			self.maxsteps=20;
			self.maxiter =30;

			%default tolerances
			self.fatol = 0;
			self.frtol = 0;
			self.gatol = 0;
			self.grtol = 0;
			self.gttol = 1e-4;

			%minimization algorithm
			PETSCMAJOR = IssmConfig('_PETSC_MAJOR_');
			PETSCMINOR = IssmConfig('_PETSC_MINOR_');
			if(PETSCMAJOR>3 | (PETSCMAJOR==3 & PETSCMINOR>=5))
				self.algorithm = 'blmvm';
			else
				self.algorithm = 'tao_blmvm';
			end

			%several responses can be used:
			self.cost_functions=101;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~self.iscontrol, return; end

			if ~IssmConfig('_HAVE_TAO_'),
				md = checkmessage(md,['TAO has not been installed, ISSM needs to be reconfigured and recompiled with TAO']);
			end

			num_controls=numel(md.inversion.control_parameters);
			num_costfunc=size(md.inversion.cost_functions,2);

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.incomplete_adjoint','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols());
			md = checkfield(md,'fieldname','inversion.maxsteps','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.maxiter','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.fatol','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.frtol','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.gatol','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.grtol','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.gttol','numel',1,'>=',0);

			PETSCMAJOR = IssmConfig('_PETSC_MAJOR_');
			PETSCMINOR = IssmConfig('_PETSC_MINOR_');
			if(PETSCMAJOR>3 | (PETSCMAJOR==3 & PETSCMINOR>=5))
				md = checkfield(md,'fieldname','inversion.algorithm','values',{'blmvm','cg','lmvm'});
			else
				md = checkfield(md,'fieldname','inversion.algorithm','values',{'tao_blmvm','tao_cg','tao_lmvm'});
			end

			md = checkfield(md,'fieldname','inversion.cost_functions','size',[1 num_costfunc],'values',supportedcostfunctions());
			md = checkfield(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices num_costfunc],'>=',0);
			md = checkfield(md,'fieldname','inversion.min_parameters','size',[NaN num_controls]);
			md = checkfield(md,'fieldname','inversion.max_parameters','size',[NaN num_controls]);

			if strcmp(solution,'BalancethicknessSolution')
				md = checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			elseif strcmp(solution,'BalancethicknessSoftSolution')
				md = checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			else
				md = checkfield(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   taoinversion parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'incomplete_adjoint','1: linear viscosity, 0: non-linear viscosity');
			fielddisplay(self,'control_parameters','ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}');
			fielddisplay(self,'maxsteps','maximum number of iterations (gradient computation)');
			fielddisplay(self,'maxiter','maximum number of Function evaluation (forward run)');
			fielddisplay(self,'fatol','convergence criterion: f(X)-f(X*) (X: current iteration, X*: "true" solution, f: cost function)');
			fielddisplay(self,'frtol','convergence criterion: |f(X)-f(X*)|/|f(X*)|');
			fielddisplay(self,'gatol','convergence criterion: ||g(X)|| (g: gradient of the cost function)');
			fielddisplay(self,'grtol','convergence criterion: ||g(X)||/|f(X)|');
			fielddisplay(self,'gttol','convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)');
			fielddisplay(self,'algorithm','minimization algorithm: ''tao_blmvm'', ''tao_cg'', ''tao_lmvm''');
			fielddisplay(self,'cost_functions','indicate the type of response for each optimization step');
			fielddisplay(self,'cost_functions_coefficients','cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter');
			fielddisplay(self,'min_parameters','absolute minimum acceptable value of the inversed parameter on each vertex');
			fielddisplay(self,'max_parameters','absolute maximum acceptable value of the inversed parameter on each vertex');
			fielddisplay(self,'vx_obs','observed velocity x component [m/yr]');
			fielddisplay(self,'vy_obs','observed velocity y component [m/yr]');
			fielddisplay(self,'vel_obs','observed velocity magnitude [m/yr]');
			fielddisplay(self,'thickness_obs','observed thickness [m]');
			fielddisplay(self,'surface_obs','observed surface elevation [m]');
			disp('Available cost functions:');
			disp('   101: SurfaceAbsVelMisfit');
			disp('   102: SurfaceRelVelMisfit');
			disp('   103: SurfaceLogVelMisfit');
			disp('   104: SurfaceLogVxVyMisfit');
			disp('   105: SurfaceAverageVelMisfit');
			disp('   201: ThicknessAbsMisfit');
			disp('   501: DragCoefficientAbsGradient');
			disp('   502: RheologyBbarAbsGradient');
			disp('   503: ThicknessAbsGradient');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','iscontrol','format','Boolean');
			WriteData(fid,prefix,'name','md.inversion.type','data',1,'format','Integer');
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','incomplete_adjoint','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxsteps','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','fatol','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','frtol','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','gatol','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','grtol','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','gttol','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','algorithm','format','String');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','cost_functions_coefficients','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','min_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','vx_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','vy_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','vz_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','thickness_obs','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','surface_obs','format','DoubleMat','mattype',1);

			%process control parameters
			num_control_parameters=numel(self.control_parameters);
			WriteData(fid,prefix,'object',self,'fieldname','control_parameters','format','StringArray');
			WriteData(fid,prefix,'data',num_control_parameters,'name','md.inversion.num_control_parameters','format','Integer');

			%process cost functions
			num_cost_functions=size(self.cost_functions,2);
			data=marshallcostfunctions(self.cost_functions);
			WriteData(fid,prefix,'data',data,'name','md.inversion.cost_functions','format','StringArray');
			WriteData(fid,prefix,'data',num_cost_functions,'name','md.inversion.num_cost_functions','format','Integer');
		end % }}}
	end
end
