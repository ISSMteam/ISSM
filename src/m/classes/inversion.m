%INVERSION class definition
%
%   Usage:
%      inversion=inversion();

classdef inversion
	properties (SetAccess=public) 
		iscontrol                   = 0
		incomplete_adjoint          = 0
		control_parameters          = NaN
		nsteps                      = 0
		maxiter_per_step            = NaN
		cost_functions              = NaN
		cost_functions_coefficients = NaN
		gradient_scaling            = NaN
		cost_function_threshold     = 0
		min_parameters              = NaN
		max_parameters              = NaN
		step_threshold              = NaN
		vx_obs                      = NaN
		vy_obs                      = NaN
		vz_obs                      = NaN
		vel_obs                     = NaN
		thickness_obs               = NaN
		surface_obs                 = NaN
	end
	methods
		function self = inversion(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self =structtoobj(inversion(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.vx_obs=project3d(md,'vector',self.vx_obs,'type','node');
			self.vy_obs=project3d(md,'vector',self.vy_obs,'type','node');
			self.vel_obs=project3d(md,'vector',self.vel_obs,'type','node');
			self.thickness_obs=project3d(md,'vector',self.thickness_obs,'type','node');
			if numel(self.cost_functions_coefficients)>1,self.cost_functions_coefficients=project3d(md,'vector',self.cost_functions_coefficients,'type','node');end;
			if numel(self.min_parameters)>1,self.min_parameters=project3d(md,'vector',self.min_parameters,'type','node');end;
			if numel(self.max_parameters)>1,self.max_parameters=project3d(md,'vector',self.max_parameters,'type','node');end;
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%default is incomplete adjoint for now
			self.incomplete_adjoint=1;

			%parameter to be inferred by control methods (only
			%drag and B are supported yet)
			self.control_parameters={'FrictionCoefficient'};

			%number of steps in the control methods
			self.nsteps=20;

			%maximum number of iteration in the optimization algorithm for
			%each step
			self.maxiter_per_step=20*ones(self.nsteps,1);

			%the inversed parameter is updated as follows:
			%new_par=old_par + gradient_scaling(n)*C*gradient with C in [0 1];
			%usually the gradient_scaling must be of the order of magnitude of the 
			%inversed parameter (10^8 for B, 50 for drag) and can be decreased
			%after the first iterations
			self.gradient_scaling=50*ones(self.nsteps,1);

			%several responses can be used:
			self.cost_functions=101;

			%step_threshold is used to speed up control method. When
			%misfit(1)/misfit(0) < self.step_threshold, we go directly to
			%the next step
			self.step_threshold=.7*ones(self.nsteps,1); %30 per cent decrement.

			%cost_function_threshold is a criteria to stop the control methods.
			%if J[n]-J[n-1]/J[n] < criteria, the control run stops
			%NaN if not applied
			self.cost_function_threshold=NaN; %not activated

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~self.iscontrol, return; end

			num_controls=numel(md.inversion.control_parameters);
			num_costfunc=size(md.inversion.cost_functions,2);

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.incomplete_adjoint','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols());
			md = checkfield(md,'fieldname','inversion.nsteps','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.maxiter_per_step','size',[md.inversion.nsteps 1],'>=',0);
			md = checkfield(md,'fieldname','inversion.step_threshold','size',[md.inversion.nsteps 1]);
			md = checkfield(md,'fieldname','inversion.cost_functions','size',[1 num_costfunc],'values',supportedcostfunctions());
			md = checkfield(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices num_costfunc],'>=',0);
			md = checkfield(md,'fieldname','inversion.gradient_scaling','size',[md.inversion.nsteps num_controls]);
			md = checkfield(md,'fieldname','inversion.min_parameters','size',[NaN num_controls]);
			md = checkfield(md,'fieldname','inversion.max_parameters','size',[NaN num_controls]);

			%Only SSA, HO and FS are supported right now
			if strcmp(solution,'StressbalanceSolution')
				if ~(md.flowequation.isSSA || md.flowequation.isMOLHO || md.flowequation.isHO || md.flowequation.isFS || md.flowequation.isL1L2),
					md = checkmessage(md,['inversion can only be performed for SSA, MOLHO, HO or FS ice flow models']);
				end
			end
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
			disp(sprintf('   inversion parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'incomplete_adjoint','1: linear viscosity, 0: non-linear viscosity');
			fielddisplay(self,'control_parameters','ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}');
			fielddisplay(self,'nsteps','number of optimization searches');
			fielddisplay(self,'cost_functions','indicate the type of response for each optimization step');
			fielddisplay(self,'cost_functions_coefficients','cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter');
			fielddisplay(self,'cost_function_threshold','misfit convergence criterion. Default is 1%, NaN if not applied');
			fielddisplay(self,'maxiter_per_step','maximum iterations during each optimization step');
			fielddisplay(self,'gradient_scaling','scaling factor on gradient direction during optimization, for each optimization step');
			fielddisplay(self,'step_threshold','decrease threshold for misfit, default is 30%');
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

			WriteData(fid,prefix,'name','md.inversion.type','data',0,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','iscontrol','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','incomplete_adjoint','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','vel_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'fieldname','nsteps','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','maxiter_per_step','format','IntMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','cost_functions_coefficients','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','gradient_scaling','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','cost_function_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','min_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','max_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','step_threshold','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','vx_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','vy_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','vz_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			if(numel(self.thickness_obs)==md.mesh.numberofelements),
				mattype=2;
			else
				mattype=1;
			end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','thickness_obs','format','DoubleMat','mattype',mattype);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','surface_obs','format','DoubleMat','mattype',mattype);

			%process control parameters
			WriteData(fid,prefix,'object',self,'fieldname','control_parameters','format','StringArray');
			WriteData(fid,prefix,'data',numel(self.control_parameters),'name','md.inversion.num_control_parameters','format','Integer');

			%process cost functions
			num_cost_functions=size(self.cost_functions,2);
			data=marshallcostfunctions(self.cost_functions);
			WriteData(fid,prefix,'data',data,'name','md.inversion.cost_functions','format','StringArray');
			WriteData(fid,prefix,'data',num_cost_functions,'name','md.inversion.num_cost_functions','format','Integer');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.inversion.iscontrol'],self.iscontrol);
			writejsdouble(fid,[modelname '.inversion.incomplete_adjoint'],self.incomplete_adjoint);
			writejscellstring(fid,[modelname '.inversion.control_parameters'],self.control_parameters);
			writejsdouble(fid,[modelname '.inversion.nsteps'],self.nsteps);
			writejs1Darray(fid,[modelname '.inversion.maxiter_per_step'],self.maxiter_per_step);
			writejs2Darray(fid,[modelname '.inversion.cost_functions'],self.cost_functions);
			writejs2Darray(fid,[modelname '.inversion.cost_functions_coefficients'],self.cost_functions_coefficients);
			writejs1Darray(fid,[modelname '.inversion.min_parameters'],self.min_parameters);
			writejs1Darray(fid,[modelname '.inversion.max_parameters'],self.max_parameters);
			writejs1Darray(fid,[modelname '.inversion.vx_obs'],self.vx_obs);
			writejs1Darray(fid,[modelname '.inversion.vy_obs'],self.vy_obs);
			writejs1Darray(fid,[modelname '.inversion.vz_obs'],self.vz_obs);
			writejs1Darray(fid,[modelname '.inversion.vel_obs'],self.vel_obs);
			writejs1Darray(fid,[modelname '.inversion.thickness_obs'],self.thickness_obs);
			writejs1Darray(fid,[modelname '.inversion.surface_obs'],self.surface_obs);
		end % }}}
	end
end
