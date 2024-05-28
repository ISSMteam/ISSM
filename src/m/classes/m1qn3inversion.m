%M1QN3INVERSION class definition
%
%   Usage:
%      m1qn3inversion=m1qn3inversion();

classdef m1qn3inversion
	properties (SetAccess=public) 
		iscontrol                   = 0
		incomplete_adjoint          = 0
		control_parameters          = NaN
		control_scaling_factors     = NaN
		maxsteps                    = 0
		maxiter                     = 0
		dxmin                       = 0
		dfmin_frac                  = 0
		gttol                       = 0
		cost_functions              = NaN
		cost_functions_coefficients = NaN
		min_parameters              = NaN
		max_parameters              = NaN
		vx_obs                      = NaN
		vy_obs                      = NaN
		vz_obs                      = NaN
		vel_obs                     = NaN
		thickness_obs               = NaN
		surface_obs                 = NaN

	end
	methods
		function self = m1qn3inversion(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(m1qn3inversion(),varargin{1});
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

			%Scaling factor for each control
			self.control_scaling_factors=1;

			%number of iterations
			self.maxsteps=20;
			self.maxiter=40;

			%several responses can be used:
			self.cost_functions=101;

			%m1qn3 parameters
			self.dxmin      = 0.1;
			self.dfmin_frac = 1;
			self.gttol      = 1e-4;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~self.iscontrol, return; end

			if ~IssmConfig('_HAVE_M1QN3_'),
				md = checkmessage(md,['M1QN3 has not been installed, ISSM needs to be reconfigured and recompiled with M1QN3']);
			end
			num_controls=numel(md.inversion.control_parameters);
			num_costfunc=size(md.inversion.cost_functions,2);

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.incomplete_adjoint','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols());
			md = checkfield(md,'fieldname','inversion.control_scaling_factors','size',[1 num_controls],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','inversion.maxsteps','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.maxiter','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.dxmin','numel',1,'>',0);
         md = checkfield(md,'fieldname','inversion.dfmin_frac','numel',1,'>=',0,'<=',1);
			md = checkfield(md,'fieldname','inversion.gttol','numel',1,'>',0);
			md = checkfield(md,'fieldname','inversion.cost_functions','size',[1 num_costfunc],'values',supportedcostfunctions());
			md = checkfield(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices num_costfunc],'>=',0);
			md = checkfield(md,'fieldname','inversion.min_parameters','size',[NaN num_controls]);
			md = checkfield(md,'fieldname','inversion.max_parameters','size',[NaN num_controls]);

			if strcmp(solution,'BalancethicknessSolution')
				md = checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','inversion.surface_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			elseif strcmp(solution,'BalancethicknessSoftSolution')
				md = checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			else
				md = checkfield(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				if ~strcmp(domaintype(md.mesh),'2Dvertical'),
					md = checkfield(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   m1qn3inversion parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'incomplete_adjoint','1: linear viscosity, 0: non-linear viscosity');
			fielddisplay(self,'control_parameters','ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}');
			fielddisplay(self,'control_scaling_factors','order of magnitude of each control (useful for multi-parameter optimization)');
			fielddisplay(self,'maxsteps','maximum number of iterations (gradient computation)');
			fielddisplay(self,'maxiter','maximum number of Function evaluation (forward run)');
			fielddisplay(self,'dxmin','convergence criterion: two points less than dxmin from eachother (sup-norm) are considered identical');
         fielddisplay(self,'dfmin_frac','expected reduction of J during the first step (e.g., 0.3=30% reduction in cost function)');
			fielddisplay(self,'gttol','convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)');
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
			WriteData(fid,prefix,'name','md.inversion.type','data',2,'format','Integer');
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','incomplete_adjoint','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','control_scaling_factors','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxsteps','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','dxmin','format','Double');
         WriteData(fid,prefix,'object',self,'class','inversion','fieldname','dfmin_frac','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','gttol','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','cost_functions_coefficients','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','min_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','vx_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','vy_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','vz_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','vel_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			if(numel(self.thickness_obs)==md.mesh.numberofelements),
				mattype=2;
			else
				mattype=1;
			end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','thickness_obs','format','DoubleMat','mattype',mattype);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','surface_obs','format','DoubleMat','mattype',mattype);

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
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.inversion.iscontrol'],self.iscontrol);
			writejsdouble(fid,[modelname '.inversion.incomplete_adjoint'],self.incomplete_adjoint);
			writejscellstring(fid,[modelname '.inversion.control_parameters'],self.control_parameters);
			writejsdouble(fid,[modelname '.inversion.control_scaling_factors'],self.control_scaling_factors);
			writejsdouble(fid,[modelname '.inversion.maxsteps'],self.maxsteps);
			writejsdouble(fid,[modelname '.inversion.maxiter'],self.maxiter);
			writejsdouble(fid,[modelname '.inversion.dxmin'],self.dxmin);
         writejsdouble(fid,[modelname '.inversion.dfmin_frac'],self.dfmin_frac);
			writejsdouble(fid,[modelname '.inversion.gttol'],self.gttol);
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
