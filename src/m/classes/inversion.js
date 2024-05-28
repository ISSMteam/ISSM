//INVERSION class definition
//
//   Usage:
//      inversion=new inversion();

function inversion (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//default is incomplete adjoint for now
		this.incomplete_adjoint=1;

		//parameter to be inferred by control methods (only
		//drag and B are supported yet)
		this.control_parameters=['FrictionCoefficient'];

		//number of steps in the control methods
		this.nsteps=20;

		//maximum number of iteration in the optimization algorithm for
		//each step
		this.maxiter_per_step=20*NewArrayFill(this.nsteps,1);

		//the inversed parameter is updated as follows:
		//new_par=old_par + gradient_scaling(n)*C*gradient with C in [0 1];
		//usually the gradient_scaling must be of the order of magnitude of the 
		//inversed parameter (10^8 for B, 50 for drag) and can be decreased
		//after the first iterations
		this.gradient_scaling=NewArrayFill(this.nsteps,50);

		//several responses can be used:
		this.cost_functions=101;

		//step_threshold is used to speed up control method. When
		//misfit(1)/misfit(0) < this.step_threshold, we go directly to
		//the next step
		this.step_threshold=NewArrayFill(this.nsteps,.7); //30 per cent decrement.

		//cost_function_threshold is a criteria to stop the control methods.
		//if J[n]-J[n-1]/J[n] < criteria, the control run stops
		//NaN if not applied
		this.cost_function_threshold=NaN; //not activated

	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   inversion parameters:'));
		fielddisplay(this,'iscontrol','is inversion activated?');
		fielddisplay(this,'incomplete_adjoint','1: linear viscosity, 0: non-linear viscosity');
		fielddisplay(this,'control_parameters',"ex: {'FrictionCoefficient'}, or {'MaterialsRheologyBbar'}");
		fielddisplay(this,'nsteps','number of optimization searches');
		fielddisplay(this,'cost_functions','indicate the type of response for each optimization step');
		fielddisplay(this,'cost_functions_coefficients','cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter');
		fielddisplay(this,'cost_function_threshold','misfit convergence criterion. Default is 1%, NaN if not applied');
		fielddisplay(this,'maxiter_per_step','maximum iterations during each optimization step');
		fielddisplay(this,'gradient_scaling','scaling factor on gradient direction during optimization, for each optimization step');
		fielddisplay(this,'step_threshold','decrease threshold for misfit, default is 30%');
		fielddisplay(this,'min_parameters','absolute minimum acceptable value of the inversed parameter on each vertex');
		fielddisplay(this,'max_parameters','absolute maximum acceptable value of the inversed parameter on each vertex');
		fielddisplay(this,'vx_obs','observed velocity x component [m/yr]');
		fielddisplay(this,'vy_obs','observed velocity y component [m/yr]');
		fielddisplay(this,'vel_obs','observed velocity magnitude [m/yr]');
		fielddisplay(this,'thickness_obs','observed thickness [m]');
		fielddisplay(this,'surface_obs','observed surface elevation [m]');
		console.log('Available cost functions:');
		console.log('   101: SurfaceAbsVelMisfit');
		console.log('   102: SurfaceRelVelMisfit');
		console.log('   103: SurfaceLogVelMisfit');
		console.log('   104: SurfaceLogVxVyMisfit');
		console.log('   105: SurfaceAverageVelMisfit');
		console.log('   201: ThicknessAbsMisfit');
		console.log('   501: DragCoefficientAbsGradient');
		console.log('   502: RheologyBbarAbsGradient');
		console.log('   503: ThicknessAbsGradient');

	}// }}}
    this.extrude = function(md) {//{{{
        this.vx_obs=project3d(md, 'vector', this.vx_obs, 'type', 'node');
        this.vy_obs=project3d(md,'vector',this.vy_obs,'type','node');
        this.vel_obs=project3d(md,'vector',this.vel_obs,'type','node');
        this.thickness_obs=project3d(md,'vector',this.thickness_obs,'type','node');

        if (this.cost_functions_coefficients.length>1) {
            this.cost_functions_coefficients=project3d(md,'vector',this.cost_functions_coefficients,'type','node');
        }			
        if (this.min_parameters.length>1) {
            this.min_parameters=project3d(md,'vector',this.min_parameters,'type','node');
        }
        if (this.max_parameters.length>1) {
            this.max_parameters=project3d(md,'vector',this.max_parameters,'type','node');
        }
        return this;
    }//}}}
	this.classname= function(){// {{{
		return "inversion";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			//Early return
			if (!this.iscontrol) return;

			num_controls=md.inversion.control_parameters.length;
			num_costfunc=md.inversion.control_parameters[0].length;

			checkfield(md,'fieldname','inversion.iscontrol','values',[0, 1]);
			checkfield(md,'fieldname','inversion.incomplete_adjoint','values',[0 ,1]);
			checkfield(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols());
			checkfield(md,'fieldname','inversion.nsteps','numel',1,'>=',0);
			checkfield(md,'fieldname','inversion.maxiter_per_step','size',[md.inversion.nsteps, 1],'>=',0);
			checkfield(md,'fieldname','inversion.step_threshold','size',[md.inversion.nsteps, 1]);
			checkfield(md,'fieldname','inversion.cost_functions','size',[1, num_costfunc],'values',supportedcostfunctions());
			checkfield(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices, num_costfunc],'>=',0);
			checkfield(md,'fieldname','inversion.gradient_scaling','size',[md.inversion.nsteps, num_controls]);
			checkfield(md,'fieldname','inversion.min_parameters','size',[md.mesh.numberofvertices , num_controls]);
			checkfield(md,'fieldname','inversion.max_parameters','size',[md.mesh.numberofvertices ,num_controls]);

			//Only SSA, HO and FS are supported right now
			if (solution=='StressbalanceSolution'){
				if (!(md.flowequation.isSSA | md.flowequation.isHO | md.flowequation.isFS | md.flowequation.isL1L2)){
					md.checkmessage('inversion can only be performed for SSA, HO or FS ice flow models');
				}
			}

			if (solution=='BalancethicknessSolution'){
				checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices ,1],'NaN',1,'Inf',1);
			}
			else if (solution=='BalancethicknessSoftSolution'){
				checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices, 1],'NaN',1,'Inf',1);
			}
			else{
				checkfield(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices ,1],'NaN',1,'Inf',1);
				checkfield(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices ,1],'NaN',1,'Inf',1);
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{

			var yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.inversion.type','data',0,'format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','iscontrol','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','incomplete_adjoint','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','vel_obs','format','DoubleMat','mattype',1,'scale',1/yts);
			if (!this.iscontrol) return;
			WriteData(fid,prefix,'object',this,'fieldname','nsteps','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','maxiter_per_step','format','IntMat','mattype',3);
			WriteData(fid,prefix,'object',this,'fieldname','cost_functions_coefficients','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','gradient_scaling','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',this,'fieldname','cost_function_threshold','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','min_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',this,'fieldname','max_parameters','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',this,'fieldname','step_threshold','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',this,'fieldname','vx_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',this,'fieldname','vy_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',this,'fieldname','vz_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			if(this.thickness_obs.length==md.mesh.numberofelements) mattype=2;
			else mattype=1;
			WriteData(fid,prefix,'object',this,'class','inversion','fieldname','thickness_obs','format','DoubleMat','mattype',mattype);
			WriteData(fid,prefix,'object',this,'class','inversion','fieldname','surface_obs','format','DoubleMat','mattype',mattype);

			//process control parameters
			num_control_parameters=this.control_parameters.length;
			WriteData(fid,prefix,'object',this,'fieldname','control_parameters','format','StringArray');
			WriteData(fid,prefix,'data',num_control_parameters,'name','md.inversion.num_control_parameters','format','Integer');

			//process cost functions
			num_cost_functions=this.cost_functions[0].length;
			data=marshallcostfunctions(this.cost_functions);
			WriteData(fid,prefix,'data',data,'name','md.inversion.cost_functions','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'data',num_cost_functions,'name','md.inversion.num_cost_functions','format','Integer');
		}//}}}
		this.fix=function() { //{{{
			this.control_parameters=NullFix(this.control_parameters,NaN);
			this.maxiter_per_step=NullFix(this.maxiter_per_step,NaN);
			this.cost_functions=NullFix(this.cost_functions,NaN);
			this.cost_functions_coefficients=NullFix(this.cost_functions_coefficients,NaN);
			this.cost_function_threshold=NullFix(this.cost_function_threshold,NaN);
			this.gradient_scaling=NullFix(this.gradient_scaling,NaN);
			this.min_parameters=NullFix(this.min_parameters,NaN);
			this.max_parameters=NullFix(this.max_parameters,NaN);
			this.step_threshold=NullFix(this.step_threshold,NaN);
			this.vx_obs=NullFix(this.vx_obs,NaN);
			this.vy_obs=NullFix(this.vy_obs,NaN);
			this.vz_obs=NullFix(this.vz_obs,NaN);
			this.vel_obs=NullFix(this.vel_obs,NaN);
			this.thickness_obs=NullFix(this.thickness_obs,NaN);
			this.surface_obs=NullFix(this.surface_obs,NaN);
		}//}}}
	//properties 
	// {{{

	this.iscontrol                   = 0;
	this.incomplete_adjoint          = 0;
	this.control_parameters          = NaN;
	this.nsteps                      = 0;
	this.maxiter_per_step            = NaN;
	this.cost_functions              = NaN;
	this.cost_functions_coefficients = NaN;
	this.gradient_scaling            = NaN;
	this.cost_function_threshold     = 0;
	this.min_parameters              = NaN;
	this.max_parameters              = NaN;
	this.step_threshold              = NaN;
	this.vx_obs                      = NaN;
	this.vy_obs                      = NaN;
	this.vz_obs                      = NaN;
	this.vel_obs                     = NaN;
	this.thickness_obs               = NaN;
	this.surface_obs                 = NaN;

	this.setdefaultparameters();
	//}}}
}
