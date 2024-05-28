%INVERSIONVALIDATION class definition
%
%   Usage:
%      inversionvalidation=inversionvalidation();

classdef inversionvalidation
	properties (SetAccess=public) 
		iscontrol                   = 0
		incomplete_adjoint          = 0
		control_parameters          = NaN
		control_scaling_factors     = NaN
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
		function self = extrude(self,md) % {{{
			self.vx_obs=project3d(md,'vector',self.vx_obs,'type','node');
			self.vy_obs=project3d(md,'vector',self.vy_obs,'type','node');
			self.vel_obs=project3d(md,'vector',self.vel_obs,'type','node');
			self.thickness_obs=project3d(md,'vector',self.thickness_obs,'type','node');
			if numel(self.cost_functions_coefficients)>1,self.cost_functions_coefficients=project3d(md,'vector',self.cost_functions_coefficients,'type','node');end;
			if numel(self.min_parameters)>1,self.min_parameters=project3d(md,'vector',self.min_parameters,'type','node');end;
			if numel(self.max_parameters)>1,self.max_parameters=project3d(md,'vector',self.max_parameters,'type','node');end;
		end % }}}
		function self = inversionvalidation(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(inversionvalidation(),varargin{1});
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

			%Scaling factor for each control
			self.control_scaling_factors=1;

			%several responses can be used:
			self.cost_functions=101;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~self.iscontrol, return; end

			num_controls=numel(md.inversion.control_parameters);
			num_costfunc=size(md.inversion.cost_functions,2);

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.incomplete_adjoint','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols());
			md = checkfield(md,'fieldname','inversion.control_scaling_factors','size',[1 num_controls],'>',0,'NaN',1,'Inf',1);
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
				if ~strcmp(domaintype(md.mesh),'2Dvertical'),
					md = checkfield(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   inversionvalidation parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'incomplete_adjoint','1: linear viscosity, 0: non-linear viscosity');
			fielddisplay(self,'control_parameters','ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}');
			fielddisplay(self,'control_scaling_factors','order of magnitude of each control (useful for multi-parameter optimization)');
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
			WriteData(fid,prefix,'name','md.inversion.type','data',3,'format','Integer');
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','incomplete_adjoint','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','control_scaling_factors','format','DoubleMat','mattype',3);
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
	end
end
