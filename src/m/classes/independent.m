%INDEPENDENT class definition
%
%   Usage:
%      independent=independent();

classdef independent
	properties (SetAccess=public) 
		name					= '';
		type					= '';
		fos_forward_index		= NaN;
		fov_forward_indices		= [];
		nods					= 0;
		min_parameters			= NaN;
		max_parameters			= NaN;
		control_scaling_factor	= 0;
		control_size			= 0;

	end
	methods
		function self = independent(varargin) % {{{

			%use provided options to change fields
			options=pairoptions(varargin{:});

			%OK get other fields
			self = setdefaultparameters(self);
			self = AssignObjectFields(pairoptions(varargin{:}), self);

			if(self.control_size == 0)
				self.control_size = 1;
			end


		end
		%}}}
		function self = setdefaultparameters(self) % {{{

			%Set default scaling factor to 1 (i.e., do not scale)
			self.control_scaling_factor = 1.0;

		end % }}}
		function md = checkconsistency(self,md,i,solution,analyses,driver) % {{{
			if ~isnan(self.fos_forward_index)
				if ~strcmpi(driver,'fos_forward'),
					error('cannot declare an independent with a fos_forward_index when the driver is not fos_forward!');
				end
				if self.nods==0,
					error('independent checkconsistency error: nods should be set to the size of the independent variable');
				end
			end

			if ~isempty(self.fov_forward_indices)
				if ~strcmpi(driver,'fov_forward'),
					error('cannot declare an independent with fov_forward_indices when the driver is not fov_forward!');
				end
				if self.nods==0,
					error('independent checkconsistency error: nods should be set to the size of the independent variable');
				end
				md = checkfield(md,'fieldname',['autodiff.independents{' num2str(i) '}.fov_forward_indices'],'>=',1,'<=',self.nods,'size',[NaN 1]);
			end

			%md = checkfield(md,'fieldname',['autodiff.independents{' num2str(i) '}.min_parameters'],'size',[md.mesh.numberofvertices 1]);
			%md = checkfield(md,'fieldname',['autodiff.independents{' num2str(i) '}.max_parameters'],'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname',['autodiff.independents{' num2str(i) '}.control_scaling_factor'],'size',[1 1],'>',0,'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   independent variable:'));

			fielddisplay(self,'name','variable name (must match corresponding String)');
			fielddisplay(self,'type','type of variable (''vertex'' or ''scalar'')');
			fielddisplay(self,'nods','size of independent variables');
			fielddisplay(self,'control_size','number of timesteps');
			fielddisplay(self,'min_parameters','absolute minimum acceptable value of the inversed parameter on each vertex');
			fielddisplay(self,'max_parameters','absolute maximum acceptable value of the inversed parameter on each vertex');
			fielddisplay(self,'control_scaling_factor','order of magnitude of each control (useful for multi-parameter optimization)');
			if ~isnan(self.fos_forward_index),
				fielddisplay(self,'fos_forward_index','index for fos_foward driver of ADOLC');
			end
			if ~isnan(self.fov_forward_indices),
				fielddisplay(self,'fov_forward_indices','indices for fov_foward driver of ADOLC');
			end
		end % }}}
	end
end
