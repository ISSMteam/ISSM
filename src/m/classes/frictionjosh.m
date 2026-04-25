%FRICTIONTEMP class definition
%
%   Usage:
%      frictionjosh=frictionjosh();

classdef frictionjosh
	properties (SetAccess=public) 
		coefficient                   = NaN;
		pressure_adjusted_temperature = NaN;
		gamma                         = 0.0;
		effective_pressure_limit      = 0.0;
		coefficient_max               = 0.0;
	end
	methods
		function self = extrude(self,md) % {{{
			self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1);
			self.pressure_adjusted_temperature=project3d(md,'vector',self.pressure_adjusted_temperature,'type','node','layer',1);
		end % }}}
		function self = frictionjosh(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictionjosh(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Default gamma: 1
			self.gamma = 1.;

			% Default max friction coefficient: 300
			self.coefficient_max = 300.0;

			%Default 0
			self.effective_pressure_limit = 0.0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end

			md = checkfield(md,'fieldname','friction.coefficient','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.pressure_adjusted_temperature','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.gamma','numel',1,'NaN',1,'Inf',1,'>',0.);
			md = checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','friction.coefficient_max','numel',1,'NaN',1,'Inf',1,'>',0.);

			%Check that temperature is provided
			md = checkfield(md,'fieldname','initialization.temperature','NaN',1,'Inf',1,'size','universal');
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters: tau_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b * 1/f(T)\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p)'));
			fielddisplay(self,'coefficient','friction coefficient [SI]');
			fielddisplay(self,'pressure_adjusted_temperature','friction pressure_adjusted_temperature (T - Tpmp) [K]');
			fielddisplay(self,'gamma','(T - Tpmp)/gamma [K]');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
			fielddisplay(self,'coefficient_max', 'effective friction C = min(coefficient_max, sqrt(exp(T_b(modern) - T_b(t))/gamma) * coefficient)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix, 'name', 'md.friction.law', 'data',9, 'format', 'Integer');
			WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'coefficient', 'format', 'DoubleMat', 'mattype',1, 'timeserieslength',md.mesh.numberofvertices+1, 'yts',md.constants.yts);
			WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'pressure_adjusted_temperature', 'format', 'DoubleMat', 'mattype',1, 'timeserieslength',md.mesh.numberofvertices+1, 'yts',md.constants.yts);
			WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'gamma', 'format', 'Double');
			WriteData(fid,prefix, 'object',self, 'class', 'friction', 'fieldname', 'effective_pressure_limit', 'format', 'Double');
			WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'coefficient_max', 'format', 'Double');
		end % }}}
	end
end
