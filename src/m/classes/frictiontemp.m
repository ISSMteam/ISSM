%FRICTIONTEMP class definition
%
%   Usage:
%      frictiontemp=frictiontemp();

classdef frictiontemp
	properties (SetAccess=public) 
		gamma       = 0;
		coefficient = NaN;
		p           = NaN;
		q           = NaN;
		coupling    = 0;
		effective_pressure = NaN;
		effective_pressure_limit = 0;
	end
	methods
		function self = extrude(self,md) % {{{
			self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1);
			self.p=project3d(md,'vector',self.p,'type','element');
			self.q=project3d(md,'vector',self.q,'type','element');
			switch self.coupling
				case 0
				case 1
				case 2
				case 3
					self.effective_pressure=project3d(md,'vector',self.effective_pressure,'type','node','layer',1);
				case 4
					error('not implemented yet');
				otherwise
					error('not supported yet');
			end
		end % }}}
		function self = frictiontemp(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictiontemp(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%By default gamma = 1;
			self.gamma = 1;
			self.coupling = 0;
			self.effective_pressure_limit = 0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end

			md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.gamma','NaN',2,'Inf',1,'numel',1,'>',0.);
			md = checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);

			%Check that temperature is provided
			md = checkfield(md,'fieldname','initialization.temperature','NaN',1,'Inf',1,'size','universal');
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters: tau_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b * 1/f(T)\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p)'));
			fielddisplay(self,'gamma','submelt sliding parameter f(T) = exp((T-Tpmp)/gamma)');
			fielddisplay(self,'coefficient','frictiontemp coefficient [SI]');
			fielddisplay(self,'p','p exponent');
			fielddisplay(self,'q','q exponent');
			fielddisplay(self,'effective_pressure','Effective Pressure for the forcing if not coupled [Pa]');
			fielddisplay(self,'coupling','Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'name','md.friction.law','data',4,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','gamma','format','Double');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','coefficient','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','p','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','q','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','coupling','format','Integer');
			WriteData(fid,prefix,'object',self,'class','friction','fieldname','effective_pressure_limit','format','Double');
			switch self.coupling
				case 0
				case 1
				case 2
				case 3
					WriteData(fid,prefix,'class','friction','object',self,'fieldname','effective_pressure','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				case 4
					error('not implemented yet');
				otherwise
					error('not supported yet');
			end
		end % }}}
	end
end
