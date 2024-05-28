%FRICTIONCOULOMB2 class definition
%
%   Usage:
%      frictioncoulomb=frictioncoulomb2();

classdef frictioncoulomb2
	properties (SetAccess=public) 
		C                        = NaN;
		m                        = NaN;
		coupling                 = 0;
		effective_pressure       = NaN;
		effective_pressure_limit = 0;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			if isstruct(self)
				disp('Recovering frictioncoulomb2 from older version');
				if isfield(self,'coefficient')
					self.C = self.coefficient;
				end
				self = structtoobj(frictioncoulomb2(),self);
			end
		end% }}}
	end
	methods
		function self = extrude(self,md) % {{{
			self.C=project3d(md,'vector',self.C,'type','node','layer',1);
			self.m=project3d(md,'vector',self.m,'type','element');
			switch self.coupling
				case 0
				case 1
					self.effective_pressure=project3d(md,'vector',self.effective_pressure,'type','node','layer',1);
				case 2
					error('not implemented yet');
				otherwise
					error('not supported yet');		
			end
		end % }}}
		function self = frictioncoulomb(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.effective_pressure_limit = 0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.m','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.coupling','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);			
			switch self.coupling
				case 0
				case 1
					md = checkfield(md,'fieldname','friction.effective_pressure','NaN',1,'Inf',1,'timeseries',1);
				case 2
					error('not implemented yet');
				otherwise
					error('not supported yet');		
			end
		end % }}}
		function disp(self) % {{{
         disp('Coulomb limited sliding law parameters:');
         disp(' ');
         disp('                     C^2 |u_b|^(m-1) * (.5*N)           ');
         disp('      tau_b = - _________________________________   u_b ');
         disp('                (C^(2/m) |u_b| + (0.5*N)^(1/m) )^m      ');
         disp(' ');
			fielddisplay(self,'C','friction coefficient [SI]');
			fielddisplay(self,'m','m exponent (Weertman would be 1/3)');
			fielddisplay(self,'effective_pressure','Effective Pressure for the forcing if not coupled [Pa]');
			fielddisplay(self,'coupling','Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled (not implemented yet)');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');	
			end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'name','md.friction.law','data',13,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','m','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','coupling','format','Integer');
			WriteData(fid,prefix,'object',self,'class','friction','fieldname','effective_pressure_limit','format','Double');
			switch self.coupling
				case 0
				case 1
					WriteData(fid,prefix,'class','friction','object',self,'fieldname','effective_pressure','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				case 2
					error('not implemented yet');
				otherwise
					error('not supported yet');		
			end

		end % }}}
	end
end
