%FRICTIONSCHOOF class definition
%
%   Usage:
%      frictionschoof=frictionschoof();

classdef frictionschoof
	properties (SetAccess=public) 
		C                        = NaN;
		Cmax                     = NaN;
		m                        = NaN;
		coupling                 = 0;
		effective_pressure       = NaN;
		effective_pressure_limit = 0;
	end
	methods
		function self = frictionschoof(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictionschoof(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.C    = project3d(md,'vector',self.C,'type','node');
			self.Cmax = project3d(md,'vector',self.Cmax,'type','node');
			self.m    = project3d(md,'vector',self.m,'type','element');
			if self.coupling==3 || self.coupling==4
				self.effective_pressure=project3d(md,'vector',self.effective_pressure,'type','node','layer',1);
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

         self.coupling = 0;
			self.effective_pressure_limit = 0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1,'>=',0.);
			md = checkfield(md,'fieldname','friction.Cmax','timeseries',1,'NaN',1,'Inf',1,'>',0.);
			md = checkfield(md,'fieldname','friction.m','NaN',1,'Inf',1,'>',0.,'size',[md.mesh.numberofelements,1]);
			md = checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);
         md = checkfield(md,'fieldname','friction.coupling','numel',[1],'values',[0:4]);
         if self.coupling==3
            md = checkfield(md,'fieldname','friction.effective_pressure','NaN',1,'Inf',1,'timeseries',1);
         end
		end % }}}
		function disp(self) % {{{
			%See Brondex et al. 2019 https://www.the-cryosphere.net/13/177/2019/
			disp('Schoof sliding law parameters:');
			disp('   Schoof''s sliding law reads:');
			disp('                         C^2 |u_b|^(m-1)                ');
			disp('      tau_b = - _____________________________   u_b   ');
			disp('               (1+(C^2/(Cmax N))^1/m |u_b| )^m          ');
			disp(' ');
			fielddisplay(self,'C','friction coefficient [SI]');
			fielddisplay(self,'Cmax','Iken''s bound (typically between 0.17 and 0.84) [SI]');
			fielddisplay(self,'m','m exponent (generally taken as m = 1/n = 1/3)');
			fielddisplay(self,'effective_pressure','Effective Pressure for the forcing if not coupled [Pa]');
         fielddisplay(self,'coupling','Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: use coupled model (not implemented yet)');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',11,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','Cmax','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','m','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','friction','fieldname','effective_pressure_limit','format','Double');
         WriteData(fid,prefix,'class','friction','object',self,'fieldname','coupling','format','Integer');
         if self.coupling==3 || self.coupling==4
            WriteData(fid,prefix,'class','friction','object',self,'fieldname','effective_pressure','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
         end
		end % }}}
	end
end
