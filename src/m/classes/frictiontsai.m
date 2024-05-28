%FRICTIONTSAI class definition
%
%   Usage:
%      frictiontsai=frictiontsai();

classdef frictiontsai
	properties (SetAccess=public) 
		C = NaN;
		f = NaN;
		m = NaN;
		effective_pressure_limit = 0;
	end
	methods
		function self = frictiontsai(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictiontsai(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			md.friction.C    = project3d(md,'vector',md.friction.C,'type','node','layer',1);
			md.friction.f = project3d(md,'vector',md.friction.f,'type','node','layer',1);
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.effective_pressure_limit = 0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1,'>',0.);
			md = checkfield(md,'fieldname','friction.f','timeseries',1,'NaN',1,'Inf',1,'>',0.);
			md = checkfield(md,'fieldname','friction.m','NaN',1,'Inf',1,'>',0.,'size',[md.mesh.numberofelements,1]);
			md = checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);
		end % }}}
		function disp(self) % {{{
			%See Brondex et al. 2017 
			disp('Tsai sliding law parameters:');
			disp('   Tsai''s sliding law reads:');
			disp('     ');
			disp('      tau_b = -  min(C |ub|^(m-1) , f N) u_b   ');
			disp('                                                   ');
			disp(' ');
			fielddisplay(self,'C','friction coefficient [SI]');
			fielddisplay(self,'f','Iken''s bound (typically between 0.17 and 0.84) [SI]');
			fielddisplay(self,'m','m exponent (generally taken as m = 1/n = 1/3)');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',12,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','f','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','m','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','friction','fieldname','effective_pressure_limit','format','Double');
		end % }}}
	end
end
