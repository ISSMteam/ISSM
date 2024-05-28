%FRICTIONREGCOULOMB2 class definition
%
%   Usage:
%      frictionregcoulomb2=frictionregcoulomb2();

classdef frictionregcoulomb2
	properties (SetAccess=public) 
		C                        = NaN;
		K                        = NaN;
		m                        = NaN;
		effective_pressure_limit = 0;
	end
	methods
		function self = frictionregcoulomb2(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictionregcoulomb2(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.C    = project3d(md,'vector',self.C,'type','node');
			self.m    = project3d(md,'vector',self.m,'type','element');
			self.K    = project3d(md,'vector',self.K,'type','node');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.effective_pressure_limit = 0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1,'>=',0.);
			md = checkfield(md,'fieldname','friction.K','NaN',1,'Inf',1,'>',0.);
			md = checkfield(md,'fieldname','friction.m','NaN',1,'Inf',1,'>',0.,'size',[md.mesh.numberofelements,1]);
			md = checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);
		end % }}}
		function disp(self) % {{{
			%See Zoet and Iverson 2020 or Choi et al., 2022 
			disp('Regularized Coulomb friction law 2 parameters:');
			disp('   Regularized Coulomb friction law reads:');
			disp('                       C N |u|^(1/m)         ');
			disp('      tau_b = -  ____________________________');
			disp('                   (|u| + (K*N)^m)^(1/m)     ');
			disp(' ');
			fielddisplay(self,'C','friction coefficient [SI]');
			fielddisplay(self,'m','m exponent');
			fielddisplay(self,'K','(K*N)^m to be velocity controlling plastic limit');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',15,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','K','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','m','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','friction','fieldname','effective_pressure_limit','format','Double');
		end % }}}
	end
end
