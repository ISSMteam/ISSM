%FRICTIONREGCOULOMB class definition
%
%   Usage:
%      frictionregcoulomb=frictionregcoulomb();

classdef frictionregcoulomb
	properties (SetAccess=public) 
		C  = NaN;
		u0 = 0.;
		m  = NaN;
	end
	methods
		function self = frictionregcoulomb(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictionregcoulomb(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.C    = project3d(md,'vector',self.C,'type','node');
			self.m    = project3d(md,'vector',self.m,'type','element');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.u0 = 1000;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1,'>=',0.);
			md = checkfield(md,'fieldname','friction.u0','NaN',1,'Inf',1,'>',0.,'numel',1);
			md = checkfield(md,'fieldname','friction.m','NaN',1,'Inf',1,'>',0.,'size',[md.mesh.numberofelements,1]);
		end % }}}
		function disp(self) % {{{
			%See Joughin et al. 2019 (equivalent form by Matt Trevers, poster at AGU 2022) https://agupubs.onlinelåibrary.wiley.com/doi/full/10.1029/2019GL082526
			disp('Regularized Coulomb friction law (Joughin et al., 2019) parameters:');
			disp('   Regularized Coulomb friction law reads:');
			disp('                       C^2 |u|^(1/m)         ');
			disp('      tau_b = -  ____________________________');
			disp('                     (|u|/u0 + 1)^(1/m)      ');
			disp(' ');
			fielddisplay(self,'C','friction coefficient [SI]');
			fielddisplay(self,'m','m exponent (set to m=3 in original paper)');
			fielddisplay(self,'u0','velocity controlling plastic limit');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',14,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','u0','format','Double','scale',1/yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','m','format','DoubleMat','mattype',2);
		end % }}}
	end
end
