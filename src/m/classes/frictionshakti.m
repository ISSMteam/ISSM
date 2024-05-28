%FRICTIONSHAKTI class definition
%
%   Usage:
%      friction=frictionshakti();

classdef frictionshakti
	properties (SetAccess=public) 
		coefficient = NaN;
	end
	methods
		function self = frictionshakti(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictionshakti(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function self = extrude(self,md) % {{{
			self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1);
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters: Sigma_b = coefficient^2 * Neff * u_b\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*(head-b))'));
			fielddisplay(self,'coefficient','friction coefficient [SI]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.friction.law','data',8,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','coefficient','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);

		end % }}}
	end
end
