%FRICTIONWATERLAYER class definition
%
%   Usage:
%      frictionwaterlayer=frictionwaterlayer();

classdef frictionwaterlayer
	properties (SetAccess=public) 
		coefficient = NaN;
		f           = NaN;
		p           = NaN;
		q           = NaN;
		water_layer = NaN;
	end
	methods
		function self = frictionwaterlayer(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictionwaterlayer(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end

			md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.f','size',[1 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','thermal.spctemperature','Inf',1,'timeseries',1,'>=',0.);

		end % }}}
		function self = extrude(self,md) % {{{
			self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1);
			self.p=project3d(md,'vector',self.p,'type','element');
			self.q=project3d(md,'vector',self.q,'type','element');
			self.water_layer=project3d(md,'vector',self.water_layer,'type','node','layer',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters: tau_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b * 1/f(T)\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*(bed+water_layer), r=q/p and s=1/p)'));
			fielddisplay(self,'coefficient','frictiontemp coefficient [SI]');
			fielddisplay(self,'f','f variable for effective pressure');
			fielddisplay(self,'p','p exponent');
			fielddisplay(self,'q','q exponent');
			fielddisplay(self,'water_layer','water thickness at the base of the ice (m)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'name','md.friction.law','data',5,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','coefficient','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','f','format','Double');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','p','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','q','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','water_layer','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		end % }}}
	end
end
