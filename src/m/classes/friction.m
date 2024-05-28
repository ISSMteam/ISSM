%FRICTION class definition
%
%	Usage:
%		friction=friction();

classdef friction
	properties (SetAccess=public) 
		coefficient              = NaN;
		p                        = NaN;
		q                        = NaN;
		coupling                 = 0;
		linearize                = 0;
		effective_pressure       = NaN;
		effective_pressure_limit = 0;
	end
	methods
		function self = extrude(self,md) % {{{
			self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1);
			self.p=project3d(md,'vector',self.p,'type','element');
			self.q=project3d(md,'vector',self.q,'type','element');
			if self.coupling==3 || self.coupling==4
				self.effective_pressure=project3d(md,'vector',self.effective_pressure,'type','node','layer',1);
			end
		end % }}}
		function self = friction(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(friction(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.linearize = 0;
			self.coupling  = 0;
			self.effective_pressure_limit = 0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			if (strcmp(solution,'TransientSolution') &  md.transient.isstressbalance ==0 & md.transient.isthermal == 0), return; end

			md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.linearize','numel',[1],'values',[0:2]);
			md = checkfield(md,'fieldname','friction.coupling','numel',[1],'values',[0:4]);
			md = checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);
         if self.coupling==3
            md = checkfield(md,'fieldname','friction.effective_pressure','NaN',1,'Inf',1,'timeseries',1);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p)'));
			fielddisplay(self,'coefficient','friction coefficient [SI]');
			fielddisplay(self,'p','p exponent');
			fielddisplay(self,'q','q exponent');
			fielddisplay(self,'coupling','Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: use coupled model (not implemented yet)');
			fielddisplay(self,'linearize','0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)');
			fielddisplay(self,'effective_pressure','Effective Pressure for the forcing if not coupled [Pa]');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.friction.law','data',1,'format','Integer');
			if(size(self.coefficient,1)==md.mesh.numberofvertices | size(self.coefficient,1)==md.mesh.numberofvertices+1),
				mattype=1;
				tsl = md.mesh.numberofvertices;
			else
				mattype=2;
				tsl = md.mesh.numberofelements;
			end
			WriteData(fid,prefix,'object',self,'fieldname','coefficient','format','DoubleMat','mattype',mattype,'timeserieslength',tsl+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','p','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'fieldname','q','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','coupling','format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','linearize','format','Integer');
			WriteData(fid,prefix,'object',self,'class','friction','fieldname','effective_pressure_limit','format','Double');
			if self.coupling==3 || self.coupling==4
				WriteData(fid,prefix,'class','friction','object',self,'fieldname','effective_pressure','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.friction.coefficient'],self.coefficient);
			writejs1Darray(fid,[modelname '.friction.p'],self.p);
			writejs1Darray(fid,[modelname '.friction.q'],self.q);
			writejs1Darray(fid,[modelname '.friction.coupling'],self.coupling);
			writejs1Darray(fid,[modelname '.friction.linearize'],self.linearize);
			writejs1Darray(fid,[modelname '.friction.effective_pressure'],self.effective_pressure);
			writejs1Darray(fid,[modelname '.friction.effective_pressure_limit'],self.effective_pressure_limit);
		end % }}}
	end
end
