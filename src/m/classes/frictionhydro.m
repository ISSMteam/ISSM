%FRICTIONWEERTMAN class definition
%
%   Usage:
%      friction=frictionhydro();

classdef frictionhydro
	properties (SetAccess=public) 
		coupling           = 0;
		q                  = NaN;
		C                  = NaN;
		As                 = NaN;
		effective_pressure = NaN;
		effective_pressure_limit = 0;
	end
	methods
		function self = frictionhydro(varargin) % {{{
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
			md = checkfield(md,'fieldname','friction.coupling','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.C','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.As','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
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
	    end
		end % }}}
		function self = extrude(self,md) % {{{
			self.q=project3d(md,'vector',self.q,'type','element');
			self.C=project3d(md,'vector',self.C,'type','element');
			self.As=project3d(md,'vector',self.As,'type','element');
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
		function disp(self) % {{{
			disp(sprintf('Effective Pressure based friction law described in Gagliardini 2007'));
			fielddisplay(self,'coupling','Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)');
			fielddisplay(self,'q','friction law exponent q>=1');
			fielddisplay(self,'C','friction law max value [SI]');
			fielddisplay(self,'As','Sliding Parameter without cavitation [m Pa^-n s^-1]');
			fielddisplay(self,'effective_pressure','Effective Pressure for the forcing if not coupled [Pa]');
			fielddisplay(self,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.friction.law','data',3,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','coupling','format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','q','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','As','format','DoubleMat','mattype',3);
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
