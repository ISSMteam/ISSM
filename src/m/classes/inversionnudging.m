%INVERSIONVALIDATION class definition
%
%   Usage:
%      inversionnudging=inversionnudging();

classdef inversionnudging
	properties (SetAccess=public) 
		iscontrol  = 0;
		maxiter    = 0;
		max_increment  = 0.0;
		min_parameters = NaN
		max_parameters = NaN
		tau_C      = 0.0;
		H0         = 0.0;
		relaxation = 0.0;
		dhdt_obs   = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.dhdt_obs =project3d(md,'vector',self.dhdt_obs, 'type', 'node');
		end % }}}
		function self = inversionnudging(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(inversionnudging(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			% maximum number of nudging steps
			self.maxiter = 100;

			%maximum step size in log space
			self.max_increment = 0.05;

			% Cadjustment timescale 100 yrs used in van den Akker et al (2025)
			self.tau_C = 100.0;

			%thickness error scale (smaller = more sensitive), 100 m used in van den Akker et al (2025)
			self.H0 = 100.0;

         %relaxation strength toward C_inv (0 = none, 1 = strong) , 0.5 used in van den Akker et al (2025)
         self.relaxation = 0.5;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if(~self.iscontrol) return; end

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.maxiter','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.max_increment','numel',1,'>',0);
			md = checkfield(md,'fieldname','inversion.min_parameters','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','inversion.max_parameters','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','inversion.H0','numel',1,'>',0);
         md = checkfield(md,'fieldname','inversion.relaxation','numel',1,'>=',0,'<=',1);
			md = checkfield(md,'fieldname','inversion.dhdt_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Nudging parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'maxiter','maximum number of nudging steps');
			fielddisplay(self,'max_increment','maximum increase *in log space* per step');
			fielddisplay(self,'min_parameters','absolute minimum acceptable value of the inversed parameter on each vertex');
			fielddisplay(self,'max_parameters','absolute maximum acceptable value of the inversed parameter on each vertex');
			fielddisplay(self,'tau_C','adjustment timescale for friction coefficient [yr]');
			fielddisplay(self,'H0','thickness error scale (smaller = more sensitive) [m]');
         fielddisplay(self,'relaxation','relaxation strength toward C_inv (0 = none, 1 = strong)');
			fielddisplay(self,'dhdt_obs','observed thickness rate of change [m/yr]');
		end % }}}
		function marshall(self, prefix, md, fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','iscontrol','format','Boolean');
			WriteData(fid,prefix,'name','md.inversion.type','data',5,'format','Integer');
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_increment','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','min_parameters','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_parameters','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','tau_C','format','Double','scale',yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','H0','format','Double');
         WriteData(fid,prefix,'object',self,'class','inversion','fieldname','relaxation','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','dhdt_obs','format','DoubleMat','mattype',1,'scale',1./yts);
		end % }}}
	end
end
