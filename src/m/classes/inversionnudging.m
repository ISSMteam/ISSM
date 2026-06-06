%INVERSIONVALIDATION class definition
%
%   Usage:
%      inversionnudging=inversionnudging();

classdef inversionnudging
	properties (SetAccess=public) 
		iscontrol          = 0;
		maxiter            = 0;
		max_increment_C    = 0.0;
		relaxation_C       = 0.0;
		tau_C              = 0.0;
		H0_C               = 0.0;
		min_C              = NaN
		max_C              = NaN
		max_increment_melt = 0.0;
		relaxation_melt    = 0.0;
		tau_melt           = 0.0;
		H0_melt            = 0.0;
		min_melt           = NaN
		max_melt           = NaN
		dhdt_obs           = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.min_C = project3d(md,'vector',self.min_C, 'type', 'node');
			self.max_C = project3d(md,'vector',self.max_C, 'type', 'node');
			self.min_melt = project3d(md,'vector',self.min_melt, 'type', 'node');
			self.max_melt = project3d(md,'vector',self.max_melt, 'type', 'node');
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

			%maximum step size in log10 space for C
			self.max_increment_C    = 0.05;
			self.max_increment_melt = 0.5;

			% Cadjustment timescale 100 yrs used in van den Akker et al (2025)
			self.tau_C    = 100.0;
			self.tau_melt = 200.0;

			%thickness error scale (smaller = more sensitive), 100 m used in van den Akker et al (2025)
			self.H0_C    = 100.0;
			self.H0_melt = 50.0;

         %relaxation strength toward C_inv (0 = none, 1 = strong) , 0.5 used in van den Akker et al (2025)
         self.relaxation_C    = 0.5;
			self.relaxation_melt = 0.3;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if(~self.iscontrol) return; end

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.maxiter','numel',1,'>=',0);
			md = checkfield(md,'fieldname','inversion.max_increment_C','numel',1,'>',0);
         md = checkfield(md,'fieldname','inversion.max_increment_melt','numel',1,'>',0);
			md = checkfield(md,'fieldname','inversion.min_C','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','inversion.max_C','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','inversion.min_melt','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','inversion.max_melt','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','inversion.H0_C','numel',1,'>',0);
         md = checkfield(md,'fieldname','inversion.H0_melt','numel',1,'>',0);
         md = checkfield(md,'fieldname','inversion.relaxation_C','numel',1,'>=',0,'<=',1);
         md = checkfield(md,'fieldname','inversion.relaxation_melt','numel',1,'>=',0,'<=',1);
			md = checkfield(md,'fieldname','inversion.dhdt_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Nudging parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'maxiter','maximum number of nudging steps');

			fielddisplay(self,'max_increment_C',   'maximum increase in C per nudging step *in log10 space*');
			fielddisplay(self,'relaxation_C','relaxation strength toward C_inv (0 = none, 1 = strong)');
			fielddisplay(self,'tau_C','adjustment timescale for friction coefficient [yr]');
			fielddisplay(self,'H0_C','thickness error scale for C (smaller = more sensitive) [m]');
			fielddisplay(self,'min_C','absolute minimum acceptable value of C');
			fielddisplay(self,'max_C','absolute maximum acceptable value of C');

			fielddisplay(self,'max_increment_melt','maximum increase in melt per nudging step');
			fielddisplay(self,'relaxation_melt','relaxation strength toward perturbation = 0 (0 = none, 1 = strong)');
         fielddisplay(self,'tau_melt','adjustment timescale for melt perturbation [yr]');
         fielddisplay(self,'H0_melt','thickness error scale for melt perturbation (smaller = more sensitive) [m]');
			fielddisplay(self,'min_melt','absolute minimum acceptable value of melt perturbation [m/yr]');
			fielddisplay(self,'max_melt','absolute maximum acceptable value of melt perturbation [m/yr]');

			fielddisplay(self,'dhdt_obs','observed thickness rate of change [m/yr]');
		end % }}}
		function marshall(self, prefix, md, fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','iscontrol','format','Boolean');
			WriteData(fid,prefix,'name','md.inversion.type','data',5,'format','Integer');
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_increment_C','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_increment_melt','format','Double','scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','min_C','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_C','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','min_melt','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','max_melt','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','tau_C','format','Double','scale',yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','tau_melt','format','Double','scale',yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','H0_C','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','H0_melt','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','relaxation_C','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','relaxation_melt','format','Double');
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','dhdt_obs','format','DoubleMat','mattype',1,'scale',1./yts);
		end % }}}
	end
end
