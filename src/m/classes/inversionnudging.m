%INVERSIONVALIDATION class definition
%
%   Usage:
%      inversionnudging=inversionnudging();

classdef inversionnudging
	properties (SetAccess=public) 
		iscontrol = 0
		tau_C     = 0.0;
		dhdt_obs  = NaN;
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

			% Cadjustment timescale 100 yrs used in van den Akker et al (2025)
			self.tau_C = 100.0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if(~self.iscontrol) return; end

			md = checkfield(md,'fieldname','inversion.iscontrol','values',[0 1]);
			md = checkfield(md,'fieldname','inversion.tau_C','numel',1,'>',0);
			md = checkfield(md,'fieldname','inversion.dhdt_obs','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Nudging parameters:'));
			fielddisplay(self,'iscontrol','is inversion activated?');
			fielddisplay(self,'tau_C','adjustment timescale for friction coefficient [yr]');
			fielddisplay(self,'dhdt_obs','observed thickness rate of change [m/yr]');
		end % }}}
		function marshall(self, prefix, md, fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','iscontrol','format','Boolean');
			WriteData(fid,prefix,'name','md.inversion.type','data',5,'format','Integer');
			if ~self.iscontrol, return; end
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','tau_C','format','Double','scale',yts);
			WriteData(fid,prefix,'object',self,'class','inversion','fieldname','dhdt_obs','format','DoubleMat','mattype',1,'scale',1./yts);
		end % }}}
	end
end
