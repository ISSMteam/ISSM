%FRICTIONWEERTMAN class definition
%
%   Usage:
%      frictionweertman=frictionweertman();

classdef frictionweertman
	properties (SetAccess=public) 
		C = NaN;
		m = NaN;
		linearize  = 0;
	end
   methods (Static)
      function self = loadobj(self) % {{{
			disp('Warning: the Weertman friciton law is updated to Sigma_b = C^2*|u_b|^(1/m-1)*u_b, since 2020-08-10');
		end
	end %}}}
	methods
		function self = frictionweertman(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			disp('-------------- file: frictionweertman.m line: 27'); 
			self.C=project3d(md,'vector',self.C,'type','node','layer',1);
			self.m=project3d(md,'vector',self.m,'type','element');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.linearize = 0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.m','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.linearize','numel',[1],'values',[0:2]);
		end % }}}
		function disp(self) % {{{
			disp('Weertman sliding law parameters:');
			disp('   Weertman''s sliding law reads:');
			disp('      v_b = C_w * Sigma_b^m');
			disp('   In ISSM, this law is rewritten as:');
			disp('      Sigma_b = C^2 * |u_b|^(1/m-1) * u_b');
			disp('   where C_w=C^(-2m)');
			disp(' ');
			fielddisplay(self,'C','friction coefficient [SI]');
			fielddisplay(self,'m','m exponent');
			fielddisplay(self,'linearize','0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',2,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','m','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','linearize','format','Integer');
		end % }}}
	end
end
