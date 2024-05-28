%BALANCETHICKNESS class definition
%
%   Usage:
%      balancethickness=balancethickness();

classdef balancethickness
	properties (SetAccess=public) 
		spcthickness      = NaN;
		thickening_rate   = NaN;
		stabilization     = 0;

		omega             = NaN;
		slopex            = NaN;
		slopey            = NaN;
	end
	methods
		function self = balancethickness(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of stabilization used
			self.stabilization=1;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if ~strcmp(solution,'BalancethicknessSolution'), return; end

			md = checkfield(md,'fieldname','balancethickness.spcthickness');
			md = checkfield(md,'fieldname','balancethickness.thickening_rate','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','balancethickness.stabilization','size',[1 1],'values',[0 1 2 3]);

			%md = checkfield(md,'fieldname','balancethickness.omega','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1,'>=',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   balance thickness solution parameters:'));

			fielddisplay(self,'spcthickness','thickness constraints (NaN means no constraint) [m]');
			fielddisplay(self,'thickening_rate','ice thickening rate used in the mass conservation (dh/dt) [m/yr]');
			fielddisplay(self,'stabilization','0: None, 1: SU, 2: SSA''s artificial diffusivity, 3:DG');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'fieldname','spcthickness','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','thickening_rate','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Integer');

			WriteData(fid,prefix,'object',self,'fieldname','slopex','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','slopey','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','omega','format','DoubleMat','mattype',1);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.balancethickness.spcthickness'],self.spcthickness);
			writejs1Darray(fid,[modelname '.balancethickness.thickening_rate'],self.thickening_rate);
			writejsdouble(fid,[modelname '.balancethickness.stabilization'],self.stabilization);
			writejs1Darray(fid,[modelname '.balancethickness.omega'],self.omega);

		end % }}}
	end
end
