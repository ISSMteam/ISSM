%LEVELSET class definition
%
%   Usage:
%      levelset=levelset();

classdef levelset
	properties (SetAccess=public) 
		stabilization		= 0;
		spclevelset			= NaN;
		reinit_frequency	= 10;
		kill_icebergs		= 0;
		migration_max		= 0.;
		fe					= 'P1';
	end
	methods
		function self = levelset(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('levelset');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Level-set parameters:'));
			fielddisplay(self,'stabilization','0: No Stabilization - No stabilization techniques applied.');
			disp('						     1: Artificial Diffusivity - Most stable, but least accurate.');
			disp('						     2: Streamline Upwinding');
			disp('						     5: SUPG - Most accurate, but may be unstable in some applications.');
			fielddisplay(self,'spclevelset','Levelset constraints (NaN means no constraint)');
			fielddisplay(self,'reinit_frequency','Amount of time steps after which the levelset function in re-initialized');
			fielddisplay(self,'kill_icebergs','remove floating icebergs to prevent rigid body motions (1: true, 0: false)');
			fielddisplay(self,'migration_max','maximum allowed migration rate (m/a)');
			fielddisplay(self,'fe','Finite Element type: ''P1'' (default), or ''P2''');
		end % }}}
		function self = extrude(self,md) % {{{

			self.spclevelset=project3d(md,'vector',self.spclevelset,'type','node');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%stabilization = 1 by default
			self.stabilization    = 1;
			self.reinit_frequency = 10;
			self.kill_icebergs    = 1;
			self.migration_max    = 1e12; % No need for general cases, unless specified

			%Linear elements by default
			self.fe='P1';

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','levelset.spclevelset','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','levelset.stabilization','values',[0 1 2 5 6]);
			md = checkfield(md,'fieldname','levelset.kill_icebergs','numel',1,'values',[0 1]);
			md = checkfield(md,'fieldname','levelset.migration_max','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','levelset.fe','values',{'P1','P2'});
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','spclevelset','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','reinit_frequency','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','kill_icebergs','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','migration_max','format','Double','scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','fe','format','String');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			writejsdouble(fid,[modelname '.levelset.stabilization'],self.stabilization);
			writejs1Darray(fid,[modelname '.levelset.spclevelset'],self.spclevelset);
			writejs1Darray(fid,[modelname '.levelset.reinit_frequency'],self.reinit_frequency);
			writejsdouble(fid,[modelname '.levelset.kill_icebergs'],self.kill_icebergs);
			writejsdouble(fid,[modelname '.levelset.migration_max'],self.migration_max);
		end % }}}
	end
end

