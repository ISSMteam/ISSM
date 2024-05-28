%FRONTALFORCINGS class definition
%
%   Usage:
%      frontalforcings=frontalforcings();

classdef frontalforcings
	properties (SetAccess=public) 
		meltingrate   = NaN;
		ablationrate   = NaN;
	end
	methods
		function self = frontalforcings(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('frontalforcings');
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
		function self = extrude(self,md) % {{{
			self.meltingrate=project3d(md,'vector',self.meltingrate,'type','node');
			self.ablationrate=project3d(md,'vector',self.ablationrate,'type','node');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			meltingrate   = NaN;
			ablationrate   = NaN;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','frontalforcings.meltingrate','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			if ~isnan(md.frontalforcings.ablationrate)
				md = checkfield(md,'fieldname','frontalforcings.ablationrate','Inf',1,'timeseries',1);
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Frontalforcings parameters:'));
			fielddisplay(self,'meltingrate','melting rate at given location [m/a]');
			fielddisplay(self,'ablationrate','frontal ablation rate at given location [m/a], it contains both calving and melting');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.frontalforcings.parameterization','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','meltingrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
			if ~isnan(md.frontalforcings.ablationrate)
				WriteData(fid,prefix,'object',self,'fieldname','ablationrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
			end
		end % }}}
	end
end
