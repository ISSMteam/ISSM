%CALVING class definition
%
%   Usage:
%      calving=calving();

classdef calving
	properties (SetAccess=public) 
		calvingrate   = NaN;
	end
	methods
		function self = calving(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calving');
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
			self.calvingrate=project3d(md,'vector',self.calvingrate,'type','node');
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.calvingrate','>=',0,'timeseries',1,'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving parameters:'));
			fielddisplay(self,'calvingrate','calving rate at given location [m/a]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','calvingrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.calving.calvingrate'],self.calvingrate);

		end % }}}
	end
end
