%CALVINGLEVERMANN class definition
%
%   Usage:
%      calvinglevermann=calvinglevermann();

classdef calvinglevermann
	properties (SetAccess=public) 
		coeff         = NaN;
	end
	methods
		function self = calvinglevermann(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvinglevermann');
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
			self.coeff=project3d(md,'vector',self.coeff,'type','node');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Proportionality coefficient in Levermann model
			self.coeff=2e13;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.coeff','>',0,'size',[md.mesh.numberofvertices 1]);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Levermann parameters:'));
			fielddisplay(self,'coeff','proportionality coefficient in Levermann model');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',3,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','coeff','format','DoubleMat','mattype',1);
		end % }}}
	end
end
