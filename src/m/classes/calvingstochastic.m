%CALVINGSTOCHASTIC class definition
%
%   Usage:
%      calvingstochastic=calvingstochastic();

classdef calvingstochastic
	properties (SetAccess=public) 
		k     = 0.0;
		d_max = 0.00
		f     = 0.0;
	end
	methods
		function self = calvingstochastic(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingstochastic');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2)
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
		function self = setdefaultparameters(self) % {{{
			
			self.k     = 22.0;
			self.d_max = 1;
         self.f     = 1e-5;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.k','numel',[1],'>',0., 'NaN', 1);
         md = checkfield(md,'fieldname','calving.d_max','numel',[1],'>',0,'<=',1,'NaN',1);
			md = checkfield(md,'fieldname','calving.f','NaN',1,'Inf',1,'>',0);
		end % }}}
		function disp(self) % {{{
			disp('   Calving Stochastic parameters:');
			disp('      T      = T0 * exp(k*(1-d))');
			disp('      max(P) = 1  - exp(1-∆t*f)');
			fielddisplay(self,'k'    ,'sensitivity of the waiting time [1/m]');
			fielddisplay(self,'d_max','max penetration threshold (between 0 and 1 for full thickness)');
			fielddisplay(self,'f',    'controls calving frequency rate [1/s]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',13,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','k','format','Double');
         WriteData(fid,prefix,'object',self,'fieldname','d_max','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','f','format','Double');
		end % }}}
	end
end
