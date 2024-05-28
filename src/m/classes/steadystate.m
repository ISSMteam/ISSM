%STEADYSTATE class definition
%
%   Usage:
%      steadystate=steadystate();

classdef steadystate
	properties (SetAccess=public) 
		reltol            = 0;
		maxiter           = 0;
		requested_outputs = {};
	end
	methods
		function self = steadystate(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			%maximum of steady state iterations
			self.maxiter=100;

			%Relative tolerance for the steadystate convertgence
			self.reltol=0.01;

			%default output
			self.requested_outputs={'default'};
		end % }}}
		function list=defaultoutputs(self,md) % {{{

			list =  [md.stressbalance.defaultoutputs(md) md.thermal.defaultoutputs(md)];

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~strcmp(solution,'SteadystateSolution'), return; end;

			if md.timestepping.time_step~=0,
				md = checkmessage(md,['for a steadystate computation, timestepping.time_step must be zero.']);
			end
			md = checkfield(md,'fieldname','steadystate.requested_outputs','stringrow',1);

			if isnan(md.stressbalance.reltol),
				md = checkmessage(md,['for a steadystate computation, stressbalance.reltol (relative convergence criterion) must be defined!']);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   steadystate solution parameters:'));

			fielddisplay(self,'reltol','relative tolerance criterion');
			fielddisplay(self,'maxiter','maximum number of iterations');
			fielddisplay(self,'requested_outputs','additional requested outputs');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','maxiter','format','Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.steadystate.requested_outputs','format','StringArray');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.steadystate.reltol'],self.reltol);
			writejsdouble(fid,[modelname '.steadystate.maxiter'],self.maxiter);
			writejscellstring(fid,[modelname '.steadystate.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
