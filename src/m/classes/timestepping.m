%TIMESTEPPING Class definition
%
%   Usage:
%      timestepping=timestepping();

classdef timestepping
	properties (SetAccess=public)
		start_time      = 0.;
		final_time      = 0.;
		time_step       = 0.;
		interp_forcing  = 1;
		average_forcing  = 0;
		cycle_forcing   = 0;
		coupling_time   = 0.;
	end
	methods
		function self = timestepping(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(timestepping(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   timestepping parameters:'));

			unit = 'yr';
			fielddisplay(self,'start_time',['simulation starting time [' unit ']']);
			fielddisplay(self,'final_time',['final time to stop the simulation [' unit ']']);
			fielddisplay(self,'time_step',['length of time steps [' unit ']']);
			fielddisplay(self,'interp_forcing','interpolate in time between requested forcing values? (0 or 1)');
			fielddisplay(self,'average_forcing','average in time if there are several forcing values between steps? (0 or 1, default is 0)');
			fielddisplay(self,'cycle_forcing','cycle through forcing? (0 or 1)');
			fielddisplay(self,'coupling_time',['length of coupling time step with ocean model [' unit ']']);

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%time between 2 time steps
			self.time_step=1./2.;

			%final time
			self.final_time=10.*self.time_step;

			%should we interpolate forcing between timesteps?
			self.interp_forcing=1;
			self.average_forcing=0;
			self.cycle_forcing=0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','timestepping.start_time','numel',[1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.final_time','numel',[1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.time_step','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.interp_forcing','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','timestepping.average_forcing','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','timestepping.cycle_forcing','numel',[1],'values',[0 1]);
			if self.final_time-self.start_time<0,
				md = checkmessage(md,'timestepping.final_time should be larger than timestepping.start_time');
			end
			if strcmp(solution,'TransientSolution'),
				md = checkfield(md,'fieldname','timestepping.time_step','numel',[1],'>',0,'NaN',1,'Inf',1);
			end
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			scale = md.constants.yts;
			WriteData(fid,prefix,'name','md.timestepping.type','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','start_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'fieldname','final_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'fieldname','time_step','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'fieldname','interp_forcing','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','average_forcing','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','cycle_forcing','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','coupling_time','format','Double','scale',scale);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejsdouble(fid,[modelname '.timestepping.start_time'],self.start_time);
			writejsdouble(fid,[modelname '.timestepping.final_time'],self.final_time);
			writejsdouble(fid,[modelname '.timestepping.time_step'],self.time_step);
			writejsdouble(fid,[modelname '.timestepping.interp_forcing'],self.interp_forcing);
			writejsdouble(fid,[modelname '.timestepping.average_forcing'],self.interp_forcing);
			writejsdouble(fid,[modelname '.timestepping.cycle_forcing'],self.cycle_forcing);

		end % }}}
	end
end
