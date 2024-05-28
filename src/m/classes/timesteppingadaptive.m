%TIMESTEPPINGADAPTIVE Class definition
%
%   Usage:
%      timesteppingadaptive=timesteppingadaptive();

classdef timesteppingadaptive
	properties (SetAccess=public)
		start_time      = 0.;
		final_time      = 0.;
		time_step_min   = 0.;
		time_step_max   = 0.;
		cfl_coefficient = 0.;
		interp_forcing  = 1;
		average_forcing = 0;
		cycle_forcing   = 0;
		coupling_time   = 0.;
	end
	methods
		function self = timesteppingadaptive(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(timesteppingadaptive(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%time between 2 time steps
			self.time_step_min = 0.01;
			self.time_step_max = 10.;

			%final time
			self.final_time=10.*self.time_step_max;

			%default CFL coefficient
			self.cfl_coefficient=0.5;

			%should we interpolate forcing between timesteps?
			self.interp_forcing=1;
			self.average_forcing=0;
			self.cycle_forcing=0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','timestepping.start_time','numel',[1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.final_time','numel',[1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.time_step_min','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.time_step_max','numel',[1],'>=',md.timestepping.time_step_min,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.cfl_coefficient','numel',[1],'>',0,'<=',1);
			md = checkfield(md,'fieldname','timestepping.interp_forcing','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','timestepping.average_forcing','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','timestepping.cycle_forcing','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','timestepping.coupling_time','numel',[1],'>=',md.timestepping.coupling_time,'NaN',1,'Inf',1);
			if self.final_time-self.start_time<0,
				md = checkmessage(md,'timestepping.final_time should be larger than timestepping.start_time');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   timesteppingadaptive parameters:'));

			unit = 'yr';
			fielddisplay(self,'start_time',['simulation starting time [' unit ']']);
			fielddisplay(self,'final_time',['final time to stop the simulation [' unit ']']);
			fielddisplay(self,'time_step_min',['minimum length of time step [' unit ']']);
			fielddisplay(self,'time_step_max',['maximum length of time step [' unit ']']);
			fielddisplay(self,'cfl_coefficient','coefficient applied to cfl condition');
			fielddisplay(self,'interp_forcing','interpolate in time between requested forcing values ? (0 or 1)');
			fielddisplay(self,'average_forcing','average in time if there are several forcing values between steps? (0 or 1, default is 0)');
			fielddisplay(self,'cycle_forcing','cycle through forcing ? (0 or 1)');
			fielddisplay(self,'coupling_time',['coupling time step with ocean model [' unit ']']);

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			scale = md.constants.yts;
			WriteData(fid,prefix,'name','md.timestepping.type','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','start_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','final_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','time_step_min','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','time_step_max','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','cfl_coefficient','format','Double');
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','interp_forcing','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','average_forcing','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','cycle_forcing','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','timestepping','fieldname','coupling_time','format','Double','scale',scale);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejsdouble(fid,[modelname '.timesteppingadaptive.start_time'],self.start_time);
			writejsdouble(fid,[modelname '.timesteppingadaptive.final_time'],self.final_time);
			writejsdouble(fid,[modelname '.timesteppingadaptive.time_step_min'],self.time_step_min);
			writejsdouble(fid,[modelname '.timesteppingadaptive.time_step_max'],self.time_step_max);
			writejsdouble(fid,[modelname '.timesteppingadaptive.cfl_coefficient'],self.cfl_coefficient);
			writejsdouble(fid,[modelname '.timesteppingadaptive.interp_forcing'],self.interp_forcing);
			writejsdouble(fid,[modelname '.timesteppingadaptive.average_forcing'],self.interp_forcing);
			writejsdouble(fid,[modelname '.timesteppingadaptive.cycle_forcing'],self.cycle_forcing);
			writejsdouble(fid,[modelname '.timesteppingadaptive.coupling_time'],self.time_step_max);

		end % }}}
	end
end
