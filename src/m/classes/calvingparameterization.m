%CALVINGPARAMETERIZATION class definition
%   For test calving laws and coefficients
%   Usage:
%      calvingparameterization=calvingparameterization();

classdef calvingparameterization
	properties (SetAccess=public) 
		min_thickness = 0.;
		use_param = 0;
		theta = 0.;
		alpha = 0;
		xoffset = 0;
		yoffset = 0;
		vel_upperbound = 8000;
		vel_threshold = 6000;
		vel_lowerbound = 0;
	end
	methods
		function self = calvingparameterization(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingparameterization');
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
		end % }}}
		function self = setdefaultparameters(self) % {{{
			%For now we turn this off by setting the threshold to 0
			self.min_thickness = 0.;

			%Parameters for the spatial temporal separation
			%The coefficient follows: gamma= f(x)
			% 0 - f(x) = y_{o} + \alpha (x+x_{o})
			% 1 - f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o}))
			self.use_param = 0;
			% the amplifier
			self.theta = 0;
			% the slope alpha
			self.alpha = 0;
			% offset in x-axis
			self.xoffset = 0;
			% offset in y-axis
			self.yoffset = 0;
			% velocity thresholds to reduce calving rate
			self.vel_upperbound = 6000; % m/a
			self.vel_lowerbound = 0; % m/a
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.min_thickness','>=',0,'NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.use_param','values',[-1, 0, 1, 2, 3, 4, 5]);
			md = checkfield(md,'fieldname','calving.theta','NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.alpha','NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.xoffset','NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.yoffset','NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.vel_lowerbound','NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.vel_threshold','NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.vel_upperbound','NaN',1,'Inf',1,'numel',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving test parameters:'));
			fielddisplay(self,'min_thickness','minimum thickness below which no ice is allowed [m]');
			fielddisplay(self,'use_param','-1 - just use frontal ablation rate, 0 - f(x) = y_{o} + \alpha (x+x_{o}), 1 - f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})), 2 - tanh(thickness), 3 - tanh(normalized vel), 4 - tanh(truncated vel), 5 - linear(truncated vel)');
			fielddisplay(self,'theta','the amplifier');
			fielddisplay(self,'alpha','the slope');
			fielddisplay(self,'xoffset','offset in x-axis');
			fielddisplay(self,'yoffset','offset in y-axis');
			fielddisplay(self,'vel_lowerbound','lowerbound of ice velocity to reduce the calving rate [m/a]');
			fielddisplay(self,'vel_threshold','threshold of ice velocity to reduce the calving rate [m/a]');
			fielddisplay(self,'vel_upperbound','upperbound of ice velocity to reduce the calving rate [m/a]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',9,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','min_thickness','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','use_param','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','theta','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','alpha','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','xoffset','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','yoffset','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','vel_lowerbound','format','Double','scale', 1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','vel_threshold','format','Double','scale', 1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','vel_upperbound','format','Double','scale', 1./yts);
		end % }}}
	end
end
