%FRICTIONPISM class definition
%
%   Usage:
%      frictionpism=frictionpism();

classdef frictionpism
	properties (SetAccess=public) 
		pseudoplasticity_exponent            = 0.;
		threshold_speed                      = 0.;
		delta                                = 0.;
		void_ratio                           = 0.;
		till_friction_angle                  = NaN;
		sediment_compressibility_coefficient = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.till_friction_angle=project3d(md,'vector',self.till_friction_angle,'type','node','layer',1);
         self.sediment_compressibility_coefficient=project3d(md,'vector',self.sediment_compressibility_coefficient,'type','node','layer',1);
		end % }}}
		function self = frictionpism(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(frictionpism(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

         self.pseudoplasticity_exponent = 0.6;
         self.threshold_speed           = 100.;
         self.delta                     = 0.02;
         self.void_ratio                = 0.69;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			if (strcmp(solution,'TransientSolution') &  md.transient.isstressbalance ==0 & md.transient.isthermal == 0), return; end

			md = checkfield(md,'fieldname','friction.pseudoplasticity_exponent','numel',[1],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.threshold_speed','numel',[1],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.delta','numel',[1],'>',0,'<',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.void_ratio','numel',[1],'>',0,'<',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.till_friction_angle','NaN',1,'Inf',1,'<',360.,'>',0.,'size',[md.mesh.numberofvertices 1]); %User should give angle in degrees, Matlab calculates in rad
			md = checkfield(md,'fieldname','friction.sediment_compressibility_coefficient','NaN',1,'Inf',1,'<',1.,'>',0.,'size',[md.mesh.numberofvertices 1]);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters for the PISM friction law (See Aschwanden et al. 2016 for more details)'));
			fielddisplay(self,'pseudoplasticity_exponent','pseudoplasticity exponent [dimensionless]');
			fielddisplay(self,'threshold_speed','threshold speed [m/yr]');
			fielddisplay(self,'delta','lower limit of the effective pressure, expressed as a fraction of overburden pressure [dimensionless]');
			fielddisplay(self,'void_ratio','void ratio at a reference effective pressure [dimensionless]');
			fielddisplay(self,'till_friction_angle','till friction angle [deg], recommended default: 30 deg');
			fielddisplay(self,'sediment_compressibility_coefficient','coefficient of compressibility of the sediment [dimensionless], recommended default: 0.12');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',10,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','pseudoplasticity_exponent','format','Double');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','threshold_speed','format','Double','scale',1./yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','delta','format','Double');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','void_ratio','format','Double');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','till_friction_angle','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','sediment_compressibility_coefficient','format','DoubleMat','mattype',1);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
         error('not implemented yet!');
		end % }}}
	end
end
