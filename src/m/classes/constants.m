%CONSTANTS class definition
%
%   Usage:
%      constants=constants();

classdef constants
	properties (SetAccess=public) 
		g                      = 0.;
		omega                  = 0.;
		yts                    = 0.;
		referencetemperature   = 0.;
		gravitational_constant = 0.;
	end
	methods
		function self = constants(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%acceleration due to gravity (m/s^2)
			self.g=9.81;

			%Earth's rotation speed 
			self.omega = 7.292*1e-5;

			%converstion from year to seconds
			self.yts=365.*24.*3600.;

			%the reference temperature for enthalpy model (cf Aschwanden)
			self.referencetemperature=223.15;
		
			%gravitational constant: 
			self.gravitational_constant = 6.67259e-11;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','constants.g','>=',0,'size',[1 1]); %We allow 0 for validation tests
			md = checkfield(md,'fieldname','constants.omega','>=',0,'size',[1 1]);
			md = checkfield(md,'fieldname','constants.yts','>',0,'size',[1 1]);
			md = checkfield(md,'fieldname','constants.referencetemperature','size',[1 1]);
			md = checkfield(md,'fieldname','constants.gravitational_constant','size',[1 1]);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   constants parameters:'));

			fielddisplay(self,'g','gravitational acceleration [m/s^2]');
			fielddisplay(self,'omega','angular velocity of Earth [rad/s]');
			fielddisplay(self,'yts','number of seconds in a year [s/yr]');
			fielddisplay(self,'referencetemperature','reference temperature used in the enthalpy model [K]');
			fielddisplay(self,'gravitational_constant','Newtonian constant of gravitation [m^3/kg/s^2]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','g','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','yts','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','referencetemperature','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','gravitational_constant','format','Double');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.constants.g'],self.g);
			writejsdouble(fid,[modelname '.constants.omega'],self.omega);
			writejsdouble(fid,[modelname '.constants.yts'],self.yts);
			writejsdouble(fid,[modelname '.constants.referencetemperature'],self.referencetemperature);
			writejsdouble(fid,[modelname '.constants.gravitational_constant'],self.gravitational_constant);

		end % }}}
	end
end
