%ROTATIONAL class definition
%
%   Usage:
%      rotational=rotational();

classdef rotational
	properties (SetAccess=public) 
		equatorialmoi         = 0;
		polarmoi              = 0;
		angularvelocity       = 0;
	end
	methods
		function self = rotational(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%moment of inertia:
			self.equatorialmoi	=8.0077*10^37; % [kg m^2]
			self.polarmoi		=8.0345*10^37; % [kg m^2]

			% mean rotational velocity of earth
			self.angularvelocity=7.2921*10^-5; % [s^-1]
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isslc==0), 
				return; 
			end

			md = checkfield(md,'fieldname','solidearth.rotational.equatorialmoi','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.rotational.polarmoi','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.rotational.angularvelocity','NaN',1,'Inf',1);

		end % }}}
		function list=defaultoutputs(self,md) % {{{
			list = {};
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   rotational parameters:'));

			fielddisplay(self,'equatorialmoi','mean equatorial moment of inertia [kg m^2]');
			fielddisplay(self,'polarmoi','polar moment of inertia [kg m^2]');
			fielddisplay(self,'angularvelocity','mean rotational velocity of earth [per second]'); 

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			
			WriteData(fid,prefix,'object',self,'fieldname','equatorialmoi','name','md.solidearth.rotational.equatorialmoi','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','polarmoi','name','md.solidearth.rotational.polarmoi','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','angularvelocity','name','md.solidearth.rotational.angularvelocity','format','Double');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.rotational.equatorialmoi'],self.equatorial_moi);
			writejsdouble(fid,[modelname '.rotational.polarmoi'],self.polar_moi);
			writejsdouble(fid,[modelname '.rotational.angularvelocity'],self.angular_velocity);
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
