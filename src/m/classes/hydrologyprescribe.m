%HYDROLOGYPRESCRIBE class definition
%
%   Usage:
%      hydrologyprescribe();

classdef hydrologyprescribe
	properties (SetAccess=public)
		head = NaN;
		requested_outputs = {};
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = hydrologyprescribe(varargin) % {{{
			switch nargin 
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(self,varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'HydrologyHead','EffectivePressure'};
		end % }}}

		function self = setdefaultparameters(self) % {{{
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologyPrescribeAnalysis',analyses)
				return;
			end

			if ~isempty(md.initialization.hydraulic_potential)
				warning('WARN: md.initialization.hydraulic_potential is defined. However, this is not used for "hydrologyprescribe" model.')
			end
			md = checkfield(md,'fieldname','hydrology.head','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('	hydrologyprescribe solution parameters:'));
			disp(sprintf('	This module is to simulate effective pressure Neff using hydraulic head from external subglacial hydrology model.'));
			disp(sprintf('	Neff = rho_i g H - Pw'));
			disp(sprintf('	Pw	= rho_w g (head - z_b))'));
			disp(sprintf('	H: ice thickness [m] / head: hydraulic head [m] / z_b: bedrock elevation'));
			fielddisplay(self,'head','subglacial hydrology water head (m)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.hydrology.model','data',10,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','head','format','DoubleMat','mattype',1);

			outputs = self.requested_outputs;
			pos = find(ismember(outputs,'default'));
			if ~isempty(pos)
				 outputs(pos) = []; %remove 'default' from outputs
				 outputs      = [outputs, defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}
	end
end
