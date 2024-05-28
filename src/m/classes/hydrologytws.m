%HYDROLOGYTWS class definition
%
%   Usage:
%      hydrologytws=hydrologytws();

classdef hydrologytws
	properties (SetAccess=public) 
		spcwatercolumn     = NaN;
		requested_outputs  = {};
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = hydrologytws(varargin) % {{{
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
			list = {''};
		end % }}}

		function self = setdefaultparameters(self) % {{{
			self.requested_outputs = {'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologyTwsAnalysis',analyses)
				return;
			end
			md = checkfield(md,'fieldname','hydrology.spcwatercolumn','Inf',1,'timeseries',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   hydrologytws solution parameters:'));
			fielddisplay(self,'spcwatercolumn','water thickness constraints (NaN means no constraint) [m]');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.hydrology.model','data',6,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','spcwatercolumn','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];  %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.hydrology.spcwatercolumn'],self.spcwatercolumn);

		end % }}}
	end
end

