%HYDROLOGYSHREVE class definition
%
%   Usage:
%      hydrologyshreve=hydrologyshreve();

classdef hydrologyshreve
	properties (SetAccess=public) 
		spcwatercolumn     = NaN;
		stabilization      = 0;
		requested_outputs  = {};
	end
	methods
		function self = hydrologyshreve(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(self,varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   hydrologyshreve solution parameters:'));
			fielddisplay(self,'spcwatercolumn','water thickness constraints (NaN means no constraint) [m]');
			fielddisplay(self,'stabilization','artificial diffusivity (default: 1). can be more than 1 to increase diffusivity.');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of stabilization to use 0:nothing 1:artificial_diffusivity
			self.stabilization     = 1;
			self.requested_outputs = {'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if ~ismember('HydrologyShreveAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.ishydrology==0), return; end

			md = checkfield(md,'fieldname','hydrology.spcwatercolumn','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.stabilization','>=',0);
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'Watercolumn','HydrologyWaterVx','HydrologyWaterVy'};
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.hydrology.model','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','spcwatercolumn','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Double');
			outputs = self.requested_outputs;
			pos = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = []; %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.hydrology.spcwatercolumn'],self.spcwatercolumn);
			writejsdouble(fid,[modelname '.hydrology.stabilization'],self.stabilization);

		end % }}}
	end
end

