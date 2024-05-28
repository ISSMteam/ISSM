%HYDROLOGYPISM class definition
%
%   Usage:
%      hydrologypism=hydrologypism();

classdef hydrologypism
	properties (SetAccess=public) 
		drainage_rate   = NaN;
		watercolumn_max = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = hydrologypism(varargin) % {{{
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
			list = {'Watercolumn'};
		end % }}}    
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologyPismAnalysis',analyses)
				return;
			end

			md = checkfield(md,'fieldname','hydrology.drainage_rate','Inf',1,'NaN',1,'>=',0,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','friction.watercolumn_max','NaN',1,'Inf',1,'>',0.,'size',[md.mesh.numberofvertices 1]);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   hydrologypism solution parameters:'));
			fielddisplay(self,'drainage_rate','fixed drainage rate [mm/yr]');
			fielddisplay(self,'watercolumn_max','maximum water column height [m], recommended default: 2 m');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.hydrology.model','data',4,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','drainage_rate','format','DoubleMat','mattype',1,'scale',1./(1000.*yts)); %from mm/yr to m/s
			WriteData(fid,prefix,'class','hydrology','object',self,'fieldname','watercolumn_max','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'data',{'Watercolumn'},'name','md.hydrology.requested_outputs','format','StringArray');
		end % }}}
	end
end

