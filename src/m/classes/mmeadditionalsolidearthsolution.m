%MMEADDITIONALSOLIDEARTHSOLUTION class definition
%
%   Usage:
%      addsol=mmeadditionalsolidearthsolution(); where the additional solid earth solutions %                               are based on a multi-model ensemble (ex: Caron et al 2017 statistics) 

classdef mmeadditionalsolidearthsolution < additionalsolidearthsolution 
	properties (SetAccess=public) 
		modelid; %index into the multi-model ensemble, each ensemble variable being defined 
	         %in the father class.
	end
	methods
		function self = mmeadditionalsolidearthsolution(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.setdefaultparameters@additionalsolidearthsolution();
			self.modelid=0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.solidearth.settings.isgrd==0), 
				error('mmeadditionalsolidearthsolution checkconsistency error message: need to run GRD solution if you are supplying a GRD additional pattern solution');
			end

			seast=length(self.displacementeast);
			snorth=length(self.displacementnorth);
			sup=length(self.displacementup);
			sgeoid=length(self.geoid);

			if (seast-snorth)~=0,
				error('mmeadditionalsolidearthsolution checkconsistency error message: displacementeast and displacementnorth should be the same size');
			end

			if (seast-sup)~=0,
				error('mmeadditionalsolidearthsolution checkconsistency error message: displacementeast and displacementup should be the same size');
			end

			if (seast-sgeoid)~=0,
				error('mmeadditionalsolidearthsolution checkconsistency error message: displacementeast and geoid should be the same size');
			end
			
			md = checkfield(md,'field',self.modelid,'NaN',1,'Inf',1,'>=',1,'<=',length(self.displacementeast));

			for i=1:seast,
				md = checkfield(md,'field',self.displacementeast{i},'NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'field',self.displacementup{i},'NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'field',self.displacementnorth{i},'NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'field',self.geoid{i},'NaN',1,'Inf',1,'timeseries',1);
			end


		end % }}}
		function disp(self) % {{{
			disp(sprintf('   external: mmeadditionalsolidearth solution:'));
			self.disp@solidearthsolution();
			fielddisplay(self,'modelid','index into the multi-model ensemble, determines which field will be used.');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'data',3,'name','md.solidearth.external.nature','format','Integer'); %code 3 for mmeadditionalsolidearthsolution class

			nummodels=length(self.displacementeast);
			WriteData(fid,prefix,'name','md.solidearth.external.nummodels','data',nummodels,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','modelid','format','Double');

			%transform our cell array of time series into cell array of time series of rates 
			for i=1:nummodels,
				displacementeast=self.displacementeast{i}; 
				displacementnorth=self.displacementnorth{i}; 
				displacementup=self.displacementup{i}; 
				geoid=self.geoid{i}; 

				time=displacementeast(end,:);
				dt=diff(time,1,2);
				
				displacementeast_rate=diff(displacementeast(1:end-1,:),1,2)./dt;
				displacementnorth_rate=diff(displacementnorth(1:end-1,:),1,2)./dt;
				displacementup_rate=diff(displacementup(1:end-1,:),1,2)./dt;
				geoid_rate=diff(geoid(1:end-1,:),1,2)./dt;

				self.displacementeast{i}=displacementeast_rate; 
				self.displacementnorth{i}=displacementnorth_rate; 
				self.displacementup{i}=displacementup_rate; 
				self.geoid{i}=geoid_rate; 
			end
			
			WriteData(fid,prefix,'object',self,'fieldname','displacementeast','format','MatArray','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1/yts);
			WriteData(fid,prefix,'object',self,'fieldname','displacementup','format','MatArray','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1/yts);
			WriteData(fid,prefix,'object',self,'fieldname','displacementnorth','format','MatArray','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1/yts);
			WriteData(fid,prefix,'object',self,'fieldname','geoid','format','MatArray','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1/yts);

		end % }}}

		function savemodeljs(self,fid,modelname) % {{{
			error('mmeadditionalsolidearthsolution error message: not implemented yet');
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
