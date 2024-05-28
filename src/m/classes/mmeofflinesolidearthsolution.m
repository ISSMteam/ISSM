%MMEOFFLINESOLIDEARTHSOLUTION class definition
%
%   Usage:
%      addsol=mmeofflinesolidearthsolution(); where the offline solid earth solutions %                               are based on a multi-model ensemble (ex: Caron et al 2017 statistics) 

classdef mmeofflinesolidearthsolution < offlinesolidearthsolution 
	properties (SetAccess=public) 
		modelid; %index into the multi-model ensemble, each ensemble variable being defined 
	         %in the father class.
	end
	methods
		function self = mmeofflinesolidearthsolution(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.setdefaultparameters@offlinesolidearthsolution();
			self.modelid=0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.solidearth.settings.isgrd==1), 
				error('mmeofflinesolidearthsolution checkconsistency error message: trying to run GRD patterns while supplying an offline solution for those patterns!'); 
			end
			
			seast=length(self.displacementeast);
			snorth=length(self.displacementnorth);
			sup=length(self.displacementup);
			sgeoid=length(self.geoid);

			if (seast-snorth)~=0,
				error('mmeofflinesolidearthsolution checkconsistency error message: displacementeast and displacementnorth should be the same size');
			end

			if (seast-sup)~=0,
				error('mmeofflinesolidearthsolution checkconsistency error message: displacementeast and displacementup should be the same size');
			end

			if (seast-sgeoid)~=0,
				error('mmeofflinesolidearthsolution checkconsistency error message: displacementeast and geoid should be the same size');
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
			disp(sprintf('   external: mmeofflinesolidearth solution:'));
			self.disp@solidearthsolution();
			fielddisplay(self,'modelid','index into the multi-model ensemble, determines which field will be used.');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'data',4,'name','md.solidearth.external.nature','format','Integer'); %code 4 for mmeofflinesolidearthsolution class
			WriteData(fid,prefix,'object',self,'fieldname','modelid','format','Double');
			nummodels=length(self.displacementeast);
			WriteData(fid,prefix,'name','md.solidearth.external.nummodels','data',nummodels,'format','Integer');

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
				barystaticsealevel_rate=diff(barystaticsealevel(1:end-1,:),1,2)./dt;

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
			error('mmeofflinesolidearthsolution error message: not implemented yet');
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
