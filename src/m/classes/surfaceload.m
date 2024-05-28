%SURFACELOAD class definition
%
%   Usage:
%      surfaceload=surfaceload();

classdef surfaceload
	properties (SetAccess=public) 
		icethicknesschange     = [];
		waterheightchange      = [];
		otherchange            = [];
	end
	methods
		function self = surfaceload(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			icethicknesschange=[];
			waterheightchange=[];
			otherchange=[];

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isslc==0), 
				return; 
			end
			if ~isempty(self.icethicknesschange),
				if isa(self.icethicknesschange,'cell'),
					for i=1:length(self.icethicknesschange),
						md = checkfield(md,'field',self.icethicknesschange{i},'NaN',0,'Inf',1,'timeserieslength',1,'Inf',1);
					end
				else
					md = checkfield(md,'field',self.icethicknesschange,'NaN',1,'Inf',1,'timeserieslength',1,'Inf',1);
				end
			end
			if ~isempty(self.waterheightchange),
				md = checkfield(md,'fieldname','solidearth.surfaceload.waterheightchange','timeseries',1,'NaN',1,'Inf',1);
			end
			if ~isempty(self.otherchange),
				md = checkfield(md,'fieldname','solidearth.surfaceload.other','timeseries',1,'NaN',1,'Inf',1);
			end

			%cross check that whereever we have an ice load, the mask is <0 on each vertex:  legacy
			%pos=find(self.deltathickness);
			%maskpos=md.mask.ice_levelset(md.mesh.elements(pos,:)); 
			%[els,vertices]=find(maskpos>0);
			%if length(els),
			%	warning('solidearth checkconsistency fail: there are elements with ice loads where some vertices are not on the ice!');
			%end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surfaceload:'));

			fielddisplay(self,'icethicknesschange','thickness change: ice height equivalent [mIce/yr]');
			fielddisplay(self,'waterheightchange','water height change: water height equivalent [mWater/yr]');
			fielddisplay(self,'otherchange','other loads (sediments) [kg/m^2/yr]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			%deal with ice thickness change: {{{
			if isempty(self.icethicknesschange),
				self.icethicknesschange=zeros(md.mesh.numberofelements+1,1);
			end

			yts=md.constants.yts;

			if isa(self.icethicknesschange,'cell'),
				%transform our cell array of time series into cell array of time series of rates 
				nummodels=length(self.icethicknesschange);
				for i=1:nummodels,
					icethicknesschange=self.icethicknesschange{i}; 
					time=icethicknesschange(end,:);
					dt=diff(time,1,2);
					icethicknesschange_rate=diff(icethicknesschange(1:end-1,:),1,2)./dt;
					self.icethicknesschange{i}=icethicknesschange_rate; 
				end
				WriteData(fid,prefix,'object',self,'fieldname','icethicknesschange','name','md.solidearth.surfaceload.icethicknesschange',...
				'format','MatArray','timeserieslength',md.mesh.numberofelements+1,'yts',yts,'scale',1/yts);
			else
				icethicknesschange=self.icethicknesschange;
				time=icethicknesschange(end,:);
				dt=diff(time,1,2);
				icethicknesschange_rate=diff(icethicknesschange(1:end-1,:),1,2)./dt;
				self.icethicknesschange=icethicknesschange_rate; 

				WriteData(fid,prefix,'object',self,'fieldname','icethicknesschange','name','md.solidearth.surfaceload.icethicknesschange',...
				'format','MatArray','timeserieslength',md.mesh.numberofelements+1,'yts',yts,'scale',1/yts);
			end
			%}}}
			%deal with water height change: {{{
			if isempty(self.waterheightchange),
				self.waterheightchange=zeros(md.mesh.numberofelements+1,1);
			end

			if isa(self.waterheightchange,'cell'),
				%transform our cell array of time series into cell array of time series of rates 
				nummodels=length(self.waterheightchange);
				for i=1:nummodels,
					waterheightchange=self.waterheightchange{i}; 
					time=waterheightchange(end,:);
					dt=diff(time,1,2);
					waterheightchange_rate=diff(waterheightchange(1:end-1,:),1,2)./dt;
					self.waterheightchange{i}=waterheightchange_rate; 
				end
				WriteData(fid,prefix,'object',self,'fieldname','waterheightchange','name','md.solidearth.surfaceload.waterheightchange',...
				'format','MatArray','timeserieslength',md.mesh.numberofelements+1,'yts',yts,'scale',1/yts);
			else
				waterheightchange=self.waterheightchange;
				time=waterheightchange(end,:);
				dt=diff(time,1,2);
				waterheightchange_rate=diff(waterheightchange(1:end-1,:),1,2)./dt;
				self.waterheightchange=waterheightchange_rate; 

				WriteData(fid,prefix,'object',self,'fieldname','waterheightchange','name','md.solidearth.surfaceload.waterheightchange',...
				'format','MatArray','timeserieslength',md.mesh.numberofelements+1,'yts',yts,'scale',1/yts);
			end
			%}}}
			%deal with other: {{{
			if isempty(self.otherchange),
				self.otherchange=zeros(md.mesh.numberofelements+1,1);
			end

			if isa(self.otherchange,'cell'),
				%transform our cell array of time series into cell array of time series of rates 
				nummodels=length(self.otherchange);
				for i=1:nummodels,
					otherchange=self.otherchange{i}; 
					time=otherchange(end,:);
					dt=diff(time,1,2);
					otherchange_rate=diff(otherchange(1:end-1,:),1,2)./dt;
					self.otherchange{i}=otherchange_rate; 
				end
				WriteData(fid,prefix,'object',self,'fieldname','otherchange','name','md.solidearth.surfaceload.otherchange',...
				'format','MatArray','timeserieslength',md.mesh.numberofelements+1,'yts',yts,'scale',1/yts);
			else
				otherchange=self.otherchange;
				time=otherchange(end,:);
				dt=diff(time,1,2);
				otherchange_rate=diff(otherchange(1:end-1,:),1,2)./dt;
				self.otherchange=otherchange_rate; 

				WriteData(fid,prefix,'object',self,'fieldname','otherchange','name','md.solidearth.surfaceload.otherchange',...
				'format','MatArray','timeserieslength',md.mesh.numberofelements+1,'yts',yts,'scale',1/yts);
			end
			%}}}

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.surfaceload.icethicknesschange'],self.icethicknesschange);
			writejs1Darray(fid,[modelname '.surfaceload.waterheightchange'],self.waterheightchange);
			writejs1Darray(fid,[modelname '.surfaceload.otherchange'],self.otherchange);
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
