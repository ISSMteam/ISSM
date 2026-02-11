%SOLIDEARTHSOLUTION class definition
%
%   Usage:
%      solidearthsolution=solidearthsolution();

classdef solidearthsolution
	properties (SetAccess=public) 
		displacementeast = [];
		displacementnorth =[];
		displacementup=[];
		geoid=[];
	end
	methods
		function self = solidearthsolution(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('%s','         units for time series is (yr)'));
			fielddisplay(self,'displacementeast','solid-Earth Eastwards bedrock displacement series (m)');
			fielddisplay(self,'displacementnorth','solid-Earth Northwards bedrock displacement time series (m)');
			fielddisplay(self,'displacementup','solid-Earth bedrock uplift time series (m)');
			fielddisplay(self,'geoid','solid-Earth geoid time series (m)');

		end % }}}
		function self = setdefaultparameters(self) % {{{
	
			self.displacementeast = [];
			self.displacementnorth =[];
			self.displacementup=[];
			self.geoid=[];
	
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','solidearth.external.displacementeast','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','solidearth.external.displacementnorth','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','solidearth.external.displacementup','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','solidearth.external.geoid','Inf',1,'timeseries',1);

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			%transform our time series into time series rates 
			if size(self.displacementeast,2)==1
				disp('External solidearthsolution warning: only one time step provided, assuming the values are rates per year');
				displacementeast_rate=[self.displacementeast;0];
				displacementnorth_rate=[self.displacementnorth;0];
				displacementup_rate=[self.displacementup;0];
				geoid_rate=[self.geoid;0];
			else
				time=self.displacementeast(end,:);
				dt=diff(time,1,2);
				displacementeast_rate=diff(self.displacementeast(1:end-1,:),1,2)./dt;
				displacementeast_rate(end+1,:)=time(2:end);
				displacementnorth_rate=diff(self.displacementnorth(1:end-1,:),1,2)./dt;
				displacementnorth_rate(end+1,:)=time(2:end);
				displacementup_rate=diff(self.displacementup(1:end-1,:),1,2)./dt;
				displacementup_rate(end+1,:)=time(2:end);
				geoid_rate=diff(self.geoid(1:end-1,:),1,2)./dt;
				geoid_rate(end+1,:)=time(2:end);
			end

			WriteData(fid, prefix, 'name', 'md.solidearth.external.nature', 'data', 0, 'format', 'Integer');
			WriteData(fid,prefix,'object',self,'fieldname','displacementeast','data',displacementeast_rate,'format','DoubleMat','name', 'md.solidearth.external.displacementeast','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','displacementup','data',displacementup_rate,'format','DoubleMat','name', 'md.solidearth.external.displacementup','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','displacementnorth','data',displacementnorth_rate,'format','DoubleMat','name', 'md.solidearth.external.displacementnorth','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','geoid','data',geoid_rate,'format','DoubleMat','name', 'md.solidearth.external.geoid','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			writejs1Darray(fid,[modelname '.solidearth.external.displacementeast'],self.displacementeast);
			writejs1Darray(fid,[modelname '.solidearth.external.displacementnorth'],self.displacementnorth);
			writejs1Darray(fid,[modelname '.solidearth.external.displacementup'],self.displacementup);
			writejs1Darray(fid,[modelname '.solidearth.external.geoid'],self.geoid);
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
