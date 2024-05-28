%GLACIERMIP class definition
%
%   Usage:
%      gmip = glaciermip('ncfile',glaciermipncfile,'version',1)
%
%

classdef glaciermip < handle
	properties (SetAccess=public) %Model fields
		ncfile  = '';
		version = NaN;

		%version 2:
		area =[];
		mass =[];
		region =[];
		time =[];
		glaciermodel =[];
		climatemodel =[];
		scenario =[];

		%from version 1:
		run=[];
		forcingmodel=[];
		realization=[];
		volume=[];
	end
	methods
		function self = glaciermip(varargin) % {{{

			if nargin==0, 
				self=setdefaultparameters(self);
			else 
				self=setdefaultparameters(self);

				options=pairoptions(varargin{:});

				self.ncfile=getfieldvalue(options,'ncfile');
				self.version=getfieldvalue(options,'version');


				%read variables:
				if self.version==1,
					self.region=ncread(self.ncfile,'region');
					self.time=ncread(self.ncfile,'time');
					self.run=ncread(self.ncfile,'run');
					self.glaciermodel=ncread(self.ncfile,'glaciermodel');
					self.forcingmodel=ncread(self.ncfile,'forcingmodel');
					self.scenario=ncread(self.ncfile,'scenario');
					self.realization=ncread(self.ncfile,'realization');
					self.area=ncread(self.ncfile,'area');
					self.volume=ncread(self.ncfile,'volume');
				elseif self.version==2,
					self.area=ncread(self.ncfile,'Area');
					self.mass=ncread(self.ncfile,'Mass');
					self.region=ncread(self.ncfile,'Region');
					self.time=ncread(self.ncfile,'Time');
					self.glaciermodel=ncread(self.ncfile,'Glacier_Model');
					self.climatemodel=ncread(self.ncfile,'Climate_Model');
					self.scenario=ncread(self.ncfile,'Scenario');
					%mass(region,time,climatemodel,glaciermodel,scenario)
				else 
					error(sprintf('glaciermipfile constructor error message: version %i for MIP not supported!'),self.version);
				end
			end
		end
		%}}}
		function inv = setdefaultparameters(inv) % {{{
		end
		%}}}
		function t=gettime(self,varargin) % {{{
			
			options=pairoptions(varargin{:});
			t=self.time;

		end % }}}
		function masses=getmass(self,varargin) % {{{
			
			options=pairoptions(varargin{:});

			rg=getfieldvalue(options,'region',1:length(self.region));
			cm=getfieldvalue(options,'climatemodel',1:length(self.climatemodel));
			gm=getfieldvalue(options,'glaciermodel',1:length(self.glaciermodel));
			sc=getfieldvalue(options,'scenario',1:length(self.scenario));
			zerostonan=getfieldvalue(options,'zerostonan',0);
			unit=getfieldvalue(options,'unit','Gt');
			sumregion=getfieldvalue(options,'sumregion',0);

			if self.version==1,
				error(sprintf('getmass not supported yet for Glacier MIP version %i',self.version));
			end

			%serialize: 
			if sumregion,
				for i=rg,
					masses_regioni=[];
					for j=cm,
						for k=gm,
							for l=sc,
								masses_regioni=[masses_regioni; squeeze(self.mass(i,:,j,k,l))];
							end
						end
					end
					if i==1, 
						masses=masses_regioni;
					else
						masses=masses+masses_regioni;
					end
				end
			else
				masses=[];
				for i=rg,
					for j=cm,
						for k=gm,
							for l=sc,
								masses=[masses; squeeze(self.mass(i,:,j,k,l))];
							end
						end
					end
				end
			end
			if zerostonan,
				masses(find(masses==0))=NaN;
			end
			if strcmpi(unit,'mmSLE'),
				masses=masses/sletogt();
			end


		end % }}}
		function massrates=getmassrates(self,varargin) % {{{
		
			options=pairoptions(varargin{:});

			unit=getfieldvalue(options,'unit','Gt/yr');

			%get mass first: 
			masses=self.getmass(varargin{:});

			%compute mass rates: 
			dt=diff(self.time);
			dm=diff(masses,1,2);

			massrates=dm;
			for i=1:size(massrates,1),
				massrates(i,:)= massrates(i,:)./dt';
			end

			if strcmpi(unit,'mmSLE/yr'),
				massrates=massrates/sletogt();
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Glacier MIP (version %i):',self.version)); 

			if self.version==1,
				fielddisplay(self,'ncfile','netcdf file for GlacierMIP results');
				fielddisplay(self,'time','time scale in yr');
				fielddisplay(self,'region','region');
				fielddisplay(self,'run','run');
				fielddisplay(self,'glaciermodel','glacier model');
				fielddisplay(self,'forcingmodel','forcing model');
				fielddisplay(self,'scenario','scenario');
				fielddisplay(self,'realization','realization');
				fielddisplay(self,'area','area');
				fielddisplay(self,'volume','volume');
			end 
			if self.version==2,
				fielddisplay(self,'ncfile','netcdf file for GlacierMIP results');
				fielddisplay(self,'time','time');
				fielddisplay(self,'region','region');
				fielddisplay(self,'glaciermodel','glaciermodel');
				fielddisplay(self,'climatemodel','climatemodel');
				fielddisplay(self,'scenario','scenario');
				fielddisplay(self,'area','area');
				fielddisplay(self,'mass','mass');
			end
		end % }}}
		function part=partition(self,md,rgi2mesh,part,value) % {{{

			for i=1:size(rgi2mesh,2),
				dh=rgi2mesh(:,i);
				pos=find(dh);
				part(pos)=value;
			end

		end  % }}}
	end 
end
