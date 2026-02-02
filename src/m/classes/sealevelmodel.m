%SEALEVELMODEL class definition
%
%   Usage:
%      slm = sealevelmodel(varargin)
%
%      where varargin is a variable list of options:
%
%   Example:
%      slm = sealevel('icecap',md_greenland,'icecap',md_antarctica,'earth',md_earth);

classdef sealevelmodel < handle
	properties (SetAccess=public) %Model fields
		% {{{
		icecaps          = {}; % list of land/ice models, name should  be change longer term.
		earth            = 0;  % model for the whole earth
		basins           = {}; % list of basins, matching icecaps, where shapefile info is held
		cluster          = 0;
		miscellaneous    = 0;
		settings         = 0;
		private          = 0;
		mergedcaps       = 0;
		transitions      = {};
		eltransitions    = {};
		planet           = '';
		%}}}
	end
	methods
		function slm = sealevelmodel(varargin) % {{{
			slm=setdefaultparameters(slm);

			if nargin==1,

				options=pairoptions(varargin{:});

				%recover all the icecap models:
				slm.icecaps=getfieldvalues(options,'ice_cap',{});
				
				%recover the earth model:
				slm.earth = getfieldvalue(options,'earth',0);

				%set planet type:
				slm.planet=getfieldvalue(options,'planet','earth');

			end
		end
		%}}}
		function checkconsistency(slm,solutiontype) % {{{

			%is the coupler turned on?
			%for i=1:length(slm.icecaps),
			%	if slm.icecaps{i}.transient.iscoupler==0,
			%		warning(sprintf('sealevelmodel.m::checkconsistency: icecap model %s should have the transient coupler option turned on!',slm.icecaps{i}.miscellaneous.name));
			%	end
			%end
				
			%if slm.earth.transient.iscoupler==0,
			%	warning('sealevelmodel.m::checkconsistency: earth model should have the transient coupler option turned on!');
			%end

			%check that the transition vectors have the right size:
			if slm.earth.mesh.numberofvertices ~= length(slm.earth.solidearth.transfercount)
				error('earth.solidearth.transfercount should be of size earth.mesh.numberofvertices') 
			end

			%check that slc is on
			if slm.earth.transient.isslc ==0
				error('earth.transient.isslc should be turned on') 
			end

			for i=1:length(slm.icecaps),

				md= slm.icecaps{i};

				%check that the transition vectors have the right size:
				if md.mesh.numberofvertices ~= length(slm.earth.solidearth.transitions{i}),
					error(['issue with size of transition vector for ice cap: ' num2str(i) ' name: ' md.miscellaneous.name]);
				end

				%check that runfrequency is the same everywhere:
				if md.solidearth.settings.runfrequency~=slm.earth.solidearth.settings.runfrequency
					error(['icecap model ' md.miscellaneous.name ' should have the same run frequency as earth!']);
				end

				%make sure steric_rate is the same everywhere:
				if ~isempty(find(md.dsl.sea_surface_height_above_geoid - slm.earth.dsl.sea_surface_height_above_geoid(slm.transitions{i})))
					error(['steric rate on ice cap ' md.miscellaneous.name ' is not the same as for the earth']);
				end

				%make sure grd is the same everywhere:
				if md.solidearth.settings.isgrd~=slm.earth.solidearth.settings.isgrd
					error(['isgrd on ice cap ' md.miscellaneous.name ' is not the same as for the earth']);
				end

				%make sure that there is no solid earth external forcing on the basins:
				if ~isempty(md.solidearth.external)
					error('cannot run external forcings on an ice sheet when running a coupling earth/ice sheet model');
				end

				%make sure that there is no solid earth external forcing on the basins:
				if ~isempty(md.solidearth.external)
					error('cannot run external forcings on an ice sheet when running a coupling earth/ice sheet model');
				end

				%make sure that we have the right grd model for computing out sealevel patterns:
				if md.solidearth.settings.grdmodel~=0
					error(['ice sheets do not run GRD module, specify solidearth.settings.grdmodel=0 on ice cap ',num2str(i)]);
				end

				%make sure sea level change solver is on
				if md.transient.isslc==0
					error(['isslc on ice cap ' md.miscellaneous.name ' is not turned off']);
				end
			end

		end
		%}}}
		function slm = setdefaultparameters(slm) % {{{

			%initialize subclasses
			slm.icecaps           = {};
			slm.earth             = {};
			slm.miscellaneous     = miscellaneous();
			slm.settings          = issmsettings();
			slm.private           = private();
			slm.cluster           = generic();
			slm.transitions       = {};
			slm.eltransitions     = {};
			slm.planet            = 'earth';
		end
		%}}}
		function disp(self) % {{{
			disp(sprintf('%19s: %-22s -- %s','icecaps'         ,['[' num2str(length(self.icecaps)) 'x1 ' class(self.icecaps) ']'],'ice caps'));
			disp(sprintf('%19s: %-22s -- %s','earth'           ,['[1x1 ' class(self.earth) ']'],'earth'));
			disp(sprintf('%19s: %-22s -- %s','settings'        ,['[1x1 ' class(self.settings) ']'],'settings properties'));
			disp(sprintf('%19s: %-22s -- %s','cluster'         ,['[1x1 ' class(self.cluster) ']'],'cluster parameters (number of cpus...)'));
			disp(sprintf('%19s: %-22s -- %s','miscellaneous'   ,['[1x1 ' class(self.miscellaneous) ']'],'miscellaneous fields'));
		end % }}}
		function self=mergeresults(self) % {{{
			champs=fieldnames(self.icecaps{1}.results.TransientSolution);
			for i=1:length(self.mergedcaps)/2,
				md=self.mergedcaps{2*(i-1)+1}; trans=self.mergedcaps{2*(i-1)+2};
				%icecaps=self.icecaps(self.range{2*(i-1)+2});
				for j=1:length(self.icecaps{1}.results.TransientSolution),
					for k=1:length(champs),
						if strcmpi(class(icecaps{1}.results.TransientSolution(j).(champs{k})),'double'),
							%vertex or element?
							if length(icecaps{1}.results.TransientSolution(j).(champs{k}))==icecaps{1}.mesh.numberofvertices,
								md.results.TransientSolution(j).(champs{k})=zeros(md.mesh.numberofvertices,1);
								for l=1:length(trans),
									resultcap=icecaps{l}.results.TransientSolution(j).(champs{k});
									md.results.TransientSolution(j).(champs{k})(trans{l})=resultcap;
								end
							else
								if strcmpi(champs{k},'IceVolume') | strcmpi(champs{k},'IceVolumeAboveFloatation') ,
									md.results.TransientSolution(j).(champs{k})=0;
									for l=1:length(trans),
										resultcap=icecaps{l}.results.TransientSolution(j).(champs{k});
										md.results.TransientSolution(j).(champs{k})= md.results.TransientSolution(j).(champs{k})+resultcap;
									end
								elseif strcmpi(champs{k},'time'),
									md.results.TransientSolution(j).(champs{k})= icecaps{1}.results.TransientSolution(j).(champs{k});
								else
									continue;
								end
							end
						else
							continue;
						end
					end
				end
				self.mergedcaps{2*(i-1)+1}=md;
			end
		end % }}}
		function listcaps(self) % {{{
			for  i=1:length(self.icecaps),
				disp(sprintf('%i: %s',i,self.icecaps{i}.miscellaneous.name));
			end
		end % }}}
		function n=ncaps(self) % {{{
			n=length(self.icecaps);
		end % }}}
		function list=continents(self) % {{{
			list={};
			for  i=1:length(self.basins),
				list{end+1}=self.basins{i}.continent;
			end
			list=unique(list);
		end % }}}
		function list=basinsfromcontinent(self,continent) % {{{
			list={};
			for  i=1:length(self.icecaps),
				if strcmpi(self.basins{i}.continent,continent),
					list{end+1}=self.basins{i}.name;
				end
			end
			list=unique(list);
		end % }}}
		function addbasin(self,bas) % {{{
			if ~strcmpi(class(bas),'basin')
				error('addbasin method only takes a ''basin'' class object as input');
			end;
			self.basins{end+1}=bas;
		end % }}}
		function intersections2d(self,varargin) % {{{

			options=pairoptions(varargin{:});
			force=getfieldvalue(options,'force',0);
			
			%initialize, to avoid issues of having more transitions than meshes.
			self.transitions={};
			self.eltransitions={};

			%for elements:
			xe=self.earth.mesh.x(self.earth.mesh.elements)*[1;1;1]/3;
			ye=self.earth.mesh.y(self.earth.mesh.elements)*[1;1;1]/3;
			
			for i=1:length(self.icecaps),
				mdi=self.icecaps{i};
		
				%for elements:
				xei=mdi.mesh.x(mdi.mesh.elements)*[1;1;1]/3;
				yei=mdi.mesh.y(mdi.mesh.elements)*[1;1;1]/3;
		
				disp(sprintf('Computing vertex intersections for basin %s',self.basins{i}.name));
			
				self.transitions{end+1}=meshintersect2d(self.earth.mesh.x,self.earth.mesh.y,mdi.mesh.x,mdi.mesh.y,'force',force);

				self.eltransitions{end+1}=meshintersect2d(xe,ye,xei,yei,'force',force);
			end
		end % }}}
		function intersections(self,varargin) % {{{

			options=pairoptions(varargin{:});
			force=getfieldvalue(options,'force',0);
			
			%initialize, to avoid issues of having more transitions than meshes.
			self.transitions={};
			self.eltransitions={};
			self.earth.solidearth.transfercount=zeros(self.earth.mesh.numberofvertices,1);

			%for elements:
			xe=self.earth.mesh.x(self.earth.mesh.elements)*[1;1;1]/3;
			ye=self.earth.mesh.y(self.earth.mesh.elements)*[1;1;1]/3;
			ze=self.earth.mesh.z(self.earth.mesh.elements)*[1;1;1]/3;
			
			for i=1:length(self.icecaps),
				mdi=self.icecaps{i};
				mdi=TwoDToThreeD(mdi,self.planet);
		
				%for elements:
				xei=mdi.mesh.x(mdi.mesh.elements)*[1;1;1]/3;
				yei=mdi.mesh.y(mdi.mesh.elements)*[1;1;1]/3;
				zei=mdi.mesh.z(mdi.mesh.elements)*[1;1;1]/3;
		
				disp(sprintf('Computing vertex intersections for basin %s',self.basins{i}.name));
			
				self.transitions{end+1}=meshintersect3d(self.earth.mesh.x,self.earth.mesh.y,self.earth.mesh.z,mdi.mesh.x,mdi.mesh.y,mdi.mesh.z,'force',force);

				self.eltransitions{end+1}=meshintersect3d(xe,ye,ze,xei,yei,zei,'force',force);

				self.earth.solidearth.transfercount(self.transitions{i})=self.earth.solidearth.transfercount(self.transitions{i})+1;
			end

			for i=1:length(self.icecaps),
				self.icecaps{i}.solidearth.transfercount=self.earth.solidearth.transfercount(self.transitions{i});
			end
		end % }}}
		function checkintersections(self) % {{{
			flags=zeros(self.earth.mesh.numberofvertices,1);
			for i=1:length(self.basins),
				flags(self.transitions{i})=i;
			end
			plotmodel(self.earth,'data',flags,'coastlines','on');

		end % }}}
		function checkbasinconsistency(self) % {{{
			for i=1:length(self.basins),
				self.basins{i}.checkconsistency();
			end

		end % }}}
		function baslist=basinindx(self,varargin) % {{{
			options=pairoptions(varargin{:});
			continent=getfieldvalue(options,'continent','all');
			bas=getfieldvalue(options,'basin','all');

			%expand continent list: {{{
			if iscell(continent),
				if length(continent)==1,
					if strcmpi(continent{1},'all'),
						%need to transform this into a list of continents:
						continent={};
						for i=1:length(self.basins),
							continent{end+1}=self.basins{i}.continent;
						end
						continent=unique(continent);
					else
						%nothing to do: assume we have a list of continents
					end
				else
					%nothing to do: assume we have a list of continents
				end
			else
				if strcmpi(continent,'all'),
					%need to transform this into a list of continents:
					continent={};
					for i=1:length(self.basins),
						continent{end+1}=self.basins{i}.continent;
					end
					continent=unique(continent);
				else
					continent={continent};
				end
			end
			%}}}
			%expand basins list using the continent list above and the extra bas discriminator: %{{{
			if iscell(bas),
				if length(bas)==1,
					if strcmpi(bas{1},'all'),
						%need to transform this into a list of basins:
						baslist=[];
						for i=1:length(self.basins),
							if self.basins{i}.iscontinentany(continent{:}),
								baslist(end+1)=i;
							end
						end
						baslist=unique(baslist);
					else
					bas=bas{1};
					baslist=[];
					for i=1:length(self.basins),
						if self.basins{i}.iscontinentany(continent{:}),
							if self.basins{i}.isnameany(bas),
								baslist(end+1)=i;
							end
						end
					end

					end
				else
					%we have a list of basin names:
					baslist=[];
					for i=1:length(bas),
						basname=bas{i};
						for j=1:length(self.basins),
							if self.basins{j}.iscontinentany(continent{:}),
								if self.basins{j}.isnameany(basname),
									baslist(end+1)=j;
								end
							end
						end
						baslist=unique(baslist);
					end
				end
			else
				if strcmpi(bas,'all'),
					baslist=[];
					for i=1:length(self.basins),
						if self.basins{i}.iscontinentany(continent{:}),
							baslist(end+1)=i;
						end
					end
					baslist=unique(baslist);
				else
					baslist=[];
					for i=1:length(self.basins),
						if self.basins{i}.iscontinentany(continent{:}),
							if self.basins{i}.isnameany(bas),
								baslist(end+1)=i;
							end
						end
					end
					baslist=unique(baslist);
				end
			end
			%}}}

		end % }}}
		function addicecap(self,md) % {{{
			if ~strcmpi(class(md),'model')
				error('addicecap method only takes a ''model'' class object as input');
			end
			self.icecaps{end+1}=md;
		end % }}}
		function basinsplot3d(self,varargin) % {{{
			for i=1:length(self.basins),
				self.basins{i}.plot3d(varargin{:});
			end
		end % }}}
		function caticecaps(self,varargin) % {{{
			
			%recover options:
			options=pairoptions(varargin{:});
			tolerance=getfieldvalue(options,'tolerance',.65);
			loneedgesdetect=getfieldvalue(options,'loneedgesdetect',0);
	
			%make 3D model:
			models=self.icecaps;
			for i=1:length(models),
				models{i}=TwoDToThreeD(models{i},self.planet);
			end
			
			%Plug all models together:
			md=models{1};
			for i=2:length(models),
				md=modelmerge3d(md,models{i},'tolerance',tolerance);
				md.private.bamg.landmask=[md.private.bamg.landmask;models{i}.private.bamg.landmask];
			end

			%Look for lone edges if asked for it: {{{
			if loneedgesdetect,
				edges=loneedges(md);
				plotmodel(md,'data',md.mask.land_levelset);
				hold on;
				for i=1:length(edges),
					ind1=edges(i,1);
					ind2=edges(i,2);
					%plot([md.mesh.x(ind1),md.mesh.x(ind2)],[md.mesh.y(ind1),md.mesh.y(ind2)],'r*-');
					plot3([md.mesh.x(ind1),md.mesh.x(ind2)],[md.mesh.y(ind1),md.mesh.y(ind2)],[md.mesh.z(ind1),md.mesh.z(ind2)],'g*-');
				end
			end %}}}
	
			%Plug into earth:
			self.earth=md;

			%Create mesh radius:
			self.earth.mesh.r=planetradius('earth')*ones(md.mesh.numberofvertices,1);

		end % }}}
		function caticecaps2d(self,varargin) % {{{

			%recover options:
			options=pairoptions(varargin{:});
			tolerance=getfieldvalue(options,'tolerance',1e-5);
			loneedgesdetect=getfieldvalue(options,'loneedgesdetect',0);
			models=self.icecaps;

			%Plug all models together:
			md=models{1};
			for i=2:length(models),
				md=modelmerge2d(md,models{i},'tolerance',tolerance);
			end

			%Look for lone edges if asked for it: {{{
			if loneedgesdetect,
				edges=loneedges(md);
				hold on;
				for i=1:length(edges),
					ind1=edges(i,1);
					ind2=edges(i,2);
					plot([md.mesh.x(ind1),md.mesh.x(ind2)],[md.mesh.y(ind1),md.mesh.y(ind2)],'g*-');
				end
			end %}}}

			%Plug into earth:
			self.earth=md;

		end % }}}
		function viscousiterations(self) % {{{
			for  i=1:length(self.icecaps),
				ic=self.icecaps{i};
				mvi=ic.results.TransientSolution(1).StressbalanceConvergenceNumSteps;
				for j=2:length(ic.results.TransientSolution)-1,
					mvi=max(mvi,ic.results.TransientSolution(j).StressbalanceConvergenceNumSteps);
				end
				disp(sprintf('%i, %s: %i',i,self.icecaps{i}.miscellaneous.name,mvi));
			end
		end % }}}
		function maxtimestep(self) % {{{
			for  i=1:length(self.icecaps),
				ic=self.icecaps{i};
				mvi=length(ic.results.TransientSolution);
				timei=ic.results.TransientSolution(end).time;
				disp(sprintf('%i, %s: %i/%g',i,self.icecaps{i}.miscellaneous.name,mvi,timei));
			end
			mvi=length(self.earth.results.TransientSolution);
			timei=self.earth.results.TransientSolution(end).time;
			disp(sprintf('Earth: %i/%g',mvi,timei));
		end % }}}
		function transfer(self,string) % {{{
			%Recover field size in one icecap:
			eval(['n=size(self.icecaps{1}.' string ',1);']);
			if n==self.icecaps{1}.mesh.numberofvertices,
				eval(['self.earth.' string '=zeros(self.earth.mesh.numberofvertices,1);']);
				for i=1:length(self.icecaps),
					eval(['self.earth.' string '(self.transitions{' num2str(i) '})=self.icecaps{' num2str(i) '}.' string ';']);
				end
			elseif n==(self.icecaps{1}.mesh.numberofvertices+1),
				%dealing with a transient dataset.
				%check that all timetags are similar between all icecaps:  %{{{
				for i=1:length(self.icecaps),
					eval(['capfieldi= self.icecaps{' num2str(i) '}.' string ';']);
					for j=(i+1):length(self.icecaps),
						eval(['capfieldj= self.icecaps{' num2str(j) '}.' string ';']);
						if ~isequal(capfieldi(end,:),capfieldj(end,:)),
							error(['Time stamps for ' string ' field are different between icecaps ' num2str(i) ' and ' num2str(j)]);
						end
					end
				end
				eval(['capfield1= self.icecaps{1}.' string ';']);
				times=capfield1(end,:);
				nsteps=length(times);
				%}}}
				%initialize:  %{{{
				eval(['field=zeros(self.earth.mesh.numberofvertices+1,' num2str(nsteps) ');']);
				eval(['field(end,:)=times;']); %transfer the times only, not the values
				%}}}
				%transfer all time fields: {{{
				for i=1:length(self.icecaps),
					eval(['capfieldi= self.icecaps{' num2str(i) '}.' string ';']);
					for j=1:nsteps,
						eval(['field(self.transitions{' num2str(i) '},' num2str(j) ')=capfieldi(1:end-1,' num2str(j) ');']); %transfer only the values, not the time.
					end
				end
				eval(['self.earth.' string '=field;']); %do not forget to plug the field variable into its final location
				%}}}
			elseif n==self.icecaps{1}.mesh.numberofelements,
				eval(['self.earth.' string '=zeros(self.earth.mesh.numberofelements,1);']);
				for i=1:length(self.icecaps),
					eval(['self.earth.' string '(self.eltransitions{' num2str(i) '})=self.icecaps{' num2str(i) '}.' string ';']);
				end
			else
				error('not supported yet');
			end
		end % }}}
		function self=homogeneize(self,noearth) % {{{
			if nargin==1,
				noearth=0;
			end
			mintimestep=Inf;
			for  i=1:length(self.icecaps),
				ic=self.icecaps{i};
				mintimestep=min(mintimestep, length(ic.results.TransientSolution));
			end
			if ~noearth,
				mintimestep=min(mintimestep, length(self.earth.results.TransientSolution));
			end
			
			for  i=1:length(self.icecaps),
				ic=self.icecaps{i};
				ic.results.TransientSolution=ic.results.TransientSolution(1:mintimestep);
				self.icecaps{i}=ic;
			end
			ic=self.earth;
			if ~noearth,
				ic.results.TransientSolution=ic.results.TransientSolution(1:mintimestep);
			end
			self.earth=ic;
		end % }}}
		function self=initializemodels(self) % {{{

			for i=1:length(self.basins),
				md=model();
				md.miscellaneous.name=self.basins{i}.name;
				self.addicecap(md);
			end
		end % }}}
	end
end
