%ISMIP6 class definition
%
%   Usage:
%      is6 = ismip6('root',rootdir,'dirs',listdirectories);
%   Example: 
%      is6 = ismip6('root','./gis/','dirs',{'UVW','ISSM'});

classdef ismip6 < handle
	properties (SetAccess=public) %Model fields
		
		root = ''; %where are the files for CMIP5
		n   = 0;   %number of files
		directories    = {};   %directories where the files are
		experiments    = {};   %names of experiments
		base           = {};   %placeholder for base
		surface        = {};   %placeholder for surface
		thickness      = {};   %placeholder for thicknesses
		deltathickness = {};   %placeholder for delta thicknesses
		deltathicknessvaf = {};   %placeholder for delta thicknesses above floatation
		deltathicknesshal = {};   %placeholder for delta thicknesses halosteric origins
		deltathicknessbar = {};   %placeholder for delta thicknesses halosteric origins
		thicknesscorrection={};  
		icemask        = {};   %placeholder for ice masks
		oceanmask      = {};   %placeholder for ocean masks
		time           = {};   %placeholder for times
		timestart      = {};   %placeholder for times
		calendar      = {};   %placeholder for times
		di            = {}; %ice densities
	end
	methods
		function self = ismip6(varargin) % {{{

			if nargin==0, 
				self=setdefaultparameters(self);
			else 
				self=setdefaultparameters(self);

				options=pairoptions(varargin{:});

				self.root=getfieldvalue(options,'root');
				self.directories=getfieldvalue(options,'directories');
				self.n=length(self.directories);

				%verify the directories exist: 
				for i=1:self.n,
					if ~exist([self.root '/' self.directories{i}],'dir'),
						error(['ismip6  constructor error: ' self.root '/' self.directories{i} ' does not exist']);
					end
				end

				%figure out names of experiments: 
				self.experiments=self.directories;
				for i=1:self.n,
					dir=self.directories{i};
					ind=findstr(dir,'exp');
					name=dir(1:ind-2);
					name=strrep(name,'/','-');
					self.experiments{i}=name;
				end

				%initialize fields: 
				self.thickness=cell(self.n,1);
				self.deltathickness=cell(self.n,1);
				self.icemask=cell(self.n,1);
				self.oceanmask=cell(self.n,1);

			end
		end
		%}}}
		function self = setdefaultparameters(self) % {{{
		end
		%}}}
		function disp(self) % {{{
			disp('   CMIP5 Ocean MIP:');

			fielddisplay(self,'n','number of files');
			fielddisplay(self,'root','where are the files for ISMIP6');
			fielddisplay(self,'directories','directories');
			fielddisplay(self,'experiments','experiments');
			fielddisplay(self,'thickness','thickness');
			fielddisplay(self,'deltathickness','deltathickness');
			fielddisplay(self,'icemask','icemask');
			fielddisplay(self,'oceanmask','oceanmask');
			fielddisplay(self,'time','time');
			fielddisplay(self,'timestart','timestart');
			fielddisplay(self,'calendar','calendar');
		end % }}}
		function listexp(self) % {{{
			disp('ISMIP6  list of experiments:');
			for i=1:self.n,
				disp(['   ' self.experiments{i}]);
			end

		end % }}}
		function [output,time,timestart,calendar]=read(self,experiment,field) % {{{

			%go through list of experiments and find the right one: 
			if strcmpi(class(experiment),'double'),
				ind=experiment;
			elseif strcmpi(class(experiment),'char'),
				for i=1:self.n,
					if strcmpi(experiment,self.experiments{i}),
						ind=i;
						break;
					end
				end
			else 
				error(['ismip6 read error message: experiment should be a string or index']);
			end

			if ind==0 
				error(['ismip6 read error message: experiment ' experiment ' could not be found!']);
			end;

			%figure out the files in this directory: 
			dir=self.directories{ind};
			currentdir=pwd;
			cd([self.root '/' dir]); 
			list=listfiles;
			cd(currentdir);

			%go through list of files and figure out which one starts with the field: 
			for i=1:length(list),
				file=list{i};
				ind=findstr(file,'_');
				file_field=file(1:ind-1);
				if strcmpi(file_field,field),
					break;
				end
			end

			%read file: 
			%output=ncread([self.root '/' dir '/' file],'file_field');
			time=ncread([self.root '/' dir '/' file],'time');
			output=ncread([self.root '/' dir '/' file ],field);

			%figure out start time: 
			info=ncinfo([self.root '/' dir '/' file]);
			attributes=[];
			for i=1:length(info.Variables),
				if strcmpi(info.Variables(i).Name,'time'),
					attributes=info.Variables(i).Attributes;
					break;
				end
			end
			for j=1:length(attributes),
				if strcmpi(attributes(j).Name,'units') | strcmpi(attributes(j).Name,'unit'),
					timestart=attributes(j).Value;
				end
				if strcmpi(attributes(j).Name,'calendar')
					calendar=attributes(j).Value;
				end
			end

			if ~exist('timestart','var'), timestart=2015; end
			if ~exist('calendar','var'), calendar=0; end

		end % }}}
		function info=readinfo(self,experiment,field) % {{{

			%go through list of experiments and find the right one: 
			if strcmpi(class(experiment),'double'),
				ind=experiment;
			elseif strcmpi(class(experiment),'char'),
				for i=1:self.n,
					if strcmpi(experiment,self.experiments{i}),
						ind=i;
						break;
					end
				end
			else 
				error(['ismip6 read error message: experiment should be a string or index']);
			end

			if ind==0 
				error(['ismip6 read error message: experiment ' experiment ' could not be found!']);
			end;

			%figure out the files in this directory: 
			dir=self.directories{ind};
			currentdir=pwd;
			cd([self.root '/' dir]); 
			list=listfiles;
			cd(currentdir);

			%go through list of files and figure out which one starts with the field: 
			for i=1:length(list),
				file=list{i};
				ind=findstr(file,'_');
				file_field=file(1:ind-1);
				if strcmpi(file_field,field),
					break;
				end
			end

			%read attributes
			info=ncinfo([self.root '/' dir '/' file]);

		end % }}}
		function interpolate(self,md,field,ismip2mesh,ismip2mesh_correction) % {{{

			for i=1:self.n,
				disp(['reading and interpolating field ' field ' for model ' self.experiments{i}]);

				%read field from disk: 
				[h,t,t0,cal]=self.read(i,field); nt=length(t);

				%map onto 1 dimension field: 
				ht=zeros(size(h,1)*size(h,2),nt);
				for j=1:size(h,3),
					hj= h(:,:,j)'; hj=hj(:); ht(:,j)=double(hj);
				end

				%map onto mesh: correct only for thicknesses
				if strcmpi(field,'lithk') | strcmpi(field,'orog') | strcmpi(field,'base'),
					hg=ismip2mesh_correction.*(ismip2mesh*ht) ;
					%hg=ismip2mesh*ht ;
				else
					hg=ismip2mesh*ht ;
				end

				%keep field:
				if strcmpi(field,'lithk'),
					pos=find(isnan(hg)); hg(pos)=0;
					self.thickness{i}=hg; 
				end
				if strcmpi(field,'orog'),
					pos=find(isnan(hg)); hg(pos)=0;
					self.surface{i}=hg; 
				end
				if strcmpi(field,'base'),
					pos=find(isnan(hg)); hg(pos)=0;
					self.base{i}=hg; 
				end
				if strcmpi(field,'sftgif'),
					hge=ones(md.mesh.numberofvertices,size(hg,2));
					for j=1:size(hg,2),
						hgj=hg(:,j);
						pos=find(hgj>0); 
						hge(md.mesh.elements(pos,:),j)=-1;
					end
					self.icemask{i}=hge; 
				end
				if strcmpi(field,'sftgrf'),
					hgv=-ones(md.mesh.numberofvertices,size(hg,2));
					for j=1:size(hg,2),
						hgj=hg(:,j);
						pos=find(hgj>.99); %we want fully grounded
						%pos=find(hgj>0); %we want slightly grounded
						hgv(md.mesh.elements(pos,:),j)=1;
					end
					self.oceanmask{i}=hgv; 
				end

				self.time{i}=t;
				self.timestart{i}=t0;
				self.calendar{i}=cal;

			end
		end  % }}}
		function part=partition(self,md,part,value) % {{{

			for i=1:self.n,
				dh=self.deltathickness{i}; 
				for j=1:size(dh,2),
					dhj=dh(:,j);
					pos=find(dhj);
					part(pos)=value;
				end
			end

		end  % }}}
	end
end
