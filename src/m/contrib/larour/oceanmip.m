%OCEANMIP class definition
%
%   Usage:
%      omip = oceanmip('root',rootdir,'files',listfiles);
%

classdef oceanmip < handle
	properties (SetAccess=public) %Model fields
		
		root = ''; %where are the files for CMIP5
		n   = 0;   %number of files
		files   = {};   %netcdf files
		zos = {}; % local sea-level anomalies with zero mean over the oceans (in mm) (lat,lon,time): 
		zostoga = {}; % global-mean thermosteric sea level anomalies (in mm) (time)
		time = {}; %time in year
		pbo = {}; %local ocean-bottom pressure changes. Also have zero mean over the oceans (lat,lon,time)
		lat = {};
		long = {};
		scenario = {}; 
		model = {}; 
		mesh_zos = {}; %interpolated (onto model mesh) zos field
		mesh_pbo = {}; %interpolated (onto model mesh) pbo field

	end
	methods
		function self = oceanmip(varargin) % {{{

			if nargin==0, 
				self=setdefaultparameters(self);
			else 
				self=setdefaultparameters(self);

				options=pairoptions(varargin{:});

				self.root=getfieldvalue(options,'root');
				self.files=getfieldvalue(options,'files');
				self.n=length(self.files);

				%read files: 
				for f=1:self.n,
					file=[self.root '/' self.files{f}];
					disp(['reading file ' file]);
					
					%figure out time interval and remove historical:
					time=ncread(file,'time');
					pos=find(diff(time)<0); 
					if isempty(pos),
						%pos=(length(time)-12*100+1):length(time);
						pos=1:length(time);
					else
						pos=[(pos+1):length(time)]; 
					end
					time2=time(pos); 
					pos2=find(time2<=2099 & time2>=2006);
					pos=pos(pos2);
					time=time(pos); nt=length(time);

					%reduce datasets: 
					time=floor(time(12:12:nt)); 
					self.time{end+1}=time;
					
					%zos: 
					zos=ncread(file,'zos'); nx=size(zos,1); ny=size(zos,2); zos=zos(:,:,pos);
					zosm=zeros(nx,ny,nt/12);
					for i=12:12:nt,
						year=i/12;
						zosm(:,:,year)=mean(zos(:,:,(i-11):i),3);
					end
					self.zos{end+1}=zosm; clear zos;

					%zostoga:
					zostoga=ncread(file,'zostoga'); zostoga=zostoga(pos);
					zostogam=zeros(nt/12,1);
					for i=12:12:nt,
						year=i/12;
						zostogam(year)=mean(zostoga(i-11:i));
					end
					%control against 2006: 
					zostogam=zostogam-zostogam(1);
					self.zostoga{end+1}=zostogam; clear zostoga;

					%pbo: 
					pbo=ncread(file,'pbo'); nx=size(pbo,1); ny=size(pbo,2); pbo=pbo(:,:,pos);
					pbom=zeros(nx,ny,nt/12);
					for i=12:12:nt,
						year=i/12;
						pbom(:,:,year)=mean(pbo(:,:,i-11:i),3);
					end
					self.pbo{end+1}=pbom; clear pbo;

					self.lat{end+1}=ncread(file,'lat');
					self.long{end+1}=ncread(file,'lon');

					%scenario and model: 
					file=self.files{f};
					ind=findstr(file,'rcp'); 
					self.scenario{end+1}=str2num(file(ind+3:ind+4));
					self.model{end+1}=file(1:ind-2);
				end

			end
		end
		%}}}
		function self = interpolate(self,md) % {{{

			%retrieve long and lat from mesh: 
			meshlong=md.mesh.long;
			pos=find(meshlong<0); meshlong(pos)=meshlong(pos)+360;
			meshlat=md.mesh.lat;

			for i=1:self.n,
				disp(['interpolating model ' self.model{i} ' onto model mesh']);

				%If we have 182 cells in long, trim by 1 on each side:
				if size(self.pbo{i},1)==182,
					self.pbo{i}=self.pbo{i}(2:181,:,:);
					self.zos{i}=self.zos{i}(2:181,:,:);
					self.long{i}=self.long{i}(2:181,:,:);
					self.lat{i}=self.lat{i}(2:181,:,:);
				end


				%triangulation: 
				long=double(self.long{i});  long=long(:); 
				pos=find(long<0); long(pos)=long(pos)+360;
				lat=double(self.lat{i}); lat=lat(:);
				[newl,uniqpos]=unique([lat,long],'rows','stable');
				long=long(uniqpos); lat=lat(uniqpos);
				index=delaunay(long,lat);

				%some checks: 
				areas=GetAreas(index,long,lat); indneg=find(areas<0);  
				index(indneg,:)=[index(indneg,1),index(indneg,3),index(indneg,2)];
				areas=GetAreas(index,long,lat); 
				ind=find(areas<1e-8); index(ind,:)=[];

				%fix if we have orphans
				[index long lat dummy newpos]=FixMesh(index,long,lat,1:length(long));

				time=self.time{i};
				%retrieve fields:
				omip_pbo=self.pbo{i};
				omip_zos=self.zos{i};
				omip_zostoga=self.zostoga{i};

				%interpolate:
				mesh_pbo=zeros(md.mesh.numberofvertices,length(time));
				mesh_zos=zeros(md.mesh.numberofvertices,length(time));

				parfor j=1:length(time),
				%for j=1:length(time),
					if mod(j,10)==0, 
						s=sprintf('   progress: %.2g ',j/length(time)*100);
						fprintf(s); pause(1e-3); fprintf(repmat('\b',1,numel(s))); 
					end

					pbo=omip_pbo(:,:,j); pbo=pbo(:); pbo=pbo(uniqpos); 
					zos=omip_zos(:,:,j); zos=zos(:); zos=zos(uniqpos);

					pboj= InterpFromMeshToMesh2d(index,long,lat,pbo(newpos),meshlong,meshlat);
					pos=find(abs(pboj)>1e4); pboj(pos)=0;
					mesh_pbo(:,j)=pboj;

					zosj=InterpFromMeshToMesh2d(index,long,lat,zos(newpos),meshlong,meshlat);
					pos=find(abs(zosj)>1e5); zosj(pos)=0;
					mesh_zos(:,j)=zosj;
				end
				self.mesh_pbo{end+1}=mesh_pbo;
				self.mesh_zos{end+1}=mesh_zos;

				%clear for memory purposes: 
				self.pbo{i}=[];
				self.zos{i}=[];
			end
		end

		%}}}
		function self = setdefaultparameters(self) % {{{
		end
		%}}}
		function self = rawclean(self) % {{{
			for i=1:self.n,
				self.zos{i}=[];
				self.pbo{i}=[];
			end
		end
		%}}}
		function [rate,time]= zostoga_mean(self) % {{{
			series=zeros(length(self.time{1}),self.n);
			for i=1:self.n,
				series(:,i)=self.zostoga{i};
			end
			rate=mean(series,2);
			time=self.time{1};
		end
		%}}}
		function [rate,time]= zostoga_std(self) % {{{
			series=zeros(length(self.time{1}),self.n);
			for i=1:self.n,
				series(:,i)=self.zostoga{i};
			end
			rate=std(series,1,2);
			time=self.time{1};
		end
		%}}}
		function [average,stddev,time]= zostoga_stats(self) % {{{
			series=zeros(length(self.time{1}),self.n);
			for i=1:self.n,
				series(:,i)=self.zostoga{i};
			end
			average=mean(series,2);
			stddev=std(series,1,2);
			time=self.time{1};
		end
		%}}}
	function array= bottompressure(self,model,gridded) % {{{
			for i=1:self.n,
				if strcmpi(model,self.model{i}),
					if gridded,
						pbo=self.pbo{i}; pbo=pbo/1000; %in meters
					else
						pbo=self.mesh_pbo{i}; pbo=pbo/1000; %in meters
					end
					array=pbo;
					break;
				end
			end
		end
		%}}}
	function arrays= dbottompressures(self,varargin) % {{{
		options=pairoptions(varargin{:});
		units=getfieldvalue(options,'units','mm/yr');
		
		arrays=cell(self.n,1);
		for i=1:self.n,
			pbo=self.mesh_pbo{i}; 
			dpbo=diff(pbo,1,2);
			if strcmpi(units,'m/yr'),
				dpbo=dpbo/1000;
			end
			time=self.time{i}; time=time(1:end-1);
			array=[dpbo;time'];
			arrays{i}=array;
		end
	end %}}}
	function arrays= dzoss(self,varargin) % {{{
		options=pairoptions(varargin{:});
		units=getfieldvalue(options,'units','mm/yr');
		
		arrays=cell(self.n,1);
		for i=1:self.n,
			zos=self.mesh_zos{i}; 
			dzos=diff(zos,1,2);
			if strcmpi(units,'m/yr'),
				dzos=dzos/1000;
			end
			time=self.time{i};time=time(1:end-1);
			array=[dzos;time'];
			arrays{i}=array;
		end
	end	%}}}
	function arrays= dzostogass(self,varargin) % {{{
		options=pairoptions(varargin{:});
		units=getfieldvalue(options,'units','mm/yr');

		arrays=cell(self.n,1);
		for i=1:self.n,
			zostoga=self.zostoga{i}; 
			dzostoga=diff(zostoga);
			if strcmpi(units,'m/yr'),
				dzostoga=dzostoga/1000;
			end
			time=self.time{i};time=time(1:end-1);
			array=[dzostoga';time'];
			arrays{i}=array;
		end
	end	%}}}
		function [lat,long]= latlong(self,model) % {{{
			for i=1:self.n,
				if strcmpi(model,self.model{i}),
					lat=self.lat{i};
					long=self.long{i};
					break;
				end
			end
		end
		%}}}
		function disp(self) % {{{
			disp('   CMIP5 Ocean MIP:');

				fielddisplay(self,'n','number of files');
				fielddisplay(self,'root','where are the files for CMIP5');
				fielddisplay(self,'files','CMIP5 ocean files');
				fielddisplay(self,'zos','local sea-level anomalies with zero mean over the oceans (in mm) (lat,lon,time)');
				fielddisplay(self,'zostoga','global-mean thermosteric sea level anomalies (in mm) (time)');
				fielddisplay(self,'time','time in years');
				fielddisplay(self,'pbo','local ocean-bottom pressure changes. Also have zero mean over the oceans (lat,lon,time)');
				fielddisplay(self,'lat','latitudes');
				fielddisplay(self,'long','longitudes');
				fielddisplay(self,'scenario','scenarios');
				fielddisplay(self,'model','model names');
		end % }}}
	end 
end
