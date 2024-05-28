%GLACIER_INVENTORY class definition
%
%   Usage:
%      inv = glacier_inventory(varargin)
%
%      where varargin is a variable list of options: 
%
%   Usage: 
%      rgi = glacier_inventory('root',shapefileroot,...
%                       'filenames',region_names,...
%                       'boxes', region_boxes);
%   Example: 
%      rgi = glacier_inventory('root','~/ModelData',...
%                       'filenames',{'01_rgi60_Alaska','02_rgi60_WesternCanadaUS'},...
%                       'boxes','00_rgi60_O2Regions.shp');
%
%   Watch out boxes do not have to match the region shapefiles. Example: O1regions vs O2 regions shapefiles. 

classdef glacier_inventory < handle
	properties (SetAccess=public) %Model fields
		root    = '';
		regions          = struct();
		boxes            = struct();
		element_connectivity = [];
		glacier_connectivity = [];
	end
	methods
		
		function self = glacier_inventory(varargin) % {{{

			options=pairoptions(varargin{:}); 

			self.root=getfieldvalue(options,'root');
			region_names=getfieldvalue(options,'filenames');
			boxes_filename=getfieldvalue(options,'boxes');

			%first read boxes regions shapefile: 
			disp('reading boxes for each region');
			self.boxes=shpread([self.root '/' boxes_filename]);

			%read the shape files and create the regions: 
			counter=0;
			self.regions=struct();
			for i=1:length(region_names),
				disp(['reading region: '  region_names{i}]);
				self.regions(i).name=region_names{i};
				self.regions(i).id=i;
				contours=shpread([self.root '/' self.regions(i).name '.shp']);

				%we can't keep all the info: build arrays of centroids instead of reading 
				%the contours.
				O1Region=zeros(length(contours),1);
				O2Region=zeros(length(contours),1);
				Area=zeros(length(contours),1);
				CenLon=zeros(length(contours),1);
				CenLat=zeros(length(contours),1);
				Connectivity=zeros(length(contours),1);
				for j=1:length(contours),
					O1Region(j)=str2num(contours(j).O1Region);
					O2Region(j)=str2num(contours(j).O2Region);
					Area(j)=contours(j).Area;
					CenLon(j)=contours(j).CenLon;
					CenLat(j)=contours(j).CenLat;
					Connectivity(j)=contours(j).Connect;
				end
				self.regions(i).Area=Area;
				self.regions(i).O1Region=O1Region;
				self.regions(i).O2Region=O2Region;
				self.regions(i).CenLat=CenLat;
				self.regions(i).CenLon=CenLon;
				self.regions(i).Connectivity=Connectivity;
				self.regions(i).lids=[1:length(contours)]';
				self.regions(i).gids=self.regions(i).lids+counter;
				counter=counter+length(contours);
			end
	end
	%}}}
		function readcontours(self) % {{{
			%we are reading the contours from each shape file for each region, reproject
			%the latlong in these contours to the local laea projection, and rewrite the shapefile
			%we a different name extension.

			for i=1:self.nregions,
				
				disp(['reading shapefile for region: '  self.regions(i).name]);
				contours=shaperead([self.root '/' self.regions(i).name '.shp']);

				disp(['concatenating contours (X,Y)']);
				count=0;
				for j=1:length(contours),
					count=count+length(contours(j).X);
				end
				xc=zeros(count,1); 
				yc=zeros(count,1);
				count=0;
				for j=1:length(contours),
					nj=length(contours(j).X);
					xc((count+j):(count+j+nj-1))=contours(j).X;
					yc((count+j):(count+j+nj-1))=contours(j).Y;
					count=count+nj;
				end

				disp(['projecting (X,Y)']);
				lat0= fix(mean(self.regions(i).CenLat)); long0=fix(mean(self.regions(i).CenLon));
				proj=laea(lat0,long0); self.regions(i).proj=proj;
				[xc,yc]=gdaltransform(xc,yc,'EPSG:4326',proj);

				disp(['plugging back into contours']);
				count=0;
				for j=1:length(contours),
					nj=length(contours(j).X);
					contours(j).X=xc((count+j):(count+j+nj-1));
					contours(j).Y=yc((count+j):(count+j+nj-1));
					count=count+nj;
				end
				disp(['saving new contours to disk']);
				shapewrite(contours,[self.root '/' self.regions(i).name '.laeaproj.shp']);
			end
		end
		%}}}
		function varargout=loadcontour(self,id) % {{{
			
			%go find the projected contours for an 'id' region: 
			disp(['reading projected shapefile for region: '  self.regions(id).name]);
			contours=shaperead([self.root '/' self.regions(id).name '.laeaproj.shp']);
			self.regions(id).contours=contours;

			if nargout==1,
				varargout{1}=contours;
			end

		end
		%}}}
		function disp(self) % {{{
			disp(sprintf('   Glacier inventory:')); 

			disp(['   number of regions: ' num2str(self.nregions())]);
			for i=1:self.nregions(),
				disp(sprintf('      region %i: ''%s'' %i glaciers ',i,self.regions(i).name,length(self.regions(i).Area)));
			end

		end % }}}
		function mesh_connectivity(self,mesh,varargin) % {{{
			
			%retrieve options: 
			options=pairoptions(varargin{:}); 

			subsetregions=getfieldvalue(options,'regions',1:length(self.boxes));
			errornotfound=getfieldvalue(options,'errornotfound',1);
			plotr=getfieldvalue(options,'plot',0);
			plotr=1;

			%The way we run this is through the O2 zones defined in boxes. We go through 
			%these  regions, figure out the centroid, figure out how many elements are close to
			%this centroid (very important to do this through vertex lat,long, and not through elemnet
			% lat,long which can be seriously warped because of the -180 to +180 longitude transition. 
			%We then project the region in lamber azimuthal equal area, along with glaciers and elements. 
			%Once projected, we figure out which glaciers belong to which element  in the local 
			%projection system. 

			%initialize glacier connectivity: 
			self.glacier_connectivity=zeros(self.nglaciers,1);

			%assume maximum width for connectivity table and initialize 
			%as sparse matrix: 
			ny=self.nglaciers(); self.element_connectivity=sparse(mesh.numberofelements,ny);
				
			disp('Building element and glacier connectivity table: ');

			%O2 regions: 
			o1=zeros(self.nglaciers(),1);
			o2=zeros(self.nglaciers(),1);
			counter=1;
			for i=1:self.nregions(),
				region=self.regions(i);
				for j=1:length(region.CenLat),
					o1(counter)=region.O1Region(j);
					o2(counter)=region.O2Region(j);
					counter=counter+1;
				end
			end

			%Go through O2 regions: 
			for i=subsetregions,
			%for i=33,
				string=self.boxes(i).RGI_CODE; 
				disp(['progressing with region ' num2str(i) ' ' string]);
				offset=findstr(string,'-'); 
				o1i=str2num(string(1:offset-1));
				o2i=str2num(string(offset+1:end));
				glaciers=find(o1==o1i & o2==o2i);

				if ~isempty(glaciers),
					%find lat,long for laea projection: 
					box=self.boxes(i).BoundingBox; 
					long0=mean(box(:,1));
					lat0=mean(box(:,2));
					proj=laea(lat0,long0);

					%find radius of box: 
					minlat=min(box(:,2)); maxlong=max(box(:,1)); 
					radius=sqrt( (lat0-minlat)^2+(long0-maxlong)^2);

					%some radius adjustment:  {{{
					switch i,
						case 4, radius=40;
						case 8, radius=60;
						case 12, radius=25;
						case 19, radius=60;
						case 32, radius=60;
						case 33, radius=10;
						case 41, radius=75;
						case 42, radius=45;
						case 61, radius=66;
						case 68, radius=10;
						case 82, radius=30;
						otherwise,
					end % }}}

					[lids,rids]=self.gidtolid(glaciers);
					%quick check on the rids, should all be the same number: 
					if min(rids)~=max(rids)error(sprintf('rids should only span on O1 region')); end;
					rid=max(rids);

					region=self.regions(rid);
					elements=flaglatlongradiuselements(mesh.elements,mesh.lat,mesh.long,lat0,long0,radius);

					if plotr, % {{{
						figure(1),clf; 
						subplot(2,1,1),hold on;
						trisurf(mesh.elements(elements,:),mesh.long,mesh.lat,mesh.lat),view(2); 
						plot3(box(1,1),box(1,2),1000,'r*','MarkerSize',10);
						plot3(box(1,1),box(2,2),1000,'r*','MarkerSize',10);

						plot3(box(2,1),box(1,2),1000,'r*','MarkerSize',10);
						plot3(box(2,1),box(2,2),1000,'r*','MarkerSize',10);

						plot3(region.CenLon(lids),region.CenLat(lids),1000*ones(length(lids),1),'k*');
					end % }}}

					%project lat,long: 
					[x,y]=gdaltransform(mesh.long,mesh.lat,'EPSG:4326',proj);
					[xlid,ylid]=gdaltransform(region.CenLon(lids),region.CenLat(lids),'EPSG:4326',proj);

					if plotr, % {{{
						figure(1),subplot(2,1,2), hold on;
						trisurf(mesh.elements(elements,:),x,y,x),view(2); 
						plot3(xlid,ylid,1000*ones(length(lids),1),'k*');
						pause(.1);
					end % }}}

					%go through lids: 
					for j=1:length(lids),
						found=0;
						x0=xlid(j); y0=ylid(j);
						for k=1:length(elements),
							el=elements(k);
							x1=x(mesh.elements(el,1)); y1=y(mesh.elements(el,1));
							x2=x(mesh.elements(el,2)); y2=y(mesh.elements(el,2));
							x3=x(mesh.elements(el,3)); y3=y(mesh.elements(el,3));

							if isintriangle(x0,x1,x2,x3,y0,y1,y2,y3),
								found=1;
								break;
							end
						end
						if ~found,
							if errornotfound,
								error(sprintf('could not find element for glacier %i with lid %i',j,lids(j)));
							end
						end
						if(found)self.glacier_connectivity(glaciers(j))=el; end;
					end
				end
			end 

			%build element connectivity table: 
			for j=1:length(self.glacier_connectivity),
				el=self.glacier_connectivity(j);
				if ~el,continue; end;
				count=self.element_connectivity(el,ny);
				if count>ny,
					error('need to enlarge connectivity table to at least');
				end
				self.element_connectivity(el,count+1)=j;
				self.element_connectivity(el,ny)=count+1;
			end

			%Reduce the number of columns (we did not initially, so we chose an arbitrary ny:
			nmax=max(self.element_connectivity(:,end));
			self.element_connectivity=self.element_connectivity(:,[1:nmax,ny]);


		end % }}}
		function mesh_connectivity2d(self,md,proj,varargin) % {{{
			
			%retrieve options: 
			options=pairoptions(varargin{:}); 

			subsetregions=getfieldvalue(options,'regions',1:self.nregions());

			%initialize glacier connectivity: 
			self.glacier_connectivity=zeros(self.nglaciers,1);

			%assume maximum width for connectivity table and initialize 
			%as sparse matrix: 
			ny=self.nglaciers(); self.element_connectivity=sparse(md.mesh.numberofelements,ny);
				
			disp('Building element and glacier connectivity table: ');
			[lat0,long0]=projlatlong(proj);
		
			[mpartition,npartition]=self.partition();
			for i=subsetregions,
				region=self.regions(i);
				disp(sprintf(' progress for region: %s',region.name));

				%project lat,long: 
				[xlid,ylid]=gdaltransform(region.CenLon,region.CenLat,'EPSG:4326',proj);

				%go through lids: 
				x0=xlid; y0=ylid;
				x1=md.mesh.x(md.mesh.elements(:,1)); y1=md.mesh.y(md.mesh.elements(:,1));
				x2=md.mesh.x(md.mesh.elements(:,2)); y2=md.mesh.y(md.mesh.elements(:,2));
				x3=md.mesh.x(md.mesh.elements(:,3)); y3=md.mesh.y(md.mesh.elements(:,3));
				in=isintrianglearraytotal(x0,x1,x2,x3,y0,y1,y2,y3);
				[els,glacs]=find(in);
				self.glacier_connectivity(mpartition(i)-1+glacs)=els;
			end

			%build element connectivity table: 
			for j=1:length(self.glacier_connectivity),
				el=self.glacier_connectivity(j);
				if ~el,continue; end;
				count=self.element_connectivity(el,ny);
				if count>ny,
					error('need to enlarge connectivity table to at least');
				end
				self.element_connectivity(el,count+1)=j;
				self.element_connectivity(el,ny)=count+1;
			end

			%Reduce the number of columns (we did not initially, so we chose an arbitrary ny:
			nmax=max(self.element_connectivity(:,end));
			self.element_connectivity=self.element_connectivity(:,[1:nmax,ny]);


		end % }}}
		function checkconnectivity(self,mesh) % {{{

			vector=find(self.element_connectivity(:,end));

			for i=1:length(vector),
				el=vector(i);

				flags=zeros(mesh.numberofelements,1);
				flags(el)=1;

				nglaciers=self.element_connectivity(el,end); 
				glaciers=self.element_connectivity(el,1:nglaciers);

				[lids,rids]=self.gidtolid(glaciers);
				lat=zeros(length(glaciers),1);
				long=zeros(length(glaciers),1);
				for j=1:nglaciers,
					lat(j)=self.regions(rids(j)).CenLat(lids(j));
					long(j)=self.regions(rids(j)).CenLon(lids(j));
				end

				proj=laea(lat(1),long(1));

				figure(1),clf,hold on;
				[x,y]=gdaltransform(mesh.long(mesh.elements(el,:)),mesh.lat(mesh.elements(el,:)),'EPSG:4326',proj);
				p=patch('XData',x,'YData',y);  set(p,'FaceColor','red')

				[x,y]=gdaltransform(long,lat,'EPSG:4326',proj);
				plot(x,y,'k*');
				pause(.1);
			end

		end % }}}
		function totalarea=area(self,varargin) % {{{
			region=-1;
			totalarea=0;
			if nargin==2,
				region=varargin{1};
			end
			if region==-1,
				%figure out the areas of everybody: 
				for i=1:self.nregions(),
					totalarea=totalarea+sum(self.regions(i).Area);
				end
			else
				totalarea=totalarea+sum(self.regions(region).Area);
			end

		end % }}}
		function [mpartition,npartition]=partition(self,varargin) % {{{
			mpartition=zeros(self.nregions(),1);
			npartition=zeros(self.nregions(),1);
			counter=0;
			for i=1:self.nregions(),
				mpartition(i)=counter+1;
				npartition(i)=counter+self.nglaciers(i);
				counter=counter+self.nglaciers(i);
			end

		end % }}}
		function n = nregions(self) % {{{
			n=length(self.regions);
		end
		%}}}
		function counter = nglaciers(self,varargin) % {{{

			if nargin==1,
				region=-1; %all regions.
			else
				region=varargin{1}; %only one.
			end

			if region==-1,
				counter=0;
				for i=1:self.nregions(),
					counter=counter+length(self.regions(i).Area);
				end
			else
				counter=length(self.regions(region).Area);
			end
		end
		%}}}
		function name=ridtoname(self,rid) % {{{
			
			fullname=self.regions(rid).name; 
			name=fullname(10:end);

		end % }}}
		function [gid]=lidtogid(self,rid,lid) % {{{
			[mpartition,npartition]=self.partition();
			gid=mpartition(rid)+lid-1;
		end % }}}
		function [lid,rid]=gidtolid(self,gid) % {{{
			[mpartition,npartition]=self.partition();
			lid=zeros(length(gid),1);
			rid=zeros(length(gid),1);
			for i=1:self.nregions(),
				pos=find(gid>=mpartition(i) & gid<=npartition(i)); 
				rid(pos)=i;
				lid(pos)=gid(pos)-mpartition(i)+1;
			end

		end % }}}
		function units(self) % {{{
			disp('Glacier inventory units: ');
			disp('   areas: km^2');
			disp('   mass: Gt');

		end
		%}}}
		function listboxes(self) % {{{

			for i=1:length(self.boxes),
				b=self.boxes(i);
				disp(sprintf('Region #%i: %s (%s)',i,b.FULL_NAME,b.RGI_CODE));
			end
			
		end
		%}}}
	end
end
