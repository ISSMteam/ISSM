%BOUNDARY class definition
%
%   Usage:
%      boundary=boundary();

classdef boundary
	properties (SetAccess=public) 
		shppath           = '';
		shpfilename       = '';
		orientation       = 'normal';  %other value is 'reverse'
		proj              = epsg2proj(4326); %Proj4.0  projection string. , default is WGS 84 (corresponds to epsg 4326)

	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here
		end% }}}
	end
	methods
		function self = boundary(varargin) % {{{
			self=setdefaultparameters(self);

			if nargin > 0,
				options=pairoptions(varargin{:});

				%recover field values: 
				self.shppath=getfieldvalue(options,'shppath','');
				self.shpfilename=getfieldvalue(options,'shpfilename','');
				self.orientation=getfieldvalue(options,'orientation','normal');

				%figure out projection string: 
				if exist(options,'epsg'),
					if exist(options,'proj'),
						error(['Error in constructor for boundary ' self.shppath '. Cannot supply epsg and proj at the same time!']);
					end
					epsg=getfieldvalue(options,'epsg');
					proj=epsg2proj(epsg); % convert to PROJ.4 format
				elseif exist(options,'proj'),
					if exist(options,'epsg'),
						error(['Error in constructor for boundary ' self.shppath '. Cannot supply proj and epsg at the same time!']);
					end
					proj=getfieldvalue(options,'proj');
				else
					proj= '+proj=longlat +datum=WGS84 +no_defs';
				end
				self.proj=proj;
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		self.shppath='';
		self.shpfilename='';
		self.orientation='normal';
		self.proj='+proj=longlat +datum=WGS84 +no_defs'; %EPSG 4326
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   boundary parameters:'));

			fielddisplay(self,'shppath','path to filename for this boundary');
			fielddisplay(self,'shpfilename','shape file name');
			fielddisplay(self,'orientation','orientation (default is ''normal'', can be ''reverse'')');
			fielddisplay(self,'proj','shape file projection string (follows PROJ.4 standard)');

		end % }}}
		function output=name(self) % {{{
			output=self.shpfilename;
		end % }}}
		function output=edges(self) % {{{

			%read domain:
			[path,name,ext]=fileparts(self.shpfilename);
			if ~strcmpi(ext,'shp'),
				ext='shp';
			end
			output=shpread([self.shppath '/' name '.' ext]);
			
			%do we reverse? 
			if strcmpi(self.orientation,'reverse'),
				output.x=flipud(output.x);
				output.y=flipud(output.y);
			end
		end % }}}
		function plot(self,varargin) % {{{
			%recover options
		
			options=pairoptions(varargin{:});

			%parse input:
			figurenumber=getfieldvalue(options,'figure',1);
			color=getfieldvalue(options,'color','r');
			linewidth=getfieldvalue(options,'linewidth',1);
			markersize=getfieldvalue(options,'markersize',1);
			unitmultiplier=getfieldvalue(options,'unit',1);
			radius=getfieldvalue(options,'radius',6371012);
			aboveground=getfieldvalue(options,'aboveground',1000);
			offset=getfieldvalue(options,'offset',.1);
			fontsize=getfieldvalue(options,'fontsize',10);
			label=getfieldvalue(options,'label','on');

			if exist(options,'epsg'),
				if exist(options,'proj'),
					error('Boundary plot error message: cannot specify epsg and proj at the same time for plot projection');
				end
				proj=epsg2proj(getfieldvalue(options,'epsg'));
			elseif exist(options,'proj'),
				proj=getfieldvalue(options,'proj');
			else
				proj=epsg2proj(4326);
			end

			%read domain:
			[path,name,ext]=fileparts(self.shpfilename);
			if ~strcmpi(ext,'shp'),
				ext='shp';
			end
			domain=shpread([self.shppath '/' name '.' ext]);

			%convert boundary to another reference frame: {{{
			for i=1:length(domain),
				try, 
					[x,y] = CoordTransform(domain(i).x,domain(i).y,self.proj,proj);
				catch me
					disp(me.message());
					self.disp();
				end
				domain(i).x=x; domain(i).y=y;
			end

			for i=1:length(domain),
				hold on;
				x=domain(i).x*unitmultiplier;
				y=domain(i).y*unitmultiplier;
				if length(x)==1,
					p=plot(x,y,'k*'); 
					set(p,'MarkerSize',markersize,'Color',color);
					if strcmpi(label,'on'),
						t=text(x,y,self.shpfilename,'FontSize',fontsize);
					end
				else
					p=plot(x,y,'k-'); 
					set(p,'MarkerSize',markersize,'Color',color);
					if strcmpi(label,'on'),
						text(sum(x)/length(x),sum(y)/length(y),self.shpfilename,'FontSize',fontsize);
					end
				end
				set(p,'Color',color);
				set(p,'LineWidth',linewidth);
			end
			%}}}
		end % }}}
		function checkconsistency(self,varargin) % {{{
			%recover options
			options=pairoptions(varargin{:});
			tolerance=getfieldvalue(options,'tolerance',1e-5);
		
			%recover edges: 
			edges=self.edges();

			if strcmpi(edges.Geometry,'Point'),
				return;
			else
				%check we don't have two identical vertices: 
				x=edges.x; y=edges.y; 
				distance=sqrt( (x(1:end-1)-x(2:end)).^2+(y(1:end-1)-y(2:end)).^2); 
				dmax=max(distance);
				isidentical=find(distance<dmax*tolerance);
				if ~isempty(isidentical),
					error(['boundary ' shpfilename ' has two vertices extremely close to one another']);
				end
			end
		end % }}}
		function plot3d(self,varargin) % {{{
			%recover options
		
			options=pairoptions(varargin{:});

			%parse input:
			figurenumber=getfieldvalue(options,'figure',1);
			color=getfieldvalue(options,'color','r');
			linewidth=getfieldvalue(options,'linewidth',1);
			markersize=getfieldvalue(options,'markersize',1);
			unitmultiplier=getfieldvalue(options,'unit',1);
			radius=getfieldvalue(options,'radius',6371012);
			aboveground=getfieldvalue(options,'aboveground',1000);
			offset=getfieldvalue(options,'offset',.1);
			fontsize=getfieldvalue(options,'fontsize',10);
			label=getfieldvalue(options,'label','on');

			if exist(options,'epsg'),
				if exist(options,'proj'),
					error('Boundary plot error message: cannot specify epsg and proj at the same time for plot projection');
				end
				proj=epsg2proj(getfieldvalue(options,'epsg'));
			elseif exist(options,'proj'),
				proj=getfieldvalue(options,'proj');
			else
				proj=epsg2proj(4326);
			end

			%read domain:
			[path,name,ext]=fileparts(self.shpfilename);
			if ~strcmpi(ext,'shp'),
				ext='shp';
			end
			domain=shpread([self.shppath '/' name '.' ext]);

			for i=1:length(domain),
				try, 
					[lat,long] = CoordTransform(domain(i).x,domain(i).y,self.proj,'EPSG:4326');
				catch me
					disp(me.message());
					self.disp();
				end
				domain(i).x=long; domain(i).y=lat;
			end

			for i=1:length(domain),

				%project on x,y,z reference frame.
				[x,y,z]=AboveGround(domain(i).x,domain(i).y,radius,aboveground);
				[xt,yt,zt]=AboveGround(domain(i).x+offset,domain(i).y+offset,radius,300000);
				hold on;
				if length(x)==1,
					p=plot3(x,y,z,'k*'); 
					set(p,'MarkerSize',markersize);
					if strcmpi(label,'on'),
						t=text(xt,yt,zt,self.shpfilename,'FontSize',fontsize);
					end
				else
					p=plot3(x,y,z,'k-'); 
					mid=floor(length(xt)/2);
					if strcmpi(label,'on'),
						text(xt(mid),yt(mid),zt(mid),self.shpfilename,'FontSize',fontsize);
					end
				end
				set(p,'Color',color);
				set(p,'LineWidth',linewidth);

			end

		end % }}}
	end
end
