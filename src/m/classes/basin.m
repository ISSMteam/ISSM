%BASIN class definition
%
%   Usage:
%      basin=basin();

classdef basin
	properties (SetAccess=public) 
		boundaries        = {};
		name              = '';
		continent         = '';
		proj              = epsg2proj(4326);
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here
		end% }}}
	end
	methods
		function self = basin(varargin) % {{{
			self=setdefaultparameters(self);

			if nargin>0,
				options=pairoptions(varargin{:}); 
		
				%recover field values: 
				self.boundaries=getfieldvalue(options,'boundaries',{});
				self.name=getfieldvalue(options,'name','');
				self.continent=getfieldvalue(options,'continent','');

				%figure out projection string: 
				if exist(options,'epsg'),
					if exist(options,'proj'),
						error(['Error in constructor for basin ' self.name '. Cannot supply epsg and proj at the same time!']);
					end
					epsg=getfieldvalue(options,'epsg');
					proj=epsg2proj(epsg); %convert to PROJ.4 format
				elseif exist(options,'proj'),
					if exist(options,'epsg'),
						error(['Error in constructor for basin ' self.name '. Cannot supply proj and epsg at the same time!']);
					end
					proj=getfieldvalue(options,'proj');
				else
					proj='+proj=longlat +datum=WGS84 +no_defs';
				end

				self.proj=proj;
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.name='';
			self.continent='';
			self.proj='+proj=longlat +datum=WGS84 +no_defs'; %EPSG 4326
			self.boundaries={};

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   basin parameters:'));
			fielddisplay(self,'continent','continent name');
			fielddisplay(self,'name','basin name');
			fielddisplay(self,'proj','basin projection string (follows PROJ.4 standard)');
			fielddisplay(self,'boundaries','list of boundary objects');
			for i=1:length(self.boundaries),
				disp(sprintf('             boundary #%i: %s',i,self.boundaries{i}.name));
			end

		end % }}}
		function boolean=isnameany(self,varargin) % {{{
			boolean=0;
			for  i=1:length(varargin),
				if strcmpi(self.name,varargin{i}), 
					boolean=1;
					break;
				end
			end
		end % }}}
		function boolean=iscontinentany(self,varargin) % {{{
			boolean=0;
			for  i=1:length(varargin),
				if strcmpi(self.continent,varargin{i}), 
					boolean=1;
					break;
				end
			end
		end % }}}
		function output=outputname(self,varargin) % {{{
		
			%recover options
			options=pairoptions(varargin{:});
			extension=getfieldvalue(options,'extension',1);

			[path,name,ext]=fileparts(self.name);
			if extension,
				output=[name ext];
			else
				output=name;
			end
		end % }}}
		function plot(self,varargin) % {{{
	
			%add option: 
			for i=1:length(self.boundaries),
				self.boundaries{i}.plot('proj',self.proj,varargin{:});
			end

		end % }}}
		function plot3d(self,varargin) % {{{
	
			%add option: 
			for i=1:length(self.boundaries),
				self.boundaries{i}.plot3d(varargin{:});
			end

		end % }}}
		function out=contour(self,varargin) % {{{
		
			%recover options
			options=pairoptions(varargin{:});
			x=[];
			y=[];

			%go through boundaries, recover edges and project them in the basin epsg reference frame: 
			for i=1:length(self.boundaries),
				boundary=self.boundaries{i};
				contour=boundary.edges();
				[contour.x,contour.y]=CoordTransform(contour.x,contour.y,boundary.proj,self.proj);
				x=[x;contour.x];
				y=[y;contour.y];
			end

			%close the contour: 
			if x(end)~=x(1) | y(end)~=y(1), 
				x(end+1)=x(1); y(end+1)=y(1);
			end

			out.x=x;
			out.y=y;
			out.nods=length(x);
		end % }}}
		function checkconsistency(self,varargin) % {{{
		
			%recover options
			options=pairoptions(varargin{:});

			%figure out if two boundaries are identical: 
			for i=1:length(self.boundaries),
				namei=self.boundaries{i}.shpfilename; 
				for j=i+1:length(self.boundaries),
					namej=self.boundaries{j}.shpfilename; 
					if strcmpi(namei,namej),
						error(['Basin ' self.name ' has two identical boundaries named ' namei ]);
					end
				end
			end

			%go through boundaries and check their consistency: 
			for i=1:length(self.boundaries),
				boundary=self.boundaries{i};
				boundary.checkconsistency();
			end
		end % }}}
		function contourplot(self,varargin) % {{{
			contour=self.contour();
			plot(contour.x,contour.y,'r*-');

		end % }}}
		function output=shapefilecrop(self,varargin) % {{{

			%recover options
			options=pairoptions(varargin{:});
			threshold=getfieldvalue(options,'threshold',.65); %.65 degrees lat,long
			inshapefile=getfieldvalue(options,'shapefile');
			outputshapefile=getfieldvalue(options,'outputshapefile','');

			if (exist(options,'epsgshapefile') & exist(options,'projshapefile')), 
				error('Basin shapefilecrop error message: cannot specify epsgshapefile and projshapefile at the same time');
			end
			if exist(options,'epsgshapefile'),
				projshapefile=epsg2proj(getfieldvalue(options,'epsgshapefile'));
			elseif exist(options,'projshapefile'),
				projshapefile=getfieldvalue(options,'projshapefile');
			else
				error('Basin shapefilecrop error message: epsgshapefile or projshapefile should be specified');
			end

			%create list of contours that have critical length > threshold (in lat,long):
			contours=shpread(inshapefile);
			llist=[];
			for i=1:length(contours),
				contour=contours(i);
				carea=polyarea(contour.x,contour.y);
				clength=sqrt(carea);
				if clength<threshold,
					llist=[llist;i];
				end
			end
			contours(llist)=[];

			%project onto reference frame:
			for i=1:length(contours),
				h=contours(i);
				[h.x,h.y]=CoordTransform(h.x,h.y,projshapefile,self.proj);
				contours(i).x=h.x;
				contours(i).y=h.y;
			end

			%only keep the contours that are inside the basin (take centroids): 
			ba=self.contour();
			flags=zeros(length(contours),1);
			for i=1:length(contours),
				h=contours(i); 
				isin=inpolygon(h.x,h.y,ba.x,ba.y);
				if ~isempty(find(isin==0)),
					flags(i)=1;
				end
			end

			pos=find(flags);
			contours(pos)=[];

			%Two options: 
			if strcmpi(outputshapefile,''),
				output=contours;
			else
				shpwrite(contours,outputshapefile);
			end

		end % }}}
	end
end
