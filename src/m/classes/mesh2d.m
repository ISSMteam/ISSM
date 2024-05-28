%MESH2D class definition
%
%   Usage:
%      mesh2d=mesh2d();

classdef mesh2d
	properties (SetAccess=public) 
		x                           = NaN;
		y                           = NaN;
		elements                    = NaN;
		numberofelements            = 0;
		numberofvertices            = 0;
		numberofedges               = 0;

		lat                         = NaN;
		long                        = NaN;
		epsg                        = 0;
		proj                        = 0;
		scale_factor                = NaN;

		vertexonboundary            = NaN;

		edges                       = NaN;
		segments                    = NaN;
		segmentmarkers              = NaN;
		vertexconnectivity          = NaN;
		elementconnectivity         = NaN;
		average_vertex_connectivity = 0;

		extractedvertices           = NaN;
		extractedelements           = NaN;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model selfect is
			% loaded. Update old properties here

			%2014 Oct. 1st
			if isstruct(self),
				oldself=self;
				%Assign property values from struct
				self=structtoobj(mesh2d(),oldself);
				if isfield(oldself,'hemisphere'),
					disp('md.mesh.hemisphere has been automatically converted to EPSG code');
					if strcmpi(oldself.hemisphere,'n'),
						self.epsg=3413;
						self.proj=epsg2proj(3413);
					else
						self.epsg=3031;
						self.proj=epsg2proj(3031);
					end
				end
			end

		end% }}}
	end
	methods
		function self = mesh2d(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=mesh2d();
					object=varargin{1};
					fields=fieldnames(object);
					for i=1:length(fields)
						field=fields{i};
						if ismember(field,properties('mesh2d')),
							self.(field)=object.(field);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%The connectivity is the average number of nodes linked to a given 
			%node through an edge. This connectivity is used to initially allocate 
			%memory to the stiffness matrix. A value of 16 seems to give a good 
			%memory/time ratio. This value can be checked in
			%test/Miscellaneous/runme.m
			self.average_vertex_connectivity=25;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			
			if strcmpi(solution,'LoveSolution'), return; end

			md = checkfield(md,'fieldname','mesh.x','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.y','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.elements','NaN',1,'Inf',1,'>',0,'values',1:md.mesh.numberofvertices);
			md = checkfield(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements 3]);
			if any(~ismember(1:md.mesh.numberofvertices,sort(unique(md.mesh.elements(:)))));
				md = checkmessage(md,'orphan nodes have been found. Check the mesh outline');
			end
			md = checkfield(md,'fieldname','mesh.numberofelements','>',0);
			md = checkfield(md,'fieldname','mesh.numberofvertices','>',0);
			md = checkfield(md,'fieldname','mesh.average_vertex_connectivity','>=',9,'message','''mesh.average_vertex_connectivity'' should be at least 9 in 2d');
			md = checkfield(md,'fieldname','mesh.segments','NaN',1,'Inf',1,'>',0,'size',[NaN 3]);
			if numel(md.mesh.scale_factor)>1,
				md = checkfield(md,'fieldname','mesh.scale_factor','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			end

			if strcmp(solution,'ThermalSolution')
					md = checkmessage(md,'thermal not supported for 2d mesh');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   2D tria Mesh (horizontal):')); 

			disp(sprintf('\n      Elements and vertices:'));
			fielddisplay(self,'numberofelements','number of elements');
			fielddisplay(self,'numberofvertices','number of vertices');
			fielddisplay(self,'elements','vertex indices of the mesh elements');
			fielddisplay(self,'x','vertices x coordinate [m]');
			fielddisplay(self,'y','vertices y coordinate [m]');
			fielddisplay(self,'edges','edges of the 2d mesh (vertex1 vertex2 element1 element2)');
			fielddisplay(self,'numberofedges','number of edges of the 2d mesh');

			disp(sprintf('\n      Properties:'));
			fielddisplay(self,'vertexonboundary','vertices on the boundary of the domain flag list');
			fielddisplay(self,'segments','edges on domain boundary (vertex1 vertex2 element)');
			fielddisplay(self,'segmentmarkers','number associated to each segment');
			fielddisplay(self,'vertexconnectivity','list of elements connected to vertex_i');
			fielddisplay(self,'elementconnectivity','list of elements adjacent to element_i');
			fielddisplay(self,'average_vertex_connectivity','average number of vertices connected to one vertex');

			disp(sprintf('\n      Extracted model:'));
			fielddisplay(self,'extractedvertices','vertices extracted from the model');
			fielddisplay(self,'extractedelements','elements extracted from the model');

			disp(sprintf('\n      Projection:'));
			fielddisplay(self,'lat','vertices latitude [degrees]');
			fielddisplay(self,'long','vertices longitude [degrees]');
			fielddisplay(self,'epsg','EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)');
			fielddisplay(self,'proj','PROJ.4 compatible projection string');
			fielddisplay(self,'scale_factor','Projection correction for volume, area, etc. computation)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.mesh.domain_type','data',['Domain' domaintype(self)],'format','String');
			WriteData(fid,prefix,'name','md.mesh.domain_dimension','data',dimension(self),'format','Integer');
			WriteData(fid,prefix,'name','md.mesh.elementtype','data',elementtype(self),'format','String');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','x','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','y','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'name','md.mesh.z','data',zeros(self.numberofvertices,1),'format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','elements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofelements','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofvertices','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','average_vertex_connectivity','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonboundary','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','segments','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','scale_factor','format','DoubleMat','mattype',1);
			if md.transient.isoceancoupling, %Need to add lat/long coordinates to couple with ocean as they rely on lat/long coordinate system
				WriteData(fid,prefix,'object',self,'class','mesh','fieldname','lat','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','mesh','fieldname','long','format','DoubleMat','mattype',1);
			end
		end % }}}
		function t = domaintype(self) % {{{
			t = '2Dhorizontal';
		end % }}}
		function d = dimension(self) % {{{
			d = 2;
		end % }}}
		function s = elementtype(self) % {{{
			s = 'Tria';
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.mesh.x'],self.x);
			writejs1Darray(fid,[modelname '.mesh.y'],self.y);
			writejs2Darray(fid,[modelname '.mesh.elements'],self.elements);
			writejsdouble(fid,[modelname '.mesh.numberofelements'],self.numberofelements);
			writejsdouble(fid,[modelname '.mesh.numberofvertices'],self.numberofvertices);
			writejsdouble(fid,[modelname '.mesh.numberofedges'],self.numberofedges);
			writejs1Darray(fid,[modelname '.mesh.lat'],self.lat);
			writejs1Darray(fid,[modelname '.mesh.long'],self.long);
			writejsdouble(fid,[modelname '.mesh.epsg'],self.epsg);
			writejsdouble(fid,[modelname '.mesh.proj'],self.proj);
			writejsdouble(fid,[modelname '.mesh.scale_factor'],self.scale_factor);
			writejs1Darray(fid,[modelname '.mesh.vertexonboundary'],self.vertexonboundary);
			writejs2Darray(fid,[modelname '.mesh.edges'],self.edges);
			writejs2Darray(fid,[modelname '.mesh.segments'],self.segments);
			writejs2Darray(fid,[modelname '.mesh.segmentmarkers'],self.segmentmarkers);
			writejs2Darray(fid,[modelname '.mesh.vertexconnectivity'],self.vertexconnectivity);
			writejs2Darray(fid,[modelname '.mesh.elementconnectivity'],self.elementconnectivity);
			writejsdouble(fid,[modelname '.mesh.average_vertex_connectivity'],self.average_vertex_connectivity);
			writejs1Darray(fid,[modelname '.mesh.extractedvertices'],self.extractedvertices);
			writejs1Darray(fid,[modelname '.mesh.extractedelements'],self.extractedelements);

		end % }}}
		function export(self,varargin) % {{{

			options=pairoptions(varargin{:});
			filename=getfieldvalue(options,'filename');
			format=getfieldvalue(options,'format','shp');
			geometry=getfieldvalue(options,'geometry','line');
			index=getfieldvalue(options,'index',[]);
			proj=getfieldvalue(options,'projection','');

			%prepare contours: 
			contours= struct([]);
			if strcmpi(geometry,'point'),
				for i=1:self.numberofvertices,
					contours(i).x = self.x(i);
					contours(i).y = self.y(i);
					contours(i).id = i;
					contours(i).Geometry = 'Point';
				end
			elseif strcmpi(geometry,'line'),
				counter=1;
				for i=1:self.numberofelements,
					el=self.elements(i,:);
					%first line:
					contours(counter).x = [self.x(el(1)) self.x(el(2))];
					contours(counter).y = [self.y(el(1)) self.y(el(2))];
					contours(counter).Geometry = 'Line';

					%second line:
					contours(counter+1).x = [self.x(el(2)) self.x(el(3))];
					contours(counter+1).y = [self.y(el(2)) self.y(el(3))];
					contours(counter+1).Geometry = 'Line';

					%third line:
					contours(counter+2).x = [self.x(el(3)) self.x(el(1))];
					contours(counter+2).y = [self.y(el(3)) self.y(el(1))];
					contours(counter+2).Geometry = 'Line';
					
					%increase counter: 
					counter=counter+3;
				end
			elseif strcmpi(geometry,'polygon'),
				if isempty(index),
					counter=1;
					for i=1:self.numberofelements,
						el=self.elements(i,:);
						contours(i).x=[self.x(el(1)) self.x(el(2)) self.x(el(3)) self.x(el(1))];
						contours(i).y=[self.y(el(1)) self.y(el(2)) self.y(el(3)) self.y(el(1))];
						contours(i).Geometry = 'Polygon';
						contours(i).Id = i;
					end
				else
					counter=1;
					for i=1:length(index),
						el=self.elements(index(i),:);
						contours(i).x=[self.x(el(1)) self.x(el(2)) self.x(el(3)) self.x(el(1))];
						contours(i).y=[self.y(el(1)) self.y(el(2)) self.y(el(3)) self.y(el(1))];
						contours(i).id = index(i);
						contours(i).Geometry = 'Polygon';
					end
				end
			else
				error(sprintf('mesh3dsurface ''export'' error message: geometry %s not supported yet (should be ''point'' or ''line''',geometry));
			end

			%write file: 
			if strcmpi(format,'shp'),
				shpwrite(contours,filename);
			elseif strcmpi(format,'exp'),
				expwrite(contours,filename);
			else
				error(sprintf('mesh3dsurface ''export'' error message: file format %s not supported yet',format));
			end

			%write projection file: 
			if ~isempty(proj),
				proj2shpprj(filename,proj);
			end
			%write style file: 
			applyqgisstyle(filename,'mesh');


		end % }}}
	end
end
