%MESH3DSURFACE class definition
%
%   Usage:
%      mesh3dsurface=mesh3dsurface();

classdef mesh3dsurface
	properties (SetAccess=public) 
		x                           = NaN;
		y                           = NaN;
		z                           = NaN;
		elements                    = NaN;
		numberofelements            = 0;
		numberofvertices            = 0;
		numberofedges               = 0;

		lat                         = NaN;
		long                        = NaN;
		r                           = NaN;
		area                        = NaN;

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
				self=structtoobj(mesh3dsurface(),oldself);
				%if isfield(oldself,'hemisphere'), %NOTE: not ever needed, epsg is not a field anymore.
				%	disp('md.mesh.hemisphere has been automatically converted to EPSG code');
				%	if strcmpi(oldself.hemisphere,'n'),
				%		self.epsg=3413;
				%	else
				%		self.epsg=3031;
				%	end
				%end
			end

		end% }}}
	end
	methods
		function self = mesh3dsurface(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=mesh3dsurface();
					object=varargin{1};
					fields=fieldnames(object);
					for i=1:length(fields)
						field=fields{i};
						if ismember(field,properties('mesh3dsurface')),
							self.(field)=object.(field);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function obj = setdefaultparameters(obj) % {{{

			%The connectivity is the average number of nodes linked to a given 
			%node through an edge. This connectivity is used to initially 
			%allocate memory to the stiffness matrix. A value of 16 seems to 
			%give a good memory/time ratio. This value can be checked in 
			%test/NightlyRun/runme.py.
			obj.average_vertex_connectivity=25;
		end % }}}
		function md = checkconsistency(obj,md,solution,analyses) % {{{

			if strcmpi(solution,'LoveSolution'), return; end

			md = checkfield(md,'fieldname','mesh.x','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.y','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.z','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.lat','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.long','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.r','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.elements','NaN',1,'Inf',1,'>',0,'values',1:md.mesh.numberofvertices);
			md = checkfield(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements 3]);
			if any(~ismember(1:md.mesh.numberofvertices,sort(unique(md.mesh.elements(:)))));
				md = checkmessage(md,'orphan nodes have been found; check the mesh outline');
			end
			md = checkfield(md,'fieldname','mesh.numberofelements','>',0);
			md = checkfield(md,'fieldname','mesh.numberofvertices','>',0);
			md = checkfield(md,'fieldname','mesh.average_vertex_connectivity','>=',9,'message','''mesh.average_vertex_connectivity'' should be at least 9 in 2d');

			if strcmp(solution,'ThermalSolution')
				md = checkmessage(md,'thermal not supported for 2d mesh');
			end
		end % }}}
		function disp(obj) % {{{
			disp(sprintf('   3D tria Mesh (surface):')); 

			disp(sprintf('\n      Elements and vertices:'));
			fielddisplay(obj,'numberofelements','number of elements');
			fielddisplay(obj,'numberofvertices','number of vertices');
			fielddisplay(obj,'elements','vertex indices of the mesh elements');
			fielddisplay(obj,'x','vertices x coordinate [m]');
			fielddisplay(obj,'y','vertices y coordinate [m]');
			fielddisplay(obj,'z','vertices z coordinate [m]');
			fielddisplay(obj,'lat','vertices latitude [degrees]');
			fielddisplay(obj,'long','vertices longitude [degrees]');
			fielddisplay(obj,'r','vertices radius [m]');
			fielddisplay(obj,'area','elemental areas [m^2]');
			
			fielddisplay(obj,'edges','edges of the 2d mesh (vertex1 vertex2 element1 element2)');
			fielddisplay(obj,'numberofedges','number of edges of the 2d mesh');

			disp(sprintf('\n      Properties:'));
			fielddisplay(obj,'vertexonboundary','vertices on the boundary of the domain flag list');
			fielddisplay(obj,'segments','edges on domain boundary (vertex1 vertex2 element)');
			fielddisplay(obj,'segmentmarkers','number associated to each segment');
			fielddisplay(obj,'vertexconnectivity','list of elements connected to vertex_i');
			fielddisplay(obj,'elementconnectivity','list of elements adjacent to element_i');
			fielddisplay(obj,'average_vertex_connectivity','average number of vertices connected to one vertex');

			disp(sprintf('\n      Extracted model:'));
			fielddisplay(obj,'extractedvertices','vertices extracted from the model');
			fielddisplay(obj,'extractedelements','elements extracted from the model'); 
		end % }}}
		function marshall(obj,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.mesh.domain_type','data',['Domain' domaintype(obj)],'format','String');
			WriteData(fid,prefix,'name','md.mesh.domain_dimension','data',dimension(obj),'format','Integer');
			WriteData(fid,prefix,'name','md.mesh.elementtype','data',elementtype(obj),'format','String');
			WriteData(fid,prefix,'object',obj,'fieldname','x','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',obj,'fieldname','y','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',obj,'fieldname','z','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',obj,'fieldname','lat','data',obj.lat,'format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',obj,'fieldname','long','data',obj.long,'format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',obj,'fieldname','r','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'name','md.mesh.z','data',zeros(obj.numberofvertices,1),'format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',obj,'fieldname','elements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',obj,'fieldname','numberofelements','format','Integer');
			WriteData(fid,prefix,'object',obj,'fieldname','numberofvertices','format','Integer');
			WriteData(fid,prefix,'object',obj,'fieldname','average_vertex_connectivity','format','Integer');
			WriteData(fid,prefix,'object',obj,'fieldname','vertexonboundary','format','DoubleMat','mattype',1);
		end % }}}
		function t = domaintype(obj) % {{{
			t = '3Dsurface';
		end % }}}
		function d = dimension(obj) % {{{
			d = 2;
		end % }}}
		function s = elementtype(obj) % {{{
			s = 'Tria';
		end % }}}
		function [x y z elements is2d isplanet] = processmesh(self,options) % {{{

			isplanet = 1;
			is2d     = 0;

			elements = self.elements;
			x        = self.x;
			y        = self.y;
			z        = self.z;
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			fprintf(fid,'%s.mesh=new mesh3dsurface();\n',modelname);
			writejs1Darray(fid,[modelname '.mesh.x'],self.x);
			writejs1Darray(fid,[modelname '.mesh.y'],self.y);
			writejs1Darray(fid,[modelname '.mesh.z'],self.z);
			writejs2Darray(fid,[modelname '.mesh.elements'],self.elements);
			writejsdouble(fid,[modelname '.mesh.numberofelements'],self.numberofelements);
			writejsdouble(fid,[modelname '.mesh.numberofvertices'],self.numberofvertices);
			writejsdouble(fid,[modelname '.mesh.numberofedges'],self.numberofedges);
			writejs1Darray(fid,[modelname '.mesh.lat'],self.lat);
			writejs1Darray(fid,[modelname '.mesh.long'],self.long);
			writejs1Darray(fid,[modelname '.mesh.r'],self.r);
			writejs1Darray(fid,[modelname '.mesh.area'],self.area);
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
					contours(i).x = self.long(i);
					contours(i).y = self.lat(i);
					contours(i).id = i;
					contours(i).Geometry = 'Point';
				end
			elseif strcmpi(geometry,'line'),
				count=1;
				for i=1:self.numberofelements,
					el=self.elements(i,:);
					%first line:
					contours(count).x = [self.long(el(1)) self.long(el(2))];
					contours(count).y = [self.lat(el(1)) self.lat(el(2))];
					contours(count).Geometry = 'Line';

					%second line:
					contours(count+1).x = [self.long(el(2)) self.long(el(3))];
					contours(count+1).y = [self.lat(el(2)) self.lat(el(3))];
					contours(count+1).Geometry = 'Line';

					%third line:
					contours(count+2).x = [self.long(el(3)) self.long(el(1))];
					contours(count+2).y = [self.lat(el(3)) self.lat(el(1))];
					contours(count+2).Geometry = 'Line';
					
					%increase count: 
					count = count+3;
				end
			elseif strcmpi(geometry,'polygon'),
				% TODO: Refactor the following to reduce repeated code, or 
				%		leave as is because it is more readable?
				if isempty(index),
					for i=1:self.numberofelements,
						el=self.elements(i,:);
						contours(i).x=[self.long(el(1)) self.long(el(2)) self.long(el(3)) self.long(el(1))];
						contours(i).y=[self.lat(el(1)) self.lat(el(2)) self.lat(el(3)) self.lat(el(1))];
						contours(i).Id = i;
						contours(i).Geometry = 'Polygon';
					end
				else
					for i=1:length(index),
						el=self.elements(index(i),:);
						contours(i).x=[self.long(el(1)) self.long(el(2)) self.long(el(3)) self.long(el(1))];
						contours(i).y=[self.lat(el(1)) self.lat(el(2)) self.lat(el(3)) self.lat(el(1))];
						contours(i).id = index(i);
						contours(i).Geometry = 'Polygon';
					end
				end
			else
				error(sprintf('mesh3dsurface ''export'' error message: geometry %s not supported yet (should be ''point'' or ''line'' or ''polygon'')',geometry));
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
