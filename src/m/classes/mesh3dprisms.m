%MESH3DPRISMS class definition
%
%   Usage:
%      mesh=mesh3dprisms();

classdef mesh3dprisms
	properties (SetAccess=public) 
		x                           = NaN;
		y                           = NaN;
		z                           = NaN;
		elements                    = NaN;
		numberoflayers              = 0;
		numberofelements            = 0;
		numberofvertices            = 0;

		lat                         = NaN;
		long                        = NaN;
		epsg                        = 0;
		proj                        = '';
		scale_factor                = NaN;

		vertexonbase                = NaN;
		vertexonsurface             = NaN;
		lowerelements               = NaN;
		lowervertex                 = NaN;
		upperelements               = NaN;
		uppervertex                 = NaN;
		vertexonboundary            = NaN;

		vertexconnectivity          = NaN;
		elementconnectivity         = NaN;
		average_vertex_connectivity = 0;

		x2d                         = NaN;
		y2d                         = NaN;
		elements2d                  = NaN;
		segments2d                  = NaN;
		numberofvertices2d          = 0;
		numberofelements2d          = 0;

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
				self=structtoobj(mesh3dprisms(),oldself);
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
		function self = mesh3dprisms(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=mesh3dprisms();
					object=varargin{1};
					fields=fieldnames(object);
					for i=1:length(fields)
						field=fields{i};
						if ismember(field,properties('mesh3dprisms')),
							self.(field)=object.(field);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%the connectivity is the averaged number of nodes linked to a
			%given node through an edge. This connectivity is used to initially
			%allocate memory to the stiffness matrix. A value of 16 seems to
			%give a good memory/time ration. This value can be checked in
			%trunk/test/Miscellaneous/runme.m
			self.average_vertex_connectivity=25;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			
			if strcmpi(solution,'LoveSolution'), return; end

			md = checkfield(md,'fieldname','mesh.x','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.y','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.z','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mesh.elements','NaN',1,'Inf',1,'>',0,'values',1:md.mesh.numberofvertices);
			md = checkfield(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements 6]);
			if any(~ismember(1:md.mesh.numberofvertices,sort(unique(md.mesh.elements(:)))));
				md = checkmessage(md,'orphan nodes have been found. Check the mesh outline');
			end
			md = checkfield(md,'fieldname','mesh.numberoflayers','>=',0);
			md = checkfield(md,'fieldname','mesh.numberofelements','>',0);
			md = checkfield(md,'fieldname','mesh.numberofvertices','>',0);
			md = checkfield(md,'fieldname','mesh.vertexonboundary','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','mesh.vertexonbase','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','mesh.vertexonsurface','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','mesh.average_vertex_connectivity','>=',24,'message','''mesh.average_vertex_connectivity'' should be at least 24 in 3d');

			%Check that mesh follows the geometry
			md = checkfield(md,'fieldname','mesh.z','>=',md.geometry.base-10^-10,'message','''mesh.z'' lower than bedrock');
			md = checkfield(md,'fieldname','mesh.z','<=',md.geometry.surface+10^-10,'message','''mesh.z'' higher than surface elevation');
			if any(max(abs(project2d(md,md.mesh.z,1)-project2d(md,md.geometry.base,1)))>1e-10),
				md = checkmessage(md,'md.mesh.z is not consistent with md.geometry.base, you changed the geometry after extrusion');
			end
			if any(max(abs(project2d(md,md.mesh.z,md.mesh.numberoflayers)-project2d(md,md.geometry.surface,1)))>1e-10),
				md = checkmessage(md,'md.mesh.z is not consistent with md.geometry.surface, you changed the geometry after extrusion !!');
			end
			if any(max(abs(project2d(md,md.mesh.z,md.mesh.numberoflayers)-project2d(md,md.mesh.z,1) - project2d(md,md.geometry.thickness,1)))>1e-10),
				md = checkmessage(md,'md.mesh.z is not consistent with md.geometry.thickness, you changed the geometry after extrusion !!');
			end
			if numel(md.mesh.scale_factor)>1,
				md = checkfield(md,'fieldname','mesh.scale_factor','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   3D prism Mesh:')); 

			disp(sprintf('\n      Elements and vertices of the original 2d mesh:'));
			fielddisplay(self,'numberofelements2d','number of elements');
			fielddisplay(self,'numberofvertices2d','number of vertices');
			fielddisplay(self,'elements2d','vertex indices of the mesh elements');
			fielddisplay(self,'segments2d','edges on 2d domain boundary (vertex1 vertex2 element)');
			fielddisplay(self,'x2d','vertices x coordinate [m]');
			fielddisplay(self,'y2d','vertices y coordinate [m]');

			disp(sprintf('\n      Elements and vertices of the extruded 3d mesh:'));
			fielddisplay(self,'numberofelements','number of elements');
			fielddisplay(self,'numberofvertices','number of vertices');
			fielddisplay(self,'elements','vertex indices of the mesh elements');
			fielddisplay(self,'x','vertices x coordinate [m]');
			fielddisplay(self,'y','vertices y coordinate [m]');
			fielddisplay(self,'z','vertices z coordinate [m]');

			disp(sprintf('\n      Properties:'));
			fielddisplay(self,'numberoflayers','number of extrusion layers');
			fielddisplay(self,'vertexonbase','lower vertices flags list');
			fielddisplay(self,'vertexonsurface','upper vertices flags list');
			fielddisplay(self,'uppervertex','upper vertex list (NaN for vertex on the upper surface)');
			fielddisplay(self,'upperelements','upper element list (NaN for element on the upper layer)');
			fielddisplay(self,'lowervertex','lower vertex list (NaN for vertex on the lower surface)');
			fielddisplay(self,'lowerelements','lower element list (NaN for element on the lower layer');
			fielddisplay(self,'vertexonboundary','vertices on the boundary of the domain flag list');

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
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','z','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','elements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberoflayers','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofelements','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofvertices','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonboundary','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonbase','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonsurface','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','lowerelements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','upperelements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','average_vertex_connectivity','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','elements2d','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','segments2d','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofvertices2d','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofelements2d','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','scale_factor','format','DoubleMat','mattype',1);
			if md.transient.isoceancoupling,
				WriteData(fid,prefix,'object',self,'class','mesh','fieldname','lat','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','mesh','fieldname','long','format','DoubleMat','mattype',1);
			end
		end % }}}
		function type = domaintype(self) % {{{
			type = '3D';
		end % }}}
		function d = dimension(self) % {{{
			d = 3;
		end % }}}
		function s = elementtype(self) % {{{
			s = 'Penta';
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			fprintf(fid,'%s.mesh=new mesh3dprisms();\n',modelname);
			writejs1Darray(fid,[modelname '.mesh.x'],self.x);
			writejs1Darray(fid,[modelname '.mesh.y'],self.y);
			writejs1Darray(fid,[modelname '.mesh.z'],self.z);
			writejs2Darray(fid,[modelname '.mesh.elements'],self.elements);
			writejsdouble(fid,[modelname '.mesh.numberoflayers'],self.numberoflayers);
			writejsdouble(fid,[modelname '.mesh.numberofelements'],self.numberofelements);
			writejsdouble(fid,[modelname '.mesh.numberofvertices'],self.numberofvertices);
			writejs1Darray(fid,[modelname '.mesh.lat'],self.lat);
			writejs1Darray(fid,[modelname '.mesh.long'],self.long);
			writejs1Darray(fid,[modelname '.mesh.epsg'],self.epsg);
			writejs1Darray(fid,[modelname '.mesh.proj'],self.proj);
			writejs1Darray(fid,[modelname '.mesh.vertexonbase'],self.vertexonbase);
			writejs1Darray(fid,[modelname '.mesh.vertexonsurface'],self.vertexonsurface);
			writejs1Darray(fid,[modelname '.mesh.lowerelements'],self.lowerelements);
			writejs1Darray(fid,[modelname '.mesh.upperelements'],self.upperelements);
			writejs1Darray(fid,[modelname '.mesh.uppervertex'],self.uppervertex);
			writejs1Darray(fid,[modelname '.mesh.vertexonboundary'],self.vertexonboundary);

			writejs2Darray(fid,[modelname '.mesh.vertexconnectivity'],self.vertexconnectivity);
			writejs2Darray(fid,[modelname '.mesh.elementconnectivity'],self.elementconnectivity);
			writejsdouble(fid,[modelname '.mesh.average_vertex_connectivity'],self.average_vertex_connectivity);
			
			writejs1Darray(fid,[modelname '.mesh.x2d'],self.x2d);
			writejs1Darray(fid,[modelname '.mesh.y2d'],self.y2d);
			writejs2Darray(fid,[modelname '.mesh.elements2d'],self.elements2d);
			writejs2Darray(fid,[modelname '.mesh.segments2d'],self.segments2d);
			writejsdouble(fid,[modelname '.mesh.numberofvertices2d'],self.numberofvertices2d);
			writejsdouble(fid,[modelname '.mesh.numberofelements2d'],self.numberofelements2d);

			writejs1Darray(fid,[modelname '.mesh.extractedvertices'],self.extractedvertices);
			writejs1Darray(fid,[modelname '.mesh.extractedelements'],self.extractedelements);

		end % }}}
	end
end
