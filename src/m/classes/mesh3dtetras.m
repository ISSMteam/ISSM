%MESH3DTETRAS class definition
%
%   Usage:
%      mesh=mesh3dtetras();

classdef mesh3dtetras
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
				self=structtoobj(mesh3dtetras(),oldself);
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
		function self = mesh3dtetras(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('mesh3dtetras');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
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
			md = checkfield(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements 4]);
			if any(~ismember(1:md.mesh.numberofvertices,sort(unique(md.mesh.elements(:)))));
				md = checkmessage(md,'orphan nodes have been found. Check the mesh outline');
			end
			md = checkfield(md,'fieldname','mesh.numberoflayers','>=',0);
			md = checkfield(md,'fieldname','mesh.numberofelements','>',0);
			md = checkfield(md,'fieldname','mesh.numberofvertices','>',0);
			md = checkfield(md,'fieldname','mesh.vertexonbase','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','mesh.vertexonsurface','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','mesh.z','>=',md.geometry.base-10^-10,'message','''mesh.z'' lower than bedrock');
			md = checkfield(md,'fieldname','mesh.z','<=',md.geometry.surface+10^-10,'message','''mesh.z'' higher than surface elevation');
			md = checkfield(md,'fieldname','mesh.average_vertex_connectivity','>=',24,'message','''mesh.average_vertex_connectivity'' should be at least 24 in 3d');
			if numel(md.mesh.scale_factor)>1,
				md = checkfield(md,'fieldname','mesh.scale_factor','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   3D tetra Mesh:')); 

			disp(sprintf('\n      Elements and vertices of the original 2d mesh:'));
			fielddisplay(self,'numberofelements2d','number of elements');
			fielddisplay(self,'numberofvertices2d','number of vertices');
			fielddisplay(self,'elements2d','vertex indices of the mesh elements');
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
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonbase','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonsurface','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','lowerelements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','upperelements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','scale_factor','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','average_vertex_connectivity','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','elements2d','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofvertices2d','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','numberofelements2d','format','Integer');
		end % }}}
		function t = domaintype(self) % {{{
			t = '3D';
		end % }}}
		function d = dimension(self) % {{{
			d = 3;
		end % }}}
		function s = elementtype(self) % {{{
			s = 'Tetra';
		end % }}}
	end
end
