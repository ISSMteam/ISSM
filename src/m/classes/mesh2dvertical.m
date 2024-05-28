%MESH2DVERTICAL class definition
%
%   Usage:
%      mesh2dvertical=mesh2dvertical();

classdef mesh2dvertical
	properties (SetAccess=public) 
		x                           = NaN;
		y                           = NaN;
		elements                    = NaN;
		numberofelements            = 0;
		numberofvertices            = 0;
		numberofedges               = 0;

		lat                         = NaN;
		long                        = NaN;
		epsg                        = NaN;
		proj                        = '';
		scale_factor                = NaN;

		vertexonboundary            = NaN;
		vertexonbase                = NaN;
		vertexonsurface             = NaN;

		edges                       = NaN;
		segments                    = NaN;
		segmentmarkers              = NaN;
		vertexconnectivity          = NaN;
		elementconnectivity         = NaN;
		average_vertex_connectivity = 0;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model selfect is
			% loaded. Update old properties here

			%2014 Oct. 1st
			if isstruct(self),
				oldself=self;
				%Assign property values from struct
				self=structtoobj(mesh2dvertical(),oldself);
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
		function self = mesh2dvertical(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('mesh2dvertical');
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
			md = checkfield(md,'fieldname','mesh.elements','NaN',1,'Inf',1,'>',0,'values',1:md.mesh.numberofvertices);
			md = checkfield(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements 3]);
			if any(~ismember(1:md.mesh.numberofvertices,sort(unique(md.mesh.elements(:)))));
				md = checkmessage(md,'orphan nodes have been found. Check the mesh outline');
			end
			md = checkfield(md,'fieldname','mesh.numberofelements','>',0);
			md = checkfield(md,'fieldname','mesh.numberofvertices','>',0);
			md = checkfield(md,'fieldname','mesh.vertexonbase','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','mesh.vertexonsurface','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','mesh.average_vertex_connectivity','>=',9,'message','''mesh.average_vertex_connectivity'' should be at least 9 in 2d');
			if numel(md.mesh.scale_factor)>1,
				md = checkfield(md,'fieldname','mesh.scale_factor','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			end

			if strcmp(solution,'ThermalSolution')
				md = checkmessage(md,'thermal not supported for 2d mesh');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   2D tria Mesh (vertical):')); 

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
			fielddisplay(self,'vertexonbase','vertices on the bed of the domain flag list');
			fielddisplay(self,'vertexonsurface','vertices on the surface of the domain flag list');
			fielddisplay(self,'segments','edges on domain boundary (vertex1 vertex2 element)');
			fielddisplay(self,'segmentmarkers','number associated to each segment');
			fielddisplay(self,'vertexconnectivity','list of elements connected to vertex_i');
			fielddisplay(self,'elementconnectivity','list of elements adjacent to element_i');
			fielddisplay(self,'average_vertex_connectivity','average number of vertices connected to one vertex');

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
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonbase','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','vertexonsurface','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','average_vertex_connectivity','format','Integer');
			WriteData(fid,prefix,'object',self,'class','mesh','fieldname','scale_factor','format','DoubleMat','mattype',1);
		end % }}}
		function t = domaintype(self) % {{{
			t = '2Dvertical';
		end % }}}
		function d = dimension(self) % {{{
			d = 2;
		end % }}}
		function s = elementtype(self) % {{{
			s = 'Tria';
		end % }}}
		function flags = vertexflags(self,value) % {{{
			flags = zeros(self.numberofvertices,1);
			pos   = self.segments(find(self.segmentmarkers==value),1:2);
			flags(pos) = 1;
		end % }}}
		function [data datatype] = processdata(self,md,data,options) % {{{

			%transpose data if necessary
			if (size(data,2) > size(data,1)),
				data=data';
			end
			datasize=size(data);

			%convert to double if necessary
			if ~isnumeric(data);
				disp('processdata info message: data is not numeric (logical?). Converted to double');
				data=double(data);
			end

			%check length
			if datasize(1)~=md.mesh.numberofvertices & datasize(1)~=md.mesh.numberofelements
				error('plotmodel error message: data not supported yet');
			end

			%quiver?
			if datasize(2)>1,
				datatype=3;
			end

			%smoothing?
			if exist(options,'smooth')
				data=averaging(md,data,getfieldvalue(options,'smooth'));
				datasize(1)=md.mesh.numberofvertices;
				%---> go to node data
			end

			%element data
			if (datasize(1)==md.mesh.numberofelements & datasize(2)==1),
				datatype=1;

				%Mask?
				if exist(options,'mask'),
					flags=getfieldvalue(options,'mask');
					pos=find(~flags);
					if length(flags)==md.mesh.numberofvertices,
						[pos2 dummy]=find(ismember(md.mesh.elements,pos));
						data(pos2,:)=NaN;
					elseif length(flags)==md.mesh.numberofelements
						data(pos,:)=NaN;
					else
						disp('plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements');
					end
				end

				%log?
				if exist(options,'log'),
					bounds=getfieldvalue(options,'caxis',[min(data(:)) max(data(:))]);
					data(find(data<bounds(1)))=bounds(1);
					if any(data<=0),
						error('Log option cannot be applied on negative values. Use caxis option (Rignot''s settings: [1.5 max(data)])');
					end
					pos=find(~isnan(data));
					data(pos)=log(data(pos))/log(getfieldvalue(options,'log'));
				end
			end

			%node data
			if (datasize(1)==md.mesh.numberofvertices & datasize(2)==1),
				datatype=2;

				%Mask?
				if exist(options,'mask'),
					flags=getfieldvalue(options,'mask');
					pos=find(~flags);
					if length(flags)==md.mesh.numberofvertices,
						data(pos,:)=NaN;
					elseif length(flags)==md.mesh.numberofelements
						data(md.mesh.elements(pos,:),:)=NaN;
					else
						disp('plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements');
					end
				end

				%log?
				if exist(options,'log'),
					%if any(data<=0),
					%	error('Log option cannot be applied on negative values. Use caxis option (Rignot''s settings: [1.5 max(data)])');
					%end
					data=log(data)/log(getfieldvalue(options,'log'));
				end
			end
		end % }}}
		function [x y z elements is2d isplanet] = processmesh(self,options) % {{{

			isplanet = 0;
			is2d     = 1;

			elements = self.elements;
			x        = self.x;
			y        = self.y;
			z        = zeros(self.numberofvertices,1);

			if exist(options,'xunit'),
				unit=getfieldvalue(options,'xunit');
				x=x*unit; % Apply to x only
			end
			if exist(options,'yunit'),
				unit=getfieldvalue(options,'yunit');
				x=x*unit; % Apply to x only
			end
		end % }}}
	end
end
