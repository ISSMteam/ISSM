%BAMGMESH class definition
%
%   Usage:
%      bamgmesh(varargin)

classdef bamgmesh
	properties (SetAccess=public) 
		% {{{
		Vertices=[];
		Edges=[];
		Triangles=[];
		IssmEdges=[];
		IssmSegments=[];
		VerticesOnGeomVertex=[];
		VerticesOnGeomEdge=[];
		EdgesOnGeomEdge=[];
		SubDomains=[];
		SubDomainsFromGeom=[];
		ElementConnectivity=[];
		NodalConnectivity=[];
		NodalElementConnectivity=[];
		CrackedVertices=[];
		CrackedEdges=[];
		PreviousNumbering=[];
		% }}}
	end
	methods
		function bg = bamgmesh(varargin)% {{{

		switch nargin
			case 0
				% if no input arguments, create a default object

			case 1

				bg=bamgmesh;
				object=varargin{1};
				fields=fieldnames(object);
				for i=1:length(fields)
					field=fields{i};
					if ismember(field,properties('bamgmesh')),
						bg.(field)=object.(field);
					end
				end

			otherwise
				error('bamgmesh constructor error message: unknown type of constructor call');
			end
		end%}}}
		function display(bm)% {{{
			disp(sprintf('\n%s = \n',inputname(1)));
			disp(struct(bm))
		end%}}}
	end
end
