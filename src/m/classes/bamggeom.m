%BAMGGEOM class definition
%
%   Usage:
%      bamggeom(varargin)

classdef bamggeom
	properties (SetAccess=public) 
		% {{{
		Vertices=[];
		Edges=[];
		TangentAtEdges=[];
		Corners=[];
		RequiredVertices=[];
		RequiredEdges=[];
		CrackedEdges=[];
		SubDomains=[];
		% }}}
	end
	methods
		function bg = bamggeom(varargin)% {{{
		%BAMGGEOM - constructor for bamggeom object
		%
		%   Usage:
		%      bamggeom = bamggeom(varargin)

		switch nargin
			case 0
				% if no input arguments, create a default object

			case 1

				bg=bamggeom;
				object=varargin{1};
				fields=fieldnames(object);
				for i=1:length(fields)
					field=fields{i};
					if ismember(field,properties('bamggeom')),
						bg.(field)=object.(field);
					end
				end

			otherwise
				error('bamggeom constructor error message: unknown type of constructor call');
			end
		end%}}}
		function display(bg)% {{{
			disp(sprintf('\n%s = \n',inputname(1)));
			disp(struct(bg))
		end%}}}
	end
end
