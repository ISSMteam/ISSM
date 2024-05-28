classdef mesh
	properties (SetAccess=public) 
		x                           = NaN;
		y                           = NaN;
		z                           = NaN
		elements                    = NaN
		dimension                   = 0;
		numberoflayers              = 0;
		numberofelements            = 0;
		numberofvertices            = 0;
		numberofedges               = 0;

		lat                         = NaN
		long                        = NaN
		hemisphere                  = NaN

		elementonbed                = NaN
		elementonsurface            = NaN
		vertexonbed                 = NaN
		vertexonsurface             = NaN
		lowerelements               = NaN
		lowervertex                 = NaN
		upperelements               = NaN
		uppervertex                 = NaN
		vertexonboundary            = NaN

		edges                       = NaN
		segments                    = NaN
		segmentmarkers              = NaN
		vertexconnectivity          = NaN
		elementconnectivity         = NaN
		average_vertex_connectivity = 0;

		x2d                         = NaN
		y2d                         = NaN
		elements2d                  = NaN
		numberofvertices2d          = 0;
		numberofelements2d          = 0;

		extractedvertices           = NaN
		extractedelements           = NaN
	end
end
