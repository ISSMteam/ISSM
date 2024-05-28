//BAMGMESH class definition
//
//   Usage:
//      bamgmesh(varargin)

function bamgmesh(){
	//methods
	this.constructor = function(args) {// {{{
		//BAMGMESH - constructor for bamgmesh object
		//
		//   Usage:
		//      bamgmesh = bamgmesh(varargin)

		//initialize list
        switch (args.length) {
            case 0:
				//if no input arguments, create a default object
                break;
            case 1:
                var object = args[0];
                for (var field in object) {
                    if (object.hasOwnProperty(field)) {
                        this[field] = object[field];
                    }
                }
                break;
            default:
				throw Error('bamggeom constructor error message: unknown type of constructor call');
        }
	}// }}}
    this.disp= function(){// {{{
        disp(sprintf('\n%s = \n', 'bamgmesh'));
        disp(this);
    }// }}}

    //properties 
    // {{{
    this.Vertices                   = [];
    this.Edges                      = [];
    this.Triangles                  = [];
    this.IssmEdges                  = [];
    this.IssmSegments               = [];
    this.VerticesOnGeomVertex       = [];
    this.VerticesOnGeomEdge         = [];
    this.EdgesOnGeomEdge            = [];
    this.SubDomains                 = [];
    this.SubDomainsFromGeom         = [];
    this.ElementConnectivity        = [];
    this.NodalConnectivity          = [];
    this.NodalElementConnectivity   = [];
    this.CrackedVertices            = [];
    this.CrackedEdges               = [];
    this.PreviousNumbering          = [];

	this.constructor(arguments);
    //}}}
}
