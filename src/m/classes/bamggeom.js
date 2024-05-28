//BAMGGEOM class definition
//
//   Usage:
//      bamggeom(varargin)

function bamggeom(){
	//methods
	this.constructor = function(args) {// {{{
		//BAMGGEOM - constructor for bamggeom object
		//
		//   Usage:
		//      bamggeom = bamggeom(varargin)

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
        disp(sprintf('\n%s = \n', 'bamggeom'));
        disp(this);
    }// }}}

    //properties 
    // {{{
    this.Vertices           = [];
    this.Edges              = [];
    this.TangentAtEdges     = [];
    this.Corners            = [];
    this.RequiredVertices   = [];
    this.RequiredEdges      = [];
    this.CrackedEdges       = [];
    this.SubDomains         = [];

	this.constructor(arguments);
    //}}}
}
