import numpy as np


class bamggeom(object):
    """
    BAMGGEOM class definition

       Usage:
          bamggeom(varargin)
    """

    def __init__(self, *args):  # {{{
        self.Vertices = np.empty((0, 3))
        self.Edges = np.empty((0, 3))
        self.TangentAtEdges = np.empty((0, 4))
        self.Corners = np.empty((0, 1))
        self.RequiredVertices = np.empty((0, 1))
        self.RequiredEdges = np.empty((0, 1))
        self.CrackedEdges = np.empty((0, 0))
        self.SubDomains = np.empty((0, 4))

        if not len(args):
            # if no input arguments, create a default object
            pass

        elif len(args) == 1:
            object = args[0]
            for field in list(object.keys()):
                if field in vars(self):
                    setattr(self, field, object[field])

        else:
            raise TypeError("bamggeom constructor error message: unknown type of constructor call")
    # }}}

    def __repr__(self):  # {{{
        s = "class '%s' object '%s'=\n" % (type(self), 'self')
        s += "    Vertices: %s\n" % str(self.Vertices)
        s += "    Edges: %s\n" % str(self.Edges)
        s += "    TangentAtEdges: %s\n" % str(self.TangentAtEdges)
        s += "    Corners: %s\n" % str(self.Corners)
        s += "    RequiredVertices: %s\n" % str(self.RequiredVertices)
        s += "    RequiredEdges: %s\n" % str(self.RequiredEdges)
        s += "    CrackedEdges: %s\n" % str(self.CrackedEdges)
        s += "    SubDomains: %s\n" % str(self.SubDomains)
        return s
    # }}}
