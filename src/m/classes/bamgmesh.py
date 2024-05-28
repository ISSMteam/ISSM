import numpy as np


class bamgmesh(object):
    """BAMGMESH class definition

    Usage:
        bamgmesh(varargin)
    """

    def __init__(self, *args):  # {{{
        self.Vertices = np.empty((0, 3))
        self.Edges = np.empty((0, 3))
        self.Triangles = np.empty((0, 0))
        self.IssmEdges = np.empty((0, 0))
        self.IssmSegments = np.empty((0, 0))
        self.VerticesOnGeomVertex = np.empty((0, 0))
        self.VerticesOnGeomEdge = np.empty((0, 0))
        self.EdgesOnGeomEdge = np.empty((0, 0))
        self.SubDomains = np.empty((0, 4))
        self.SubDomainsFromGeom = np.empty((0, 0))
        self.ElementConnectivity = np.empty((0, 0))
        self.NodalConnectivity = np.empty((0, 0))
        self.NodalElementConnectivity = np.empty((0, 0))
        self.CrackedVertices = np.empty((0, 0))
        self.CrackedEdges = np.empty((0, 0))

        if not len(args):
            # if no input arguments, create a default object
            pass

        elif len(args) == 1:
            object = args[0]
            for field in list(object.keys()):
                if field in vars(self):
                    setattr(self, field, object[field])

        else:
            raise TypeError("bamgmesh constructor error message: unknown type of constructor call")
    # }}}

    def __repr__(self):  # {{{
        s = "class '%s' object '%s' = \n" % (type(self), 'self')
        s += "    Vertices: %s\n" % str(self.Vertices)
        s += "    Edges: %s\n" % str(self.Edges)
        s += "    Triangles: %s\n" % str(self.Triangles)
        s += "    IssmEdges: %s\n" % str(self.IssmEdges)
        s += "    IssmSegments: %s\n" % str(self.IssmSegments)
        s += "    VerticesOnGeomVertex: %s\n" % str(self.VerticesOnGeomVertex)
        s += "    VerticesOnGeomEdge: %s\n" % str(self.VerticesOnGeomEdge)
        s += "    EdgesOnGeomEdge: %s\n" % str(self.EdgesOnGeomEdge)
        s += "    SubDomains: %s\n" % str(self.SubDomains)
        s += "    SubDomainsFromGeom: %s\n" % str(self.SubDomainsFromGeom)
        s += "    ElementConnectivity: %s\n" % str(self.ElementConnectivity)
        s += "    NodalConnectivity: %s\n" % str(self.NodalConnectivity)
        s += "    NodalElementConnectivity: %s\n" % str(self.NodalElementConnectivity)
        s += "    CrackedVertices: %s\n" % str(self.CrackedVertices)
        s += "    CrackedEdges: %s\n" % str(self.CrackedEdges)
        return s
    # }}}
