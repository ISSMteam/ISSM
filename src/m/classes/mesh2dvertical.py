import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
import MatlabFuncs as m
from WriteData import WriteData


class mesh2dvertical(object):
    """
    MESH2DVERTICAL class definition

       Usage:
          mesh2dvertical = mesh2dvertical()
    """

    def __init__(self):  # {{{
        self.x = float('NaN')
        self.y = float('NaN')
        self.elements = float('NaN')
        self.numberofelements = 0
        self.numberofvertices = 0
        self.numberofedges = 0

        self.lat = float('NaN')
        self.long = float('NaN')
        self.epsg = float('NaN')
        self.scale_factor = float('NaN')

        self.vertexonboundary = float('NaN')
        self.vertexonbase = float('NaN')
        self.vertexonsurface = float('NaN')

        self.edges = float('NaN')
        self.segments = float('NaN')
        self.segmentmarkers = float('NaN')
        self.vertexconnectivity = float('NaN')
        self.elementconnectivity = float('NaN')
        self.average_vertex_connectivity = 0

    #set defaults
        self.setdefaultparameters()

    # }}}
    def __repr__(self):  # {{{
        string = "   2D tria Mesh (vertical):"

        string = "%s\n%s" % (string, "\n      Elements and vertices:")
        string = "%s\n%s" % (string, fielddisplay(self, "numberofelements", "number of elements"))
        string = "%s\n%s" % (string, fielddisplay(self, "numberofvertices", "number of vertices"))
        string = "%s\n%s" % (string, fielddisplay(self, "elements", "vertex indices of the mesh elements"))
        string = "%s\n%s" % (string, fielddisplay(self, "x", "vertices x coordinate [m]"))
        string = "%s\n%s" % (string, fielddisplay(self, "y", "vertices y coordinate [m]"))
        string = "%s\n%s" % (string, fielddisplay(self, "edges", "edges of the 2d mesh (vertex1 vertex2 element1 element2)"))
        string = "%s\n%s" % (string, fielddisplay(self, "numberofedges", "number of edges of the 2d mesh"))

        string = "%s%s" % (string, "\n\n      Properties:")
        string = "%s\n%s" % (string, fielddisplay(self, "vertexonboundary", "vertices on the boundary of the domain flag list"))
        string = "%s\n%s" % (string, fielddisplay(self, 'vertexonbase', 'vertices on the bed of the domain flag list'))
        string = "%s\n%s" % (string, fielddisplay(self, 'vertexonsurface', 'vertices on the surface of the domain flag list'))
        string = "%s\n%s" % (string, fielddisplay(self, "segments", "edges on domain boundary (vertex1 vertex2 element)"))
        string = "%s\n%s" % (string, fielddisplay(self, "segmentmarkers", "number associated to each segment"))
        string = "%s\n%s" % (string, fielddisplay(self, "vertexconnectivity", "list of elements connected to vertex_i"))
        string = "%s\n%s" % (string, fielddisplay(self, "elementconnectivity", "list of elements adjacent to element_i"))
        string = "%s\n%s" % (string, fielddisplay(self, "average_vertex_connectivity", "average number of vertices connected to one vertex"))

        string = "%s%s" % (string, "\n\n      Projection:")
        string = "%s\n%s" % (string, fielddisplay(self, "lat", "vertices latitude [degrees]"))
        string = "%s\n%s" % (string, fielddisplay(self, "long", "vertices longitude [degrees]"))
        string = "%s\n%s" % (string, fielddisplay(self, "epsg", "EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)"))
        string = "%s\n%s" % (string, fielddisplay(self, "scale_factor", "Projection correction for volume, area, etc. computation"))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        #the connectivity is the averaged number of nodes linked to a
        #given node through an edge. This connectivity is used to initially
        #allocate memory to the stiffness matrix. A value of 16 seems to
        #give a good memory / time ration. This value can be checked in
        #trunk / test / Miscellaneous / runme.m
        self.average_vertex_connectivity = 25.

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if(solution == 'LoveSolution'):
            return

        md = checkfield(md, 'fieldname', 'mesh.x', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.y', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.elements', 'NaN', 1, 'Inf', 1, '>', 0, 'values', np.arange(1, md.mesh.numberofvertices + 1))
        md = checkfield(md, 'fieldname', 'mesh.elements', 'size', [md.mesh.numberofelements, 3])
        if np.any(np.logical_not(m.ismember(np.arange(1, md.mesh.numberofvertices + 1), md.mesh.elements))):
            md.checkmessage("orphan nodes have been found. Check the mesh outline")
        md = checkfield(md, 'fieldname', 'mesh.numberofelements', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.numberofvertices', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.vertexonbase', 'size', [md.mesh.numberofvertices], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'mesh.vertexonsurface', 'size', [md.mesh.numberofvertices], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'mesh.average_vertex_connectivity', '>=', 9, 'message', "'mesh.average_vertex_connectivity' should be at least 9 in 2d")
        if(np.size(self.scale_factor) > 1):
            md = checkfield(md, 'fieldname', 'mesh.scale_factor', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])

        if solution == 'ThermalSolution':
            md.checkmessage("thermal not supported for 2d mesh")

        return md
    # }}}

    def domaintype(self):  # {{{
        return "2Dvertical"
    # }}}

    def dimension(self):  # {{{
        return 2
    # }}}

    def elementtype(self):  # {{{
        return "Tria"
    # }}}

    def vertexflags(self, value):  # {{{
        flags = np.zeros((self.numberofvertices, ))
        pos = self.segments[np.where(self.segmentmarkers == value), 0:2] - 1
        flags[pos] = 1
        return flags
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.mesh.domain_type', 'data', "Domain" + self.domaintype(), 'format', 'String')
        WriteData(fid, prefix, 'name', 'md.mesh.domain_dimension', 'data', self.dimension(), 'format', 'Integer')
        WriteData(fid, prefix, 'name', 'md.mesh.elementtype', 'data', self.elementtype(), 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'x', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'y', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'name', 'md.mesh.z', 'data', np.zeros(self.numberofvertices), 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'elements', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofelements', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofvertices', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'vertexonbase', 'format', 'BooleanMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'vertexonsurface', 'format', 'BooleanMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'average_vertex_connectivity', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'scale_factor', 'format', 'DoubleMat', 'mattype', 1)
    # }}}
