import numpy as np
from fielddisplay import fielddisplay
from checkfield import *
import MatlabFuncs as m
from WriteData import WriteData


class mesh3dprisms(object):
    """
    MESH3DPRISMS class definition

       Usage:
          mesh3d = mesh3dprisms()
    """

    def __init__(self):  # {{{
        self.x = float('NaN')
        self.y = float('NaN')
        self.z = float('NaN')
        self.elements = float('NaN')
        self.numberoflayers = 0
        self.numberofelements = 0
        self.numberofvertices = 0

        self.lat = float('NaN')
        self.long = float('NaN')
        self.epsg = 0
        self.scale_factor = float('NaN')

        self.vertexonbase = float('NaN')
        self.vertexonsurface = float('NaN')
        self.lowerelements = float('NaN')
        self.lowervertex = float('NaN')
        self.upperelements = float('NaN')
        self.uppervertex = float('NaN')
        self.vertexonboundary = float('NaN')

        self.vertexconnectivity = float('NaN')
        self.elementconnectivity = float('NaN')
        self.average_vertex_connectivity = 0

        self.x2d = float('NaN')
        self.y2d = float('NaN')
        self.elements2d = float('NaN')
        self.numberofvertices2d = 0
        self.numberofelements2d = 0

        self.extractedvertices = float('NaN')
        self.extractedelements = float('NaN')

    #set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        string = "   3D prism Mesh:"

        string = "%s\n%s" % (string, "\n      Elements and vertices of the original 2d mesh3dprisms:")

        string = "%s\n%s" % (string, fielddisplay(self, "numberofelements2d", "number of elements"))
        string = "%s\n%s" % (string, fielddisplay(self, "numberofvertices2d", "number of vertices"))
        string = "%s\n%s" % (string, fielddisplay(self, "elements2d", "vertex indices of the mesh3dprisms elements"))
        string = "%s\n%s" % (string, fielddisplay(self, "x2d", "vertices x coordinate [m]"))
        string = "%s\n%s" % (string, fielddisplay(self, "y2d", "vertices y coordinate [m]"))

        string = "%s\n%s" % (string, "\n\n      Elements and vertices of the extruded 3d mesh3dprisms:")
        string = "%s\n%s" % (string, fielddisplay(self, "numberofelements", "number of elements"))
        string = "%s\n%s" % (string, fielddisplay(self, "numberofvertices", "number of vertices"))
        string = "%s\n%s" % (string, fielddisplay(self, "elements", "vertex indices of the mesh3dprisms elements"))
        string = "%s\n%s" % (string, fielddisplay(self, "x", "vertices x coordinate [m]"))
        string = "%s\n%s" % (string, fielddisplay(self, "y", "vertices y coordinate [m]"))
        string = "%s\n%s" % (string, fielddisplay(self, "z", "vertices z coordinate [m]"))

        string = "%s%s" % (string, "\n\n      Properties:")
        string = "%s\n%s" % (string, fielddisplay(self, "numberoflayers", "number of extrusion layers"))
        string = "%s\n%s" % (string, fielddisplay(self, "vertexonbase", "lower vertices flags list"))
        string = "%s\n%s" % (string, fielddisplay(self, "vertexonsurface", "upper vertices flags list"))
        string = "%s\n%s" % (string, fielddisplay(self, "uppervertex", "upper vertex list (NaN for vertex on the upper surface)"))
        string = "%s\n%s" % (string, fielddisplay(self, "upperelements", "upper element list (NaN for element on the upper layer)"))
        string = "%s\n%s" % (string, fielddisplay(self, "lowervertex", "lower vertex list (NaN for vertex on the lower surface)"))
        string = "%s\n%s" % (string, fielddisplay(self, "lowerelements", "lower element list (NaN for element on the lower layer)"))
        string = "%s\n%s" % (string, fielddisplay(self, "vertexonboundary", "vertices on the boundary of the domain flag list"))
        string = "%s\n%s" % (string, fielddisplay(self, "vertexconnectivity", "list of elements connected to vertex_i"))
        string = "%s\n%s" % (string, fielddisplay(self, "elementconnectivity", "list of elements adjacent to element_i"))
        string = "%s\n%s" % (string, fielddisplay(self, "average_vertex_connectivity", "average number of vertices connected to one vertex"))

        string = "%s%s" % (string, "\n\n      Extracted model:")
        string = "%s\n%s" % (string, fielddisplay(self, "extractedvertices", "vertices extracted from the model"))
        string = "%s\n%s" % (string, fielddisplay(self, "extractedelements", "elements extracted from the model"))

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
        self.average_vertex_connectivity = 25

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'mesh.x', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.y', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.z', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.elements', 'NaN', 1, 'Inf', 1, '>', 0, 'values', np.arange(1, md.mesh.numberofvertices + 1))
        md = checkfield(md, 'fieldname', 'mesh.elements', 'size', [md.mesh.numberofelements, 6])
        if np.any(np.logical_not(m.ismember(np.arange(1, md.mesh.numberofvertices + 1), md.mesh.elements))):
            md.checkmessage("orphan nodes have been found. Check the mesh3dprisms outline")
        md = checkfield(md, 'fieldname', 'mesh.numberoflayers', '>=', 0)
        md = checkfield(md, 'fieldname', 'mesh.numberofelements', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.numberofvertices', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.vertexonbase', 'size', [md.mesh.numberofvertices], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'mesh.vertexonsurface', 'size', [md.mesh.numberofvertices], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'mesh.average_vertex_connectivity', '>=', 24, 'message', "'mesh.average_vertex_connectivity' should be at least 24 in 3d")
        if(np.size(self.scale_factor) > 1):
            md = checkfield(md, 'fieldname', 'mesh.scale_factor', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])

        return md
    # }}}

    def domaintype(self):  # {{{
        return "3D"
    # }}}

    def dimension(self):  # {{{
        return 3
    # }}}

    def elementtype(self):  # {{{
        return "Penta"
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.mesh.domain_type', 'data', "Domain" + self.domaintype(), 'format', 'String')
        WriteData(fid, prefix, 'name', 'md.mesh.domain_dimension', 'data', self.dimension(), 'format', 'Integer')
        WriteData(fid, prefix, 'name', 'md.mesh.elementtype', 'data', self.elementtype(), 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'x', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'y', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'z', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'elements', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberoflayers', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofelements', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofvertices', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'vertexonbase', 'format', 'BooleanMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'vertexonsurface', 'format', 'BooleanMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'lowerelements', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'upperelements', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'average_vertex_connectivity', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'elements2d', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofvertices2d', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofelements2d', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'scale_factor', 'format', 'DoubleMat', 'mattype', 1)
        if md.transient.isoceancoupling:
            WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'lat', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'long', 'format', 'DoubleMat', 'mattype', 1)
    # }}}
