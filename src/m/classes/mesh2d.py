import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
import MatlabFuncs as m
from WriteData import WriteData


class mesh2d(object):
    """mesh2d class definition

    Usage:
        mesh2d = mesh2d()
    """

    def __init__(self):  # {{{
        self.x = np.nan
        self.y = np.nan
        self.elements = np.nan
        self.numberofelements = 0
        self.numberofvertices = 0
        self.numberofedges = 0

        self.lat = np.nan
        self.long = np.nan
        self.epsg = 0
        self.scale_factor = np.nan

        self.vertexonboundary = np.nan
        self.edges = np.nan
        self.segments = np.nan
        self.segmentmarkers = np.nan
        self.vertexconnectivity = np.nan
        self.elementconnectivity = np.nan
        self.average_vertex_connectivity = 0

        self.extractedvertices = np.nan
        self.extractedelements = np.nan

        # Set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        s = '   2D tria Mesh (horizontal):\n'

        s += '{}\n'.format('      Elements and vertices:')
        s += '{}\n'.format(fielddisplay(self, 'numberofelements', 'number of elements'))
        s += '{}\n'.format(fielddisplay(self, 'numberofvertices', 'number of vertices'))
        s += '{}\n'.format(fielddisplay(self, 'elements', 'vertex indices of the mesh elements'))
        s += '{}\n'.format(fielddisplay(self, 'x', 'vertices x coordinate [m]'))
        s += '{}\n'.format(fielddisplay(self, 'y', 'vertices y coordinate [m]'))
        s += '{}\n'.format(fielddisplay(self, 'edges', 'edges of the 2d mesh (vertex1 vertex2 element1 element2)'))
        s += '{}\n'.format(fielddisplay(self, 'numberofedges', 'number of edges of the 2d mesh'))
        s += '\n'
        s += '{}\n'.format('      Properties:')
        s += '{}\n'.format(fielddisplay(self, 'vertexonboundary', 'vertices on the boundary of the domain flag list'))
        s += '{}\n'.format(fielddisplay(self, 'segments', 'edges on domain boundary (vertex1 vertex2 element)'))
        s += '{}\n'.format(fielddisplay(self, 'segmentmarkers', 'number associated to each segment'))
        s += '{}\n'.format(fielddisplay(self, 'vertexconnectivity', 'list of elements connected to vertex_i'))
        s += '{}\n'.format(fielddisplay(self, 'elementconnectivity', 'list of elements adjacent to element_i'))
        s += '{}\n'.format(fielddisplay(self, 'average_vertex_connectivity', 'average number of vertices connected to one vertex'))
        s += '\n'
        s += '{}\n'.format('      Extracted model:')
        s += '{}\n'.format(fielddisplay(self, 'extractedvertices', 'vertices extracted from the model'))
        s += '{}\n'.format(fielddisplay(self, 'extractedelements', 'elements extracted from the model'))
        s += '\n'
        s += '{}\n'.format('      Projection:')
        s += '{}\n'.format(fielddisplay(self, 'lat', 'vertices latitude [degrees]'))
        s += '{}\n'.format(fielddisplay(self, 'long', 'vertices longitude [degrees]'))
        s += '{}\n'.format(fielddisplay(self, 'epsg', 'EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)'))
        s += '{}\n'.format(fielddisplay(self, 'scale_factor', 'Projection correction for volume, area, etc. computation'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # The connectivity is the average number of nodes linked to a given 
        # node through an edge. This connectivity is used to initially allocate 
        # memory to the stiffness matrix. A value of 16 seems to give a good 
        # memory/time ratio. This value can be checked in 
        # test/Miscellaneous/runme.m
        self.average_vertex_connectivity = 25

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
            md.checkmessage('orphan nodes have been found. Check the mesh outline')
        md = checkfield(md, 'fieldname', 'mesh.numberofelements', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.numberofvertices', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.average_vertex_connectivity', '>=', 9, 'message', '\'mesh.average_vertex_connectivity\' should be at least 9 in 2d')
        md = checkfield(md, 'fieldname', 'mesh.segments', 'NaN', 1, 'Inf', 1, '>', 0, 'size', [np.nan, 3])
        if(np.size(self.scale_factor) > 1):
            md = checkfield(md, 'fieldname', 'mesh.scale_factor', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])

        if solution == 'ThermalSolution':
            md.checkmessage('thermal not supported for 2d mesh')

        return md
    # }}}

    def domaintype(self):  # {{{
        return '2Dhorizontal'
    # }}}

    def dimension(self):  # {{{
        return 2
    # }}}

    def elementtype(self):  # {{{
        return 'Tria'
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.mesh.domain_type', 'data', 'Domain' + self.domaintype(), 'format', 'String')
        WriteData(fid, prefix, 'name', 'md.mesh.domain_dimension', 'data', self.dimension(), 'format', 'Integer')
        WriteData(fid, prefix, 'name', 'md.mesh.elementtype', 'data', self.elementtype(), 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'x', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'y', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'name', 'md.mesh.z', 'data', np.zeros(self.numberofvertices), 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'elements', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofelements', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'numberofvertices', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'average_vertex_connectivity', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'vertexonboundary', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'segments', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'scale_factor', 'format', 'DoubleMat', 'mattype', 1)
        if md.transient.isoceancoupling:
            WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'lat', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'mesh', 'fieldname', 'long', 'format', 'DoubleMat', 'mattype', 1)
    # }}}
