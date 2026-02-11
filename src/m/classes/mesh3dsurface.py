import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from MatlabFuncs import *
from model import *
from WriteData import WriteData


class mesh3dsurface(object):
    """MESH3DSURFACE class definition

    Usage:
        mesh3dsurface = mesh3dsurface()
    """

    def __init__(self, *args):  # {{{
        self.x = np.nan
        self.y = np.nan
        self.z = np.nan
        self.elements = np.nan
        self.numberofelements = 0
        self.numberofvertices = 0
        self.numberofedges = 0

        self.lat = np.nan
        self.long = np.nan
        self.r = np.nan

        self.vertexonboundary = np.nan
        self.edges = np.nan
        self.segments = np.nan
        self.segmentmarkers = np.nan
        self.vertexconnectivity = np.nan
        self.elementconnectivity = np.nan
        self.average_vertex_connectivity = 0

        self.extractedvertices = np.nan
        self.extractedelements = np.nan

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            self = mesh3dsurface()
            arg = args[1]
            fields = fieldnames(arg)
            for i in range(len(fields)):
                field = fields[i]
                if ismember(field, properties('mesh3dsurface')):
                    self.field = arg.field
        else:
            raise RuntimeError('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   2D tria Mesh (surface):'

        s += '\n      Elements and vertices:'
        s = "%s\n%s" % (s, fielddisplay(self, 'numberofelements', 'number of elements'))
        s = "%s\n%s" % (s, fielddisplay(self, 'numberofvertices', 'number of vertices'))
        s = "%s\n%s" % (s, fielddisplay(self, 'elements', 'vertex indices of the mesh elements'))
        s = "%s\n%s" % (s, fielddisplay(self, 'x', 'vertices x coordinate [m]'))
        s = "%s\n%s" % (s, fielddisplay(self, 'y', 'vertices y coordinate [m]'))
        s = "%s\n%s" % (s, fielddisplay(self, 'z', 'vertices z coordinate [m]'))
        s = "%s\n%s" % (s, fielddisplay(self, 'lat', 'vertices latitude [degrees]'))
        s = "%s\n%s" % (s, fielddisplay(self, 'long', 'vertices longitude [degrees]'))
        s = "%s\n%s" % (s, fielddisplay(self, 'r', 'vertices radius [m]'))

        s = "%s\n%s" % (s, fielddisplay(self, 'edges', 'edges of the 2d mesh (vertex1 vertex2 element1 element2)'))
        s = "%s\n%s" % (s, fielddisplay(self, 'numberofedges', 'number of edges of the 2d mesh'))

        s += '\n      Properties:'
        s = "%s\n%s" % (s, fielddisplay(self, 'vertexonboundary', 'vertices on the boundary of the domain flag list'))
        s = "%s\n%s" % (s, fielddisplay(self, 'segments', 'edges on domain boundary (vertex1 vertex2 element)'))
        s = "%s\n%s" % (s, fielddisplay(self, 'segmentmarkers', 'number associated to each segment'))
        s = "%s\n%s" % (s, fielddisplay(self, 'vertexconnectivity', 'list of elements connected to vertex_i'))
        s = "%s\n%s" % (s, fielddisplay(self, 'elementconnectivity', 'list of elements adjacent to element_i'))
        s = "%s\n%s" % (s, fielddisplay(self, 'average_vertex_connectivity', 'average number of vertices connected to one vertex'))

        s += '\n      Extracted model():'
        s = "%s\n%s" % (s, fielddisplay(self, 'extractedvertices', 'vertices extracted from the model()'))
        s = "%s\n%s" % (s, fielddisplay(self, 'extractedelements', 'elements extracted from the model()'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        #The connectivity is the average number of nodes linked to a given node 
        #through an edge. This connectivity is used to initially allocate 
        #memory to the stiffness matrix. A value of 16 seems to give a good 
        #memory/time ratio. This value can be checked in 
        #test/NightlyRun/runme.py.
        self.average_vertex_connectivity = 25
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'mesh.x', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.y', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.z', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.lat', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.long', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.r', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'mesh.elements', 'NaN', 1, 'Inf', 1, '>', 0, 'values', np.arange(1, md.mesh.numberofvertices + 1))
        md = checkfield(md, 'fieldname', 'mesh.elements', 'size', [md.mesh.numberofelements, 3])
        if hasattr(np, 'isin'): #Numpy 2017+
            tmp = np.isin(np.arange(1, md.mesh.numberofvertices + 1), md.mesh.elements.flat)
        else: #For backward compatibility
            tmp = np.in1d(np.arange(1, md.mesh.numberofvertices + 1), md.mesh.elements.flat)
        if np.any(np.logical_not(tmp)):
            md = md.checkmessage('orphan nodes have been found; check the mesh outline')

        md = checkfield(md, 'fieldname', 'mesh.numberofelements', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.numberofvertices', '>', 0)
        md = checkfield(md, 'fieldname', 'mesh.average_vertex_connectivity', '>=', 9, 'message', '"mesh.average_vertex_connectivity" should be at least 9 in 2d')

        if solution == 'ThermalSolution':
            md = md.checkmessage('thermal not supported for 2d mesh')

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.mesh.domain_type', 'data', 'Domain' + self.domaintype(), 'format', 'String')
        WriteData(fid, prefix, 'name', 'md.mesh.domain_dimension', 'data', self.dimension(), 'format', 'Integer')
        WriteData(fid, prefix, 'name', 'md.mesh.elementtype', 'data', self.elementtype(), 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'x', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'y', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'z', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'lat', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'long', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'r', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'name', 'md.mesh.z', 'data', np.zeros(md.mesh.numberofvertices), 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'elements', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'numberofelements', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'numberofvertices', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'average_vertex_connectivity', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vertexonboundary', 'format', 'DoubleMat', 'mattype', 1)
    # }}}

    def domaintype(self):  # {{{
        return '3Dsurface'
    # }}}

    def dimension(self):  # {{{
        return 2
    # }}}

    def elementtype(self):  # {{{
        return 'Tria'
    # }}}

    def processmesh(self, options):  # {{{
        isplanet = 1
        is2d = 0

        elements = self.elements
        x = self.x
        y = self.y
        z = self.z

        return x, y, z, elements, is2d, isplanet
    # }}}

    def savemodeljs(self, fid, modelname):  # {{{
        fid.write('  #s.mesh = new mesh3dsurface()\n' % modelname)
        writejs1Darray(fid, [modelname, '.mesh.x'], self.x)
        writejs1Darray(fid, [modelname, '.mesh.y'], self.y)
        writejs1Darray(fid, [modelname, '.mesh.z'], self.z)
        writejs2Darray(fid, [modelname, '.mesh.elements'], self.elements)
        writejsdouble(fid, [modelname, '.mesh.numberofelements'], self.numberofelements)
        writejsdouble(fid, [modelname, '.mesh.numberofvertices'], self.numberofvertices)
        writejsdouble(fid, [modelname, '.mesh.numberofedges'], self.numberofedges)
        writejs1Darray(fid, [modelname, '.mesh.lat'], self.lat)
        writejs1Darray(fid, [modelname, '.mesh.long'], self.long)
        writejs1Darray(fid, [modelname, '.mesh.r'], self.r)
        writejs1Darray(fid, [modelname, '.mesh.vertexonboundary'], self.vertexonboundary)
        writejs2Darray(fid, [modelname, '.mesh.edges'], self.edges)
        writejs2Darray(fid, [modelname, '.mesh.segments'], self.segments)
        writejs2Darray(fid, [modelname, '.mesh.segmentmarkers'], self.segmentmarkers)
        writejs2Darray(fid, [modelname, '.mesh.vertexconnectivity'], self.vertexconnectivity)
        writejs2Darray(fid, [modelname, '.mesh.elementconnectivity'], self.elementconnectivity)
        writejsdouble(fid, [modelname, '.mesh.average_vertex_connectivity'], self.average_vertex_connectivity)
        writejs1Darray(fid, [modelname, '.mesh.extractedvertices'], self.extractedvertices)
        writejs1Darray(fid, [modelname, '.mesh.extractedelements'], self.extractedelements)
    # }}}

    def export(self, *args):  # {{{
        options = pairoptions(*args)

        filename    = options.getfieldvalue('filename')
        format      = options.getfieldvalue('format', 'shp')
        geometry    = options.getfieldvalue('geometry', 'line')
        index       = options.getfieldvalue('index', [])
        proj        = options.getfieldvalue('projection', '')

        #prepare contours:
        contours = []
        if geometry == 'point':
            for i in range(len(self.numberofvertices)):
                contour             = OrderedStruct()
                contour.x           = self.long[i]
                contour.y           = self.lat[i]
                contour.id          = i
                contour.Geometry    = 'Point'
                contours[i]         = contour
        elif geometry == 'line':
            count = 0
            for i in range(len(self.numberofelements)): 
                el = self.elements[i]

                #first line:
                contour             = OrderedStruct()
                contour.x           = [self.long[el[0]], self.long[el[1]]]
                contour.y           = [self.lat[el[0]], self.lat[el[1]]]
                contour.Geometry    = 'Line' 
                contours[count]     = contour

                #second line:
                contour             = OrderedStruct()
                contour.x           = [self.long[el[1]], self.long[el[2]]]
                contour.y           = [self.lat[el[1]], self.lat[el[2]]]
                contour.Geometry    = 'Line'
                contours[count + 1] = contour

                #second line:
                contour             = OrderedStruct()
                contour.x           = [self.long[el[2]], self.long[el[0]]]
                contour.y           = [self.lat[el[2]], self.lat[el[0]]]
                contour.Geometry    = 'Line'
                contours[count + 2] = contour

                #increase count:
                count += 3
        elif geometry == 'polygon':
            # TODO: Refactor the following to reduce repeated code, or 
            #       leave as is because it is more readable?
            if index == []:
                for i in range(len(self.numberofelements)):
                    el = self.elements[i]

                    contour             = OrderedStruct()
                    contour.x           = [self.long[el[0]], self.long[el[1]], self.long[el[2]]]
                    contour.y           = [self.lat[el[0]], self.lat[el[1]], self.lat[el[2]]]
                    contour.Id          = i
                    contour.Geometry    = 'Polygon'
                    contours[i]         = contour
            else:
                for i in range(len(index)):
                    el = self.elements[index[i]]

                    contour             = OrderedStruct()
                    contour.x           = [self.long[el[0]], self.long[el[1]], self.long[el[2]]]
                    contour.y           = [self.lat[el[0]], self.lat[el[1]], self.lat[el[2]]]
                    contour.Id          = index[i]
                    contour.Geometry    = 'Polygon'
                    contours[i]         = contour
        else:
            raise RuntimeError("mesh3dsurface 'export' error message: geometry %s not supported yet (should be 'point' or 'line' or 'polygon')" % geometry)

        #write file:
        if format == 'shp':
            shpwrite(contours, filename)
        elif format == 'exp':
            expwrite(contours, filename)
        else:
            raise RuntimeError("mesh3dsurface 'export' error message: file format %s not supported yet" % format)

        #write projection file:
        if proj != '':
            proj2shpprj(filename, proj)

        #write style file:
        applyqgisstyle(filename, 'mesh')
    # }}}
