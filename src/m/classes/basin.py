from collections import OrderedDict
import math

from bamg import bamg
import matplotlib.path as path
import numpy as np

from boundary import boundary
from epsg2proj import epsg2proj
from fielddisplay import fielddisplay
from helpers import fileparts
from inpolygon import inpolygon
from laea import laea
from pairoptions import pairoptions
from polyarea import polyarea
from shpread import shpread


class basin(object):  # {{{
    """basin class definition

    Usage:
        basin = basin()
    """

    def __init__(self, *args):  # {{{
        self.boundaries = []
        self.name       = ''
        self.continent  = ''
        self.proj       = epsg2proj(4326)

        self.setdefaultparameters()

        if len(args):
            options = pairoptions(*args)

            #recover field values
            self.boundaries = options.getfieldvalue('boundaries', [])
            self.name = options.getfieldvalue('name', '')
            self.continent = options.getfieldvalue('continent', '')

            # figure out projection string
            if options.exist('epsg'):
                if options.exist('proj'):
                    raise RuntimeError('Error in constructor for basin %s. Cannot supply epsg and proj at the same time!'.format(self.name))
                epsg = options.getfieldvalue('epsg', '')
                proj = epsg2proj(epsg)#convert to PROJ.4 format
            elif options.exist('proj'):
                if options.exist('epsg'):
                    raise RuntimeError('Error in constructor for basin {}. Cannot supply proj and epsg at the same time!'.format(self.name))
                proj = options.getfieldvalue('proj', '')
            else:
                proj = '+proj=longlat +datum=WGS84 +no_defs'

            self.proj = proj
    # }}}

    def __repr__(self):  # {{{
        s = '   basin parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'continent', 'continent name'))
        s += '{}\n'.format(fielddisplay(self, 'name', 'basin name'))
        s += '{}\n'.format(fielddisplay(self, 'proj', 'basin projection string (follows PROJ.4 standard'))
        s += '{}\n'.format(fielddisplay(self, 'boundaries','list of boundary objects'))
        for i in range(len(self.boundaries)):
            s += '             boundary #{}: {}\n'.format((i + 1), self.boundaries[i])

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.name       = ''
        self.continent  = ''
        self.proj       = '+proj=longlat +datum=WGS84 +no_defs' #EPSG 4326
        self.boundaries = []

        return self
    # }}}

    def isnameany(self, *args):  # {{{
        boolean = 0
        for arg in args:
            if type(arg) in [np.ndarray, list]:
                for name in arg:
                    if name == self.name:
                        boolean = 1
                        break
            elif arg == self.name:
                boolean = 1
                break
        return boolean
    # }}}

    def iscontinentany(self, *args):  # {{{
        boolean = 0
        for arg in args:
            if type(arg) in [np.ndarray, list]:
                for continent in arg:
                    if continent == self.continent:
                        boolean = 1
                        break
            elif arg == self.continent:
                boolean = 1
                break
        return boolean
    # }}}

    def outputname(self, *args):  # {{{
        #recover options
        options = pairoptions(*args)
        extension = options.getfieldvalue('extension', 1)

        path, name, ext = fileparts(*args)
        if extension:
            output = '{}{}'.format(name, ext)
        else:
            output = name

        return output
    # }}}

    def plot(self, *args):  # {{{
        #add option
        for i in range(len(self.boundaries)):
            self.boundaries[i].plot('proj', self.proj, *args)
    # }}}

    def plot3d(self, *args):  # {{{
        #add option
        for i in range(len(self.boundaries)):
            self.boundaries[i].plot3d(*args)
    # }}}

    def contour(self, *args):  # {{{
        #recover options
        options = pairoptions(*args)
        x = []
        y = []

        #go through boundaries, recover edges and project them in the basin epsg reference frame
        for i in range(len(self.boundaries)):
            boundary = self.boundaries[i]
            contours = boundary.edges() # NOTE: Differs from MATLAB in that shpread returns a list
            for j in range(len(contours)):
                contours[j]['x'], contours[j]['y'] = CoordTransform(contours[j]['x'], contours[j]['y'], boundary.proj, self.proj)
                if x == []:
                    x = np.array(contours[j]['x'])
                    y = np.array(contours[j]['y'])
                else:
                    x = np.concatenate((x, contours[j]['x']), axis=0)
                    y = np.concatenate((y, contours[j]['y']), axis=0)

        #close the contour
        if x[-1] != x[0] or y[-1] != y[0]:
            x = np.append(x, x[0])
            y = np.append(y, y[0])

        out = OrderedDict()
        out['x'] = x
        out['y'] = y
        out['nods'] = len(x)

        return out
    # }}}

    def checkconsistency(self, *args):  # {{{
        #recover options
        options = pairoptions(*args)

        #figure out if two boundaries are identical
        for i in range(len(self.boundaries)):
            namei = self.boundaries[i].shpfilename
            for j in range(1, len(self.boundaries)):
                namej = self.boundaries[j].shpfilename
                if namei == namej:
                    raise RuntimeError('Basin {} has two identical boundaries named {}'.format(self.name, namei))

        #go through boundaries and check their consistency
        for i in range(len(self.boundaries)):
            boundary == self.boundaries[i]
            boundary.checkconsistency()
    # }}}

    def contourplot(self, *args):  # {{{
        contour = self.contour()
        plot(contour.x, contour.y, 'r*-')
    # }}}

    def shapefilecrop(self, *args):  # {{{
        #recover options
        options = pairoptions(*args)
        threshold = options.getfieldvalue('threshold', .65) #.65 degrees lat, long
        inshapefile = options.getfieldvalue('shapefile')
        outputshapefile = options.getfieldvalue('outputshapefile', '')

        if options.exist('epsgshapefile') and options.exist('projshapefile'):
            raise RuntimeError('Basin shapefilecrop error messgae: cannot specify epsgshapefile and projshapefile at the same time')

        if options.exist('epsgshapefile'):
            projshapefile = epsg2proj(options.getfieldvalue('epsgshapefile'))
        elif options.exist('projshapefile'):
            projshapefile = options.getfieldvalue('projshapefile')
        else:
            raise RuntimeError('Basin shapefilecrop error message: epsgshapefile or projshapefile should be specified')

        #create list of contours that have critical length > threshold (in lat, long)
        contours = shpread(inshapefile)
        llist = []
        for i in range(len(contours)):
            contour = contours[i]
            carea = polyarea(contour['x'], contour['y'])
            clength = math.sqrt(carea)
            if clength < threshold:
                llist.append(i)

        contours_above_threshold = []
        for i in range(len(contours)):
            if i not in llist:
                contours_above_threshold.append(contours[i])
        contours = contours_above_threshold

        #project onto reference frame
        for i in range(len(contours)):
            h = contours[i]
            h['x'], h['y'] = CoordTransform(h['x'], h['y'], projshapefile, self.proj)
            contours[i]['x'] = h['x']
            contours[i]['y'] = h['y']

        #only keep the contours that are inside the basin (take centroids)
        ba = self.contour()
        flags = np.zeros(len(contours))
        for i in range(len(contours)):
            h = contours[i]
            isin = inpolygon(h['x'], h['y'], ba['x'], ba['y'])
            if len(np.where(isin == 0)[0]):
                flags[i] = 1

        pos = flags.nonzero()[0]

        contours_in_basin = []
        for i in range(len(contours)):
            if i not in pos:
                contours_in_basin.append(contours[i])
        contours = contours_in_basin

        #Two options
        output = None
        if outputshapefile == '':
            output = contours
        else:
            shpwrite(contours, outputshapefile)

        return output
    # }}}
# }}}
