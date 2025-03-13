import math

import numpy as np

from epsg2proj import epsg2proj
from fielddisplay import fielddisplay
from helpers import *
from pairoptions import pairoptions
from shpread import shpread


class boundary(object):  # {{{
    """boundary class definition

    Usage:
        boundary = boundary()
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
            self.shppath = options.getfieldvalue('shppath', '')
            self.shpfilename = options.getfieldvalue('shpfilename', '')
            self.orientation = options.getfieldvalue('orientation', 'normal')

            # If shppath is missing trailing slash, add it
            if self.shppath[-1] != '/':
                self.shppath += '/'

            #figure out projection string:
            if options.exist('epsg'):
                if options.exist('proj'):
                    raise RuntimeError('Error in constructor for boundary {}. Cannot supply epsg and proj at the same time!'.format(self.shppath))
                epsg = options.getfieldvalue('epsg')
                proj = epsg2proj(epsg) # convert to PROJ.4 format
            elif options.exist('proj'):
                if options.exist('epsg'):
                    raise RuntimeError('Error in constructor for boundary {}. Cannot supply proj and epsg at the same time!'.format(self.shppath))
                proj = options.getfieldvalue('proj')
            else:
                proj = '+proj=longlat +datum=WGS84 +no_defs'

            self.proj = proj
    # }}}

    def __repr__(self):  # {{{
        s = '   boundary parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'shppath', 'path to filename for this boundary'))
        s += '{}\n'.format(fielddisplay(self, 'shpfilename', 'shape file name'))
        s += '{}\n'.format(fielddisplay(self, 'orientation', 'orientation (default is \'normal\', can be \'reverse\')'))
        s += '{}\n'.format(fielddisplay(self, 'proj', 'shape file projection string (follows PROJ.4 standard)'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.shppath = ''
        self.shpfilename = ''
        self.orientation = 'normal'
        self.proj = '+proj=longlat +datum=WGS84 +no_defs' #EPSG 4326

        return self
    # }}}

    def name(self):  # {{{
        output = self.shpfilename

        return output
    # }}}

    def edges(self):  # {{{
        #read domain
        path, name, ext = fileparts(self.shpfilename)
        if ext != '.shp':
            ext = '.shp'
        output = shpread('{}{}{}'.format(self.shppath, name, ext))

        #do we reverse?
        if self.orientation == 'reverse':
            for i in range(len(output)):
                output[i]['x'] = np.flipud(output[i]['x'])
                output[i]['y'] = np.flipud(output[i]['y'])

        return output
    # }}}

    def plot(self, *args):  # {{{
        #recover options
        options = pairoptions(*args)

        #parse input
        # TODO: Sort out what options we should set (see 
        # src/m/classes/boundary.m)
        
        if options.exist('epsg'):
            if options.exist('proj'):
                raise RuntimeError('Boundary plot error message: cannot specify epsg and proj at the same time for plot projection')
            proj = epsg2proj(options.getfieldvalue('epsg'))
        elif options.exist('proj'):
            proj = options.getfieldvalue('proj')
        else:
            proj = epsg2proj(4326)

        #read domain
        path, name, ext = fileparts(self.shpfilename)
        if ext != '.shp':
            ext = '.shp'
        domain = shpread('{}{}{}'.format(self.shppath, name, ext))

        #convert boundary to another reference frame #{{{
        for i in range(len(domain)):
            try:
                x, y = CoordTransform(domain[i]['x'], domain[i]['y'], self.proj, proj)
            except error as e:
                print(e)
                print(self)

        # TODO: Figure out how to recover figure here: do we pass 'fig' and 
        # 'ax' in args?
        #for i in range(len(domain)):
            # x = domain[i].x * unitmultiplier
            # y = domain[i].y * unitmultiplier
            # if len(x) == 1:
        # }}}

        #TODO: Finish translating from MATLAB after test2004.py runs without plot
    # }}}

    def checkconsistency(self, *args):  # {{{
        #recover options
        options = pairoptions(*args)
        tolerance = getfieldvalue(options, 'tolerance', 1e-5)

        #recover edges
        edges = self.edges()

        if edges.Geometry == 'Point':
            return
        else:
            #check that we don't have two identical vertices
            x = edges.x
            y = edges.y
            distance = math.sqrt((x[:-2] - x[1:-1]) ** 2 + (y[:-2] - y[1:-1]) ** 2)
            dmax = distance.max()
            isidentical = np.where(np.asarray(distance) < dmax * tolerance)
            for elem in isidentical: # distance should be a 1D array, but return from np.where is a tuple of arrays
                if len(elem) != 0:
                    raise Exception('boundary {} has two vertices extermely close to one another'.format(shpfilename))

    def plot3d(self, *args):  # {{{
        #recover options
        options = pairoptions(*args)

        #parse input
        # TODO: Sort out what options we should set (see 
        # src/m/classes/boundary.m)
        
        if options.exist('epsg'):
            if options.exist('proj'):
                raise RuntimeError('Boundary plot error message: cannot specify epsg and proj at the same time for plot projection')
            proj = epsg2proj(options.getfieldvalue('epsg'))
        elif options.exist('proj'):
            proj = options.getfieldvalue('proj')
        else:
            proj = epsg2proj(4326)

        #read domain
        path, name, ext = fileparts(self.shpfilename)
        if ext != '.shp':
            ext = '.shp'
        domain = shpread('{}{}{}'.format(self.shppath, name, ext))

        #TODO: Finish translating from MATLAB after test2004.py runs without plot
    # }}}
# }}}
