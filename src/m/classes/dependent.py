import os.path
import numpy as np
from pairoptions import pairoptions
from fielddisplay import fielddisplay
from MatlabFuncs import *
from MeshProfileIntersection import MeshProfileIntersection


class dependent(object):
    """
    DEPENDENT class definition

       Usage:
          dependent = dependent()
    """

    def __init__(self, *args):  # {{{
        self.name = ''
        self.fos_reverse_index = float('NaN')
        self.exp = ''
        self.segments = []
        self.index = -1
        self.nods = 0

    #set defaults
        self.setdefaultparameters()

    #use provided options to change fields
        options = pairoptions(*args)

        self.name = options.getfieldvalue('name', '')
        self.exp = options.getfieldvalue('exp', '')
        self.segments = options.getfieldvalue('segments', [])
        self.index = options.getfieldvalue('index', -1)
        self.nods = options.getfieldvalue('nods', 0)

        #if name is mass flux:
        if strcmpi(self.name, 'MassFlux'):
            #make sure that we supplied a file and that it exists!
            if not os.path.exists(self.exp):
                raise IOError("dependent checkconsistency: specified 'exp' file does not exist!")
            #process the file and retrieve segments
            mesh = options.getfieldvalue('mesh')
            self.segments = MeshProfileIntersection(mesh.elements, mesh.x, mesh.y, self.exp)[0]
    # }}}
    def __repr__(self):  # {{{
        s = "   dependent variable:\n"

        s += "%s\n" % fielddisplay(self, 'name', "variable name (must match corresponding String)")

        if not np.isnan(self.fos_reverse_index):
            s += "%s\n" % fielddisplay(self, 'fos_reverse_index', "index for fos_reverse driver of ADOLC")
        if self.exp:
            s += "%s\n" % fielddisplay(self, 'exp', "file needed to compute dependent variable")
            s += "%s\n" % fielddisplay(self, 'segments', "mass flux segments")

        return s
    # }}}
    def setdefaultparameters(self):  # {{{
        #do nothing
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        if strcmpi(self.name, 'MassFlux'):
            if not self.segments:
                raise RuntimeError("dependent checkconsistency error: need segments to compute this dependent response")
            if self.index < 0:
                raise RuntimeError("dependent checkconsistency error: index for segments should be >= 0")

        if not np.isnan(self.fos_reverse_index):
            if not strcmpi(driver, 'fos_reverse'):
                raise TypeError("cannot declare a dependent with a fos_reverse_index when the driver is not fos_reverse!")
            if self.nods == 0:
                raise TypeError("dependent checkconsistency error: nods should be set to the size of the independent variable")

        return md
    # }}}
