import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class calving(object):
    """CALVING class definition

    Usage:
        calving = calving()
    """

    def __init__(self):  # {{{
        self.calvingrate = np.nan
        #self.setdefaultparameters() # Uncomment if/when setdefaultparameters is used
    # }}}

    def __repr__(self):  # {{{
        s = '   Calving parameters:'
        s += '{}\n'.format(fielddisplay(self, 'calvingrate', 'calving rate at given location [m/a]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.calvingrate = project3d(md, 'vector', self.calvingrate, 'type', 'node')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if solution != 'TransientSolution' or not md.transient.ismovingfront:
            return md

        md = checkfield(md, 'fieldname', 'calving.calvingrate', '>=', 0, 'timeseries', 1, 'NaN', 1, 'Inf', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 1, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'calvingrate', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts, 'scale', 1. / yts)
    # }}}
