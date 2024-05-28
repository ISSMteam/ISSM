import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class frontalforcings(object):
    """frontalforcings class definition

    Usage:
        frontalforcings = frontalforcings()
    """

    def __init__(self, *args):  # {{{
        self.meltingrate = np.nan
        self.ablationrate = np.nan

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # TODO: Replace the following with constructor
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   Frontalforcings parameters:'
        s += '{}\n'.format(fielddisplay(self, 'meltingrate', 'melting rate at given location [m/a]'))
        s += '{}\n'.format(fielddisplay(self, 'ablationrate', 'frontal ablation rate at given location [m/a], it contains both calving and melting'))

        return s
    # }}}

    def extrude(self, md):  # {{{
        self.meltingrate = project3d(md, 'vector', self.meltingrate, 'type', 'node')
        self.ablationrate = project3d(md, 'vector', self.ablationrate, 'type', 'node')
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        self.meltingrate = np.nan
        self.ablationrate = np.nan
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if (solution != 'TransientSolution') or (not md.transient.ismovingfront):
            return md

        md = checkfield(md, 'fieldname', 'frontalforcings.meltingrate', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)
        if not np.isnan(md.frontalforcings.ablationrate):
            md = checkfield(md, 'fieldname', 'frontalforcings.ablationrate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.frontalforcings.parameterization', 'data', 1, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'meltingrate', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts, 'scale', 1. / yts)
        if not np.isnan(md.frontalforcings.ablationrate):
            WriteData(fid, prefix, 'object', self, 'fieldname', 'ablationrate', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts, 'scale', 1. / yts)
    # }}}
