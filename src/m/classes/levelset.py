import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class levelset(object):
    """LEVELSET class definition

    Usage:
        levelset = levelset()
    """

    def __init__(self):  # {{{
        self.stabilization = 0
        self.spclevelset = np.nan
        self.reinit_frequency = 5
        self.kill_icebergs = 0
        self.migration_max = 0
        self.fe = 'P1'

        # Set defaults
        self.setdefaultparameters()
    # }}}
    def __repr__(self):  # {{{
        s = '   Level-set parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'stabilization', '0: no, 1: artificial_diffusivity, 2: streamline upwinding'))
        s += '{}\n'.format(fielddisplay(self, 'spclevelset', 'Levelset constraints (NaN means no constraint)'))
        s += '{}\n'.format(fielddisplay(self, 'reinit_frequency', 'Amount of time steps after which the levelset function in re-initialized'))
        s += '{}\n'.format(fielddisplay(self, 'kill_icebergs', 'remove floating icebergs to prevent rigid body motions (1: true, 0: false)'))
        s += '{}\n'.format(fielddisplay(self, 'migration_max', 'maximum allowed migration rate (m/a)'))
        s += '{}\n'.format(fielddisplay(self, 'fe', 'Finite Element type: \'P1\' (default), or \'P2\''))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Stabilization = 1 by default
        self.stabilization = 1
        self.reinit_frequency = 5
        self.kill_icebergs = 1
        self.migration_max = 1e12 # No need for general cases, unless specified

        # Linear elements by default
        self.fe = 'P1'

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if (solution != 'TransientSolution') or (not md.transient.ismovingfront):
            return md

        md = checkfield(md, 'fieldname', 'levelset.spclevelset', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'levelset.stabilization', 'numel', [1], 'values', [0, 1, 2, 5, 6])
        md = checkfield(md, 'fieldname', 'levelset.kill_icebergs', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'levelset.migration_max', 'numel', [1], 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'levelset.fe', 'values', ['P1', 'P2'])

        return md
    # }}}

    def extrude(self, md):  # {{{
        self.spclevelset = project3d(md, 'vector', self.spclevelset, 'type', 'node')
        return self
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'object', self, 'fieldname', 'stabilization', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spclevelset', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'reinit_frequency', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'kill_icebergs', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'migration_max', 'format', 'Double', 'scale', 1 / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'fe', 'format', 'String')
    # }}}
