import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class frictionregcoulomb(*args):
    """frictionregcoulomb class definition

    Usage:
        frictionregcoulomb = frictionregcoulomb()
    """

    def __init__(self, *args):  # {{{
        self.C = np.nan
        self.u0 = 0
        self.m = np.nan

        nargs = len(args)

        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # See ./frictionregcoulomb.m
            raise Exception('not implemented')
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        # See Joughin et al. 2019 (equivalent form by Matt Trevers, poster at AGU 2022) https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL082526
        s = 'Regularized Coulomb friction law (Joughin et al., 2019) parameters:\n'
        s += '   Regularized Coulomb friction law reads:\n'
        s += '                       C^2 |u|^(1/m)         \n'
        s += '      tau_b = -  ____________________________\n'
        s += '                     (|u|/u0 + 1)^(1/m)      \n'
        s += '\n'
        s += '{}\n'.format(fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(fielddisplay(self, 'm', 'm exponent (set to m = 3 in original paper)'))
        s += '{}\n'.format(fielddisplay(self, 'u0', 'velocity controlling plastic limit'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.u0 = 1000
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.C = project3d(md, 'vector', self.C, 'type', 'node')
        self.m = project3d(md, 'vector', self.m, 'type', 'element')
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        md = checkfield(md, 'fieldname', 'friction.C', 'timeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0.)
        md = checkfield(md, 'fieldname', 'friction.u0', 'NaN', 1, 'Inf', 1, '>', 0., 'numel', 1)
        md = checkfield(md, 'fieldname', 'friction.m', 'NaN', 1, 'Inf', 1, '>', 0., 'size', [md.mesh.numberofelements, 1]);
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 14, 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'C', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'u0', 'format', 'Double', 'scale', 1 / yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'm', 'format', 'DoubleMat', 'mattype', 2)
    # }}}
