import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class frictionregcoulomb2(*args):
    """frictionregcoulomb2 class definition

    Usage:
        frictionregcoulomb2 = frictionregcoulomb2()
    """

    def __init__(self, *args):  # {{{
        self.C = np.nan
        self.K = np.nan
        self.m = np.nan
        self.effective_pressure_limit = 0

        nargs = len(args)

        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # See ./frictionregcoulomb2.m
            raise Exception('not implemented')
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        # See Zoet and Iverson 2020 or Choi et al., 2022
        s = 'Regularized Coulomb friction law 2 parameters:\n'
        s += '   Regularized Coulomb friction law reads:\n'
        s += '                       C N |u|^(1/m)         \n'
        s += '      tau_b = -  ____________________________\n'
        s += '                   (|u| + (K*N)^m)^(1/m)     \n'
        s += '\n'
        s += '{}\n'.format(fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += '{}\n'.format(fielddisplay(self, 'm', 'm exponent'))
        s += '{}\n'.format(fielddisplay(self, 'K', '(K * N) ^ m to be velocity controlling plastic limit'))
        s += '{}\n'.format(fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.effective_pressure_limit = 0
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.C = project3d(md, 'vector', self.C, 'type', 'node')
        self.m = project3d(md, 'vector', self.m, 'type', 'element')
        self.K = project3d(md, 'vector', self.K, 'type', 'node')
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        md = checkfield(md, 'fieldname', 'friction.C', 'timeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0.)
        md = checkfield(md, 'fieldname', 'friction.K', 'NaN', 1, 'Inf', 1, '>', 0.)
        md = checkfield(md, 'fieldname', 'friction.m', 'NaN', 1, 'Inf', 1, '>', 0., 'size', [md.mesh.numberofelements, 1]);
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 15, 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'C', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'K', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'm', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'friction', 'fieldname', 'effective_pressure_limit', 'format', 'Double')
    # }}}
