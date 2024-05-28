import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from structtoobj import structtoobj
from WriteData import WriteData


class frictionschoof(object):
    """FRICTIONSCHOOF class definition

    Usage:
        friction = frictionschoof()
    """

    def __init__(self, *args):  # {{{
        self.C                        = np.nan
        self.Cmax                     = np.nan
        self.m                        = np.nan
        self.coupling                 = 0
        self.effective_pressure       = np.nan
        self.effective_pressure_limit = 0

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            self = structtoobj(self, args[0])
        else:
            raise Exception('constructor not supported')
    # }}}
    def __repr__(self):  # {{{
        # See Brondex et al. 2019 https://www.the-cryosphere.net/13/177/2019/
        s = 'Schoof sliding law parameters:\n'
        s += '   Schoof\'s sliding law reads:\n'
        s += '                         C^2 |u_b|^(m-1)                \n'
        s += '      tau_b = - _____________________________   u_b   \n'
        s += '               (1+(C^2/(Cmax N))^1/m |u_b| )^m          \n'
        s += '\n'
        s += "{}\n".format(fielddisplay(self, 'C', 'friction coefficient [SI]'))
        s += "{}\n".format(fielddisplay(self, 'Cmax', 'Iken\'s bound (typically between 0.17 and 0.84) [SI]'))
        s += "{}\n".format(fielddisplay(self, 'm', 'm exponent (generally taken as m = 1/n = 1/3)'))
        s += '{}\n'.format(fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += "{}\n".format(fielddisplay(self, 'effective_pressure_limit', 'fNeff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)'))
        return s
    # }}}
    def setdefaultparameters(self):  # {{{
        self.effective_pressure_limit = 0
        return self
    # }}}
    def extrude(self, md):  # {{{
        self.C = project3d(md, 'vector', self.C, 'type', 'node')
        self.Cmax = project3d(md, 'vector', self.Cmax, 'type', 'node')
        self.m = project3d(md, 'vector', self.m, 'type', 'element')
        if self.coupling in [3, 4]:
            self.effective_pressure = project3d(md, 'vector', self.effective_pressure, 'type', 'node', 'layer', 1)
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        md = checkfield(md, 'fieldname', 'friction.C', 'timeseries', 1, 'NaN', 1, 'Inf', 1, '>',0.)
        md = checkfield(md, 'fieldname', 'friction.Cmax', 'timeseries', 1, 'NaN', 1, 'Inf', 1, '>', 0.)
        md = checkfield(md, 'fieldname', 'friction.m', 'NaN', 1, 'Inf', 1, '>', 0., 'size', [md.mesh.numberofelements, 1])
        md = checkfield(md, 'fieldname', 'friction.effective_pressure_limit', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'friction.coupling', 'numel', [1], 'values', [0, 1, 2, 3, 4])
        if self.coupling == 3:
            md = checkfield(md, 'fieldname', 'friction.effective_pressure', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 11, 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'C', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'Cmax', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'm', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'coupling', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'friction', 'fieldname', 'effective_pressure_limit', 'format', 'Double')
        if self.coupling == 3 or self.coupling == 4:
            WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'effective_pressure', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
