import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class friction(object):
    """friction class definition

    Usage:
        friction = friction()
    """

    def __init__(self):  # {{{
        self.coefficient = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.linearize = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0
        self.setdefaultparameters()
    # }}}
    def __repr__(self):  # {{{
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * base, r = q / p and s = 1 / p)\n'
        s += '{}\n'.format(fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        s += '{}\n'.format(fielddisplay(self, 'p', 'p exponent'))
        s += '{}\n'.format(fielddisplay(self, 'q', 'q exponent'))
        s += '{}\n'.format(fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        s += '{}\n'.format(fielddisplay(self, 'linearize', '0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s
    # }}}
    def setdefaultparameters(self):  # {{{
        self.linearize = 0
        self.coupling = 0
        self.effective_pressure_limit = 0
        return self
    # }}}
    def extrude(self, md):  # {{{
        self.coefficient = project3d(md, 'vector', self.coefficient, 'type', 'node', 'layer', 1)
        self.p = project3d(md, 'vector', self.p, 'type', 'element')
        self.q = project3d(md, 'vector', self.q, 'type', 'element')
        if self.coupling in[3, 4]:
            self.effective_pressure = project3d(md, 'vector', self.effective_pressure, 'type', 'node', 'layer', 1)
        return self
    # }}}
    def defaultoutputs(self, md):  # {{{
        list = []
        return list
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        if solution == 'TransientSolution' and not md.transient.isstressbalance and not md.transient.isthermal:
            return md
        md = checkfield(md, 'fieldname', 'friction.coefficient', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.q', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.p', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.linearize', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'friction.coupling', 'numel', [1], 'values', [0, 1, 2, 3, 4])
        md = checkfield(md, 'fieldname', 'friction.effective_pressure_limit', 'numel', [1], '>=', 0)
        if self.coupling == 3:
            md = checkfield(md, 'fieldname', 'friction.effective_pressure', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 1, 'format', 'Integer')
        if type(self.coefficient) in [np.ndarray] and (self.coefficient.shape[0] == md.mesh.numberofvertices or self.coefficient.shape[0] == (md.mesh.numberofvertices + 1)):
            mattype = 1
            tsl = md.mesh.numberofvertices
        else:
            mattype = 2
            tsl = md.mesh.numberofelements
        WriteData(fid, prefix, 'object', self, 'fieldname', 'coefficient', 'format', 'DoubleMat', 'mattype', mattype, 'timeserieslength', tsl + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'p', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'q', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'coupling', 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'linearize', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'friction', 'fieldname', 'effective_pressure_limit', 'format', 'Double')
        if self.coupling == 3 or self.coupling == 4:
            WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'effective_pressure', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
