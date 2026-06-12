import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class inversionnudging(object):
    """INVERSIONNUDGING class definition

    Usage:
        inversionnudging = inversionnudging()
    """

    def __init__(self):  # {{{
        self.iscontrol          = 0
        self.maxiter            = 0
        self.C0                 = 0.0
        self.max_increment_C    = 0.0
        self.relaxation_C       = 0.0
        self.tau_C              = 0.0
        self.H0_C               = 0.0
        self.min_C              = np.nan
        self.max_C              = np.nan
        self.melt0              = 0.0
        self.max_increment_melt = 0.0
        self.relaxation_melt    = 0.0
        self.tau_melt           = 0.0
        self.H0_melt            = 0.0
        self.min_melt           = np.nan
        self.max_melt           = np.nan
        self.dhdt_obs           = np.nan

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   Nudging parameters:\n\n'
        s += '{}\n'.format(fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter', 'maximum number of nudging steps'))
        s += '{}\n'.format(fielddisplay(self, 'dhdt_obs', 'observed thickness rate of change [m/yr]'))
        s += '\n     Friction parameters:\n\n'
        s += '         1   dC     H-Hobs    1  dH    rC\n'
        s += '         --  -- = - ------- - - --- - ---  (C - Ci)\n'
        s += '         C0  dt     tauC H0   H0 dt   tauC\n\n'
        s += '{}\n'.format(fielddisplay(self, 'C0', 'Friction scaling factor'))
        s += '{}\n'.format(fielddisplay(self, 'tau_C', 'adjustment timescale for friction coefficient [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'H0_C', 'thickness error scale for C (smaller = more sensitive) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'relaxation_C', 'relaxation strength toward C_inv (0 = none, 1 = strong)'))
        s += '{}\n'.format(fielddisplay(self, 'max_increment_C', 'maximum increase in C per nudging step *in log10 space*'))
        s += '{}\n'.format(fielddisplay(self, 'min_C', 'absolute minimum acceptable value of C'))
        s += '{}\n'.format(fielddisplay(self, 'max_C', 'absolute maximum acceptable value of C'))
        s += '\n     Melt perturbation parameters:\n\n'
        s += '         1   dP     H-Hobs    1  dH   r_m\n'
        s += '         --  -- = + ------- + - --- - ---  P\n'
        s += '       melt0 dt     tauM H0   H0 dt  tau_m\n\n'
        s += '{}\n'.format(fielddisplay(self, 'melt0', 'melt scaling factor [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'tau_melt', 'adjustment timescale for melt perturbation [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'H0_melt', 'thickness error scale for melt perturbation (smaller = more sensitive) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'relaxation_melt', 'relaxation strength toward perturbation = 0 (0 = none, 1 = strong)'))
        s += '{}\n'.format(fielddisplay(self, 'max_increment_melt', 'maximum increase in melt per nudging step'))
        s += '{}\n'.format(fielddisplay(self, 'min_melt', 'absolute minimum acceptable value of melt perturbation [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'max_melt', 'absolute maximum acceptable value of melt perturbation [m/yr]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # maximum number of nudging steps
        self.maxiter = 100

        # maximum step size in log10 space for C
        self.max_increment_C    = 0.05
        self.max_increment_melt = 0.5

        # C adjustment timescale, 100 yrs used in van den Akker et al (2025)
        self.tau_C    = 100.0
        self.tau_melt = 200.0

        # thickness error scale (smaller = more sensitive), 100 m used in van den Akker et al (2025)
        self.H0_C    = 100.0
        self.H0_melt = 50.0

        # relaxation strength toward C_inv (0 = none, 1 = strong), 0.5 used in van den Akker et al (2025)
        self.relaxation_C    = 0.5
        self.relaxation_melt = 0.3

        # default orders of magnitude for melt and C
        self.C0    = 1.0  # SI
        self.melt0 = 1.0  # m/yr

        return self
    # }}}

    def extrude(self, md):  # {{{
        self.min_C    = project3d(md, 'vector', self.min_C,    'type', 'node')
        self.max_C    = project3d(md, 'vector', self.max_C,    'type', 'node')
        self.min_melt = project3d(md, 'vector', self.min_melt, 'type', 'node')
        self.max_melt = project3d(md, 'vector', self.max_melt, 'type', 'node')
        self.dhdt_obs = project3d(md, 'vector', self.dhdt_obs, 'type', 'node')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if not self.iscontrol:
            return md

        md = checkfield(md, 'fieldname', 'inversion.iscontrol', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'inversion.maxiter', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.C0', 'numel', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'inversion.max_increment_C', 'numel', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'inversion.max_increment_melt', 'numel', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'inversion.min_C', 'size', [md.mesh.numberofvertices, 1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'inversion.max_C', 'size', [md.mesh.numberofvertices, 1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'inversion.melt0', 'numel', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'inversion.min_melt', 'size', [md.mesh.numberofvertices, 1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'inversion.max_melt', 'size', [md.mesh.numberofvertices, 1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'inversion.H0_C', 'numel', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'inversion.H0_melt', 'numel', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'inversion.relaxation_C', 'numel', [1], '>=', 0, '<=', 1)
        md = checkfield(md, 'fieldname', 'inversion.relaxation_melt', 'numel', [1], '>=', 0, '<=', 1)
        md = checkfield(md, 'fieldname', 'inversion.dhdt_obs', 'size', [md.mesh.numberofvertices, 1], 'NaN', 1, 'Inf', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'iscontrol', 'format', 'Boolean')
        WriteData(fid, prefix, 'name', 'md.inversion.type', 'data', 5, 'format', 'Integer')
        if not self.iscontrol:
            return
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'maxiter', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'C0', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'melt0', 'format', 'Double', 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'max_increment_C', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'max_increment_melt', 'format', 'Double', 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'min_C', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'max_C', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'min_melt', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'max_melt', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'tau_C', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'tau_melt', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'H0_C', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'H0_melt', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'relaxation_C', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'relaxation_melt', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'dhdt_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
    # }}}
