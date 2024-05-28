import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class frictionpism(object):
    """FRICTIONPISM class definition

    Usage:
        frictionpism = frictionpism()
    """

    def __init__(self):
        self.pseudoplasticity_exponent = 0.
        self.threshold_speed = 0.
        self.delta = 0.
        self.void_ratio = 0.
        self.till_friction_angle = np.nan
        self.sediment_compressibility_coefficient = np.nan

        self.setdefaultparameters()
        self.requested_outputs = []
    # }}}

    def extrude(self, md):  # {{{
        self.till_friction_angle = project3d(md, 'vector', self.till_friction_angle, 'type', 'node', 'layer', 1)
        self.sediment_compressibility_coefficient = project3d(md, 'vector', self.sediment_compressibility_coefficient, 'type', 'node', 'layer', 1)
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        self.pseudoplasticity_exponent = 0.6
        self.threshold_speed = 100.
        self.delta = 0.02
        self.void_ratio = 0.69
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        if solution == 'TransientSolution' and not md.transient.isstressbalance and not md.transient.isthermal:
            return md
        md = checkfield(md, 'fieldname', 'friction.pseudoplasticity_exponent', 'numel', [1], '>', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.threshold_speed', 'numel', [1], '>', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.delta', 'numel', [1], '>', 0, '<', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.void_ratio', 'numel', [1], '>', 0, '<', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.till_friction_angle', 'NaN', 1, 'Inf', 1, '<', 360., '>', 0., 'size', [md.mesh.numberofvertices])  #User should give angle in degrees, Matlab calculates in rad
        md = checkfield(md, 'fieldname', 'friction.sediment_compressibility_coefficient', 'NaN', 1, 'Inf', 1, '<', 1., '>', 0., 'size', [md.mesh.numberofvertices])
        return md
    # }}}

    def __repr__(self):  # {{{
        s = 'Basal shear stress parameters for the PISM friction law (See Aschwanden et al. 2016 for more details)\n'
        s += "{}\n".format(fielddisplay(self, 'pseudoplasticity_exponent', 'pseudoplasticity exponent [dimensionless]'))
        s += "{}\n".format(fielddisplay(self, 'threshold_speed', 'threshold speed [m / yr]'))
        s += "{}\n".format(fielddisplay(self, 'delta', 'lower limit of the effective pressure, expressed as a fraction of overburden pressure [dimensionless]'))
        s += "{}\n".format(fielddisplay(self, 'void_ratio', 'void ratio at a reference effective pressure [dimensionless]'))
        s += "{}\n".format(fielddisplay(self, 'till_friction_angle', 'till friction angle [deg], recommended default: 30 deg'))
        s += "{}\n".format(fielddisplay(self, 'sediment_compressibility_coefficient', 'coefficient of compressibility of the sediment [dimensionless], recommended default: 0.12'))
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 10, 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'pseudoplasticity_exponent', 'format', 'Double')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'threshold_speed', 'format', 'Double', 'scale', 1. / yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'delta', 'format', 'Double')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'void_ratio', 'format', 'Double')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'till_friction_angle', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'sediment_compressibility_coefficient', 'format', 'DoubleMat', 'mattype', 1)
    # }}}
