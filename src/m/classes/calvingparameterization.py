from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class calvingparameterization(object):
    """calvingparameterization class definition
    For test calving laws and coefficients

    Usage:
        calvingparameterization = calvingparameterization()
    """

    def __init__(self, *args):  # {{{
        self.min_thickness = 0
        self.use_param = 0
        self.theta = 0
        self.alpha = 0
        self.xoffset = 0
        self.yoffset = 0
        self.vel_upperbound = 0
        self.vel_threshold = 0
        self.vel_lowerbound = 0

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
        s = '   Calving test parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        s += '{}\n'.format(fielddisplay(self, 'use_param', '-1 - just use frontal ablation rate, 0 - f(x) = y_{o} + \alpha (x+x_{o}), 1 - f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})), 2 - tanh(thickness), 3 - tanh(normalized vel), 4 - tanh(truncated vel), 5 - linear(truncated vel)'))
        s += '{}\n'.format(fielddisplay(self, 'theta', 'the amplifier'))
        s += '{}\n'.format(fielddisplay(self, 'alpha', 'the slope'))
        s += '{}\n'.format(fielddisplay(self, 'xoffset', 'offset in x-axis'))
        s += '{}\n'.format(fielddisplay(self, 'yoffset', 'offset in y-axis'))
        s += '{}\n'.format(fielddisplay(self, 'vel_lowerbound', 'lowerbound of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(fielddisplay(self, 'vel_threshold', 'threshold of ice velocity to reduce the calving rate [m/a]'))
        s += '{}\n'.format(fielddisplay(self, 'vel_upperbound', 'upperbound of ice velocity to reduce the calving rate [m/a]'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # For now we turn this off by setting the threshold to 0
        self.min_thickness = 0.

        # Parameters for the spatial temporal separation
        # The coefficient follows: gamma= f(x)
        # 0 - f(x) = y_{o} + \alpha (x+x_{o})
        # 1 - f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o}))
        self.use_param = 0

        # The amplifier
        self.theta = 0

        # The slope alpha
        self.alpha = 0

        # Offset in x-axis
        self.xoffset

        # Offset in y-axis
        self.yoffset

        # Velocity thresholds to reduce calving rate
        self.vel_upperbound = 6000 # m/a
        self.vel_lowerbound = 0 # m/a
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if not solution == 'TransientSolution' or not md.transient.ismovingfront:
            return

        md = checkfield(md, 'fieldname', 'calving.min_thickness', '>=', 0, 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'calving.use_param', 'values', [-1, 0, 1, 2, 3, 4, 5])
        md = checkfield(md, 'fieldname', 'calving.theta', 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'calving.alpha', 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'calving.xoffset', 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'calving.yoffset', 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'calving.vel_lowerbound', 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'calving.vel_threshold', 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'calving.vel_upperbound', 'NaN', 1, 'Inf', 1, 'numel', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 9, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_thickness', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'use_param', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'theta', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'alpha', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'xoffset', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'yoffset', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vel_lowerbound', 'format', 'Double', 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vel_threshold','format', 'Double', 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vel_upperbound', 'format', 'Double', 'scale', 1. / yts)
    # }}}
