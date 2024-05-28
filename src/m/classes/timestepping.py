from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class timestepping(object):
    """TIMESTEPPING Class definition

    Usage:
        timestepping = timestepping()
    """

    def __init__(self, *args):  #{{{
        self.start_time = 0
        self.final_time = 0
        self.time_step = 0
        self.interp_forcing = 1
        self.average_forcing = 0
        self.cycle_forcing = 0
        self.coupling_time = 0

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise RuntimeError('constructor not supported')
    # }}}

    def __repr__(self):  #{{{
        s = '   timestepping parameters:\n'
        unit = 'yr'
        s += '{}\n'.format(fielddisplay(self, 'start_time', 'simulation starting time [' + unit + ']'))
        s += '{}\n'.format(fielddisplay(self, 'final_time', 'final time to stop the simulation [' + unit + ']'))
        s += '{}\n'.format(fielddisplay(self, 'time_step', 'length of time steps [' + unit + ']'))
        s += '{}\n'.format(fielddisplay(self, 'interp_forcing', 'interpolate in time between requested forcing values? (0 or 1)'))
        s += '{}\n'.format(fielddisplay(self, 'average_forcing', 'average in time if there are several forcing values between steps? (0 or 1, default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'cycle_forcing', 'cycle through forcing? (0 or 1)'))
        s += '{}\n'.format(fielddisplay(self, 'coupling_time', 'length of coupling time steps with ocean model [' + unit + ']'))
        return s
    # }}}

    def setdefaultparameters(self):  #{{{
        # Time between 2 time steps
        self.time_step = 1 / 2

        # Final time
        self.final_time = 10 * self.time_step

        # Should we interpolate forcing between timesteps?
        self.interp_forcing = 1
        self.average_forcing = 0
        self.cycle_forcing = 0

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        md = checkfield(md, 'fieldname', 'timestepping.start_time', 'numel', [1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'timestepping.final_time', 'numel', [1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'timestepping.time_step', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'timestepping.interp_forcing', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'timestepping.average_forcing', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'timestepping.cycle_forcing', 'numel', [1], 'values', [0, 1])
        if (self.final_time - self.start_time) < 0:
            md.checkmessage('timestepping.final_time should be larger than timestepping.start_time')
        if solution == 'TransientSolution':
            md = checkfield(md, 'fieldname', 'timestepping.time_step', 'numel', [1], '>', 0, 'NaN', 1, 'Inf', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        scale = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.timestepping.type', 'data', 1, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'start_time', 'format', 'Double', 'scale', scale)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'final_time', 'format', 'Double', 'scale', scale)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'time_step', 'format', 'Double', 'scale', scale)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'interp_forcing', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'average_forcing', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'cycle_forcing', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'coupling_time', 'format', 'Double', 'scale', scale)
    # }}}
