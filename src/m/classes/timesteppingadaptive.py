from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class timesteppingadaptive(object):
    """
    TIMESTEPPINGADAPTIVE Class definition

       Usage:
          timesteppingadaptive = timesteppingadaptive()
    """

    def __init__(self, *args):  # {{{
        if not len(args):
            self.start_time = 0.
            self.final_time = 0.
            self.time_step_min = 0.
            self.time_step_max = 0.
            self.cfl_coefficient = 0.
            self.interp_forcing = 1
            self.average_forcing = 0
            self.cycle_forcing = 0
            self.coupling_time = 0.

            #set defaults
            self.setdefaultparameters()

        elif len(args) == 1 and args[0].__module__ == 'timestepping':
            old = args[0]
            #first call setdefaultparameters:
            self.setdefaultparameters()
            self.start_time = old.start_time
            self.final_time = old.final_time
            self.interp_forcing = old.interp_forcing
            self.average_forcing = old.average_forcing
            self.cycle_forcing = old.cycle_forcing
            self.coupling_time = old.coupling_time

        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        string = "   timesteppingadaptive parameters:"
        string = '{}\n{}'.format (string, fielddisplay(self, "start_time", "simulation starting time [yr]"))
        string = '{}\n{}'.format(string, fielddisplay(self, "final_time", "final time to stop the simulation [yr]"))
        string = '{}\n{}'.format (string, fielddisplay(self, "time_step_min", "minimum length of time steps [yr]"))
        string = '{}\n{}'.format (string, fielddisplay(self, "time_step_max", "maximum length of time steps [yr]"))
        string = '{}\n{}'.format (string, fielddisplay(self, "cfl_coefficient", "coefficient applied to cfl condition"))
        string = '{}\n{}'.format (string, fielddisplay(self, "interp_forcing", "interpolate in time between requested forcing values ? (0 or 1)"))
        string = '{}\n{}'.format(string, fielddisplay(self, 'average_forcing', 'average in time if there are several forcing values between steps? (0 or 1, default is 0)'))
        string = '{}\n{}'.format(string, fielddisplay(self, "cycle_forcing", "cycle through forcing ? (0 or 1)"))
        string = '{}\n{}'.format(string, fielddisplay(self, "coupling_time", "coupling time steps with ocean model [yr]"))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        #time between 2 time steps
        self.time_step_min = 0.01
        self.time_step_max = 10.
        #final time
        self.final_time = 10. * self.time_step_max
        #time adaptation?
        self.cfl_coefficient = 0.5
        #should we interpolate forcing between timesteps?
        self.interp_forcing = 1
        self.average_forcing = 0
        self.cycle_forcing   = 0
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'timestepping.start_time', 'numel', [1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'timestepping.final_time', 'numel', [1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'timestepping.time_step_min', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'timestepping.time_step_max', 'numel', [1], '>=', md.timestepping.time_step_min, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'timestepping.cfl_coefficient', 'numel', [1], '>', 0, '<=', 1)
        if self.final_time - self.start_time < 0:
            md.checkmessage("timestepping.final_time should be larger than timestepping.start_time")
        md = checkfield(md, 'fieldname', 'timestepping.interp_forcing', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'timestepping.average_forcing', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'timestepping.cycle_forcing', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'timestepping.coupling_time', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.timestepping.type', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'start_time', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'final_time', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'time_step_min', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'time_step_max', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'cfl_coefficient', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'interp_forcing', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'average_forcing', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'cycle_forcing', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'timestepping', 'fieldname', 'coupling_time', 'format', 'Double', 'scale', yts)
    # }}}
