import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class steadystate(object):
    """
    STEADYSTATE class definition

       Usage:
          steadystate = steadystate()
    """

    def __init__(self):  # {{{
        self.reltol = 0
        self.maxiter = 0
        self.requested_outputs = []

    #set defaults
        self.setdefaultparameters()

    # }}}
    def __repr__(self):  # {{{
        string = '   steadystate solution parameters:'
        string = "%s\n%s" % (string, fielddisplay(self, 'reltol', 'relative tolerance criterion'))
        string = "%s\n%s" % (string, fielddisplay(self, 'maxiter', 'maximum number of iterations'))
        string = "%s\n%s" % (string, fielddisplay(self, 'requested_outputs', 'additional requested outputs'))
        return string
    # }}}

    def defaultoutputs(self, md):  # {{{
        return md.stressbalance.defaultoutputs(md) + md.thermal.defaultoutputs(md)

    # }}}
    def setdefaultparameters(self):  # {{{
        #maximum of steady state iterations
        self.maxiter = 100
        #Relative tolerance for the steadystate convertgence
        self.reltol = 0.01
        #default output
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if not solution == 'SteadystateSolution':
            return md

        if not md.timestepping.time_step == 0:
            md.checkmessage("for a steadystate computation, timestepping.time_step must be zero.")

        if np.isnan(md.stressbalance.reltol):
            md.checkmessage("for a steadystate computation, stressbalance.reltol (relative convergence criterion) must be defined!")

        md = checkfield(md, 'fieldname', 'steadystate.requested_outputs', 'stringrow', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'reltol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'maxiter', 'format', 'Integer')

    #process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.steadystate.requested_outputs', 'format', 'StringArray')
    # }}}
