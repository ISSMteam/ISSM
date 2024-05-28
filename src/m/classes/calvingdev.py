from fielddisplay import fielddisplay
from project3d import project3d
from checkfield import checkfield
from WriteData import WriteData


class calvingdev(object):
    """
    CALVINGDEV class definition

       Usage:
          calvingdev = calvingdev()
    """

    def __init__(self):  # {{{

        self.stress_threshold_groundedice = 0.
        self.stress_threshold_floatingice = 0.

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = '   Calving Pi parameters:'
        string = "%s\n%s" % (string, fielddisplay(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))

        return string
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        #Default sigma max
        self.stress_threshold_groundedice = 1e6
        self.stress_threshold_floatingice = 150e3
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if solution == 'TransientSolution' or md.transient.ismovingfront == 0:
            return

        md = checkfield(md, 'fieldname', 'calving.stress_threshold_groundedice', '>', 0, 'nan', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'calving.stress_threshold_floatingice', '>', 0, 'nan', 1, 'Inf', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stress_threshold_groundedice', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stress_threshold_floatingice', 'format', 'DoubleMat', 'mattype', 1)
    # }}}
