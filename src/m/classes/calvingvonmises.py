from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class calvingvonmises(object):
    """CALVINGVONMISES class definition

    Usage:
        calvingvonmises = calvingvonmises()
    """

    def __init__(self):  # {{{
        self.stress_threshold_groundedice = 0
        self.stress_threshold_floatingice = 0
        self.min_thickness = 0

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        s = '   Calving VonMises parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'stress_threshold_groundedice', 'sigma_max applied to grounded ice only [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'stress_threshold_floatingice', 'sigma_max applied to floating ice only [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed [m]'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Default sigma max
        self.stress_threshold_groundedice = 1e6
        self.stress_threshold_floatingice = 150e3

        # Turn off min_thickness by default
        self.min_thickness = 0.
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if solution == 'TransientSolution' or not md.transient.ismovingfront:
            return

        md = checkfield(md, 'fieldname', 'calving.stress_threshold_groundedice', '>', 0, 'nan', 1, 'Inf', 1, 'size', 'universal')
        md = checkfield(md, 'fieldname', 'calving.stress_threshold_floatingice', '>', 0, 'nan', 1, 'Inf', 1, 'size', 'universal')
        md = checkfield(md, 'fieldname', 'calving.min_thickness', '>=', 0, 'NaN', 1, 'Inf', 1, 'numel', [1])

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stress_threshold_groundedice', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts )
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stress_threshold_floatingice', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_thickness', 'format', 'Double')
    # }}}
