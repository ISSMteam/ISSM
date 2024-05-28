from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class calvingminthickness(object):
    """
    CALVINGMINTHICKNESS class definition

       Usage:
          calvingminthickness = calvingminthickness()
    """

    def __init__(self):  # {{{

        self.min_thickness = 0.

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = '   Calving Minimum thickness:'
        string = "%s\n%s" % (string, fielddisplay(self, 'min_thickness', 'minimum thickness below which no ice is allowed'))
        return string
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        #minimum thickness is 100 m by default
        self.min_thickness = 100.
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if solution == 'TransientSolution' or md.transient.ismovingfront == 0:
            return

        md = checkfield(md, 'fieldname', 'calving.min_thickness', '>', 0, 'NaN', 1, 'Inf', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 4, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_thickness', 'format', 'Double')
    # }}}
