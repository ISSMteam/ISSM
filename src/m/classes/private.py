from collections import OrderedDict
from fielddisplay import fielddisplay


class private(object):
    """
    PRIVATE class definition

       Usage:
          private = private()
    """

    def __init__(self):  # {{{
        self.isconsistent = True
        self.runtimename = ''
        self.bamg = OrderedDict()
        self.solution = ''

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = '   private parameters: do not change'

        string = "%s\n%s" % (string, fielddisplay(self, 'isconsistent', 'is model self consistent'))
        string = "%s\n%s" % (string, fielddisplay(self, 'runtimename', 'name of the run launched'))
        string = "%s\n%s" % (string, fielddisplay(self, 'bamg', 'structure with mesh properties constructed if bamg is used to mesh the domain'))
        string = "%s\n%s" % (string, fielddisplay(self, 'solution', 'type of solution launched'))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        return md
    # }}}
