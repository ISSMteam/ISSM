from collections import OrderedDict
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class miscellaneous(object):
    """
    MISCELLANEOUS class definition

       Usage:
          miscellaneous = miscellaneous()
    """

    def __init__(self):  # {{{
        self.notes = ''
        self.name = ''
        self.dummy = OrderedDict()

    #set defaults
        self.setdefaultparameters()

    # }}}
    def __repr__(self):  # {{{
        string = '   miscellaneous parameters:'

        string = "%s\n%s" % (string, fielddisplay(self, 'notes', 'notes in a cell of strings'))
        string = "%s\n%s" % (string, fielddisplay(self, 'name', 'model name'))
        string = "%s\n%s" % (string, fielddisplay(self, 'dummy', 'empty field to store some data'))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'miscellaneous.name', 'empty', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  #  {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'name', 'format', 'String')
    # }}}
