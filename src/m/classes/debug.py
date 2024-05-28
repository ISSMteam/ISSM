from fielddisplay import fielddisplay
from WriteData import *


class debug(object):
    """
    DEBUG class definition

       Usage:
          debug = debug()
    """

    def __init__(self):  # {{{
        self.valgrind = False
        self.gprof = False
        self.profiling = False

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = "   debug parameters:"

        string = "%s\n%s" % (string, fielddisplay(self, "valgrind", "use Valgrind to debug (0 or 1)"))
        string = "%s\n%s" % (string, fielddisplay(self, "gprof", "use gnu - profiler to find out where the time is spent"))
        string = "%s\n%s" % (string, fielddisplay(self, 'profiling', 'enables profiling (memory, flops, time)'))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'profiling', 'format', 'Boolean')
    # }}}
