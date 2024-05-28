from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData
import numpy as np


class outputdefinition(object):
    """
    OUTPUTDEFINITION class definition

       Usage:
          outputdefinition = outputdefinition()
    """

    def __init__(self):  # {{{
        self.definitions = []
    # }}}

    def __repr__(self):  # {{{
        string = "   Outputdefinitions:"

        string = "%s\n%s" % (string, fielddisplay(self, "definitions", "list of potential outputs that can be requested, but which need additional data to be defined"))

        return string
    # }}}

    def extrude(self, md):  # {{{
        for definition in self.definitions:
            definition.extrude(md)

        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{

        md = checkfield(md, 'fieldname', 'outputdefinition.definitions', 'cell', 1)
        for definition in self.definitions:
            definition.checkconsistency(md, solution, analyses)

    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        data = []
        for i in range(len(self.definitions)):
            self.definitions[i].marshall(prefix, md, fid)
            classdefinition = self.definitions[i].__class__.__name__
            classdefinition = classdefinition[0].upper() + classdefinition[1:]
            data.append(classdefinition)

        data = np.unique(data)
        WriteData(fid, prefix, 'data', data, 'name', 'md.outputdefinition.list', 'format', 'StringArray')
    # }}}
