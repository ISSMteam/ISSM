import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from pairoptions import pairoptions
from WriteData import WriteData


class nodalvalue(object):
    """NODALVALUE class definition

    Usage:
        nodalvalue=nodalvalue()
        nodalvalue=nodalvalue(
            'name', 'SealevelchangeSNodalValue',
            'definitionstring', 'Outputdefinition1',
            'model_string', 'SealevelchangeS',
            'node', 1
        )
    """

    def __init__(self, *args):  #{{{
        self.name = ''
        self.definitionstring = ''  # string that identifies this output definition uniquely, from 'Outputdefinition[1-10]'
        self.model_string = ''  # string for field that is being retrieved
        self.node = np.nan  #for which node are we retrieving the value?

        #use provided options to change fields
        options = pairoptions(*args)

        # Get name
        self.name = options.getfieldvalue('name', '')
        self.definitionstring = options.getfieldvalue('definitionstring', '')
        self.model_string = options.getfieldvalue('model_string', '')
        self.node = options.getfieldvalue('node', '')
    # }}}

    def __repr__(self):  # {{{
        s = '   Nodalvalue:\n'
        s += '{}\n'.format(fielddisplay(self, 'name', 'identifier for this nodalvalue response'))
        s += '{}\n'.format(fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from \'Outputdefinition[1-10]\''))
        s += '{}\n'.format(fielddisplay(self, 'model_string', 'string for field that is being retrieved'))
        s += '{}\n'.format(fielddisplay(self, 'node', 'vertex index at which we retrieve the value'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if not isinstance(self.name, str):
            raise Exception('nodalvalue error message: \'name\' field should be a string!')
        OutputdefinitionStringArray = []
        for i in range(100):
            OutputdefinitionStringArray.append('Outputdefinition{}'.format(i))
        md = checkfield(md, 'fieldname', 'self.definitionstring', 'field', self.definitionstring, 'values', OutputdefinitionStringArray)
        md = checkfield(md, 'fieldname', 'self.node', 'field', self.node, 'values', range(md.mesh.numberofvertices))
        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        WriteData(fid, prefix, 'data', self.name, 'name', 'md.nodalvalue.name', 'format', 'String')
        WriteData(fid, prefix, 'data', self.definitionstring, 'name', 'md.nodalvalue.definitionenum', 'format', 'String')
        WriteData(fid, prefix, 'data', self.model_string, 'name', 'md.nodalvalue.model_enum', 'format', 'String')
        WriteData(fid, prefix, 'data', self.node, 'name', 'md.nodalvalue.node', 'format', 'Integer')
    # }}}
