from project3d import project3d
from fielddisplay import fielddisplay
from pairoptions import pairoptions
from checkfield import checkfield
from WriteData import WriteData
from ContourToMesh import ContourToMesh
import numpy as np
import os


class regionaloutput(object):
    """
    REGIONALOUTPUT class definition

       Usage:
          regionaloutput = regionaloutput()
          regionaloutput = regionaloutput('name', 'Volume1', 'definitionstring', 'Outputdefinition1', 'outputnamestring', 'IceVolume', 'mask', mask)
          regionaloutput = regionaloutput('name', 'Volume1', 'definitionstring', 'Outputdefinition1', 'outputnamestring', 'IceVolume', 'maskexpstring', 'Exp/Mask.exp', 'model', md)

       where mask is a vectorial field of size md.mesh.numberofvertices, 1 : where vertices with values > 1 are to be included in the calculated region.
       Alternatively, the user can pass in an Argus file and model object instead of a mask, and mask will be calculated for the user
    """

    def __init__(self, *args):  # {{{

        self.name = ''
        self.model= ''
        self.definitionstring = ''
        self.outputnamestring = ''
        self.mask = np.nan
        self.maskexpstring = ''

    #set defaults
        self.setdefaultparameters()

    #use provided options to change fields
        options = pairoptions(*args)

    #OK get other fields
        self = options.AssignObjectFields(self)

    #get name
        if options.getfieldvalue('model', 0):
            if options.getfieldvalue('maskexpstring', 0):
                modelname = options.getfieldvalue('model')
                self.maskexpstring = options.getfieldvalue('maskexpstring')
                self.setmaskfromexp(modelname)

        # if (len(self.mask) <= 1 & np.any(np.isnan(self.mask))):
        #     raise IOError('regionaloutput error message: ''mask'' field or ''maskexpstring'' and ''model'' fields should be defined!')

    # }}}

    def __repr__(self):  # {{{
        string = "   Regionaloutput:"
        string = "%s\n%s" % (string, fielddisplay(self, 'name', 'identifier for this regional response'))
        string = "%s\n%s" % (string, fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from Outputdefinition[1 - 100]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'outputnamestring', 'string that identifies the type of output you want, eg. IceVolume, TotalSmb, GroudedArea'))
        string = "%s\n%s" % (string, fielddisplay(self, 'mask', 'mask vectorial field which identifies the region of interest (value > 0 will be included)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'maskexpstring', 'name of Argus file that can be passed in to define the regional mask'))
        return string
    # }}}

    def extrude(self, md):  # {{{
        self.mask = project3d(md, 'vector', self.mask, 'type', 'node')
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def setmaskfromexp(self, md):  # {{{
        if len(self.maskexpstring) > 0:
            self.mask = ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y, self.maskexpstring, 'node', 1)

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{

        if not isinstance(self.name, str):
            raise RuntimeError("regionaloutput error message: 'name' field should be a string!")

        if not isinstance(self.outputnamestring, str):
            raise RuntimeError("regionaloutput error message: 'outputnamestring' field should be a string!")

        if len(self.maskexpstring) > 0:
            if not os.path.isfile(self.maskexpstring):
                raise RuntimeError("regionaloutput error message: file name for mask exp does not point to a legitimate file on disk!")
            else:
                self.setmaskfromexp(md)

        OutputdefinitionStringArray = []
        for i in range(1, 100):
            x = 'Outputdefinition' + str(i)
            OutputdefinitionStringArray.append(x)

        md = checkfield(md, 'field', self.definitionstring, 'values', OutputdefinitionStringArray)
        md = checkfield(md, 'field', self.mask, 'size', [md.mesh.numberofvertices], 'NaN', 1, 'Inf', 1)
        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{

        #before marshalling, make sure mask is set:
        self.setmaskfromexp(md)

    #ok, marshall strings and mask:
        WriteData(fid, prefix, 'data', self.name, 'name', 'md.regionaloutput.name', 'format', 'String')
        WriteData(fid, prefix, 'data', self.definitionstring, 'name', 'md.regionaloutput.definitionstring', 'format', 'String')
        WriteData(fid, prefix, 'data', self.outputnamestring, 'name', 'md.regionaloutput.outputnamestring', 'format', 'String')
        WriteData(fid, prefix, 'data', self.mask, 'name', 'md.regionaloutput.mask', 'format', 'DoubleMat', 'mattype', 1)

    # }}}
