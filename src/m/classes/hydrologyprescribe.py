#!/usr/bin/env python3
import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData

class hydrologyprescribe(object):
    """
    HydrologyPrescribe class definition

    Usage:
      hydrologyprescribe = hydrologyprescribe()
    """

    def __init__(self): # {{{
        self.head = np.nan

        #set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        string = '   hydrologypism solution parameters:\n'
        string+= '   This module is to simulate effective pressure Neff using hydraulic head from external subglacial hydrology model\n'
        string+= '   Neff = rho_i g H - Pw\n'
        string+= '   Pw   = rho_w g (head - z_b)\n'
        string+= '   H: ice thickness (m) / head: hydrology head (m) / z_b: bedrock elevation'

        string+= "{}\n".format(fielddisplay(self, 'head', 'subglacial hydrology water head (m)'
        return string
        # }}}

    def extrude(self,md): # {{{
        return self
    # }}}

    def defaultoutputs(self,md): # {{{
        lists=['HydrologyHead','EffectivePressure']
        return lists
    # }}}

    def setdefaultparameters(self): # {{{
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{

        #Early return
        if 'HydrologyPrescribe' not in analyses:
            return

        if not np.isempty(md.initialization.hydraulic_potential):
            warnings.warn('WARN: md.initialization.hydrology is defined. However, this is not used for "hydrologyprescribe" model.')

        md = checkfield(md, 'fieldname', 'hydrology.head','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1)
    # }}}

    def marshall(self, prefix, md, fid): # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 10, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'head', 'format', 'DoubleMat', 'mattype', 1)

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray')
    # }}}

