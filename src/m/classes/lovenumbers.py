import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from getlovenumbers import getlovenumbers
from pairoptions import pairoptions
from WriteData import WriteData


class lovenumbers(object):  #{{{
    """lovenumbers class definition

    Usage:
        lovenumbers = lovenumbers()
        lovenumbers = lovenumbers('maxdeg', 10000, 'referenceframe', 'CF');

    Choose numbers of degrees required (1000 by default) and reference frame 
    (between CF and CM; CM by default)
    """

    def __init__(self, *args):  #{{{
        # Loading love numbers
        self.h = [] # Provided by PREM model
        self.k = [] # idem
        self.l = [] # idem

        # Tidal love numbers for computing rotational feedback
        self.th = []
        self.tk = []
        self.tl = []
        self.tk2secular = 0 # deg 2 secular number
        self.pmtf_colinear = []
        self.pmtf_ortho = []

        # Time/frequency for visco-elastic love numbers
        self.timefreq = []
        self.istime = 1

        options = pairoptions(*args)
        maxdeg = options.getfieldvalue('maxdeg', 1000)
        referenceframe = options.getfieldvalue('referenceframe', 'CM')
        self.setdefaultparameters(maxdeg, referenceframe)
    # }}}

    def __repr__(self):  #{{{
        s = '   lovenumbers parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'h', 'load Love number for radial displacement'))
        s += '{}\n'.format(fielddisplay(self, 'k', 'load Love number for gravitational potential perturbation'))
        s += '{}\n'.format(fielddisplay(self, 'l', 'load Love number for horizontal displacements'))
        s += '{}\n'.format(fielddisplay(self, 'th', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(fielddisplay(self, 'tk', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(fielddisplay(self, 'tl', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(fielddisplay(self, 'tk2secular', 'secular fluid Love number'))
        s += '{}\n'.format(fielddisplay(self, 'pmtf_colinear', 'Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)'))
        s += '{}\n'.format(fielddisplay(self, 'pmtf_ortho', 'Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)'))
        s += '{}\n'.format(fielddisplay(self, 'istime', 'time (default: 1) or frequency love numbers (0)'))
        s += '{}\n'.format(fielddisplay(self, 'timefreq', 'time/frequency vector (yr or 1/yr)'))
        s += '{}\n'.format(fielddisplay(self, 'pmtf_colinear', 'Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)'))
        s += '{}\n'.format(fielddisplay(self, 'pmtf_ortho', 'Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)'))
        return s
    # }}}

    def setdefaultparameters(self, maxdeg, referenceframe):  #{{{
        # Initialize love numbers
        self.h = getlovenumbers('type', 'loadingverticaldisplacement', 'referenceframe', referenceframe, 'maxdeg', maxdeg).reshape(-1,1)
        self.k = getlovenumbers('type', 'loadinggravitationalpotential', 'referenceframe', referenceframe, 'maxdeg', maxdeg).reshape(-1,1)
        self.l = getlovenumbers('type', 'loadinghorizontaldisplacement', 'referenceframe', referenceframe, 'maxdeg', maxdeg).reshape(-1,1)
        self.th = getlovenumbers('type', 'tidalverticaldisplacement', 'referenceframe', referenceframe, 'maxdeg', maxdeg).reshape(-1,1)
        self.tk = getlovenumbers('type', 'tidalgravitationalpotential', 'referenceframe', referenceframe, 'maxdeg', maxdeg).reshape(-1,1)
        self.tl = getlovenumbers('type', 'tidalhorizontaldisplacement', 'referenceframe', referenceframe, 'maxdeg', maxdeg).reshape(-1,1)

        # Secular fluid love number
        self.tk2secular = 0.942

        self.pmtf_colinear = np.array([0.0]).reshape(-1, 1)
        self.pmtf_ortho = np.array([0.0]).reshape(-1, 1)
        if maxdeg >= 2:
            self.pmtf_colinear = ((1.0 + self.k[2, :]) / (1.0 - self.tk[2, :] / self.tk2secular)).reshape(-1, 1) # Valid only for elastic regime, not viscous. Also neglects chandler wobble.
            self.pmtf_ortho = np.array([0.0]).reshape(-1, 1)
        # Time
        self.istime = 1 # Temporal love numbers by default
        self.timefreq = np.zeros(1) # Elastic case by default
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.h', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.k', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.l', 'NaN', 1, 'Inf', 1)

        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.th', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.tk', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.tl', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.tk2secular', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.pmtf_colinear', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.pmtf_ortho', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.timefreq', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.lovenumbers.istime', 'NaN', 1, 'Inf', 1, 'values', [0, 1])

        # Check that love numbers are provided at the same level of accuracy
        if (self.h.shape[0] != self.k.shape[0]) or (self.h.shape[0] != self.l.shape[0]):
            raise ValueError('lovenumbers error message: love numbers should be provided at the same level of accuracy')

        ntf = len(self.timefreq)
        if (np.shape(self.h)[1] != ntf or np.shape(self.k)[1] != ntf or np.shape(self.l)[1] != ntf or np.shape(self.th)[1] != ntf or np.shape(self.tk)[1] != ntf or np.shape(self.tl)[1] != ntf or np.shape(self.pmtf_colinear)[1] != ntf or np.shape(self.pmtf_ortho)[1] != ntf):
            raise ValueError('lovenumbers error message: love numbers should have as many time/frequency steps as the time/frequency vector')

        if self.istime and self.timefreq[0] != 0:
            raise ValueError('temporal love numbers must start with elastic response, i.e. timefreq[0] = 0')
        return md
    # }}}

    def defaultoutputs(self, md):  #{{{
        return[]
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'h', 'name', 'md.solidearth.lovenumbers.h', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'k', 'name', 'md.solidearth.lovenumbers.k', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'l', 'name', 'md.solidearth.lovenumbers.l', 'format', 'DoubleMat', 'mattype', 1)

        WriteData(fid, prefix, 'object', self, 'fieldname', 'th', 'name', 'md.solidearth.lovenumbers.th', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'tk', 'name', 'md.solidearth.lovenumbers.tk', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'tl', 'name', 'md.solidearth.lovenumbers.tl', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'pmtf_colinear','name','md.solidearth.lovenumbers.pmtf_colinear','format','DoubleMat','mattype',1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'pmtf_ortho','name','md.solidearth.lovenumbers.pmtf_ortho','format','DoubleMat','mattype',1)
        WriteData(fid, prefix, 'object', self, 'data', self.tk2secular, 'fieldname', 'lovenumbers.tk2secular', 'format', 'Double')

        if (self.istime):
            scale = md.constants.yts
        else:
            scale = 1.0 / md.constants.yts
        WriteData(fid, prefix, 'object', self, 'fieldname', 'istime', 'name', 'md.solidearth.lovenumbers.istime', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'timefreq', 'name', 'md.solidearth.lovenumbers.timefreq', 'format', 'DoubleMat', 'mattype', 1, 'scale', scale);
    # }}}

    def extrude(self, md):  #{{{
        return
    # }}}
