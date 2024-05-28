import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from lovenumbers import lovenumbers
from MatlabFuncs import *
from planetradius import planetradius
from rotational import rotational
from solidearthsettings import solidearthsettings
from solidearthsolution import solidearthsolution
from WriteData import WriteData


class solidearth(object):
    """SOLIDEARTH class definition

    Usage:
        solidearth = solidearth()
        solidearth = solidearth('earth')

    TODO:
    - Update translation from solidearth.m
    """

    def __init__(self, *args):  # {{{
        self.settings          = solidearthsettings()
        self.external          = None
        self.lovenumbers       = lovenumbers()
        self.rotational        = rotational()
        self.planetradius      = planetradius('earth')
        self.requested_outputs = []
        self.transfercount     = []
        self.transitions       = []
        self.partitionice      = []
        self.partitionhydro    = []
        self.partitionocean    = []

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters('earth')
        elif nargs == 1:
            self.setdefaultparameters(args[0])
        else:
            raise Exception('solidearth constructor error message: zero or one argument only!')
    # }}}

    def __repr__(self):  # {{{
        s = '   solidearthinputs, forcings and settings:\n'
        s += '{}\n'.format(fielddisplay(self, 'planetradius', 'planet radius [m]'))
        s += '{}\n'.format(fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(fielddisplay(self, 'transfercount', 'number of icecaps vertices are part of'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(fielddisplay(self, 'partitionice', 'ice partition vector for barystatic contribution'))
        s += '{}\n'.format(fielddisplay(self, 'partitionhydro', 'hydro partition vector for barystatic contribution'))
        s += '{}\n'.format(fielddisplay(self, 'partitionocean', 'ocean partition vector for barystatic contribution'))
        if not self.external:
            s += '{}\n'.format(fielddisplay(self, 'external', 'external solution, of the type solidearthsolution'))
        print(self.settings)
        print(self.lovenumbers)
        print(self.rotational)
        try:
            print(self.external)
        except TypeError:
            pass
        return s
    # }}}

    def setdefaultparameters(self, planet):  # {{{
        # Output default
        self.requested_outputs = ['default']

        # Transitions should be a list
        self.transitions = []
        self.transfercount = [0]

        # No partitions requested for barystatic contribution
        self.partitionice = []
        self.partitionhydro = []
        self.partitionocean = []

        # No external solutions by default
        self.external = None

        # Planet radius
        self.planetradius = planetradius(planet)
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return md

        md = checkfield(md, 'fieldname', 'solidearth.requested_outputs', 'stringrow', 1)

        self.settings.checkconsistency(md, solution, analyses)
        self.lovenumbers.checkconsistency(md, solution, analyses)
        self.rotational.checkconsistency(md, solution, analyses)
        if self.external:
            if not isa(self.external, solidearthsolution):
                raise Exception('solidearth consistency check: external field should be a solidearthsolution')
            self.external.checkconsistency(md, solution, analyses)
        return md
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['Sealevel']
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'planetradius', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'transitions', 'format', 'MatArray')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'transfercount', 'format', 'DoubleMat', 'mattype', 1)

        if len(self.partitionice):
            npartice = np.max(self.partitionice) + 2
        else:
            npartice = 0

        if len(self.partitionhydro):
            nparthydro = np.max(self.partitionhydro) + 2
        else:
            nparthydro = 0

        if len(self.partitionocean):
            npartocean = np.max(self.partitionocean) + 2
        else:
            npartocean = 0

        WriteData(fid, prefix, 'object', self, 'fieldname', 'partitionice', 'mattype', 1, 'format', 'DoubleMat');
        WriteData(fid, prefix, 'data', npartice, 'format', 'Integer', 'name', 'md.solidearth.npartice');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'partitionhydro', 'mattype', 1, 'format', 'DoubleMat');
        WriteData(fid, prefix, 'data', nparthydro,'format', 'Integer', 'name','md.solidearth.nparthydro');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'partitionocean', 'mattype', 1, 'format', 'DoubleMat');
        WriteData(fid, prefix, 'data', npartocean,'format', 'Integer', 'name','md.solidearth.npartocean');

        self.settings.marshall(prefix, md, fid)
        self.lovenumbers.marshall(prefix, md, fid)
        self.rotational.marshall(prefix, md, fid)
        if self.external:
            WriteData(fid, prefix, 'data', 1, 'format', 'Integer', 'name', 'md.solidearth.isexternal')
            self.external.marshall(prefix, md, fid)
        else:
            WriteData(fid, prefix, 'data', 0, 'format', 'Integer', 'name', 'md.solidearth.isexternal')

        # Process requested outputs
        outputs = self.requested_outputs
        pos = np.where(np.asarray(outputs) == 'default')[0]
        if len(pos):
            outputs = np.delete(outputs, pos) # remove 'default' from outputs
            outputs = np.append(outputs, self.defaultoutputs(md)) # add defaults
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.solidearth.requested_outputs', 'format', 'StringArray')
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}
