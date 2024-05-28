from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class rotational(object):
    """ROTATIONAL class definition

    Usage:
        rotational = rotational()
    """

    def __init__(self, *args):  #{{{
        self.equatorialmoi = 0
        self.polarmoi = 0
        self.langularvelocity = 0

        nargin = len(args)
        if nargin == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   rotational parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'equatorialmoi', 'mean equatorial moment of inertia [kg m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'polarmoi', 'polar moment of inertia [kg m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'angularvelocity', 'mean rotational velocity of earth [per second]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Moment of inertia
        self.equatorialmoi = 8.0077 * pow(10, 37) # [kg m^2]
        self.polarmoi = 8.0345 * pow(10, 37) # [kg m^2]

        # Mean rotational velocity of earth
        self.angularvelocity = 7.2921 * pow(10, -5) # [s^-1]
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return md
        md = checkfield(md, 'fieldname', 'solidearth.rotational.equatorialmoi', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.rotational.polarmoi', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'solidearth.rotational.angularvelocity', 'NaN', 1, 'Inf', 1)
        return md
    # }}}

    def defaultoutputs(self, md):  #{{{
        return []
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'equatorialmoi', 'name', 'md.solidearth.rotational.equatorialmoi', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'polarmoi', 'name', 'md.solidearth.rotational.polarmoi', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'angularvelocity', 'name', 'md.solidearth.rotational.angularvelocity', 'format', 'Double')
    # }}}

    def extrude(self, md):  #{{{
        return self
    # }}}
