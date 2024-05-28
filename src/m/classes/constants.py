from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class constants(object):
    """CONSTANTS class definition

       Usage:
          constants = constants()
    """

    def __init__(self):  # {{{
        self.g                      = 0.
        self.omega                  = 0.
        self.yts                    = 0.
        self.referencetemperature   = 0.
        self.gravitational_constant = 0.;

        #set defaults
        self.setdefaultparameters()
    # }}}
    def __repr__(self):  # {{{
        s = '   constants parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'g', 'gravitational acceleration [m/s^2]'))
        s += '{}\n'.format(fielddisplay(self, 'omega', 'angular velocity of Earth [rad/s]'))
        s += '{}\n'.format(fielddisplay(self, 'yts', 'number of seconds in a year [s/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'referencetemperature', 'reference temperature used in the enthalpy model [K]'))
        s += '{}\n'.format(fielddisplay(self, 'gravitational_constant', 'Newtonian constant of gravitation [m^3/kg/s^2]'))
        return s
    # }}}
    def setdefaultparameters(self):  # {{{
        # Acceleration due to gravity (m / s^2)
        self.g = 9.81

        # Earth's rotation speed
        self.omega = 7.292 * 1e-5

        # Converstion from year to seconds
        self.yts = 365.0 * 24.0 * 3600.0

        # The reference temperature for enthalpy model (cf Aschwanden)
        self.referencetemperature = 223.15

        # Gravitational constant:
        self.gravitational_constant = 6.67259e-11

        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'constants.g', '>=', 0, 'size', [1]) # We allow 0 for validation tests
        md = checkfield(md, 'fieldname', 'constants.omega', '>=', 0, 'size', [1])
        md = checkfield(md, 'fieldname', 'constants.yts', '>', 0, 'size', [1])
        md = checkfield(md, 'fieldname', 'constants.referencetemperature', 'size', [1])
        md = checkfield(md, 'fieldname', 'constants.gravitational_constant','size',[1]);

        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'g', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'yts', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'referencetemperature', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'gravitational_constant', 'format', 'Double');
    # }}}
