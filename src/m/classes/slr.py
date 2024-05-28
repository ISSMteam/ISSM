import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from MatlabFuncs import *
from model import *
from planetradius import planetradius
from WriteData import WriteData
from project3d import project3d


class slr(object):
    """SLR class definition

    Usage:
        slr = slr()
    """
    def __init__(self):  # {{{
        self.deltathickness = np.nan
        self.sealevel = np.nan
        self.spcthickness = np.nan
        self.maxiter = 0
        self.reltol = 0
        self.abstol = 0
        self.love_h = 0  #provided by PREM model()
        self.love_k = 0  #ideam
        self.love_l = 0  #ideam
        self.tide_love_k = 0  #ideam
        self.tide_love_h = 0  #ideam
        self.fluid_love = 0
        self.equatorial_moi = 0
        self.polar_moi = 0
        self.angular_velocity = 0
        self.rigid = 0
        self.elastic = 0
        self.rotation = 0
        self.ocean_area_scaling = 0
        self.hydro_rate = 0  #rate of steric expansion from hydrological effects.
        self.geodetic_run_frequency = 1  #how many time steps we skip before we run the geodetic part of the solver during transient
        self.geodetic = 0  #compute geodetic SLR? (in addition to steric?)
        self.degacc = 0
        self.horiz = 0
        self.planetradius = planetradius('earth')
        self.requested_outputs = []
        self.transitions = []

        #set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   slr parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'deltathickness', 'thickness change: ice height equivalent [m]'))
        s += '{}\n'.format(fielddisplay(self, 'sealevel', 'current sea level (prior to computation) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'spcthickness', 'thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'reltol', 'sea level rise relative convergence criterion, (NaN: not applied)'))
        s += '{}\n'.format(fielddisplay(self, 'abstol', 'sea level rise absolute convergence criterion, (default, NaN: not applied)'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter', 'maximum number of nonlinear iterations'))
        s += '{}\n'.format(fielddisplay(self, 'love_h', 'load Love number for radial displacement'))
        s += '{}\n'.format(fielddisplay(self, 'love_k', 'load Love number for gravitational potential perturbation'))
        s += '{}\n'.format(fielddisplay(self, 'love_l', 'load Love number for horizontal displaements'))
        s += '{}\n'.format(fielddisplay(self, 'tide_love_k', 'tidal load Love number (degree 2)'))
        s += '{}\n'.format(fielddisplay(self, 'tide_love_h', 'tidal load Love number (degree 2)'))
        s += '{}\n'.format(fielddisplay(self, 'fluid_love', 'secular fluid Love number'))
        s += '{}\n'.format(fielddisplay(self, 'equatorial_moi', 'mean equatorial moment of inertia [kg m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'polar_moi', 'polar moment of inertia [kg m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'angular_velocity', 'mean rotational velocity of earth [per second]'))
        s += '{}\n'.format(fielddisplay(self, 'ocean_area_scaling', 'correction for model representation of ocean area [default: No correction]'))
        s += '{}\n'.format(fielddisplay(self, 'hydro_rate', 'rate of hydrological expansion [mm / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'geodetic', 'compute geodetic SLR? (in addition to steric?) default 0'))
        s += '{}\n'.format(fielddisplay(self, 'geodetic_run_frequency', 'how many time steps we skip before we run SLR solver during transient (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'rigid', 'rigid earth graviational potential perturbation'))
        s += '{}\n'.format(fielddisplay(self, 'elastic', 'elastic earth graviational potential perturbation'))
        s += '{}\n'.format(fielddisplay(self, 'rotation', 'earth rotational potential perturbation'))
        s += '{}\n'.format(fielddisplay(self, 'degacc', 'accuracy (default .01 deg) for numerical discretization of the Green''s functions'))
        s += '{}\n'.format(fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Convergence criterion: absolute, relative and residual
        self.reltol = 0.01 # default
        self.abstol = np.nan # 1 mm of sea level rise
        # Maximum number of non-linear iterations
        self.maxiter = 5
        # Computational flags
        self.geodetic = 0
        self.rigid = 1
        self.elastic = 1
        self.ocean_area_scaling = 0
        self.rotation = 1
        # Tidal love numbers
        self.tide_love_h = 0.6149 # degree 2
        self.tide_love_k = 0.3055 # degree 2
        # Secular fluid love number
        self.fluid_love = 0.942
        # Moment of inertia
        self.equatorial_moi = 8.0077e37  # [kg m^2]
        self.polar_moi = 8.0345e37  # [kg m^2]
        # Mean rotational velocity of earth
        self.angular_velocity = 7.2921e-5  # [s^-1]
        # Numerical discretization accuracy
        self.degacc = 0.01
        # Hydro
        self.hydro_rate = 0
        # How many time steps we skip before we run SLR solver during transient
        self.geodetic_run_frequency = 1
        # Output default
        self.requested_outputs = ['default']
        # Transitions should be a list of vectors
        self.transitions = []
        # Horizontal displacement? (not on by default)
        self.horiz = 0
        # Earth area
        self.planetradius = planetradius('earth')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'SealevelchangeAnalysis' not in analyses or (solution == 'TransientSolution' and not md.transient.isslr):
            return md

        md = checkfield(md, 'fieldname', 'slr.deltathickness', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'slr.sealevel', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'slr.spcthickness', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'slr.love_h', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.love_k', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.love_l', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.tide_love_h', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.tide_love_k', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.fluid_love', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.equatorial_moi', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.polar_moi', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.angular_velocity', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'slr.reltol', 'size', [1, 1])
        md = checkfield(md, 'fieldname', 'slr.abstol', 'size', [1, 1])
        md = checkfield(md, 'fieldname', 'slr.maxiter', 'size', [1, 1], '>=', 1)
        md = checkfield(md, 'fieldname', 'slr.geodetic_run_frequency', 'size', [1, 1], '>=', 1)
        md = checkfield(md, 'fieldname', 'slr.hydro_rate', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'slr.degacc', 'size', [1, 1], '>=', 1e-10)
        md = checkfield(md, 'fieldname', 'slr.requested_outputs', 'stringrow', 1)
        md = checkfield(md, 'fieldname', 'slr.horiz', 'NaN', 1, 'Inf', 1, 'values', [0, 1])

        # Check that love numbers are provided at the same level of accuracy:
        if (size(self.love_h, 0) != size(self.love_k, 0)) or (size(self.love_h, 0) != size(self.love_l, 0)):
            raise Exception('slr error message: love numbers should be provided at the same level of accuracy')

        # Cross check that whereever we have an ice load, the mask is < 0 on each vertex:
        pos = np.where(self.deltathickness)
        maskpos = md.mask.ice_levelset[md.mesh.elements[pos, :]]
        els = np.where(maskpos > 0)
        if len(els[0]) > 0:
            print('Warning: slr.py::checkconsistency: there are elements with ice loads where some vertices are not on the ice!')

        # Check that if geodetic is requested, we are a mesh3dsurface model (planet), or if we are not, a coupler to a planet model is provided
        if self.geodetic:
            if md.transient.iscoupler:
                # We are good
                pass
            else:
                if md.mesh.__class__.__name__ == 'mesh3dsurface':
                    # We are good
                    pass
                else:
                    raise Exception('model is requesting geodetic computations without being a mesh3dsurface, or being coupled to one!')
        return md
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['Sealevel']
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deltathickness', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sealevel', 'mattype', 1, 'format', 'DoubleMat', 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spcthickness', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'reltol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'abstol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'maxiter', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'love_h', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'love_k', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'love_l', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'tide_love_h', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'tide_love_k', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'fluid_love', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'equatorial_moi', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'polar_moi', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'angular_velocity', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'rigid', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'elastic', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'rotation', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ocean_area_scaling', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'geodetic_run_frequency', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hydro_rate', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1e-3 / md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'degacc', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'transitions', 'format', 'MatArray')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'horiz', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'geodetic', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'planetradius', 'format', 'Double')

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.slr.requested_outputs', 'format', 'StringArray')
    # }}}
