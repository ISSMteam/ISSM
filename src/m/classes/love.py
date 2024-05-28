import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData

class love(object):
    """LOVE - class definition

    Usage:
        md.love = love()
    """

    def __init__(self):  # {{{
        self.nfreq = 0
        self.frequencies = 0
        self.sh_nmax = 0
        self.sh_nmin = 0
        self.g0 = 0
        self.r0 = 0
        self.mu0 = 0
        self.Gravitational_Constant = 0
        self.chandler_wobble = 0
        self.allow_layer_deletion = 0
        self.underflow_tol = 0
        self.pw_threshold = 0
        self.min_integration_steps = 0
        self.max_integration_dr = 0
        self.integration_scheme = 0
        self.istemporal = 0
        self.n_temporal_iterations = 0
        self.time = 0
        self.love_kernels = 0
        self.forcing_type = 0
        self.inner_core_boundary = 0
        self.core_mantle_boundary = 0
        self.complex_computation = 0
        self.quad_precision = 0
        self.hypergeom_table1 = 0
        self.hypergeom_table2 = 0
        self.hypergeom_table1 = 0
        self.hypergeom_nz = 0
        self.hypergeom_z = 0

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '{}\n'.format(fielddisplay(self, 'nfreq', 'number of frequencies sampled (default: 1, elastic) [Hz]'))
        s += '{}\n'.format(fielddisplay(self, 'frequencies', 'frequencies sampled (convention defaults to 0 for the elastic case) [Hz]'))
        s += '{}\n'.format(fielddisplay(self, 'sh_nmax', 'maximum spherical harmonic degree (default: 256, .35 deg, or 40 km at equator)'))
        s += '{}\n'.format(fielddisplay(self, 'sh_nmin', 'minimum spherical harmonic degree (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'g0', 'adimensioning constant for gravity (default: 10) [m/s^2]'))
        s += '{}\n'.format(fielddisplay(self, 'r0', 'adimensioning constant for radius (default: 6371*10^3) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'mu0', 'adimensioning constant for stress (default: 10^11) [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'Gravitational_Constant', 'Newtonian constant of gravitation (default: 6.67259e-11 [m^3 kg^-1 s^-2])'))
        s += '{}\n'.format(fielddisplay(self, 'chandler_wobble', 'includes the inertial terms for the chandler wobble in the rotational feedback love numbers, only for forcing_type=11 (default: 0) (/!\\ 1 has not been validated yet)'))
        s += '{}\n'.format(fielddisplay(self, 'allow_layer_deletion', 'allow for migration of the integration boundary with increasing spherical harmonics degree (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'underflow_tol', 'threshold of deep to surface love number ratio to trigger the deletion of layers (default: 1e-16)'))
        s += '{}\n'.format(fielddisplay(self, 'pw_threshold', 'if relative variation across frequencies is smaller than this ratio, the post-widder transform for time-dependent love numbers is bypassed (default (1e-3)'))
        s += '{}\n'.format(fielddisplay(self, 'min_integration_steps', 'minimum number of radial steps to propagate the yi system from the bottom to the top of each layer (default: 50)'))
        s += '{}\n'.format(fielddisplay(self, 'max_integration_dr', 'maximum length of radial steps to propagate the yi system from the bottom to the top of each layer (default: 10e3) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'istemporal', ['1 for time-dependent love numbers, 0 for frequency-dependent or elastic love numbers (default: 0)', 'If 1: use fourierlove function build_frequencies_from_time to meet consistency']))
        s += '{}\n'.format(fielddisplay(self, 'n_temporal_iterations', 'max number of iterations in the inverse Laplace transform. Also the number of spectral samples per time step requested (default: 8)'))
        s += '{}\n'.format(fielddisplay(self, 'time', 'time vector for deformation if istemporal (default: 0) [s]'))
        s += '{}\n'.format(fielddisplay(self, 'love_kernels', 'compute love numbers at depth? (default: 0)'))
        s += '{}\n'.format(fielddisplay(self, 'forcing_type', 'integer indicating the nature and depth of the forcing for the Love number calculation (default: 11):'))
        s += '{}\n'.format('                                                     1:  Inner core boundary -- Volumic Potential')
        s += '{}\n'.format('                                                     2:  Inner core boundary -- Pressure')
        s += '{}\n'.format('                                                     3:  Inner core boundary -- Loading')
        s += '{}\n'.format('                                                     4:  Inner core boundary -- Tangential traction')
        s += '{}\n'.format('                                                     5:  Core mantle boundary -- Volumic Potential')
        s += '{}\n'.format('                                                     6:  Core mantle boundary -- Pressure')
        s += '{}\n'.format('                                                     7:  Core mantle boundary -- Loading')
        s += '{}\n'.format('                                                     8:  Core mantle boundary -- Tangential traction')
        s += '{}\n'.format('                                                     9:  Surface -- Volumic Potential')
        s += '{}\n'.format('                                                     10: Surface -- Pressure')
        s += '{}\n'.format('                                                     11: Surface -- Loading')
        s += '{}\n'.format('                                                     12: Surface -- Tangential traction ')
        s += '{}\n'.format(fielddisplay(self, 'inner_core_boundary', 'interface index in materials.radius locating forcing. Only used for forcing_type 1--4 (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'core_mantle_boundary', 'interface index in materials.radius locating forcing. Only used for forcing_type 5--8 (default: 2)'))
        s += '{}\n'.format(fielddisplay(self, 'complex_computation', 'return love numbers as 0: real (useful for elastic or temporal forms), 1: complex numbers (useful for Fourier spectral form) (default: 0)'))
        s += '{}\n'.format(fielddisplay(self, 'quad_precision', 'toggle computation love numbers and post-widder transform with 32 digit precision, useful for temporal form (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'debug', 'outputs yi system matrix prior to solving (default: 0)'))
        s += '{}\n'.format(fielddisplay(self, 'hypergeom_table1', 'table 1 for hypergeometric function, only for EBM rheology (default: [1])'))
        s += '{}\n'.format(fielddisplay(self, 'hypergeom_table2', 'table 2 for hypergeometric function, only for EBM rheology (default: [1])'))
        s += '{}\n'.format(fielddisplay(self, 'hypergeom_nalpha', 'length of hypergeometric table, only for EBM rheology (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'hypergeom_nz', 'width of hypergeometric table, only for EBM rheology (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'hypergeom_z', 'abscissa for hypergeometric table, only for EBM rheology (default: [0])'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # We setup an elastic love number computation by default
        self.nfreq = 1
        self.frequencies = [0] # Hz
        self.sh_nmax = 256 # .35 degree, 40 km at the equator
        self.sh_nmin = 1
        # Work on Python script for computing g0 for given Earth's structure
        self.g0 = 9.81 # m/s^2
        self.r0 = 6371 * 1e3 # m
        self.mu0 = 1e11 # Pa
        self.Gravitational_Constant = 6.67259e-11 # m^3 kg^-1 s^-2
        self.chandler_wobble = 0
        self.allow_layer_deletion = 1
        self.underflow_tol = 1e-16 # Threshold of deep to surface love number ratio to trigger the deletion of layer
        self.pw_threshold = 1e-3 # If relative variation across frequencies is smaller than this ratio, the post-widder transform for time-dependent love numbers is bypassed 
        self.min_integration_steps=50
        self.max_integration_dr=1e4
        self.integration_scheme=1
        self.istemporal = 0
        self.n_temporal_iterations = 8
        self.time = [0] # s
        self.love_kernels = 0
        self.forcing_type = 11 # Surface loading
        self.inner_core_boundary = 1
        self.core_mantle_boundary = 2
        self.complex_computation = 0
        self.quad_precision = 0
        self.debug = 0
        self.hypergeom_table1 = 1
        self.hypergeom_table2 = 1
        self.hypergeom_nalpha = 1
        self.hypergeom_nz = 1
        self.hypergeom_z = 0
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'LoveAnalysis' not in analyses:
            return md

        md = checkfield(md, 'fieldname', 'love.nfreq', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.frequencies', 'NaN', 1, 'Inf', 1, 'numel', md.love.nfreq)
        md = checkfield(md, 'fieldname', 'love.sh_nmax', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.sh_nmin', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.g0', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.r0', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.mu0', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.Gravitational_Constant', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.chandler_wobble', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'love.allow_layer_deletion', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'love.underflow_tol', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.pw_threshold', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.min_integration_steps', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.max_integration_dr', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.integration_scheme', 'NaN', 1, 'Inf', 1, 'numel', 1, '>=', 0, '<=', 2)
        md = checkfield(md, 'fieldname', 'love.love_kernels', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'love.forcing_type', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0, '<=', 12)
        md = checkfield(md, 'fieldname', 'love.complex_computation', 'NaN', 1, 'Inf', 1, 'numel', 1, 'values', [0, 1])

        md = checkfield(md, 'fieldname', 'love.istemporal', 'values', [0, 1])

        md = checkfield(md, 'fieldname', 'love.hypergeom_nalpha', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.hypergeom_nz', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'love.hypergeom_z', 'NaN', 1, 'Inf', 1, 'numel', md.love.hypergeom_nz)
        md = checkfield(md, 'fieldname', 'love.hypergeom_table1', 'NaN', 1, 'Inf', 1, 'numel', md.love.hypergeom_nz * md.love.hypergeom_nalpha)
        md = checkfield(md, 'fieldname', 'love.hypergeom_table2', 'NaN', 1, 'Inf', 1, 'numel', md.love.hypergeom_nz * md.love.hypergeom_nalpha)

        if md.love.istemporal:
            md = checkfield(md, 'fieldname', 'love.n_temporal_iterations', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
            md = checkfield(md, 'fieldname', 'love.time', 'NaN', 1, 'Inf', 1, 'numel', md.love.nfreq / 2 / md.love.n_temporal_iterations)
        if md.love.sh_nmin <= 1 and (md.love.forcing_type == 1 or md.love.forcing_type == 5 or md.love.forcing_type == 9):
            raise RuntimeError('Degree 1 not supported for forcing type {}. Use sh_min >= 2 for this kind of calculation.'.format(md.love.forcing_type))

        if md.love.chandler_wobble  == 1:
            print('Warning: Chandler wobble in Love number calculator has not been validated yet')

        # Need 'litho' material
        if md.materials.__class__.__name__ != 'materials' or 'litho' not in md.materials.nature:
            raise RuntimeError('Need a \'litho\' material to run a Fourier Love number analysis')

        mat = np.where(np.array(md.materials.nature) == 'litho')[0]
        if md.love.forcing_type <= 4:
            md = checkfield(md, 'fieldname', 'love.inner_core_boundary', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0, '<=', md.materials[mat].numlayers)
        elif md.love.forcing_type <= 8:
            md = checkfield(md, 'fieldname', 'love.core_mantle_boundary', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0, '<=', md.materials[mat].numlayers)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'nfreq', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'frequencies', 'format', 'DoubleMat', 'mattype',3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sh_nmax', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sh_nmin', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'g0', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'r0', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'mu0', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'Gravitational_Constant', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'chandler_wobble', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'allow_layer_deletion', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'underflow_tol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'pw_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_integration_steps', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'max_integration_dr', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'integration_scheme', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'istemporal', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'n_temporal_iterations', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'complex_computation', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'quad_precision', 'format', 'Boolean')
        # Note: no need to marshall the time vector, we have frequencies
        WriteData(fid, prefix, 'object', self, 'fieldname', 'love_kernels', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'forcing_type', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'inner_core_boundary', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'core_mantle_boundary', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'debug', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hypergeom_table1', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hypergeom_table2', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hypergeom_nalpha', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hypergeom_nz', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hypergeom_z', 'format', 'DoubleMat', 'mattype', 3)
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def build_frequencies_from_time(self):  # {{{
        if not self.istemporal:
            raise RuntimeError('cannot build frequencies for temporal love numbers if love.istemporal==0')
        print('Temporal love numbers: Overriding md.love.nfreq and md.love.frequencies')
        self.nfreq = len(self.time) * 2 * self.n_temporal_iterations
        self.frequencies = np.zeros((self.nfreq,))
        for i in range(len(self.time)):
            for j in range(2 * self.n_temporal_iterations):
                if self.time[i] == 0:
                    self.frequencies[(i - 1) * 2 * self.n_temporal_iterations + j] = 0 # Convention to avoid marshalling infinite numbers
                else:
                    self.frequencies[(i - 1) * 2 * self.n_temporal_iterations + j] = j * np.log(2) / self.time[i] / 2 / np.pi
        return self
    # }}}
