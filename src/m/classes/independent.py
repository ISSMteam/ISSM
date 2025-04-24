import numpy as np
from pairoptions import pairoptions
from fielddisplay import fielddisplay
from checkfield import checkfield
from MatlabFuncs import *


class independent(object):
    """independent class definition

    Usage:
        independent = independent()
    """

    def __init__(self, *args):  # {{{
        self.name = ''
        self.type = ''
        self.fos_forward_index = np.nan
        self.fov_forward_indices = np.array([])
        self.nods = 0
        self.min_parameters = np.nan
        self.max_parameters = np.nan
        self.control_scaling_factor = np.nan
        self.control_size = 0

        # Set defaults
        self.setdefaultparameters()

        # Use provided options to change fields
        options = pairoptions(*args)

        # Get other fields
        self = options.AssignObjectFields(self)

        if self.control_size == 0:
            self.control_size = 1
    # }}}

    def __repr__(self):  # {{{
        s = '   independent variable:\n'

        s += '{}\n'.format(fielddisplay(self, 'name', 'variable name (must match corresponding String)'))
        s += '{}\n'.format(fielddisplay(self, 'type', 'type of variable (\'vertex\' or \'scalar\')'))
        s += '{}\n'.format(fielddisplay(self, 'nods', 'size of independent variables'))
        s += '{}\n'.format(fielddisplay(self, 'control_size', 'number of timesteps'))
        s += '{}\n'.format(fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'control_scaling_factor', 'order of magnitude of each control (useful for multi-parameter optimization)'))
        if not np.isnan(self.fos_forward_index):
            s += '{}\n'.format(fielddisplay(self, 'fos_forward_index', 'index for fos_foward driver of ADOLC'))
        if np.any(np.logical_not(np.isnan(self.fov_forward_indices))):
            s += '{}\n'.format(fielddisplay(self, 'fov_forward_indices', 'indices for fov_foward driver of ADOLC'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Do nothing
        return self
    # }}}

    def checkconsistency(self, md, i, solution, analyses, driver):  # {{{
        if not np.isnan(self.fos_forward_index):
            if not strcmpi(driver, 'fos_forward'):
                raise TypeError('cannot declare an independent with a fos_forward_index when the driver is not fos_forward!')
            if self.nods == 0:
                raise TypeError('independent checkconsistency error: nods should be set to the size of the independent variable')

        if len(self.fov_forward_indices) > 0:
            if not strcmpi(driver, 'fov_forward'):
                raise TypeError('cannot declare an independent with fov_forward_indices when the driver is not fov_forward!')
            if self.nods == 0:
                raise TypeError('independent checkconsistency error: nods should be set to the size of the independent variable')
            md = checkfield(md, 'fieldname', 'autodiff.independents[%d].fov_forward_indices' % i, '>=', 1, '<=', self.nods)

        return md
    # }}}
