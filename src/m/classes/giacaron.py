import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class giacaron(object):
    """
    GIA class definition for Caron model (Caron et al, Geophysical Journal 
    International, 2017)

       Usage:
          giacaron = giacaron()
    """

    def __init__(self, *args):  # {{{
        #Physical constants
        self.gravitational_constant         = np.nan
        self.surface_radius                 = np.nan
        self.core_mantle_boudary_radius     = np.nan
        self.inner_core_boudary_radius      = np.nan
        self.stress_norm                    = np.nan
        self.gravity_norm                   = np.nan
        self.radius_norm                    = np.nan

        #Numerical parameters
        self.allow_layer_deletion           = np.nan
        self.verbose_mode                   = np.nan

        #GIA problem setup
        self.forcing_type                   = np.nan
        self.isincompressible               = np.nan
        self.benchmark_mode                 = np.nan
        self.calculate_sea_level            = np.nan
        self.calculate_rotational_feedback  = np.nan
        self.subtract_present_day           = np.nan
        self.ntime                          = np.nan
        self.nphi                           = np.nan

        #Earth model
        



        nargin = len(args)

        if nargin == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}
