import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class solidearthsettings(object):
    """solidearthsettings class definition

    Usage:
        solidearthsettings = solidearthsettings()
    """

    def __init__(self, *args):  # {{{
        self.reltol                 = 0
        self.abstol                 = 0
        self.maxiter                = 0
        self.selfattraction         = 1
        self.elastic                = 1
        self.viscous                = 1
        self.rotation               = 1
        self.grdocean               = 1
        self.ocean_area_scaling     = 0
        self.runfrequency           = 1 # How many time steps we skip before we run grd_core
        self.sealevelloading        = 1 # Will sea-level loads be computed?
        self.isgrd                  = 0 # Will GRD patterns be computed?
        self.compute_bp_grd         = 0 # Will GRD patterns for bottom pressures be computed?
        self.degacc                 = 0 # Degree increment for resolution of Green tables
        self.timeacc                = 0 # Time step accuracy required to compute Green tables
        self.horiz                  = 0 # Compute horizontal deformation?
        self.grdmodel               = 0 # GRD model (0 by default, 1 for (visco-)elastic, 2 for Ivins)
        self.cross_section_shape    = 0 # Cross section only used when GRD model is Ivins

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   solidearth settings:\n'
        s += '{}\n'.format(fielddisplay(self, 'reltol', 'sea level change relative convergence criterion (default, NaN: not applied)'))
        s += '{}\n'.format(fielddisplay(self, 'abstol', 'sea level change absolute convergence criterion (default, NaN: not applied)'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter', 'maximum number of nonlinear iterations'))
        s += '{}\n'.format(fielddisplay(self, 'grdocean', 'does this planet have an ocean, if set to 1: global water mass is conserved in GRD module (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'ocean_area_scaling', 'correction for model representation of ocean area (default: No correction)'))
        s += '{}\n'.format(fielddisplay(self, 'sealevelloading', 'enables surface loading from sea-level change (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'isgrd', 'compute GRD patterns (default: 1'))
        s += '{}\n'.format(fielddisplay(self, 'compute_bp_grd', 'compute GRD patterns for bottom pressure loads (default 1)'))
        s += '{}\n'.format(fielddisplay(self, 'runfrequency', 'how many time steps we skip before we run solidearthsettings solver during transient (default: 1)'))
        s += '{}\n'.format(fielddisplay(self, 'selfattraction', 'enables surface mass load to perturb the gravity field'))
        s += '{}\n'.format(fielddisplay(self, 'elastic', 'enables elastic deformation from surface loading'))
        s += '{}\n'.format(fielddisplay(self, 'viscous', 'enables viscous deformation from surface loading'))
        s += '{}\n'.format(fielddisplay(self, 'rotation', 'enables polar motion to feedback on the GRD fields'))
        s += '{}\n'.format(fielddisplay(self, 'degacc', 'accuracy (default: .01 deg) for numerical discretization of the Green\'s functions'))
        s += '{}\n'.format(fielddisplay(self, 'timeacc', 'time accuracy (default: 1 year) for numerical discretization of the Green\'s functions'))
        s += '{}\n'.format(fielddisplay(self, 'grdmodel', 'type of deformation model, 0 for no GRD, 1 for spherical GRD model (SESAW model), 2 for half-space planar GRD (visco-elastic model from Ivins)'))
        s += '{}\n'.format(fielddisplay(self, 'cross_section_shape', '1: square-edged (default). 2: elliptical. See iedge in GiaDeflectionCore'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Convergence criterion: absolute, relative, and residual
        self.reltol = 0.01  # 1 percent
        self.abstol = np.nan  # default

        # Maximum of non-linear iterations
        self.maxiter = 5

        # Computational flags
        self.selfattraction = 1
        self.elastic = 1
        self.viscous = 1
        self.rotation = 1
        self.grdocean = 1
        self.ocean_area_scaling = 0
        self.compute_bp_grd = 0
        self.isgrd = 0
        self.sealevelloading = 1

        # Numerical discretization accuracy
        self.degacc = 0.01
        self.timeacc = 1

        # How many time steps we skip before we run solidearthsettings solver during transient
        self.runfrequency = 1

        # Horizontal displacement? (not on by default)
        self.horiz = 0

        # Cross section for Ivins model
        self.cross_section_shape = 1 # Square as default (see iedge in GiaDeflectionCorex)

        # GRD model by default
        self.grdmodel = 1
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return md

        md = checkfield(md, 'fieldname', 'solidearth.settings.reltol', 'size', [1])
        md = checkfield(md, 'fieldname', 'solidearth.settings.abstol', 'size', [1])
        md = checkfield(md, 'fieldname', 'solidearth.settings.maxiter', 'size', [1], '>=', 1)
        md = checkfield(md, 'fieldname', 'solidearth.settings.runfrequency', 'size', [1], '>=', 1)
        md = checkfield(md, 'fieldname', 'solidearth.settings.degacc', 'size', [1], '>=', 1e-10)
        md = checkfield(md, 'fieldname', 'solidearth.settings.timeacc', 'size', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'solidearth.settings.horiz', 'NaN', 1, 'Inf', 1, 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'solidearth.settings.grdmodel', '>=', 0, '<=', 2)
        md = checkfield(md, 'fieldname', 'solidearth.settings.cross_section_shape', 'numel', [1], 'values', [1, 2])

        if self.elastic and not self.selfattraction:
            raise Exception('solidearthsettings checkconsistency error message: need selfattraction on if elastic flag is set')
        if self.viscous and not self.elastic:
            raise Exception('solidearthsettings checkconsistency error message: need elastic on if viscous flag is set')

        # A GRD computation has been requested, make some checks on the nature of the meshes provided
        if self.isgrd:
            if md.mesh.__class__.__name__ == 'mesh3dsurface':
                if self.grdmodel == 2:
                    raise Exception('model requires a 2D mesh to run gia Ivins computations (change mesh from mesh3dsurface to mesh2d)')
            else:
                if self.grdmodel == 1:
                    raise Exception('model requires a 3D surface mesh to run GRD computations (change mesh from mesh2d to mesh3dsurface)')
            if self.sealevelloading and not self.grdocean:
                raise Exception('solidearthsettings checkconsistency error message: need grdocean on if sealevelloading flag is set')

        if self.compute_bp_grd and not md.solidearth.settings.isgrd:
            raise Exception('solidearthsettings checkconsistency error message; if bottom pressure grd patterns are requested, solidearth settings class should have isgrd flag on')

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'reltol', 'name', 'md.solidearth.settings.reltol', 'format', 'Double');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'abstol', 'name', 'md.solidearth.settings.abstol', 'format', 'Double');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'maxiter', 'name', 'md.solidearth.settings.maxiter', 'format', 'Integer');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'selfattraction', 'name', 'md.solidearth.settings.selfattraction', 'format', 'Boolean');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'elastic', 'name', 'md.solidearth.settings.elastic', 'format', 'Boolean');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'viscous', 'name', 'md.solidearth.settings.viscous', 'format', 'Boolean');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'rotation', 'name', 'md.solidearth.settings.rotation', 'format', 'Boolean');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'grdocean', 'name', 'md.solidearth.settings.grdocean', 'format', 'Boolean');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ocean_area_scaling', 'name', 'md.solidearth.settings.ocean_area_scaling', 'format', 'Boolean');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'runfrequency', 'name', 'md.solidearth.settings.runfrequency', 'format', 'Integer');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'degacc', 'name', 'md.solidearth.settings.degacc', 'format', 'Double');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'timeacc', 'name', 'md.solidearth.settings.timeacc', 'format', 'Double', 'scale', md.constants.yts);
        WriteData(fid, prefix, 'object', self, 'fieldname', 'horiz', 'name', 'md.solidearth.settings.horiz', 'format', 'Integer');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sealevelloading', 'name', 'md.solidearth.settings.sealevelloading', 'format', 'Integer');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isgrd', 'name', 'md.solidearth.settings.isgrd', 'format', 'Integer');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'compute_bp_grd', 'name', 'md.solidearth.settings.compute_bp_grd', 'format', 'Integer');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'grdmodel', 'name', 'md.solidearth.settings.grdmodel', 'format', 'Integer');
        WriteData(fid, prefix, 'object', self, 'fieldname', 'cross_section_shape', 'name', 'md.solidearth.settings.cross_section_shape', 'format', 'Integer');
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}
