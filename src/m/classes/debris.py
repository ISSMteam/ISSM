import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class debris(object):
    """debris class definition

    Usage:
        debris = debris()
    """

    def __init__(self, *args):  # {{{
        self.spcthickness = np.nan
        self.min_thickness = 0
        self.stabilization = 0
        self.packingfraction = 0
        self.removalmodel = 0
        self.displacementmodel = 0
        self.max_displacementvelocity = 0
        self.removal_slope_threshold = 0
        self.removal_stress_threshold = 0
        self.vertex_pairing = np.nan
        self.requested_outputs = []

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # TODO: Replace the following with constructor
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   debris solution parameters:\n'
        s += '{}\n'.format(fielddisplay(self,'spcthickness','debris thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(fielddisplay(self,'min_thickness','minimum debris thickness allowed [m]'))
        s += '{}\n'.format(fielddisplay(self,'packingfraction','fraction of debris covered in the ice'))
        s += '{}\n'.format(fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(fielddisplay(self,'removalmodel','frontal removal of debris. 0: no removal, 1: Slope-triggered debris removal, 2: driving-stress triggered debris removal'))
        s += '{}\n'.format(fielddisplay(self,'displacementmodel','debris displacement. 0: no displacement, 1: ...'))
        s += '{}\n'.format(fielddisplay(self,'max_displacementvelocity','maximum velocity of debris transport (v_ice + v_displacement) (m/a)'))
        s += '{}\n'.format(fielddisplay(self,'removal_slope_threshold','critical slope (degrees) for removalmodel (1)'))
        s += '{}\n'.format(fielddisplay(self,'removal_stress_threshold','critical stress (Pa) for removalmodel (2)'))

        s += '\n      {}\n'.format('Penalty options:')
        s += '{}\n'.format(fielddisplay(self,'vertex_pairing','pairs of vertices that are penalized'))
        s += '{}\n'.format(fielddisplay(self,'requested_outputs','additional outputs requested'))
        return s
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['DebrisThickness', 'DebrisMaskNodeActivation', 'VxDebris', 'VyDebris']
    # }}}

    def setdefaultparameters(self):  # {{{
        # Type of stabilization to use 0:nothing 1:artificial_diffusivity 3:Discontinuous Galerkin
        self.stabilization = 2

        # Minimum debris thickness that can be used
        self.min_thickness = 0

        # Fraction of debris covered in the ice
        self.packingfraction = 0.01

        # Type of frontal debris removal
        self.removalmodel = 0

        # Type of debris displacement
        self.displacementmodel = 0

        # Slope threshold for removalmodel (1)
        self.removal_slope_threshold = 0

        # Stress threshold for removalmodel (2)
        self.removal_stress_threshold = 0

        # Max velocity for displacementmodel (1)
        self.max_displacementvelocity = 0

        # Default output
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if not 'MasstransportAnalysis' in analyses or solution == 'TransientSolution' and not md.transient.isdebris:
            return md

        md = checkfield(md, 'fieldname', 'debris.spcthickness')
        md = checkfield(md, 'fieldname', 'debris.stabilization', 'values', [0, 1, 2, 3, 4, 5])
        md = checkfield(md, 'fieldname', 'debris.min_thickness', '>=', 0)
        md = checkfield(md, 'fieldname', 'debris.packingfraction', '>=', 0)
        md = checkfield(md, 'fieldname', 'debris.removalmodel', 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'debris.displacementmodel', 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'debris.max_displacementvelocity', '>=', 0)
        md = checkfield(md, 'fieldname', 'debris.removal_slope_threshold', '>=', 0)
        md = checkfield(md, 'fieldname', 'debris.removal_stress_threshold', '>=', 0)
        md = checkfield(md, 'fieldname', 'debris.requested_outputs', 'stringrow', 1)

        if not np.any(np.isnan(md.stressbalance.vertex_pairing)) and len(md.stressbalance.vertex_pairing) > 0:
            md = checkfield(md, 'fieldname', 'stressbalance.vertex_pairing', '>', 0)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spcthickness', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_thickness', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stabilization', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'removalmodel', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'displacementmodel', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'max_displacementvelocity', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'removal_slope_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'removal_stress_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'packingfraction', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vertex_pairing', 'format', 'DoubleMat', 'mattype', 3)

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.debris.requested_outputs', 'format', 'StringArray')
    # }}}

    def extrude(self, md):  #{{{
        self.spcthickness = project3d(md, 'vector', self.spcthickness, 'type', 'node')
        return
    # }}}
