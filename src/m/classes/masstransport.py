import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class masstransport(object):
    """MASSTRANSPORT class definition

    Usage:
        masstransport = masstransport()
    """

    def __init__(self):  # {{{
        self.spcthickness = float('NaN')
        self.isfreesurface = 0
        self.min_thickness = 0.
        self.hydrostatic_adjustment = 0
        self.stabilization = 0
        self.vertex_pairing = float('NaN')
        self.penalty_factor = 0
        self.requested_outputs = []

        # Set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        s = '   Masstransport solution parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'spcthickness', 'thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'isfreesurface', 'do we use free surfaces (FS only) or mass conservation'))
        s += '{}\n'.format(fielddisplay(self, 'min_thickness', 'minimum ice thickness allowed [m]'))
        s += '{}\n'.format(fielddisplay(self, 'hydrostatic_adjustment', 'adjustment of ice shelves surface and bed elevations: ''Incremental'' or ''Absolute'' '))
        s += '{}\n'.format(fielddisplay(self, 'stabilization', '0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: discontinuous Galerkin, 4: flux corrected transport, 5: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.spcthickness = project3d(md, 'vector', self.spcthickness, 'type', 'node')
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['Thickness', 'Surface', 'Base']

    # }}}

    def setdefaultparameters(self):  # {{{
        # Type of stabilization to use 0:nothing 1:artificial_diffusivity 3:Discontinuous Galerkin
        self.stabilization = 1
        # Factor applied to compute the penalties kappa = max(stiffness matrix) * 1.0**penalty_factor
        self.penalty_factor = 3
        # Minimum ice thickness that can be used
        self.min_thickness = 1
        # Hydrostatic adjustment
        self.hydrostatic_adjustment = 'Absolute'
        # Default output
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if ('MasstransportAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.ismasstransport):
            return md

        md = checkfield(md, 'fieldname', 'masstransport.spcthickness', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'masstransport.isfreesurface', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'masstransport.hydrostatic_adjustment', 'values', ['Absolute', 'Incremental'])
        md = checkfield(md, 'fieldname', 'masstransport.stabilization', 'values', [0, 1, 2, 3, 4, 5])
        md = checkfield(md, 'fieldname', 'masstransport.min_thickness', '>', 0)
        md = checkfield(md, 'fieldname', 'masstransport.requested_outputs', 'stringrow', 1)
        if not np.any(np.isnan(self.vertex_pairing)) and len(self.vertex_pairing) > 0:
            md = checkfield(md, 'fieldname', 'stressbalance.vertex_pairing', '>', 0)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spcthickness', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isfreesurface', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_thickness', 'format', 'Double')
        WriteData(fid, prefix, 'data', self.hydrostatic_adjustment, 'format', 'String', 'name', 'md.masstransport.hydrostatic_adjustment')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stabilization', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vertex_pairing', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'penalty_factor', 'format', 'Double')

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.masstransport.requested_outputs', 'format', 'StringArray')
    # }}}
