import sys

import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
import MatlabFuncs as m
from project3d import project3d
from WriteData import WriteData


class stressbalance(object):
    """STRESSBALANCE class definition

    Usage:
        stressbalance = stressbalance()
    """

    def __init__(self):  # {{{
        self.spcvx = np.nan
        self.spcvy = np.nan
        self.spcvx_base = np.nan
        self.spcvy_base = np.nan
        self.spcvx_shear = np.nan
        self.spcvy_shear = np.nan
        self.spcvz = np.nan
        self.restol = 0
        self.reltol = 0
        self.abstol = 0
        self.ishydrologylayer = 0
        self.isnewton = 0
        self.FSreconditioning = 0
        #self.icefront = np.nan -- no longer in use
        self.maxiter = 0
        self.shelf_dampening = 0
        self.vertex_pairing = np.nan
        self.penalty_factor = np.nan
        self.rift_penalty_lock = np.nan
        self.rift_penalty_threshold = 0
        self.referential = np.nan
        self.loadingforce = np.nan
        self.requested_outputs = []

        # Set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        s = '   StressBalance solution parameters:\n'
        s += '      Convergence criteria:\n'
        s += '{}\n'.format(fielddisplay(self, 'restol', 'mechanical equilibrium residual convergence criterion'))
        s += '{}\n'.format(fielddisplay(self, 'reltol', 'velocity relative convergence criterion, NaN: not applied'))
        s += '{}\n'.format(fielddisplay(self, 'abstol', 'velocity absolute convergence criterion, NaN: not applied'))
        s += '{}\n'.format(fielddisplay(self, 'isnewton', '0: Picard\'s fixed point, 1: Newton\'s method, 2: hybrid'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter', 'maximum number of nonlinear iterations'))
        s += '      boundary conditions:\n'
        s += '{}\n'.format(fielddisplay(self, 'spcvx', 'x-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'spcvy', 'y-axis velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'spcvz', 'z-axis velocity constraint (NaN means no constraint) [m / yr]'))
        #s += '{}\n'.format(fielddisplay(self, 'icefront', 'segments on ice front list (last column 0: Air, 1: Water, 2: Ice'))
        s += '      MOLHO boundary conditions:\n'
        s += '{}\n'.format(fielddisplay(self, 'spcvx_base', 'x-axis basal velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'spcvy_base', 'y-axis basal velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'spcvx_shear', 'x-axis shear velocity constraint (NaN means no constraint) [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'spcvy_shear', 'y-axis shear velocity constraint (NaN means no constraint) [m / yr]'))
        s += '      Rift options:\n'
        s += '{}\n'.format(fielddisplay(self, 'rift_penalty_threshold', 'threshold for instability of mechanical constraints'))
        s += '{}\n'.format(fielddisplay(self, 'rift_penalty_lock', 'number of iterations before rift penalties are locked'))
        s += '      Penalty options:\n'
        s += '{}\n'.format(fielddisplay(self, 'penalty_factor', 'offset used by penalties: penalty = Kmax * 10^offset'))
        s += '{}\n'.format(fielddisplay(self, 'vertex_pairing', 'pairs of vertices that are penalized'))
        s += '      Hydrology layer:\n'
        s += '{}\n'.format(fielddisplay(self, 'ishydrologylayer', '(SSA only) 0: no subglacial hydrology layer in driving stress, 1: hydrology layer in driving stress'));
        s += '      Other:\n'
        s += '{}\n'.format(fielddisplay(self, 'shelf_dampening', 'use dampening for floating ice ? Only for FS model'))
        s += '{}\n'.format(fielddisplay(self, 'FSreconditioning', 'multiplier for incompressibility equation. Only for FS model'))
        s += '{}\n'.format(fielddisplay(self, 'referential', 'local referential'))
        s += '{}\n'.format(fielddisplay(self, 'loadingforce', 'loading force applied on each point [N/m^3]'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.spcvx = project3d(md, 'vector', self.spcvx, 'type', 'node')
        self.spcvy = project3d(md, 'vector', self.spcvy, 'type', 'node')
        self.spcvz = project3d(md, 'vector', self.spcvz, 'type', 'node')
        self.referential = project3d(md, 'vector', self.referential, 'type', 'node')
        self.loadingforce = project3d(md, 'vector', self.loadingforce, 'type', 'node')

        if md.flowequation.isMOLHO:
            self.spcvx_base = project3d(md, 'vector', self.spcvx_base, 'type', 'node')
            self.spcvy_base = project3d(md, 'vector', self.spcvy_base, 'type', 'node')
            self.spcvx_shear = project3d(md, 'vector', self.spcvx_shear, 'type', 'poly', 'degree', 4)
            self.spcvy_shear = project3d(md, 'vector', self.spcvy_shear, 'type', 'poly', 'degree', 4)

        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Maximum of non-linear iterations
        self.maxiter = 100

        # Convergence criterion: absolute, relative and residual
        self.restol = pow(10, -4)
        self.reltol = 0.01
        self.abstol = 10

        self.FSreconditioning = pow(10, 13)
        self.shelf_dampening = 0

        # Penalty factor applied kappa = max(stiffness matrix) * 1.0^penalty_factor
        self.penalty_factor = 3

        # Stop the iterations of rift if below a threshold
        self.rift_penalty_threshold = 0

        # In some solutions, it might be needed to stop a run when only a few 
        # constraints remain unstable. For thermal computation, this parameter 
        # is often used.
        self.rift_penalty_lock = 10

        # Output default
        self.requested_outputs = ['default']
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        if md.mesh.dimension() == 3:
            list = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
        else:
            list = ['Vx', 'Vy', 'Vel', 'Pressure']
        return list
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses:
            return md
        if solution == 'TransientSolution' and not md.transient.isstressbalance:
            return md

        md = checkfield(md, 'fieldname', 'stressbalance.spcvx', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'stressbalance.spcvy', 'Inf', 1, 'timeseries', 1)
        if m.strcmp(md.mesh.domaintype(), '3D'):
            md = checkfield(md, 'fieldname', 'stressbalance.spcvz', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'stressbalance.restol', 'size', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'stressbalance.reltol', 'size', [1])
        md = checkfield(md, 'fieldname', 'stressbalance.abstol', 'size', [1])
        md = checkfield(md, 'fieldname', 'stressbalance.ishydrologylayer', 'numel', [1], 'values', [0, 1]);
        md = checkfield(md, 'fieldname', 'stressbalance.isnewton', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'stressbalance.FSreconditioning', 'size', [1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'stressbalance.maxiter', 'size', [1], '>=', 1)
        md = checkfield(md, 'fieldname', 'stressbalance.referential', 'size', [md.mesh.numberofvertices, 6])
        md = checkfield(md, 'fieldname', 'stressbalance.loadingforce', 'size', [md.mesh.numberofvertices, 3])
        md = checkfield(md, 'fieldname', 'stressbalance.requested_outputs', 'stringrow', 1)
        if not np.any(np.isnan(self.vertex_pairing)) and len(self.vertex_pairing) > 0:
            md = checkfield(md, 'fieldname', 'stressbalance.vertex_pairing', '>', 0)
        # Singular solution
        if (not np.any(np.logical_or(np.logical_not(np.isnan(md.stressbalance.spcvx)), np.logical_not(np.isnan(md.stressbalance.spcvy))))) & (not np.any(md.mask.ocean_levelset>0)):
            print('\n !!! Warning: no spc applied, model might not be well posed if no basal friction is applied, check for solution crash\n')
        # CHECK THAT EACH LINES CONTAIN ONLY NAN VALUES OR NO NAN VALUES
        if np.any(np.logical_and(np.sum(np.isnan(md.stressbalance.referential), axis=1) != 0, np.sum(np.isnan(md.stressbalance.referential), axis=1) != 6)):
            md.checkmessage('Each line of stressbalance.referential should contain either only NaN values or no NaN values')
        # CHECK THAT THE TWO VECTORS PROVIDED ARE ORTHOGONAL
        if np.any(np.sum(np.isnan(md.stressbalance.referential), axis=1) == 0):
            pos = [i for i, item in enumerate(np.sum(np.isnan(md.stressbalance.referential), axis=1)) if item == 0]
            for item in md.stressbalance.referential[pos, :]:
                if np.abs(np.inner(item[0:2], item[3:5])) > sys.float_info.epsilon:
                    md.checkmessage('Vectors in stressbalance.referential (columns 1 to 3 and 4 to 6) must be orthogonal')
        # CHECK THAT NO rotation specified for FS Grounded ice at base
        if m.strcmp(md.mesh.domaintype(), '3D') and md.flowequation.isFS:
            pos = np.nonzero(np.logical_and(md.mask.ocean_levelset, md.mesh.vertexonbase))
            if np.any(np.logical_not(np.isnan(md.stressbalance.referential[pos, :]))):
                md.checkmessage('no referential should be specified for basal vertices of grounded ice')
        if md.flowequation.isMOLHO:
            md = checkfield(md, 'fieldname', 'stressbalance.spcvx_base', 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'stressbalance.spcvy_base', 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'stressbalance.spcvx_shear', 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'stressbalance.spcvy_shear', 'Inf', 1, 'timeseries', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{

        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'vertex_pairing', 'format', 'DoubleMat', 'mattype', 3)

        yts = md.constants.yts

        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'spcvx', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'spcvy', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'spcvz', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'restol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'reltol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'abstol', 'format', 'Double', 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'ishydrologylayer', 'format', 'Boolean');
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'isnewton', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'FSreconditioning', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'maxiter', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'shelf_dampening', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'penalty_factor', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'rift_penalty_lock', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'rift_penalty_threshold', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'referential', 'format', 'DoubleMat', 'mattype', 1)
        if isinstance(self.loadingforce, (list, tuple, np.ndarray)) and np.size(self.loadingforce, 1) == 3:
            WriteData(fid, prefix, 'data', self.loadingforce[:, 0], 'format', 'DoubleMat', 'mattype', 1, 'name', 'md.stressbalance.loadingforcex')
            WriteData(fid, prefix, 'data', self.loadingforce[:, 1], 'format', 'DoubleMat', 'mattype', 1, 'name', 'md.stressbalance.loadingforcey')
            WriteData(fid, prefix, 'data', self.loadingforce[:, 2], 'format', 'DoubleMat', 'mattype', 1, 'name', 'md.stressbalance.loadingforcez')
        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.stressbalance.requested_outputs', 'format', 'StringArray')
        # MOLHO
        if md.flowequation.isMOLHO:
            WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'spcvx_base', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'spcvy_base', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'spcvx_shear', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'stressbalance', 'fieldname', 'spcvy_shear', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
    # }}}
