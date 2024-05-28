import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class thermal(object):
    """THERMAL class definition

       Usage:
          thermal = thermal()
    """

    def __init__(self):  # {{{
        self.spctemperature = np.nan
        self.penalty_threshold = 0
        self.stabilization = 0
        self.reltol = 0
        self.maxiter = 0
        self.penalty_lock = 0
        self.penalty_factor = 0
        self.isenthalpy = 0
        self.isdynamicbasalspc = 0
        self.isdrainicecolumn = 0
        self.watercolumn_upperlimit = 0
        self.fe = 'P1'
        self.requested_outputs = []
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   Thermal solution parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'spctemperature', 'temperature constraints (NaN means no constraint) [K]'))
        s += '{}\n'.format(fielddisplay(self, 'stabilization', '0: no, 1: artificial_diffusivity, 2: SUPG'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter', 'maximum number of non linear iterations'))
        s += '{}\n'.format(fielddisplay(self, 'reltol', 'relative tolerance criterion'))
        s += '{}\n'.format(fielddisplay(self, 'penalty_lock', 'stabilize unstable thermal constraints that keep zigzagging after n iteration (default is 0, no stabilization)'))
        s += '{}\n'.format(fielddisplay(self, 'penalty_threshold', 'threshold to declare convergence of thermal solution (default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'isenthalpy', 'use an enthalpy formulation to include temperate ice (default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'isdynamicbasalspc', 'enable dynamic setting of basal forcing. required for enthalpy formulation (default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'isdrainicecolumn', 'wether waterfraction drainage is enabled for enthalpy formulation (default is 1)'))
        s += '{}\n'.format(fielddisplay(self, 'watercolumn_upperlimit', 'upper limit of basal watercolumn for enthalpy formulation (default is 1000m)'))
        s += '{}\n'.format(fielddisplay(self, 'fe', 'Finite Element type: ''P1'' (default), ''P1xP2'''))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.spctemperature = project3d(md, 'vector', self.spctemperature, 'type', 'node', 'layer', md.mesh.numberoflayers, 'padding', np.nan)
        if isinstance(md.initialization.temperature, np.ndarray) and np.size(md.initialization.temperature, axis=0) == md.mesh.numberofvertices:
            self.spctemperature = float('NaN') * np.ones((md.mesh.numberofvertices))
            pos = np.where(md.mesh.vertexonsurface)[0]
            self.spctemperature[pos] = md.initialization.temperature[pos]  #impose observed temperature on surface
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        if self.isenthalpy:
            return ['Enthalpy', 'Temperature', 'Waterfraction', 'Watercolumn', 'BasalforcingsGroundediceMeltingRate']
        else:
            return ['Temperature', 'BasalforcingsGroundediceMeltingRate']
    # }}}

    def setdefaultparameters(self):  # {{{
        #Number of unstable constraints acceptable
        self.penalty_threshold = 0
        #Type of stabilization used
        self.stabilization = 1
        #Relative tolerance for the enthalpy convergence
        self.reltol = 0.01
        #Maximum number of iterations
        self.maxiter = 100
        #factor used to compute the values of the penalties: kappa = max(stiffness matrix) * 1.0**penalty_factor
        self.penalty_factor = 3
        #Should we use cold ice (default) or enthalpy formulation
        self.isenthalpy = 0
        #will basal boundary conditions be set dynamically
        self.isdynamicbasalspc = 0
        #wether waterfraction drainage is enabled
        self.isdrainicecolumn = 1
        #set an upper limit for local stored watercolumn
        self.watercolumn_upperlimit = 1000
        #Finite element interpolation
        self.fe = 'P1'
        #default output
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if ('ThermalAnalysis' not in analyses and 'EnthalpyAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isthermal):
            return md
        md = checkfield(md, 'fieldname', 'thermal.stabilization', 'numel', [1], 'values', [0, 1, 2, 3])
        md = checkfield(md, 'fieldname', 'thermal.spctemperature', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md,'fieldname', 'thermal.fe', 'values', ['P1', 'P1xP2', 'P1xP3'])
        md = checkfield(md, 'fieldname', 'thermal.requested_outputs', 'stringrow', 1)
        if 'EnthalpyAnalysis' in analyses and md.thermal.isenthalpy and md.mesh.dimension() == 3:
            md = checkfield(md, 'fieldname', 'thermal.isdrainicecolumn', 'numel', [1], 'values', [0, 1])
            md = checkfield(md, 'fieldname', 'thermal.watercolumn_upperlimit', '>=', 0)

            TEMP = md.thermal.spctemperature[:-1].flatten()
            pos = np.where(~np.isnan(TEMP))
            try:
                spccol = np.size(md.thermal.spctemperature, 1)
            except IndexError:
                spccol = 1

            replicate = np.tile(md.geometry.surface - md.mesh.z, (spccol)).flatten()
            control = md.materials.meltingpoint - md.materials.beta * md.materials.rho_ice * md.constants.g * replicate + 1.0e-5
            md = checkfield(md, 'fieldname', 'thermal.spctemperature', 'field', md.thermal.spctemperature.flatten()[pos], '<=', control[pos], 'message', "spctemperature should be below the adjusted melting point")
            md = checkfield(md, 'fieldname', 'thermal.isenthalpy', 'numel', [1], 'values', [0, 1])
            md = checkfield(md, 'fieldname', 'thermal.isdynamicbasalspc', 'numel', [1], 'values', [0, 1])
            if(md.thermal.isenthalpy):
                if np.isnan(md.stressbalance.reltol):
                    md.checkmessage("for a steadystate computation, thermal.reltol (relative convergence criterion) must be defined!")
                md = checkfield(md, 'fieldname', 'thermal.reltol', '>', 0., 'message', "reltol must be larger than zero")
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spctemperature', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'penalty_threshold', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stabilization', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'reltol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'maxiter', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'penalty_lock', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'penalty_factor', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isenthalpy', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isdrainicecolumn', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'watercolumn_upperlimit', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'fe', 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isdynamicbasalspc', 'format', 'Boolean')

    #process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.thermal.requested_outputs', 'format', 'StringArray')
    # }}}
