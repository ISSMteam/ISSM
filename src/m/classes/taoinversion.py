import numpy as np
from checkfield import checkfield
from IssmConfig import IssmConfig
from marshallcostfunctions import marshallcostfunctions
from fielddisplay import fielddisplay
from project3d import project3d
from supportedcontrols import *
from supportedcostfunctions import *
from WriteData import WriteData


class taoinversion(object):  # {{{
    """TAOINVERSION class definition

    Usage:
        inversion = taoinversion()
    """

    def __init__(self):
        self.iscontrol = 0
        self.incomplete_adjoint = 0
        self.control_parameters = np.nan
        self.maxsteps = 0
        self.maxiter = 0
        self.fatol = 0
        self.frtol = 0
        self.gatol = 0
        self.grtol = 0
        self.gttol = 0
        self.algorithm = ''
        self.cost_functions = np.nan
        self.cost_functions_coefficients = np.nan
        self.min_parameters = np.nan
        self.max_parameters = np.nan
        self.vx_obs = np.nan
        self.vy_obs = np.nan
        self.vz_obs = np.nan
        self.vel_obs = np.nan
        self.thickness_obs = np.nan
        self.surface_obs = np.nan

        self.setdefaultparameters()
    # }}}

    def __repr__(self):
        s = '   taoinversion parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        # s += '{}\n'.format(fielddisplay(self, 'mantle_viscosity', 'mantle viscosity constraints (NaN means no constraint) (Pa s)'))
        # s += '{}\n'.format(fielddisplay(self, 'lithosphere_thickness', 'lithosphere thickness constraints (NaN means no constraint) (m)'))
        # s += '{}\n'.format(fielddisplay(self, 'cross_section_shape', "1: square-edged, 2: elliptical - edged surface"))
        s += '{}\n'.format(fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(fielddisplay(self, 'fatol', 'convergence criterion: f(X) - f(X * ) (X: current iteration, X * : "true" solution, f: cost function)'))
        s += '{}\n'.format(fielddisplay(self, 'frtol', 'convergence criterion: |f(X) - f(X * )| / |f(X * )|'))
        s += '{}\n'.format(fielddisplay(self, 'gatol', 'convergence criterion: ||g(X)|| (g: gradient of the cost function)'))
        s += '{}\n'.format(fielddisplay(self, 'grtol', 'convergence criterion: ||g(X)|| / |f(X)|'))
        s += '{}\n'.format(fielddisplay(self, 'gttol', 'convergence criterion: ||g(X)|| / ||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(fielddisplay(self, 'algorithm', 'minimization algorithm: ''tao_blmvm'', ''tao_cg'', ''tao_lmvm'''))
        s += '{}\n'.format(fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
        s += '{}\n'.format(fielddisplay(self, 'surface_obs', 'observed surface elevation [m]'))
        s += '{}\n'.format('Available cost functions:')
        s += '{}\n'.format('   101: SurfaceAbsVelMisfit')
        s += '{}\n'.format('   102: SurfaceRelVelMisfit')
        s += '{}\n'.format('   103: SurfaceLogVelMisfit')
        s += '{}\n'.format('   104: SurfaceLogVxVyMisfit')
        s += '{}\n'.format('   105: SurfaceAverageVelMisfit')
        s += '{}\n'.format('   201: ThicknessAbsMisfit')
        s += '{}\n'.format('   501: DragCoefficientAbsGradient')
        s += '{}\n'.format('   502: RheologyBbarAbsGradient')
        s += '{}\n'.format('   503: ThicknessAbsGradient')
        return s

    def setdefaultparameters(self):
        #default is incomplete adjoint for now
        self.incomplete_adjoint = 1
        #parameter to be inferred by control methods (only
        #drag and B are supported yet)
        self.control_parameters = ['FrictionCoefficient']
        #number of iterations and steps
        self.maxsteps = 20
        self.maxiter = 30
        #default tolerances
        self.fatol = 0
        self.frtol = 0
        self.gatol = 0
        self.grtol = 0
        self.gttol = 1e-4
        #minimization algorithm
        PETSCMAJOR = IssmConfig('_PETSC_MAJOR_')[0]
        PETSCMINOR = IssmConfig('_PETSC_MINOR_')[0]
        if(PETSCMAJOR > 3 or (PETSCMAJOR == 3 and PETSCMINOR >= 5)):
            self.algorithm = 'blmvm'
        else:
            self.algorithm = 'tao_blmvm'
        #several responses can be used:
        self.cost_functions = 101
        return self

    def extrude(self, md):
        self.vx_obs = project3d(md, 'vector', self.vx_obs, 'type', 'node')
        self.vy_obs = project3d(md, 'vector', self.vy_obs, 'type', 'node')
        self.vel_obs = project3d(md, 'vector', self.vel_obs, 'type', 'node')
        self.thickness_obs = project3d(md, 'vector', self.thickness_obs, 'type', 'node')

        if numel(self.cost_functions_coefficients) > 1:
            self.cost_functions_coefficients = project3d(md, 'vector', self.cost_functions_coefficients, 'type', 'node')

        if numel(self.min_parameters) > 1:
            self.min_parameters = project3d(md, 'vector', self.min_parameters, 'type', 'node')

        if numel(self.max_parameters) > 1:
            self.max_parameters = project3d(md, 'vector', self.max_parameters, 'type', 'node')

        return self

    def checkconsistency(self, md, solution, analyses):
        if not self.iscontrol:
            return md
        if not IssmConfig('_HAVE_TAO_')[0]:
            md = md.checkmessage('TAO has not been installed, ISSM needs to be reconfigured and recompiled with TAO')

        num_controls = np.size(md.inversion.control_parameters)
        num_costfunc = np.size(md.inversion.cost_functions)

        md = checkfield(md, 'fieldname', 'inversion.iscontrol', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'inversion.incomplete_adjoint', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'inversion.control_parameters', 'cell', 1, 'values', supportedcontrols())
        md = checkfield(md, 'fieldname', 'inversion.maxsteps', 'numel', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.maxiter', 'numel', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.fatol', 'numel', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.frtol', 'numel', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.gatol', 'numel', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.grtol', 'numel', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.gttol', 'numel', 1, '>=', 0)

        PETSCMAJOR = IssmConfig('_PETSC_MAJOR_')[0]
        PETSCMINOR = IssmConfig('_PETSC_MINOR_')[0]
        if(PETSCMAJOR > 3 or (PETSCMAJOR == 3 and PETSCMINOR >= 5)):
            md = checkfield(md, 'fieldname', 'inversion.algorithm', 'values', ['blmvm', 'cg', 'lmvm'])
        else:
            md = checkfield(md, 'fieldname', 'inversion.algorithm', 'values', ['tao_blmvm', 'tao_cg', 'tao_lmvm'])

        md = checkfield(md, 'fieldname', 'inversion.cost_functions', 'size', [num_costfunc], 'values', supportedcostfunctions())
        md = checkfield(md, 'fieldname', 'inversion.cost_functions_coefficients', 'size', [md.mesh.numberofvertices, num_costfunc], '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.min_parameters', 'size', [md.mesh.numberofvertices, num_controls])
        md = checkfield(md, 'fieldname', 'inversion.max_parameters', 'size', [md.mesh.numberofvertices, num_controls])

        if solution == 'BalancethicknessSolution':
            md = checkfield(md, 'fieldname', 'inversion.thickness_obs', 'size', [md.mesh.numberofvertices], 'NaN', 1, 'Inf', 1)
        elif solution == 'BalancethicknessSoftSolution':
            md = checkfield(md, 'fieldname', 'inversion.thickness_obs', 'size', [md.mesh.numberofvertices], 'NaN', 1, 'Inf', 1)
        else:
            md = checkfield(md, 'fieldname', 'inversion.vx_obs', 'size', [md.mesh.numberofvertices], 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'inversion.vy_obs', 'size', [md.mesh.numberofvertices], 'NaN', 1, 'Inf', 1)

    def marshall(self, prefix, md, fid):

        yts = md.constants.yts
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'iscontrol', 'format', 'Boolean')
        WriteData(fid, prefix, 'name', 'md.inversion.type', 'data', 1, 'format', 'Integer')
        if not self.iscontrol:
            return
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'incomplete_adjoint', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'maxsteps', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'maxiter', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'fatol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'frtol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'gatol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'grtol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'gttol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'algorithm', 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'cost_functions_coefficients', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'min_parameters', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'max_parameters', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'vx_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'vy_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'vz_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'thickness_obs', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'surface_obs', 'format', 'DoubleMat', 'mattype', 1)

    #process control parameters
        num_control_parameters = np.size(self.control_parameters)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'control_parameters', 'format', 'StringArray')
        WriteData(fid, prefix, 'data', num_control_parameters, 'name', 'md.inversion.num_control_parameters', 'format', 'Integer')

    #process cost functions
        num_cost_functions = np.size(self.cost_functions)
        data = marshallcostfunctions(self.cost_functions)
        WriteData(fid, prefix, 'data', data, 'name', 'md.inversion.cost_functions', 'format', 'StringArray')
        WriteData(fid, prefix, 'data', num_cost_functions, 'name', 'md.inversion.num_cost_functions', 'format', 'Integer')
