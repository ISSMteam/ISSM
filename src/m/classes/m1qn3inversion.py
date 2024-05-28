import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from marshallcostfunctions import marshallcostfunctions
from project3d import project3d
from supportedcontrols import supportedcontrols
from supportedcostfunctions import supportedcostfunctions
from WriteData import WriteData


class m1qn3inversion(object):
    """M1QN3 class definition

    Usage:
        m1qn3inversion = m1qn3inversion()
    """

    def __init__(self, *args):  # {{{
        if not len(args):
            print('empty init')
            self.iscontrol                   = 0
            self.incomplete_adjoint          = 0
            self.control_parameters          = np.nan
            self.control_scaling_factors     = np.nan
            self.maxsteps                    = 0
            self.maxiter                     = 0
            self.dxmin                       = 0.
            self.dfmin_frac                  = 0.
            self.gttol                       = 0.
            self.cost_functions              = np.nan
            self.cost_functions_coefficients = np.nan
            self.min_parameters              = np.nan
            self.max_parameters              = np.nan
            self.vx_obs                      = np.nan
            self.vy_obs                      = np.nan
            self.vz_obs                      = np.nan
            self.vel_obs                     = np.nan
            self.thickness_obs               = np.nan

            self.setdefaultparameters()
        elif len(args) == 1 and args[0].__module__ == 'inversion':
            print('converting inversion to m1qn3inversion')
            inv = args[0]
            #first call setdefaultparameters:
            self.setdefaultparameters()

            #then go fish whatever is available in the inversion object provided to the constructor
            self.iscontrol = inv.iscontrol
            self.incomplete_adjoint = inv.incomplete_adjoint
            self.control_parameters = inv.control_parameters
            self.maxsteps = inv.nsteps
            self.cost_functions = inv.cost_functions
            self.cost_functions_coefficients = inv.cost_functions_coefficients
            self.min_parameters = inv.min_parameters
            self.max_parameters = inv.max_parameters
            self.vx_obs = inv.vx_obs
            self.vy_obs = inv.vy_obs
            self.vz_obs = inv.vz_obs
            self.vel_obs = inv.vel_obs
            self.thickness_obs = inv.thickness_obs
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   m1qn3inversion parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(fielddisplay(self, 'control_parameters', 'ex: [''FrictionCoefficient''], or [''MaterialsRheologyBbar'']'))
        s += '{}\n'.format(fielddisplay(self, 'control_scaling_factors', 'order of magnitude of each control (useful for multi - parameter optimization)'))
        s += '{}\n'.format(fielddisplay(self, 'maxsteps', 'maximum number of iterations (gradient computation)'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter', 'maximum number of Function evaluation (forward run)'))
        s += '{}\n'.format(fielddisplay(self, 'dxmin', 'convergence criterion: two points less than dxmin from eachother (sup - norm) are considered identical'))
        s += '{}\n'.format(fielddisplay(self, 'dfmin_frac', 'expected reduction of J during the first step (e.g., 0.3=30% reduction in cost function)'))
        s += '{}\n'.format(fielddisplay(self, 'gttol', '||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)'))
        s += '{}\n'.format(fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'vx_obs', 'observed velocity x component [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vy_obs', 'observed velocity y component [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m / yr]'))
        s += '{}\n'.format(fielddisplay(self, 'thickness_obs', 'observed thickness [m]'))
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
    # }}}

    def setdefaultparameters(self):  # {{{
        #default is incomplete adjoint for now
        self.incomplete_adjoint = 1
        #parameter to be inferred by control methods (only
        #drag and B are supported yet)
        self.control_parameters = 'FrictionCoefficient'
        #Scaling factor for each control
        self.control_scaling_factors = 1
        #number of iterations
        self.maxsteps = 20
        self.maxiter = 40
        #several responses can be used:
        self.cost_functions = 101
        #m1qn3 parameters
        self.dxmin = 0.1
        self.dfmin_frac = 1.
        self.gttol = 1e-4

        return self
    # }}}

    def extrude(self, md):  # {{{
        self.vx_obs = project3d(md, 'vector', self.vx_obs, 'type', 'node')
        self.vy_obs = project3d(md, 'vector', self.vy_obs, 'type', 'node')
        self.vel_obs = project3d(md, 'vector', self.vel_obs, 'type', 'node')
        self.thickness_obs = project3d(md, 'vector', self.thickness_obs, 'type', 'node')
        if not np.any(np.isnan(self.cost_functions_coefficients)):
            self.cost_functions_coefficients = project3d(md, 'vector', self.cost_functions_coefficients, 'type', 'node')
        if not np.any(np.isnan(self.min_parameters)):
            self.min_parameters = project3d(md, 'vector', self.min_parameters, 'type', 'node')
        if not np.any(np.isnan(self.max_parameters)):
            self.max_parameters = project3d(md, 'vector', self.max_parameters, 'type', 'node')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if not self.iscontrol:
            return md

        num_controls = np.size(md.inversion.control_parameters)
        num_costfunc = np.size(md.inversion.cost_functions)

        md = checkfield(md, 'fieldname', 'inversion.iscontrol', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'inversion.incomplete_adjoint', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'inversion.control_parameters', 'cell', 1, 'values', supportedcontrols())
        md = checkfield(md, 'fieldname', 'inversion.control_scaling_factors', 'size', [num_controls], '>', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'inversion.maxsteps', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.maxiter', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.dxmin', 'numel', [1], '>', 0.)
        md = checkfield(md, 'fieldname', 'inversion.dfmin_frac', 'numel', [1], '>=', 0., '<=', 1.)
        md = checkfield(md, 'fieldname', 'inversion.gttol', 'numel', [1], '>', 0.)
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
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'iscontrol', 'format', 'Boolean')
        WriteData(fid, prefix, 'name', 'md.inversion.type', 'data', 2, 'format', 'Integer')
        if not self.iscontrol:
            return
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'incomplete_adjoint', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'control_scaling_factors', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'maxsteps', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'maxiter', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'dxmin', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'dfmin_frac', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'gttol', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'cost_functions_coefficients', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'min_parameters', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'max_parameters', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'vx_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'vy_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'vz_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'inversion', 'fieldname', 'thickness_obs', 'format', 'DoubleMat', 'mattype', 1)

        # Process control parameters
        num_control_parameters = len(self.control_parameters)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'control_parameters', 'format', 'StringArray')
        WriteData(fid, prefix, 'data', num_control_parameters, 'name', 'md.inversion.num_control_parameters', 'format', 'Integer')

        # Process cost functions
        num_cost_functions = np.size(self.cost_functions)
        data = marshallcostfunctions(self.cost_functions)
        WriteData(fid, prefix, 'data', data, 'name', 'md.inversion.cost_functions', 'format', 'StringArray')
        WriteData(fid, prefix, 'data', num_cost_functions, 'name', 'md.inversion.num_cost_functions', 'format', 'Integer')
    # }}}
