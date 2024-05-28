import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from marshallcostfunctions import marshallcostfunctions
from project3d import project3d
from supportedcontrols import supportedcontrols
from supportedcostfunctions import supportedcostfunctions
from WriteData import WriteData


class inversion(object):
    """INVERSION class definition

    Usage:
        inversion = inversion()
    """

    def __init__(self):  # {{{
        self.iscontrol = 0
        self.incomplete_adjoint = 0
        self.control_parameters = np.nan
        self.nsteps = 0
        self.maxiter_per_step = np.nan
        self.cost_functions = ''
        self.cost_functions_coefficients = np.nan
        self.gradient_scaling = np.nan
        self.cost_function_threshold = 0
        self.min_parameters = np.nan
        self.max_parameters = np.nan
        self.step_threshold = np.nan
        self.vx_obs = np.nan
        self.vy_obs = np.nan
        self.vz_obs = np.nan
        self.vel_obs = np.nan
        self.thickness_obs = np.nan
        self.surface_obs = np.nan

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   inversion parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'iscontrol', 'is inversion activated?'))
        s += '{}\n'.format(fielddisplay(self, 'incomplete_adjoint', '1: linear viscosity, 0: non - linear viscosity'))
        s += '{}\n'.format(fielddisplay(self, 'control_parameters', 'ex: {''FrictionCoefficient''}, or {''MaterialsRheologyBbar''}'))
        s += '{}\n'.format(fielddisplay(self, 'nsteps', 'number of optimization searches'))
        s += '{}\n'.format(fielddisplay(self, 'cost_functions', 'indicate the type of response for each optimization step'))
        s += '{}\n'.format(fielddisplay(self, 'cost_functions_coefficients', 'cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
        s += '{}\n'.format(fielddisplay(self, 'cost_function_threshold', 'misfit convergence criterion. Default is 1%, NaN if not applied'))
        s += '{}\n'.format(fielddisplay(self, 'maxiter_per_step', 'maximum iterations during each optimization step'))
        s += '{}\n'.format(fielddisplay(self, 'gradient_scaling', 'scaling factor on gradient direction during optimization, for each optimization step'))
        s += '{}\n'.format(fielddisplay(self, 'step_threshold', 'decrease threshold for misfit, default is 30%'))
        s += '{}\n'.format(fielddisplay(self, 'min_parameters', 'absolute minimum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'max_parameters', 'absolute maximum acceptable value of the inversed parameter on each vertex'))
        s += '{}\n'.format(fielddisplay(self, 'vx_obs', 'observed velocity x component [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vy_obs', 'observed velocity y component [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vel_obs', 'observed velocity magnitude [m/yr]'))
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
    # }}}

    def setdefaultparameters(self):  # {{{
        #default is incomplete adjoint for now
        self.incomplete_adjoint = 1
        #parameter to be inferred by control methods (only
        #drag and B are supported yet)
        self.control_parameters = 'FrictionCoefficient'
        #number of steps in the control methods
        self.nsteps = 20
        #maximum number of iteration in the optimization algorithm for
        #each step
        self.maxiter_per_step = 20 * np.ones(self.nsteps)
        #the inversed parameter is updated as follows:
        #new_par = old_par + gradient_scaling(n) * C * gradient with C in [0 1]
        #usually the gradient_scaling must be of the order of magnitude of the
        #inversed parameter (1.0e8 for B, 50 for drag) and can be decreased
        #after the first iterations
        self.gradient_scaling = 50 * np.ones((self.nsteps, 1))
        #several responses can be used:
        self.cost_functions = [101, ]
        #step_threshold is used to speed up control method. When
        #misfit(1) / misfit(0) < self.step_threshold, we go directly to
        #the next step
        self.step_threshold = 0.7 * np.ones(self.nsteps)  #30 per cent decrement
        #cost_function_threshold is a criteria to stop the control methods.
        #if J[n] - J[n - 1] / J[n] < criteria, the control run stops
        #NaN if not applied
        self.cost_function_threshold = np.nan  #not activated
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
        md = checkfield(md, 'fieldname', 'inversion.nsteps', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.maxiter_per_step', 'size', [md.inversion.nsteps], '>=', 0)
        md = checkfield(md, 'fieldname', 'inversion.step_threshold', 'size', [md.inversion.nsteps])
        md = checkfield(md, 'fieldname', 'inversion.cost_functions', 'size', [num_costfunc], 'values', supportedcostfunctions())
        if num_costfunc == 1:
            md.inversion.cost_functions_coefficients = np.squeeze(md.inversion.cost_functions_coefficients)
            md = checkfield(md, 'fieldname', 'inversion.cost_functions_coefficients', 'size', [md.mesh.numberofvertices], '>=', 0)
        else:
            md = checkfield(md, 'fieldname', 'inversion.cost_functions_coefficients', 'size', [md.mesh.numberofvertices, num_costfunc], '>=', 0)

        if num_controls == 1:
            md.inversion.gradient_scaling = np.squeeze(md.inversion.gradient_scaling)
            md.inversion.min_parameters = np.squeeze(md.inversion.min_parameters)
            md.inversion.max_parameters = np.squeeze(md.inversion.max_parameters)
            md = checkfield(md, 'fieldname', 'inversion.gradient_scaling', 'size', [md.inversion.nsteps])
            md = checkfield(md, 'fieldname', 'inversion.min_parameters', 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'inversion.max_parameters', 'size', [md.mesh.numberofvertices])
        else:
            md = checkfield(md, 'fieldname', 'inversion.gradient_scaling', 'size', [md.inversion.nsteps, num_controls])
            md = checkfield(md, 'fieldname', 'inversion.min_parameters', 'size', [md.mesh.numberofvertices, num_controls])
            md = checkfield(md, 'fieldname', 'inversion.max_parameters', 'size', [md.mesh.numberofvertices, num_controls])

        # Only SSA, HO and FS are supported right now
        if solution == 'StressbalanceSolution':
            if not (md.flowequation.isSSA or md.flowequation.isMOLHO or md.flowequation.isHO or md.flowequation.isFS or md.flowequation.isL1L2):
                md.checkmessage("'inversion can only be performed for SSA, MOLHO, HO or FS ice flow models")
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

        WriteData(fid, prefix, 'name', 'md.inversion.type', 'data', 0, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'iscontrol', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'incomplete_adjoint', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vel_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        if not self.iscontrol:
            return
        WriteData(fid, prefix, 'object', self, 'fieldname', 'nsteps', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'maxiter_per_step', 'format', 'IntMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'cost_functions_coefficients', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'gradient_scaling', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'cost_function_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_parameters', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'max_parameters', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'step_threshold', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vx_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vy_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vz_obs', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'thickness_obs', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'surface_obs', 'format', 'DoubleMat', 'mattype', 1)

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
