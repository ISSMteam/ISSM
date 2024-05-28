import numpy as np
from pairoptions import pairoptions


def parametercontroldrag(md, *args):
    """
    PARAMETERCONTROLDRAG - parameterization for control method on drag

    It is possible to specify the number of steps, values for the
    minimum and maximum values of the drag, the
    kind of cm_responses to use or the the optscal.

    Usage:
       md = parametercontroldrag(md, varargin)

    Example:
      md = parametercontroldrag(md)
      md = parametercontroldrag(md, 'nsteps', 20, 'cm_responses', 0)
      md = parametercontroldrag(md, 'cm_min', 1, 'cm_max', 150, 'cm_jump', 0.99, 'maxiter', 20)
      md = parametercontroldrag(md, eps_cm', 1.0e-4, 'optscal', [1.0e7 1.0e8])

    See also PARAMETERCONTROLB
    """

    #process options
    options = pairoptions(*args)

    #control type
    md.inversion.control_parameters = 'FrictionCoefficient'

    #weights
    weights = options.getfieldvalue('weights', np.ones(md.mesh.numberofvertices))
    if np.size(weights) != md.mesh.numberofvertices:
        md.inversion.cost_functions_coefficients = np.ones(md.mesh.numberofvertices)
    else:
        md.inversion.cost_functions_coefficients = weights

    #nsteps
    nsteps = options.getfieldvalue('nsteps', 100)
    if (np.size(nsteps) != 1) | (nsteps <= 0) | (np.floor(nsteps) != nsteps):
        md.inversion.nsteps = 100
    else:
        md.inversion.nsteps = nsteps

    #cm_min
    cm_min = options.getfieldvalue('cm_min', np.ones(md.mesh.numberofvertices))
    if (np.size(cm_min) == 1):
        md.inversion.min_parameters = cm_min * np.ones(md.mesh.numberofvertices)
    elif (np.size(cm_min) == md.mesh.numberofvertices):
        md.inversion.min_parameters = cm_min
    else:
        md.inversion.min_parameters = cm_min

    #cm_max
    cm_max = options.getfieldvalue('cm_max', 250 * np.ones(md.mesh.numberofvertices))
    if (np.size(cm_max) == 1):
        md.inversion.max_parameters = cm_max * np.ones(md.mesh.numberofvertices)
    elif (np.size(cm_max) == md.mesh.numberofvertices):
        md.inversion.max_parameters = cm_max
    else:
        md.inversion.max_parameters = cm_max

    #eps_cm
    eps_cm = options.getfieldvalue('eps_cm', float('nan'))
    if (np.size(eps_cm) != 1 | eps_cm < 0):
        md.inversion.cost_function_threshold = float('nan')
    else:
        md.inversion.cost_function_threshold = eps_cm

    #maxiter
    maxiter = options.getfieldvalue('maxiter', 10 * np.ones(md.inversion.nsteps))
    if (np.any(maxiter < 0) | np.any(np.floor(maxiter) != maxiter)):
        md.inversion.maxiter_per_step = 10 * np.ones(md.inversion.nsteps)
    else:
        raise RuntimeError("not implemented yet, see below matlab lines")
    #md.inversion.maxiter_per_step = repmat(maxiter(:), md.inversion.nsteps, 1)
    #md.inversion.maxiter_per_step(md.inversion.nsteps + 1:end) = []

    #cm_jump
    cm_jump = options.getfieldvalue('cm_jump', 0.8 * np.ones(md.inversion.nsteps))
    if not np.isreal(cm_jump):
        md.inversion.step_threshold = 0.8 * np.ones(md.inversion.nsteps)
    else:
        raise RuntimeError("not implemented yet, see below matlab lines")
    #md.inversion.step_threshold = repmat(cm_jump(:), md.inversion.nsteps, 1)
    #md.inversion.step_threshold(md.inversion.nsteps + 1:end) = []

    #cm_responses
    found = 0
    if options.exist('cm_responses'):
        cm_responses = options.getfieldvalue('cm_responses')
        if not any(cm_responses not in [101, 105]):
            md.inversion.cost_functions = np.tile(cm_responses[:], (md.inversion.nsteps, 1))
            md.inversion.cost_functions[md.inversion.nsteps:] = []
            found = 1
    if not found:
        third = np.ceil(md.inversion.nsteps / 3)
        md.inversion.cost_functions = [103 * np.ones(third, 1),
                                       101 * np.ones(third, 1),
                                       np.tile([101, 101, 103, 101], (third, 1))]
        md.inversion.cost_functions[md.inversion.nsteps:] = []

    #optscal
    found = 0
    if options.exist('optscal'):
        optscal = options.getfieldvalue('optscal')
        if not any(optscal < 0):
            md.inversion.gradient_scaling = np.tile(optscal[:], (md.inversion.nsteps, 1))
            md.inversion.gradient_scaling[md.inversion.nsteps:] = []
            found = 1

    if not found:
        third = np.ceil(md.inversion.nsteps / 3)
        md.inversion.gradient_scaling = [50 * np.ones(3, 1),
                                         15 * np.ones(third - 3, 1),
                                         10 * np.ones(third, 1),
                                         np.tile([10, 10, 20, 10], (third, 1))]
        md.inversion.gradient_scaling[md.inversion.nsteps:] = []

    return md
