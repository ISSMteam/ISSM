#Test Name: SquareShelfCMBMOLHO
from model import *
from socket import gethostname
import numpy as np
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from generic import generic
from SetMOLHOBC import SetMOLHOBC

md = triangle(model(), '../Exp/Square.exp', 200000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'MOLHO', 'all')

# control parameters

md.inversion.iscontrol = 1
md.inversion.control_parameters = ['MaterialsRheologyBbar']
md.inversion.min_parameters = 1.0e6 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = 2.0e9 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.nsteps = 2
md.inversion.cost_functions = [101]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.gradient_scaling = 1.0e7 * np.ones((md.inversion.nsteps, len(md.inversion.control_parameters)))
md.inversion.maxiter_per_step = 2. * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.3 * np.ones((md.inversion.nsteps))
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy


md.cluster = generic('name', gethostname(), 'np', 3)
md = SetMOLHOBC(md)
md = solve(md, 'Stressbalance')


# Fields and tolerances to track changes

field_names = ['Gradient', 'Misfits', 'MaterialsRheologyBbar', 'Pressure', 'Vel', 'Vx', 'Vy']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Gradient1,
                md.results.StressbalanceSolution.J,
                md.results.StressbalanceSolution.MaterialsRheologyBbar,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy]
