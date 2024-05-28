#Test Name: SquareShelfCMBHO

from model import *
from socket import gethostname
import numpy as np
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 200000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md.extrude(3, 1.)
md = setflowequation(md, 'HO', 'all')


# control parameters
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['MaterialsRheologyBbar']
md.inversion.min_parameters = 1e6 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = 2. * 1e9 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.nsteps = 2
md.inversion.cost_functions = [101]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.gradient_scaling = 1e7 * np.ones((md.inversion.nsteps, len(md.inversion.control_parameters)))
md.inversion.maxiter_per_step = 2. * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.3 * np.ones((md.inversion.nsteps))
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')


# Fields and tolerances to track changes

field_names = ['Gradient', 'Misfits', 'MaterialsRheologyBbar', 'Pressure', 'Vel', 'Vx', 'Vy']
field_tolerances = [1e-07, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08]
field_values = [md.results.StressbalanceSolution.Gradient1,
                md.results.StressbalanceSolution.J,
                md.results.StressbalanceSolution.MaterialsRheologyBbar,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy]
