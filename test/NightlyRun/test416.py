#Test Name: SquareSheetShelfCMDragSteaHO
import numpy as np
from model import *
from socket import gethostname

from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 170000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.extrude(3, 1.)
md = setflowequation(md, 'HO', 'all')

#control parameters
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['FrictionCoefficient']
md.inversion.min_parameters = 1. * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = 200. * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.nsteps = 2
md.inversion.cost_functions = [102, 501]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, 2))
md.inversion.cost_functions_coefficients[:, 1] = 2. * 10**-7
md.inversion.gradient_scaling = 3. * np.ones((md.inversion.nsteps, len(md.inversion.control_parameters)))
md.inversion.maxiter_per_step = 2 * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.3 * np.ones((md.inversion.nsteps))
md.timestepping.time_step = 0.
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy


md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Steadystate')

#Fields and tolerances to track changes
field_names = ['Gradient', 'Misfits', 'FrictionCoefficient', 'Pressure', 'Vel', 'Vx', 'Vy', 'Vz', 'Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-08, 1e-07, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-07, 1e-08, 1e-05]
field_values = [md.results.SteadystateSolution.Gradient1,
                md.results.SteadystateSolution.J,
                md.results.SteadystateSolution.FrictionCoefficient,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate]
