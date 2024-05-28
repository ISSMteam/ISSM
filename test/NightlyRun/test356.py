#Test Name: TransientFrictionSchoof
import numpy as np

from frictionschoof import frictionschoof
from socket import gethostname
from model import *
from parameterize import parameterize
from setflowequation import setflowequation
from setmask import setmask
from solve import solve
from transient import transient
from triangle import triangle


md = triangle(model(), '../Exp/Square.exp', 200000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')
md = setflowequation(md, 'SSA', 'all')

# Use Schoof's law
Cmax = 0.8
md.friction = frictionschoof()
md.friction.m    = 1.0 / 3.0 * np.ones((md.mesh.numberofelements, 1))
md.friction.Cmax = Cmax * np.ones((md.mesh.numberofvertices, 1))
md.friction.C = 200 * np.ones((md.mesh.numberofvertices, 1))

# Control parameters
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['FrictionC']
md.inversion.min_parameters = 1. * np.ones((md.mesh.numberofvertices, 1))
md.inversion.max_parameters = 10000. * np.ones((md.mesh.numberofvertices, 1))
md.inversion.nsteps = 2
md.inversion.cost_functions = [102, 501]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, 2))
md.inversion.cost_functions_coefficients[:, 1] = 2e-7
md.inversion.gradient_scaling = 3. * np.ones((md.inversion.nsteps, 1))
md.inversion.maxiter_per_step = 2 * np.ones((md.inversion.nsteps, 1))
md.inversion.step_threshold = 0.3 * np.ones((md.inversion.nsteps, 1))
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs= md.initialization.vy

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Gradient', 'Misfits', 'FrictionC', 'Pressure', 'Vel', 'Vx', 'Vy']
field_tolerances = [1e-12, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [
    md.results.StressbalanceSolution.Gradient1,
    md.results.StressbalanceSolution.J,
    md.results.StressbalanceSolution.FrictionC,
    md.results.StressbalanceSolution.Pressure,
    md.results.StressbalanceSolution.Vel,
    md.results.StressbalanceSolution.Vx,
    md.results.StressbalanceSolution.Vy
]
