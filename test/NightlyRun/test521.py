#Test Name: PigNudgingSSA2d
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from generic import generic
from inversionnudging import inversionnudging

md = triangle(model(), '../Exp/Pig.exp', 20000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

# Time stepping for nudging
md.timestepping.time_step = .5
md.timestepping.final_time = 1.
md.inversion = inversionnudging(md.inversion)
md.inversion.dhdt_obs = np.zeros((md.mesh.numberofvertices, 1))
md.inversion.maxiter = 3
# friction
md.inversion.H0_C = 10.  # default is 100 m
md.inversion.tau_C = 10
md.inversion.relaxation_C = .4
md.inversion.max_increment_C = 0.05  # in log10 space
md.inversion.min_C = 1 * np.ones((md.mesh.numberofvertices, 1))
md.inversion.max_C = 10000 * np.ones((md.mesh.numberofvertices, 1))
# melt
md.inversion.tau_melt = 200
md.inversion.H0_melt = 50
md.inversion.relaxation_melt = 0.3
md.inversion.max_increment_melt = 0.5
md.inversion.min_melt = -30 * np.ones((md.mesh.numberofvertices, 1))
md.inversion.max_melt = 50 * np.ones((md.mesh.numberofvertices, 1))
md.inversion.iscontrol = 1
md.verbose.solution = 1

md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['PerturbationMeltingRate', 'FrictionCoefficient', 'J']
field_tolerances = [1e-12, 2e-12, 1e-12]
field_values = [md.results.TransientSolution[0].BasalforcingsPerturbationMeltingRate,
                md.results.TransientSolution[0].FrictionCoefficient,
                md.results.TransientSolution[0].J]
