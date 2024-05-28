#Test Name: PigCMBFS
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Pig.exp', 11000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')

#impose hydrostatic equilibrium (required by Stokes)
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md.extrude(3, 1.)
md = setflowequation(md, 'FS', 'all')
md = md.extract(md.mask.ocean_levelset < 0.)

#control parameters
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['MaterialsRheologyBbar']
md.inversion.min_parameters = 10.**6 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = 2. * 10**9 * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.nsteps = 2
md.inversion.cost_functions = [101]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.gradient_scaling = 10.**8 * np.ones((md.inversion.nsteps, len(md.inversion.control_parameters)))
md.inversion.maxiter_per_step = 2. * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy

md.cluster = generic('name', gethostname(), 'np', 1)
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Gradient', 'Misfits', 'MaterialsRheologyB', 'Pressure', 'Vel', 'Vx', 'Vy']
field_tolerances = [5e-11, 5e-11, 5e-11, 1e-09, 1e-11, 5e-11, 1e-11]
field_values = [md.results.StressbalanceSolution.Gradient1,
                md.results.StressbalanceSolution.J,
                md.results.StressbalanceSolution.MaterialsRheologyBbar,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy]
