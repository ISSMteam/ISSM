#Test Name: 79NorthCMBalThic2dCG
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/79North.exp', 10000.)
md = setmask(md, '../Exp/79NorthShelf.exp', '')
md = parameterize(md, '../Par/79North.py')
md = setflowequation(md, 'SSA', 'all')

#control parameters
md.inversion.nsteps = 2
md.masstransport.stabilization = 1
md.inversion.iscontrol = 1
md.inversion.control_parameters = ['BalancethicknessThickeningRate']
md.inversion.thickness_obs = md.geometry.thickness
md.inversion.min_parameters = -50. * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.max_parameters = 50. * np.ones((md.mesh.numberofvertices, len(md.inversion.control_parameters)))
md.inversion.cost_functions = [201]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.gradient_scaling = 10. / md.constants.yts * np.ones((md.inversion.nsteps, len(md.inversion.control_parameters)))
md.inversion.maxiter_per_step = 4 * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Balancethickness')

#Fields and tolerances to track changes
field_names = ['Gradient', 'Misfits', 'BalancethicknessThickeningRate', 'Thickness']
field_tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12]
field_values = [md.results.BalancethicknessSolution.Gradient1,
                md.results.BalancethicknessSolution.J,
                md.results.BalancethicknessSolution.BalancethicknessThickeningRate,
                md.results.BalancethicknessSolution.Thickness]
