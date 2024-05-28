#Test Name: 79NorthCMBalThicVxVy
import numpy as np
import copy
from model import *
from socket import gethostname
from triangle import *
from meshconvert import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/79North.exp', 10000.)
md = meshconvert(md)
md = setmask(md, '../Exp/79NorthShelf.exp', '')
md = parameterize(md, '../Par/79North.py')
md = setflowequation(md, 'SSA', 'all')

#Ice sheet only
md = model.extract(md, md.mask.ocean_levelset > 0.)
pos = np.nonzero(md.mesh.vertexonboundary)
md.balancethickness.spcthickness[pos] = md.geometry.thickness[pos]

#control parameters
md.inversion.thickness_obs = copy.deepcopy(md.geometry.thickness)
md.inversion.iscontrol = 1
md.inversion.nsteps = 2
md.inversion.control_parameters = ['Vx', 'Vy']
md.balancethickness.stabilization = 1
md.inversion.gradient_scaling = np.vstack((10. / md.constants.yts * np.ones((md.inversion.nsteps)), 10. / md.constants.yts * np.ones((md.inversion.nsteps)))).T
md.inversion.min_parameters = np.vstack((-2000. * np.ones((md.mesh.numberofvertices)), -2000. * np.ones((md.mesh.numberofvertices)))).T
md.inversion.max_parameters = np.vstack((+ 2000. * np.ones((md.mesh.numberofvertices)), +2000. * np.ones((md.mesh.numberofvertices)))).T
md.inversion.cost_functions = [201]
md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, len(md.inversion.cost_functions)))
md.inversion.maxiter_per_step = 4 * np.ones((md.inversion.nsteps))
md.inversion.step_threshold = 0.99 * np.ones((md.inversion.nsteps))

md.verbose.control = 1
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Balancethickness')

#Fields and tolerances to track changes
field_names = ['Gradient1', 'Gradient2', 'Misfits', 'Vx', 'Vy', 'Thickness']
field_tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12]
field_values = [md.results.BalancethicknessSolution.Gradient1,
                md.results.BalancethicknessSolution.Gradient2,
                md.results.BalancethicknessSolution.J,
                md.results.BalancethicknessSolution.Vx,
                md.results.BalancethicknessSolution.Vy,
                md.results.BalancethicknessSolution.Thickness]
