#Test Name: SquareShelfConstrained
from model import *
from socket import gethostname
import numpy as np
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from paterson import *
from solve import *
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 150000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

# redo the parameter file for this special shelf.
# constant thickness, constrained (vy = 0) flow into an icefront,
# from 0 m / yr at the grounding line.

# tighten
md.stressbalance.restol = 1e-4

# needed later
ymin = min(md.mesh.y)
ymax = max(md.mesh.y)
xmin = min(md.mesh.x)
xmax = max(md.mesh.x)

di = md.materials.rho_ice / md.materials.rho_water

h = 1000.
md.geometry.thickness = h * np.ones((md.mesh.numberofvertices))
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness

# Initial velocity and pressure
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

# Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

# Boundary conditions:
md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))

# constrain flanks to 0 normal velocity
pos = np.where(np.logical_or(md.mesh.x == xmin, md.mesh.x == xmax))
md.stressbalance.spcvx[pos] = 0
md.stressbalance.spcvz[pos] = np.nan

# constrain grounding line to 0 velocity
pos = np.where(md.mesh.y == ymin)
md.stressbalance.spcvx[pos] = 0
md.stressbalance.spcvy[pos] = 0

# icefront
nodeonicefront = np.zeros(md.mesh.numberofvertices)
pos = np.where(md.mesh.y == ymax)
nodeonicefront[pos] = 1
md.mask.ice_levelset = -1 + nodeonicefront

md = solve(md, 'Stressbalance')

# create analytical solution: strain rate is constant = ((rho_ice * g * h) / 4B)^3 (Paterson, 4th Edition, page 292.
# ey_c=(md.materials.rho_ice * md.constants.g * (1 - di) * md.geometry.thickness. / (4 * md.materials.rheology_B)).^3
# vy_c = ey_c. * md.mesh.y * md.constants.yts

# Fields and tolerances to track changes
field_names = ['Vy']
field_tolerances = [1e-13]
field_values = [md.results.StressbalanceSolution.Vy]
