#Test Name: SquareSheetShelfGroundingLine2dSoft
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.initialization.vx[:] = 0.
md.initialization.vy[:] = 0.
md.geometry.base = -700. - (md.mesh.y - 500000.) / 1000.
md.geometry.bed = -700. - (md.mesh.y - 500000.) / 1000.
md.geometry.thickness[:] = 1300.
md.geometry.surface = md.geometry.base + md.geometry.thickness
md.transient.isstressbalance = 1
md.transient.isgroundingline = 1
md.groundingline.migration = 'AggressiveMigration'

md.timestepping.time_step = .1
md.timestepping.final_time = 1

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')
vel1 = md.results.TransientSolution[-1].Vel

#get same results with offset in bed and sea level:
md.geometry.base = -700. - (md.mesh.y - 500000.) / 1000.
md.geometry.bed = -700. - (md.mesh.y - 500000.) / 1000.
md.geometry.thickness[:] = 1300.
md.geometry.surface = md.geometry.base + md.geometry.thickness

md.geometry.base = md.geometry.base + 1000
md.geometry.bed = md.geometry.bed + 1000
md.geometry.surface = md.geometry.surface + 1000
md.solidearth.initialsealevel = 1000 * np.ones((md.mesh.numberofvertices, ))

md = solve(md, 'Transient', 'checkconsistency', 'no')
vel2 = md.results.TransientSolution[-1].Vel

#Fields and tolerances to track changes
field_names = ['Vel', 'Veloffset']
field_tolerances = [1e-13, 1e-13]
field_values = [vel1, vel2]
