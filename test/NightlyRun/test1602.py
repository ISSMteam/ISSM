#Test Name: SquareSheetShelfHORotation
import numpy as np
import sys
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
md.extrude(5, 1.)
md = setflowequation(md, 'HO', 'all')
md.stressbalance.spcvx[np.nonzero(md.mesh.y > 0.)] = np.nan
md.initialization.vx[:] = 0.
md.initialization.vy[:] = 0.
md.initialization.vel = np.zeros_like(md.initialization.vx)

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')
vel0 = md.results.StressbalanceSolution.Vel

theta = 30. * np.pi / 180.
x = md.mesh.x
y = md.mesh.y
md.mesh.x = np.cos(theta) * x - np.sin(theta) * y
md.mesh.y = np.sin(theta) * x + np.cos(theta) * y

md.stressbalance.referential[:, 0:3] = np.tile([np.cos(theta), np.sin(theta), 0], (md.mesh.numberofvertices, 1))
md.stressbalance.referential[:, 3:] = np.tile([0, 0, 1], (md.mesh.numberofvertices, 1))
md = solve(md, 'Stressbalance')
vel1 = md.results.StressbalanceSolution.Vel

#plotmodel(md, 'data', vel0, 'data', vel1, 'data', vel1 - vel0, 'title', 'Cartesian CS', 'title', 'Rotated CS', 'title', 'difference', 'view  #all', 2)
print("Error between Cartesian and rotated CS: {}".format(np.max(np.abs(vel0 - vel1)) / (np.max(np.abs(vel0)) + sys.float_info.epsilon)))

#Fields and tolerances to track changes
field_names = ['vel1']
field_tolerances = [1e-9]
field_values = [vel1]
