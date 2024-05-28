#Test Name: SquareShelfStressSSA2dDamageUpdate
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from matdamageice import matdamageice
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md.materials = matdamageice()
md = parameterize(md, '../Par/SquareShelf.py')
md.damage.isdamage = 1
md.damage.D = 0. * np.ones(md.mesh.numberofvertices)
md.damage.spcdamage = np.nan * np.ones(md.mesh.numberofvertices)
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

md.stressbalance.requested_outputs = ['default', 'NewDamage']
md.damage.stress_threshold = 1.3e5
md.damage.kappa = 2.8

md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure', 'NewDamage']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.NewDamage]
