#Test Name: SquareShelfStressSSAHOTiling

from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md.extrude(3, 2.)
md = setflowequation(md, 'HO', '../Exp/SquareHalfRight.exp', 'fill', 'SSA')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')

# Fields and tolerances to track changes

field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [5e-09, 5e-09, 5e-09, 5e-09, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
