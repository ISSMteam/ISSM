#Test Name: SquareShelfStressHOFS3dTiling

from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 120000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf2.py')
md.extrude(2, 1.)
md = setflowequation(md, 'FS', '../Exp/SquareHalfRight.exp', 'fill', 'HO')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')

# Fields and tolerances to track changes

field_names = ['Vx', 'Vy', 'Vz', 'Vel']
field_tolerances = [1e-09, 1e-09, 5e-06, 1e-09]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel
]
