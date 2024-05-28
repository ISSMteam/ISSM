#Test Name: SquareShelfStressSSA2dHigherOrder
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

from ContourToMesh import *

md = triangle(model(), '../Exp/Square.exp', 150000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

field_names = []
field_tolerances = []
field_values = []
for i in ['P1bubble', 'P1bubblecondensed', 'P2']:
    md.flowequation.fe_SSA = i
    md = solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [1e-12, 1e-13, 1e-13, 1e-13]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]
