#Test Name: SquareSheetShelfStressHOHigherOrder
from model import *
from socket import gethostname

from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.extrude(5, 1.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

field_names = []
field_tolerances = []
field_values = []
for i in ['P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P1xP3', 'P2xP4']:
    md.flowequation.fe_HO = i
    md = solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vz' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [1e-07, 6e-08, 6e-08, 6e-08, 3e-13]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vz,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]
