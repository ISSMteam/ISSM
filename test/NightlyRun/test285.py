#Test Name: SquareShelfStressHOHigherOrder
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

from ContourToMesh import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

field_names = []
field_tolerances = []
field_values = []
for i in ['P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P1xP3', 'P2xP4']:
    md.flowequation.fe_HO = i
    md = solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vz' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [6.7e-08, 5e-08, 2e-08, 5e-08, 1e-13]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vz,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]
