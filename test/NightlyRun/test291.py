#Test Name: SquareShelfStressFSP4z
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
md.extrude(2, 1.)
md = setflowequation(md, 'FS', 'all')
md.flowequation.fe_FS = 'OneLayerP4z'
md.cluster = generic('name', gethostname(), 'np', 1)
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [5e-5, 5e-5, 8e-5, 5e-5, 1e-7]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
