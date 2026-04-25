#Test Name: SquareSheetConstrainedStressSIASSA3d
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')
md.extrude(3, 1.)
md = setflowequation(md, 'SIA', '../Exp/SquareHalfRight.exp', 'fill', 'SSA')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
