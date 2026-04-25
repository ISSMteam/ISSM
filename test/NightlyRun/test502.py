#Test Name: PigStressHO
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Pig.exp', 20000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md.extrude(4, 0.9)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure']
field_tolerances = [1e-08, 1e-08, 1e-08, 1e-08, 1e-08]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
