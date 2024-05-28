#Test Name: PigBalVel1
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
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Balancevelocity')

# Fields and tolerances to track changes
field_names = ['DrivingStressX', 'DrivingStressY', 'Vel']
field_tolerances = [1e-13, 1e-13, 1e-13]
field_values = [md.results.BalancevelocitySolution.DrivingStressX,
                md.results.BalancevelocitySolution.DrivingStressY,
                md.results.BalancevelocitySolution.Vel]
