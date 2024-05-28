#Test Name: SquareShelfStressSSA2dEnhanced
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from matenhancedice import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md.materials = matenhancedice()
md.materials.rheology_B = 3.15e8 * np.ones(md.mesh.numberofvertices, )
md.materials.rheology_n = 3 * np.ones(md.mesh.numberofelements, )
md.materials.rheology_E = np.ones(md.mesh.numberofvertices, )
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure]
