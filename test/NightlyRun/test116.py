#Test Name: SquareShelfConstrainedBalThic2d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 150000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
#Add boundary conditions on thickness on the border
pos = np.nonzero(md.mesh.vertexonboundary)
md.balancethickness.spcthickness[pos] = md.geometry.thickness[pos]
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Balancethickness')

#Fields and tolerances to track changes
field_names = ['Thickness']
field_tolerances = [1e-13]
field_values = [md.results.BalancethicknessSolution.Thickness]
