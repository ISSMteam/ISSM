#Test Name: 79NorthBedSlop2d
from model import *
from socket import gethostname

from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/79North.exp', 10000.)
md = setmask(md, '../Exp/79NorthShelf.exp', '')
md = parameterize(md, '../Par/79North.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'BedSlope')

#Fields and tolerances to track changes
field_names = ['BedSlopeX', 'BedSlopeY']
field_tolerances = [1e-13, 1e-13]
field_values = [md.results.BedSlopeSolution.BedSlopeX,
                md.results.BedSlopeSolution.BedSlopeY]
