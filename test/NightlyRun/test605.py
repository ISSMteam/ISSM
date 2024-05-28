#Test Name: 79NorthSurfSlop3d
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
md.extrude(5, 1.5)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'SurfaceSlope')

#Fields and tolerances to track changes
field_names = ['SurfaceSlopeX', 'SurfaceSlopeY']
field_tolerances = [1e-13, 1e-13]
field_values = [md.results.SurfaceSlopeSolution.SurfaceSlopeX,
                md.results.SurfaceSlopeSolution.SurfaceSlopeY]
