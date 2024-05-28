#Test Name: SquareSheetShelfTherStea
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.extrude(4, 1.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.timestepping.time_step = 0.
md = solve(md, 'Thermal')

#Fields and tolerances to track changes
field_names = ['Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-13, 1e-5]
field_values = [md.results.ThermalSolution.Temperature,
                md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate]
