#Test Name: SquareShelfTherStea

from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 180000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md.extrude(3, 1.)
md = setflowequation(md, 'SSA', 'all')
md.timestepping.time_step = 0
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Thermal')


# Fields and tolerances to track changes

field_names = ['Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-13, 6e-6]
field_values = [md.results.ThermalSolution.Temperature, md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate]
