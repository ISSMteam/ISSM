#Test Name: SquareSheetConstrainedEnthalpyStea
import numpy as np
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
md = setflowequation(md, 'SSA', 'all')
md.timestepping.time_step = 0.
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md.thermal.isenthalpy = 1
md.thermal.isdynamicbasalspc = 1

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Thermal')

#Fields and tolerances to track changes
field_names = ['Enthalpy', 'Waterfraction', 'Temperature']
field_tolerances = [1e-13, 5e-13, 1e-13]
field_values = [md.results.ThermalSolution.Enthalpy,
                md.results.ThermalSolution.Waterfraction,
                md.results.ThermalSolution.Temperature]
