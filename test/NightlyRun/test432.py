#Test Name: SquareSheetShelfSteaEnthalpyHO
import numpy as np
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
md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.timestepping.time_step = 0.
md.thermal.isenthalpy = 1
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md = solve(md, 'Steadystate')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Temperature', 'Waterfraction', 'Enthalpy']
field_tolerances = [3e-09, 1e-09, 1e-09, 1e-09, 1e-13, 2e-10, 5e-10, 3e-9]
field_values = [md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.Waterfraction,
                md.results.SteadystateSolution.Enthalpy]
