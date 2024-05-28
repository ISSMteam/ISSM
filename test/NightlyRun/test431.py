#Test Name: SquareSheetShelfSteaEnthalpySSA3d
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
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.timestepping.time_step = 0.
md.thermal.isenthalpy = 1
md.thermal.isdynamicbasalspc = 1
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md = solve(md, 'Steadystate')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Temperature', 'Waterfraction', 'Enthalpy']
field_tolerances = [8e-4, 5e-4, 5e-4, 5e-4, 1e-13, 1e-4, 5e-4, 5e-4]
field_values = [md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.Waterfraction,
                md.results.SteadystateSolution.Enthalpy]
