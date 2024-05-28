#Test Name: SquareSheetShelfSteaEnthalpyRheologiesHO
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.timestepping.time_step = 0.
md.thermal.isenthalpy = 1
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices, ))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices, ))

#Go solve
field_names = []
field_tolerances = []
field_values = []
for i in ['LliboutryDuval', 'CuffeyTemperate']:
    print(' ')
    print(' ====== Testing rheology law: ' + i + ' = ')

    md.materials.rheology_law = i
    md = solve(md, 'Steadystate')
    field_names += ['Vx' + i, 'Vy' + i, 'Vz' + i, 'Vel' + i, 'Pressure' + i,
                    'Temperature' + i, 'Waterfraction' + i, 'Enthalpy' + i]
    field_tolerances += [2e-09, 1e-09, 1e-09, 1e-09, 1e-13, 2e-10, 6e-10, 1e-9]
    field_values += [md.results.SteadystateSolution.Vx,
                     md.results.SteadystateSolution.Vy,
                     md.results.SteadystateSolution.Vz,
                     md.results.SteadystateSolution.Vel,
                     md.results.SteadystateSolution.Pressure,
                     md.results.SteadystateSolution.Temperature,
                     md.results.SteadystateSolution.Waterfraction,
                     md.results.SteadystateSolution.Enthalpy]
