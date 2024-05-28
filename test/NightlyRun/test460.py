#Test Name: SquareSheetShelfStressFSEstar
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from matestar import *

md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = md.extrude(3, 1.)
md.materials = matestar()
md.materials.rheology_B = 3.15e8 * np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_Ec = np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_Es = 3 * np.ones((md.mesh.numberofvertices, ))
md.cluster = generic('name', gethostname(), 'np', 3)

#Go solve
field_names = []
field_tolerances = []
field_values = []
#md.initialization.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface - md.mesh.y)
for i in ['SSA', 'HO', 'FS']:
    md = setflowequation(md, i, 'all')
    md = solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vz' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [7e-06, 2e-05, 2e-06, 5e-06, 8e-07]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vz,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]
