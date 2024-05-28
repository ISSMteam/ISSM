#Test Name: SquareSheetShelfThermalFSEstar
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
md.materials.rheology_Es = 3. * np.ones((md.mesh.numberofvertices, ))

md = setflowequation(md, 'FS', 'all')
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices, ))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices, ))
md.transient.isstressbalance = 0
md.transient.ismasstransport = 0
md.transient.issmb = 1
md.transient.isthermal = 1
md.transient.isgroundingline = 0
md.thermal.isenthalpy = 1
md.thermal.isdynamicbasalspc = 1
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Enthalpy1', 'Waterfraction1', 'Temperature1',
               'Enthalpy2', 'Waterfraction2', 'Temperature2',
               'Enthalpy3', 'Waterfraction3', 'Temperature3']
field_tolerances = [1e-12, 1e-11, 1e-12,
                    1e-12, 1e-10, 1e-12,
                    1e-12, 1e-9, 1e-12]
field_values = [md.results.TransientSolution[0].Enthalpy,
                md.results.TransientSolution[0].Waterfraction,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[1].Enthalpy,
                md.results.TransientSolution[1].Waterfraction,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[2].Enthalpy,
                md.results.TransientSolution[2].Waterfraction,
                md.results.TransientSolution[2].Temperature]
