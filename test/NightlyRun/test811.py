#Test Name: ValleyGlacierLevelsetCalvingMOLHO2d
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from SetMOLHOBC import SetMOLHOBC
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 50000)
md = setmask(md, '', '')
md = parameterize(md, '../Par/ValleyGlacierShelf.py')
md = setflowequation(md, 'MOLHO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#Transient
md.transient.isstressbalance = True
md.transient.ismasstransport = True
md.transient.issmb = True
md.transient.isthermal = False
md.transient.isgroundingline = True
md.transient.ismovingfront = True

md.calving.calvingrate = 1000. * np.ones((md.mesh.numberofvertices))
md.frontalforcings.meltingrate = np.zeros((md.mesh.numberofvertices))
md.levelset.migration_max = 1e10

md = SetMOLHOBC(md)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Thickness1', 'Surface1', 'MaskIceLevelset1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Thickness2', 'Surface2', 'MaskIceLevelset2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Thickness3', 'Surface3', 'MaskIceLevelset3']
field_tolerances = [2e-13, 4e-13, 2e-13, 3e-13, 3e-13, 2e-13, 2e-13,
                    2e-13, 2e-13, 2e-13, 3e-13, 3e-13, 2e-13, 2e-13,
                    2e-13, 6e-13, 2e-13, 3e-13, 3e-13, 2e-13, 2e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].MaskIceLevelset,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].MaskIceLevelset,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].MaskIceLevelset]
