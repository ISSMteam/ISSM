#Test Name: ValleyGlacierLevelsetCalvingSSA2dCrevassedepth
import numpy as np
from calvingcrevassedepth import *
from socket import gethostname
from model import *
from parameterize import *
from setflowequation import *
from setmask import *
from solve import *
from triangle import *

md = triangle(model(), '../Exp/Square.exp', 50000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/ValleyGlacierShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.timestepping.time_step = 2
md.timestepping.final_time = 50
md.cluster = generic('name', gethostname(), 'np', 3)

#Transient
md.transient.isstressbalance = 1
md.transient.ismovingfront = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 1

md.calving = calvingcrevassedepth()
md.calving.crevasse_opening_stress=1
md.calving.water_height = 50 * np.ones((md.mesh.numberofvertices, ))
md.frontalforcings.meltingrate = np.zeros((md.mesh.numberofvertices, ))
md.levelset.spclevelset = np.nan * np.ones((md.mesh.numberofvertices, ))
md.levelset.reinit_frequency = 1
md.levelset.migration_max = 1e10

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Thickness1', 'Surface1', 'MaskIceLevelset1'
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Thickness2', 'Surface2', 'MaskIceLevelset2'
               'Vx10', 'Vy10', 'Vel10', 'Pressure10', 'Thickness10', 'Surface10', 'MaskIceLevelset10']
field_tolerances = [1e-8, 1e-8, 1e-8, 1e-9, 1e-9, 1e-9, 3e-9,
                    1e-8, 1e-8, 1e-8, 1e-9, 1e-9, 1e-9, 3e-9,
                    1e-8, 1e-8, 1e-8, 1e-9, 1e-9, 1e-9, 3e-9]
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
                md.results.TransientSolution[9].Vx,
                md.results.TransientSolution[9].Vy,
                md.results.TransientSolution[9].Vel,
                md.results.TransientSolution[9].Pressure,
                md.results.TransientSolution[9].Thickness,
                md.results.TransientSolution[9].Surface,
                md.results.TransientSolution[9].MaskIceLevelset]
