#Test Name: SquareShelfLevelsetCalvingSSA2dMinThickness
import numpy as np
from calvingminthickness import *
from socket import gethostname
from model import *
from parameterize import *
from setflowequation import *
from setmask import *
from solve import *
from triangle import *

md = triangle(model(), '../Exp/Square.exp', 30000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#Do not kill ice bergs as all is floating
md.levelset.kill_icebergs = 0

x = md.mesh.x
xmin = min(x)
xmax = max(x)
Lx = (xmax - xmin)
alpha = 2. / 3.
md.mask.ice_levelset = -1 + 2 * (md.mesh.y > 9e5)

md.timestepping.time_step = 1
md.timestepping.final_time = 3

#Transient
md.transient.isstressbalance = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 0
md.transient.ismovingfront = 1

md.calving = calvingminthickness()
md.calving.min_thickness = 400
md.frontalforcings.meltingrate = np.zeros((md.mesh.numberofvertices, ))
md.levelset.spclevelset = np.nan * np.ones((md.mesh.numberofvertices, ))
md.levelset.reinit_frequency = 1
md.levelset.migration_max = 1e10

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Thickness1', 'Surface1', 'MaskIceLevelset1'
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Thickness2', 'Surface2', 'MaskIceLevelset2'
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Thickness3', 'Surface3', 'MaskIceLevelset3']
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
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].MaskIceLevelset]
