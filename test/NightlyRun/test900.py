#Test Name: SquareNoDynUnConfinedHydroDC
import numpy as np
from model import *
from setmask import *
from triangle import triangle
from parameterize import parameterize
from solve import solve
from socket import gethostname
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 100000.)
md = setmask(md, '', '')
#reduced square (20m long)
md.mesh.x = md.mesh.x / 5.0e4
md.mesh.y = md.mesh.y / 5.0e4
md = parameterize(md, '../Par/SquareNoDyn.py')
md.cluster = generic('name', gethostname(), 'np', 1)

md.transient.ishydrology = True
md.hydrology = hydrologydc()
md.hydrology = md.hydrology.initialize(md)


#Hydro Model Parameters
md.hydrology.isefficientlayer = 0
md.hydrology.sedimentlimit_flag = 0
md.hydrology.mask_thawed_node = np.ones((md.mesh.numberofvertices))
md.hydrology.rel_tol = 1.0e-6
md.hydrology.penalty_lock = 0
md.hydrology.max_iter = 200
md.hydrology.transfer_flag = 0
md.hydrology.unconfined_flag = 1
#Sediment
md.hydrology.sediment_porosity = 0.1
md.hydrology.sediment_thickness = 10.0
md.hydrology.sediment_transmitivity = (1.0e-3 * md.hydrology.sediment_thickness) * np.ones((md.mesh.numberofvertices))
#init
md.initialization.sediment_head = -5.0 * np.ones((md.mesh.numberofvertices))
#BC
md.hydrology.spcsediment_head = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.spcsediment_head[np.where(md.mesh.x == 0)] = 0.5

md.timestepping.time_step = 5 / md.constants.yts  #5s steppin
md.settings.output_frequency = 2
md.timestepping.final_time = 300 / md.constants.yts  #500s run

md = solve(md, 'Transient')

#fields to track, results can also be found in
#Wang 2009 Fig 6b (journal of Hydrology)
field_names = ['SedimentWaterHead1', 'SedimentWaterHead2']
field_tolerances = [1e-13, 1e-13]
field_values = [md.results.TransientSolution[10].SedimentHead, md.results.TransientSolution[30].SedimentHead]
