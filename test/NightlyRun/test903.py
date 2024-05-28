#Test Name: SquareNoDynExtrudedHydrologyDCOneLayer
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from solve import solve
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 100000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareNoDyn.py')
md.cluster = generic('name', gethostname(), 'np', 1)

md.transient.ishydrology = True
md.hydrology = hydrologydc()
md.hydrology = md.hydrology.initialize(md)

md.hydrology.isefficientlayer = 0
md.hydrology.sedimentlimit_flag = 1
md.hydrology.sedimentlimit = 8000.0
md.hydrology.mask_thawed_node = np.ones((md.mesh.numberofvertices))
md.initialization.sediment_head = np.zeros((md.mesh.numberofvertices))
md.hydrology.spcsediment_head = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.spcsediment_head[np.where(md.mesh.y == 0)] = 0.0

md.basalforcings.groundedice_melting_rate = 2.0 * np.ones((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = 0.0 * np.ones((md.mesh.numberofvertices))
md.hydrology.sediment_transmitivity = 3.0 * np.ones((md.mesh.numberofvertices))

md.timestepping.time_step = 0
md.timestepping.final_time = 1.0
md.extrude(3, 1.)
md = solve(md, 'Hydrology')

#Fields and tolerances to track changes
#you can also compare with an analitic solution, but it is exact
#only if no limits are applied
#analitic=(md.mesh.y.^2 - 2 * md.mesh.y * 1.0e6) * (-2.0 / (2 * md.constants.yts * md.hydrology.sediment_transmitivity))
field_names = ['SedimentWaterHead', 'SedimentHeadResidual']
field_tolerances = [1e-13, 3e-10]
field_values = [md.results.HydrologySolution.SedimentHead, md.results.HydrologySolution.SedimentHeadResidual]
