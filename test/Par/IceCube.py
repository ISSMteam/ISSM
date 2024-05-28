import os.path
import numpy as np
import inspect
from verbose import verbose
from paterson import paterson
from SetIceSheetBC import SetIceSheetBC
from arch import *

#Start defining model parameters here

#Geometry
md.geometry.thickness = 1000.0 * np.ones((md.mesh.numberofvertices))
md.geometry.base = np.zeros((md.mesh.numberofvertices))
md.geometry.surface = md.geometry.base + md.geometry.thickness

md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#Calving
md.calving.calvingrate = np.zeros((md.mesh.numberofvertices))
md.levelset.spclevelset = np.nan * np.ones((md.mesh.numberofvertices))

#Friction
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.where(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

#Numerical parameters
md.masstransport.stabilization = 1.
md.thermal.stabilization = 1.
md.verbose = verbose(0)
md.settings.waitonlock = 30
md.stressbalance.restol = 0.05
md.steadystate.reltol = 0.05
md.stressbalance.reltol = 0.05
md.stressbalance.abstol = float('NaN')
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.
md.groundingline.migration = 'None'

#Boundary conditions:
md = SetIceSheetBC(md)

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
