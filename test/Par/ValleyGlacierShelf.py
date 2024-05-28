import os.path
from arch import *
import numpy as np
import inspect
from verbose import verbose
from paterson import paterson
from SetIceShelfBC import SetIceShelfBC

#Start defining model parameters here
x = md.mesh.x
y = md.mesh.y
xmin, xmax = min(x), max(x)
ymin, ymax = min(y), max(y)
Lx = (xmax - xmin)
Ly = (ymax - ymin)
xm, ym = (xmin + xmax) / 2., (ymin + ymax) / 2.

#Geometry: U - shaped valley in y direction
thk_center = 1000.
thk_margin = 0.5 * thk_center
bmax = 0.
bmin = -thk_center * md.materials.rho_ice / md.materials.rho_water

alpha = 2. / 3.
slope = 0.9 * (bmin - bmax) * (x - xmin) / (Lx * alpha) + 0.1 * (bmin - bmax) * (y - ymin) / (Ly) + bmax
md.geometry.surface = (thk_center + bmax) + slope
md.geometry.base = bmax + slope + 4. / Ly**2 * (thk_center - thk_margin) * (np.power(y - ym, 2))
md.geometry.thickness = md.geometry.surface - md.geometry.base
md.geometry.bed = md.geometry.base

#Mask
md.mask.ice_levelset = x - alpha * Lx
md.mask.ocean_levelset = np.ones((md.mesh.numberofvertices))

#Initial velocity
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature = (273.15 - 5.) * np.ones((md.mesh.numberofvertices))
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#Thermal
md.thermal.isenthalpy = False
md.thermal.spctemperature = float('nan') * np.ones((md.mesh.numberofvertices))

#Groundingline
md.groundingline.migration = 'SubelementMigration'

#Surface mass balance and basal melting
md.smb.mass_balance = 0.3 * np.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = md.smb.mass_balance
md.basalforcings.floatingice_melting_rate = md.smb.mass_balance

#Friction
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

#Transient
md.transient.isstressbalance = True
md.transient.ismovingfront = True
md.transient.ismasstransport = False
md.transient.isthermal = False
md.transient.isgroundingline = True

#Stressbalance
md.stressbalance.maxiter = 100
md.stressbalance.restol = 0.05
md.stressbalance.reltol = 0.05
md.stressbalance.abstol = float('nan')

#Masstransport
md.calving.calvingrate = 0. * np.ones((md.mesh.numberofvertices))
md.frontalforcings.meltingrate = 0. * np.ones((md.mesh.numberofvertices))
md.levelset.spclevelset = float('NaN') * np.ones((md.mesh.numberofvertices))
md.masstransport.stabilization = 1.

#Numerical parameters
md.thermal.stabilization = 1.
md.settings.waitonlock = 30
md.steadystate.reltol = 0.05
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.

#Verbose
md.verbose = verbose(0)

#Deal with boundary conditions:
md = SetIceShelfBC(md)

#Change name so that no tests have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
