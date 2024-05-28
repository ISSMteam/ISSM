import os.path
from arch import archread
import numpy as np
import inspect
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from nye import nye
from SetIceSheetBC import SetIceSheetBC

#Start defining model parameters here
CO2_temp = 150.
CO2_n = 7.
CO2_meltingPoint = 195.
CO2_latentHeat = 199000.
CO2_rhoIce = 1562.
CO2_heatCapacity = 700.
CO2_thermalCond = 0.5
CO2_dynViscosity = 13.72 * 10**-6
CO2_rhoLiquidZeroDeg = 929.
md.materials.rho_ice = CO2_rhoIce
md.materials.rho_freshwater = CO2_rhoLiquidZeroDeg
md.materials.thermalconductivity = CO2_thermalCond
md.materials.heatcapacity = CO2_heatCapacity
md.materials.meltingpoint = CO2_meltingPoint
md.materials.latentheat = CO2_latentHeat
md.materials.mu_water = CO2_dynViscosity

#Geometry
hmin = 300.
hmax = 1000.
ymin = np.min(md.mesh.y)
ymax = np.max(md.mesh.y)
xmin = min(md.mesh.x)
xmax = max(md.mesh.x)
md.geometry.thickness = hmax + (hmin - hmax) * (md.mesh.y - ymin) / (ymax - ymin) + 0.1 * (hmin - hmax) * (md.mesh.x - xmin) / (xmax - xmin)
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness + 20.
md.geometry.surface = md.geometry.base + md.geometry.thickness

#Initial velocity
x = np.array(archread('../Data/SquareSheetConstrained.arch', 'x'))
y = np.array(archread('../Data/SquareSheetConstrained.arch', 'y'))
vx = np.array(archread('../Data/SquareSheetConstrained.arch', 'vx'))
vy = np.array(archread('../Data/SquareSheetConstrained.arch', 'vy'))
index = archread('../Data/SquareSheetConstrained.arch', 'index').astype(int)

md.initialization.vx = InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)
md.initialization.vy = InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature = CO2_temp * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = nye(md.initialization.temperature, 1)
md.materials.rheology_n = CO2_n * np.ones((md.mesh.numberofelements))

#Surface mass balance and basal melting
md.smb.mass_balance = 10. * np.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = 5. * np.ones((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = 5. * np.ones((md.mesh.numberofvertices))

#Friction
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

#Numerical parameters
md.masstransport.stabilization = 1
md.thermal.stabilization = 1
md.verbose = verbose(0)
md.settings.waitonlock = 30
md.stressbalance.restol = 0.05
md.stressbalance.reltol = 0.05
md.steadystate.reltol = 0.05
md.stressbalance.abstol = float('nan')
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.
md.groundingline.migration = 'None'

#Deal with boundary conditions:
md = SetIceSheetBC(md)

#Change name so that no tests have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
