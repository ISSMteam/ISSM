import os.path
import numpy as np
import copy
import inspect
from paterson import paterson
from verbose import verbose

#Start defining model parameters here

di = md.materials.rho_ice / md.materials.rho_water
rad = 1.e6
shelfextent = 2.e5
#Geometry
hmin = 300.
hmax = 1000.
radius = np.sqrt(md.mesh.x * md.mesh.x + md.mesh.y * md.mesh.y.reshape(-1))
ymin = np.min(radius)
ymax = np.max(radius)
md.geometry.thickness = hmax + (hmin - hmax) * (radius - ymin) / (ymax - ymin)
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness

pos = np.nonzero(md.mask.ocean_levelset > 0.)[0]
md.geometry.base[pos] = md.geometry.base[pos] - 300. * (radius[pos] - (rad - shelfextent)) / (rad - shelfextent)
md.geometry.surface = md.geometry.base + md.geometry.thickness

pos = np.nonzero(radius < 200000.)
md.geometry.thickness[pos] = 100.
md.geometry.base[pos] = -di * md.geometry.thickness[pos] - 20.
md.geometry.surface[pos] = md.geometry.base[pos] + md.geometry.thickness[pos]

pos = np.nonzero(np.logical_and(np.logical_and(md.mesh.x < 0.2 * 1.e6, md.mesh.x > -0.2 * 1.e6), md.mesh.y > 0.))
md.geometry.thickness[pos] = 100.
md.geometry.base[pos] = -di * md.geometry.thickness[pos] - 20.
md.geometry.surface[pos] = md.geometry.base[pos] + md.geometry.thickness[pos]

pos = np.nonzero(np.logical_and(np.logical_and(md.mesh.x < 0.1 * 1.e6, md.mesh.x > -0.1 * 1.e6), np.logical_and(md.mesh.y < -0.5 * 1.e6, md.mesh.y > -0.6 * 1.e6)))
md.geometry.thickness[pos] = 100.
md.geometry.base[pos] = -di * md.geometry.thickness[pos] - 20.
md.geometry.surface[pos] = md.geometry.base[pos] + md.geometry.thickness[pos]

#plug holes into the ice sheet, to test for grounding line migration.
di = md.materials.rho_ice / md.materials.rho_water
rad = np.sqrt(md.mesh.x**2 + md.mesh.y**2)
pos = np.nonzero(rad < 200000.)
md.geometry.thickness[pos] = 100.
md.geometry.base[pos] = -di * md.geometry.thickness[pos] - 20.
md.geometry.surface[pos] = md.geometry.base[pos] + md.geometry.thickness[pos]

pos = np.nonzero(np.logical_and(np.logical_and(md.mesh.x < 0.2 * 1.e6, md.mesh.x > -0.2 * 1.e6), md.mesh.y > 0.))
md.geometry.thickness[pos] = 100.
md.geometry.base[pos] = -di * md.geometry.thickness[pos] - 20.
md.geometry.surface[pos] = md.geometry.base[pos] + md.geometry.thickness[pos]

pos = np.nonzero(np.logical_and(np.logical_and(md.mesh.x < 0.1 * 1.e6, md.mesh.x > -0.1 * 1.e6), np.logical_and(md.mesh.y < -0.5 * 1.e6, md.mesh.y > -0.6 * 1.e6)))
md.geometry.thickness[pos] = 100.
md.geometry.base[pos] = -di * md.geometry.thickness[pos] - 20.
md.geometry.surface[pos] = md.geometry.base[pos] + md.geometry.thickness[pos]

#Initial velocity
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#Surface mass balance and basal melting
md.smb.mass_balance = -10. * np.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
pos = np.nonzero(md.mask.ocean_levelset > 0.)[0]
md.basalforcings.groundedice_melting_rate[pos] = 10.
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.geothermalflux = np.ones((md.mesh.numberofvertices))

#Friction
radius = 1.e6
shelfextent = 2.e5
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices))
xelem = np.mean(md.mesh.x[md.mesh.elements.astype(int) - 1], axis=1)
yelem = np.mean(md.mesh.y[md.mesh.elements.astype(int) - 1], axis=1)
rad = np.sqrt(xelem**2 + yelem**2)
flags = np.zeros(md.mesh.numberofelements)
pos = np.nonzero(rad >= (radius - shelfextent))
md.friction.coefficient[md.mesh.elements[pos, :] - 1] = 0.
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
md.timestepping.time_step = 5.
md.timestepping.final_time = 5.

#bathymetry and grounding line migration:
md.groundingline.migration = 'AggressiveMigration'
md.geometry.bed = copy.deepcopy(md.geometry.base)
pos = np.nonzero(md.mask.ocean_levelset < 0.)[0]
md.geometry.bed[pos] = md.geometry.base[pos] - 900.

#Deal with boundary conditions:
md.stressbalance.spcvx = float('nan') * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvy = float('nan') * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvz = float('nan') * np.ones((md.mesh.numberofvertices))

pos = np.nonzero(np.logical_and(md.mesh.x == 0, md.mesh.y == 0))
md.stressbalance.spcvx[pos] = 0
md.stressbalance.spcvy[pos] = 0

pos = np.nonzero(md.mesh.vertexonboundary)
md.mask.ice_levelset[pos] = 0
md.balancethickness.spcthickness = float('nan') * np.ones((md.mesh.numberofvertices))
md.masstransport.spcthickness = float('nan') * np.ones((md.mesh.numberofvertices))
md.stressbalance.referential = float('nan') * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.loadingforce = 0 * np.ones((md.mesh.numberofvertices, 3))
md.thermal.spctemperature = 737. * np.ones((md.mesh.numberofvertices))

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
