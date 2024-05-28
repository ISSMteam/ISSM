import numpy as np
from paterson import paterson
from SetMarineIceSheetBC import SetMarineIceSheetBC

#Ok, start defining model parameters here

md.timestepping.time_step = 0
md.groundingline.migration = 'None'

print("      creating thickness")
h = 1000.
md.geometry.thickness = h * np.ones((md.mesh.numberofvertices))
md.geometry.base = -1000. * np.ones((md.mesh.numberofvertices))
md.geometry.surface = md.geometry.base + md.geometry.thickness

print("      creating velocities")
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))

print("      creating drag")
md.friction.coefficient = 200. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

print("      creating temperatures")
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices, ))
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices, ))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices, ))

print("      creating flow law parameter")
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

print("      creating surface mass balance")
md.smb.mass_balance = np.ones((md.mesh.numberofvertices)) / md.constants.yts  #1m / a
#md.basalforcings.melting_rate = 0. * np.ones((md.mesh.numberofvertices)) / md.constants.yts  #1m / a
md.basalforcings.groundedice_melting_rate = 0. * np.ones((md.mesh.numberofvertices)) / md.constants.yts  #1m / a
md.basalforcings.floatingice_melting_rate = 0. * np.ones((md.mesh.numberofvertices)) / md.constants.yts  #1m / a

#Deal with boundary conditions:

print("      boundary conditions for stressbalance model")
md = SetMarineIceSheetBC(md, '../Exp/SquareFront.exp')

print("      boundary conditions for thermal model")
md.thermal.spctemperature[:] = md.initialization.temperature
md.basalforcings.geothermalflux = np.zeros((md.mesh.numberofvertices))
md.basalforcings.geothermalflux[np.nonzero(md.mask.ocean_levelset > 0.)[0]] = 1. * 10**-3  #1 mW / m^2
