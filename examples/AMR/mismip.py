import numpy as np
from SetIceShelfBC import SetIceShelfBC

# creating thickness
md.geometry.bed = -100 - np.abs(md.mesh.x) / 1000
md.geometry.base = -90 * np.ones((md.mesh.numberofvertices))
md.geometry.surface = 10 * np.ones((md.mesh.numberofvertices))
md.geometry.thickness = md.geometry.surface - md.geometry.base
md.mask.ocean_levelset = -1 * np.ones((md.mesh.numberofvertices))

# creating basal drag
md.friction.coefficient = np.sqrt(10**7) * np.ones((md.mesh.numberofvertices))  #q = 1.
md.friction.p = 3 * np.ones((md.mesh.numberofelements))
md.friction.q = np.zeros((md.mesh.numberofelements))

# creating flow law paramter
md.materials.rheology_B = 1 / ((1e-25)**(1 / 3)) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements))
md.materials.rheology_law = 'None'

# creating boundary conditions
md = SetIceShelfBC(md, './front.exp')
md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))
pos = np.nonzero(np.logical_or((np.abs(md.mesh.y - 50000) < 0.1), np.abs(md.mesh.y) < 0.1))
md.stressbalance.spcvy[pos] = 0
pos2 = np.nonzero(np.abs(md.mesh.x) < 0.1)
md.stressbalance.spcvx[pos2] = 0
md.stressbalance.spcvz[pos] = np.nan
md.stressbalance.spcvz[pos2] = np.nan

# creating forcing conditions
md.smb.mass_balance = 0.5 * np.ones((md.mesh.numberofvertices))
md.basalforcings.geothermalflux = 0.5 * np.ones((md.mesh.numberofvertices))
md.thermal.spctemperature = np.nan * np.ones((md.mesh.numberofvertices))
md.groundingline.migration = 'SubelementMigration'

# setting parameters
md.materials.rho_ice = 900
md.materials.rho_water = 1000
md.constants.g = 9.8
md.constants.yts = 3600 * 24 * 365
md.transient.isthermal = 0
md.transient.isgroundingline = 1
md.stressbalance.isnewton = 0

# setting inital condition
md.initialization.vx = np.ones((md.mesh.numberofvertices))
md.initialization.vy = np.ones((md.mesh.numberofvertices))
md.initialization.vz = np.ones((md.mesh.numberofvertices))
md.initialization.vel = np.sqrt(2) * np.ones((md.mesh.numberofvertices))
md.initialization.pressure = md.constants.g * md.materials.rho_ice * md.geometry.thickness
md.initialization.temperature = 273 * np.ones((md.mesh.numberofvertices))
