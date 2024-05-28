import numpy as np
from paterson import paterson
from SetIceShelfBC import SetIceShelfBC
#Start defining model parameters here

print('      creating thickness')
hmin = 300
hmax = 1000
ymin = np.nanmin(md.mesh.y)
ymax = np.nanmax(md.mesh.y)
md.geometry.thickness = hmax + (hmin - hmax) * (md.mesh.y - ymin) / (ymax - ymin)
md.geometry.base = - md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness

print('      creating drag')
md.friction.coefficient = np.where(md.mask.ocean_levelset < 0., 0, 200)
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

print('      initial velocity')
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.vel = np.zeros((md.mesh.numberofvertices))

print('      creating flow law parameter')
md.materials.rheology_B = paterson((273 - 20) * np.ones((md.mesh.numberofvertices)))
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements))
md.damage.D = np.zeros((md.mesh.numberofvertices))

print('      creating boundary conditions')
md = SetIceShelfBC(md, 'Front.exp')
