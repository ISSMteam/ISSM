import numpy as np
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
md.friction.coefficient = np.where(md.mask.ocean_levelset < 0., 0., 200)
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

print('      creating flow law paramter')
md.materials.rheology_B = np.where(md.mesh.x < md.mesh.y, 1.4 * 1e8, 1.8 * 1e8)
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements))

print('      creating boundary conditions')
md = SetIceShelfBC(md, 'Front.exp')
