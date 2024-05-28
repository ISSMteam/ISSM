import os.path
import inspect
from arch import *
import numpy as np
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson
from SetMarineIceSheetBC import SetMarineIceSheetBC

#Start defining model parameters here

#Geometry
hmin = 300.
hmax = 1000.
ymin = min(md.mesh.y)
ymax = max(md.mesh.y)
xmin = min(md.mesh.x)
xmax = max(md.mesh.x)
md.geometry.thickness = hmax + (hmin - hmax) * (md.mesh.y - ymin) / (ymax - ymin) + 0.1 * (hmin - hmax) * (md.mesh.x - xmin) / (xmax - xmin)
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
bed_sheet = -md.materials.rho_ice / md.materials.rho_water * (hmax + (hmin - hmax) * (ymax / 2 - ymin) / (ymax - ymin))
pos = np.nonzero(md.mesh.y <= ymax / 2.)
md.geometry.base[pos] = bed_sheet
md.geometry.surface = md.geometry.base + md.geometry.thickness

#Initial velocity
x = np.array(archread('../Data/SquareSheetShelf.arch', 'x'))
y = np.array(archread('../Data/SquareSheetShelf.arch', 'y'))
vx = np.array(archread('../Data/SquareSheetShelf.arch', 'vx'))
vy = np.array(archread('../Data/SquareSheetShelf.arch', 'vy'))
index = np.array(archread('../Data/SquareSheetShelf.arch', 'index')).astype(int)

md.initialization.vx = InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)
md.initialization.vy = InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#Accumulation and melting
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
md.steadystate.reltol = 0.05
md.stressbalance.reltol = 0.05
md.stressbalance.abstol = float('NaN')
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.
md.groundingline.migration = 'None'

#Deal with boundary conditions:
md = SetMarineIceSheetBC(md, '../Exp/SquareFront.exp')

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
