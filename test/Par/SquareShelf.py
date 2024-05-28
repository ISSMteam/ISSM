from arch import *
import inspect
import os.path

import numpy as np

from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson
from SetIceShelfBC import SetIceShelfBC
from verbose import verbose

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
md.geometry.surface = md.geometry.base + md.geometry.thickness
md.geometry.bed = md.geometry.base - 500.

#Initial velocity and pressure
x = np.array(archread('../Data/SquareShelf.arch', 'x'))
y = np.array(archread('../Data/SquareShelf.arch', 'y'))
vx = np.array(archread('../Data/SquareShelf.arch', 'vx'))
vy = np.array(archread('../Data/SquareShelf.arch', 'vy'))
index = archread('../Data/SquareShelf.arch', 'index').astype(int)

#dbg - begin
#  #print 'vars in SquareShelf.nc:'
#  #for v in iVelF.variables:
#  #    print v
#dbg - end

# x = x[0][:]
# y = y[0][:]
# vx = vx[0][:]
# vy = vy[0][:]
# index = index[0][:]

md.initialization.vx = InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)
md.initialization.vy = InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)
x = y = vx = vy = index = None
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

#dbg - begin
#print '...vx:'
#print md.initialization.vx
#print '...vy:'
#print md.initialization.vy
#  #print '...vz:'
#  #print md.initialization.vz
#  #print '...pressure:'
#  #print md.initialization.pressure
#dbg - end

# Materials
md.initialization.temperature = (273 - 20) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements))

# Friction
md.friction.coefficient = 20 * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.where(md.mask.ocean_levelset < 0)[0]] = 0
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# Numerical parameters
md.masstransport.stabilization = 1
md.thermal.stabilization = 1
md.settings.waitonlock = 30
md.verbose = verbose(0)
md.stressbalance.restol = 0.10
md.steadystate.reltol = 0.02
md.stressbalance.reltol = 0.02
md.stressbalance.abstol = np.nan
md.timestepping.time_step = 1
md.timestepping.final_time = 3
md.groundingline.migration = 'None'

# Boundary conditions:
md = SetIceShelfBC(md, '../Exp/SquareFront.exp')

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
