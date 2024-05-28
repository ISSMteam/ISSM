#Geometry specific to Experiments A and B
import os.path
import inspect

import numpy as np

from arch import *
from InterpFromMeshToMesh2d import *
from materials import *
from paterson import *
from SetIceSheetBC import *
from verbose import *


rad = 600000.
nv = md.mesh.numberofvertices
if (np.isnan(md.geometry.thickness)):
    md.geometry.thickness = np.zeros((md.mesh.numberofvertices, 1))
for i in range(nv):
    dist = np.sqrt(md.mesh.x[i]**2 + md.mesh.y[i]**2)
    if (dist <= rad):
        md.geometry.thickness[i] = 2000.0
    else:
        md.geometry.thickness[i] = 0

md.geometry.thickness = md.geometry.thickness.reshape(-1, 1)
md.geometry.base = np.zeros((md.mesh.numberofvertices, 1))
md.geometry.surface = md.geometry.thickness + md.geometry.base.reshape(-1, 1)  #would otherwise create a 91x91 matrix

# GIA parameters specific to Experiments A and B
md.materials = materials('litho', 'ice')
md.materials.radius = [10, 6271e3, 6371e3] # 100km lithosphere, the rest is mantle material
md.materials.density = [3.38e3, 3.36e3]
md.materials.lame_u = [1.45e11, 6.7e10]
md.materials.viscosity = [1e21, 1e40]

# Ice density used for benchmarking, not 917 kg/m^3
md.materials.rho_ice = 1000 # kg/m^3

# Initial velocity
x = archread('../Data/SquareSheetConstrained.arch', 'x')
y = archread('../Data/SquareSheetConstrained.arch', 'y')
vx = archread('../Data/SquareSheetConstrained.arch', 'vx')
vy = archread('../Data/SquareSheetConstrained.arch', 'vy')
index = archread('../Data/SquareSheetConstrained.arch', 'index').astype(int)

md.initialization.vx = np.array(InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)).reshape(-1, 1)
md.initialization.vy = np.array(InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)).reshape(-1, 1)
vx = None
vy = None
x = None
y = None
index = None
md.initialization.vz = np.zeros((md.mesh.numberofvertices, 1))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices, 1))

# Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices, 1))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements, 1))

# Friction
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices, 1))
md.friction.coefficient[np.where(md.mask.ocean_levelset < 0.)] = 0.
md.friction.p = np.ones((md.mesh.numberofelements, 1))
md.friction.q = np.ones((md.mesh.numberofelements, 1))

# Numerical parameters
md.groundingline.migration = 'None'
md.masstransport.stabilization = 1
md.thermal.stabilization = 1.
md.verbose = verbose(0)
md.settings.waitonlock = 30.
md.stressbalance.restol = 0.05
md.steadystate.reltol = 0.05
md.stressbalance.reltol = 0.05
md.stressbalance.abstol = np.nan
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.

# Boundary conditions:
md = SetIceSheetBC(md)

# Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
