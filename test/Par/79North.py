import os.path
import inspect
from arch import *
import numpy
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson
from SetMarineIceSheetBC import SetMarineIceSheetBC

#Start defining model parameters here

#Geometry and observation
x = numpy.array(archread('../Data/79North.arch', 'x'))
y = numpy.array(archread('../Data/79North.arch', 'y'))
vx = numpy.array(archread('../Data/79North.arch', 'vx'))
vy = numpy.array(archread('../Data/79North.arch', 'vy'))
index = numpy.array(archread('../Data/79North.arch', 'index')).astype(int)
surface = numpy.array(archread('../Data/79North.arch', 'surface'))
thickness = numpy.array(archread('../Data/79North.arch', 'thickness'))

md.initialization.vx = InterpFromMeshToMesh2d(index, x, y, vx, md.mesh.x, md.mesh.y)[:, 0]
md.initialization.vy = InterpFromMeshToMesh2d(index, x, y, vy, md.mesh.x, md.mesh.y)[:, 0]
md.geometry.surface = InterpFromMeshToMesh2d(index, x, y, surface, md.mesh.x, md.mesh.y)[:, 0]
md.geometry.thickness = InterpFromMeshToMesh2d(index, x, y, thickness, md.mesh.x, md.mesh.y)[:, 0]
md.geometry.base = md.geometry.surface - md.geometry.thickness

#Materials
md.initialization.temperature = (273. - 20.) * numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * numpy.ones((md.mesh.numberofelements))
md.initialization.temperature = md.initialization.temperature

#Friction
md.friction.coefficient = 50. * numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = numpy.ones((md.mesh.numberofelements))
md.friction.q = numpy.ones((md.mesh.numberofelements))

#Ice shelf melting and surface mass balance
md.basalforcings.floatingice_melting_rate = numpy.zeros((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate[numpy.nonzero(md.mask.ocean_levelset < 0.)[0]] = 0.
md.basalforcings.groundedice_melting_rate = numpy.zeros((md.mesh.numberofvertices))
md.smb.mass_balance = 15 * numpy.ones((md.mesh.numberofvertices))

#Numerical parameters
md.masstransport.stabilization = 1
md.thermal.stabilization = 1
md.verbose = verbose(0)
md.settings.waitonlock = 30
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.
md.stressbalance.restol = 0.05
md.stressbalance.reltol = 0.005
md.steadystate.reltol = 0.005
md.stressbalance.abstol = float('NaN')
md.groundingline.migration = 'None'

#Boundary conditions:
md = SetMarineIceSheetBC(md)
pos = numpy.nonzero(md.mesh.vertexonboundary)
md.balancethickness.spcthickness[pos] = md.geometry.thickness[pos]
md.masstransport.spcthickness[pos] = md.geometry.thickness[pos]

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
