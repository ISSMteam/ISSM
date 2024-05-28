#Test Name: PigBamgMesh
import numpy as np
import copy
from model import *
from bamg import *
from setmask import *
from parameterize import *
from ComputeHessian import *
from ComputeMetric import *


#Simple mesh 1
hVertices = 10000. * np.ones((27))
hVertices[0:5] = 1000.
md = bamg(model(), 'domain', '../Exp/Pig.exp', 'hmax', 20000., 'hVertices', hVertices, 'gradation', 3.)
x1 = md.mesh.x
y1 = md.mesh.y

#Simple mesh 2
md = bamg(model(), 'domain', '../Exp/Pig.exp', 'hmax', 10000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
x2 = md.mesh.x
y2 = md.mesh.y

#refine existing mesh 1
hessian = ComputeHessian(md.mesh.elements, md.mesh.x, md.mesh.y, md.inversion.vy_obs, 'node')
metric = ComputeMetric(hessian, 2. / 9., 1., 1000., 25. * 10.**3, [])
md.miscellaneous.dummy = metric
md2 = bamg(copy.deepcopy(md), 'metric', md.miscellaneous.dummy, 'hmin', 1000., 'hmax', 20000., 'gradation', 3.)
x3 = md2.mesh.x
y3 = md2.mesh.y

#refine existing mesh 2
md2 = bamg(copy.deepcopy(md), 'metric', md.miscellaneous.dummy, 'hmin', 1000., 'hmax', 20000., 'gradation', 3., 'anisomax', 1.)
x4 = md2.mesh.x
y4 = md2.mesh.y

#refine existing mesh 3
hVertices = np.nan * np.ones((md.mesh.numberofvertices))
hVertices[np.nonzero(md.mesh.vertexonboundary)] = 500.
md2 = bamg(copy.deepcopy(md), 'metric', md.miscellaneous.dummy, 'hmin', 1000., 'hmax', 20000., 'gradation', 3., 'anisomax', 1., 'hVertices', hVertices)
x5 = md2.mesh.x
y5 = md2.mesh.y

#refine existing mesh 4
md2 = bamg(copy.deepcopy(md), 'field', md.inversion.vy_obs, 'hmin', 1000., 'hmax', 20000., 'gradation', 3., 'Hessiantype', 0, 'err', np.array([1.]))
x6 = md2.mesh.x
y6 = md2.mesh.y

#refine existing mesh 5
md2 = bamg(copy.deepcopy(md), 'field', np.vstack((md.inversion.vy_obs, md.geometry.thickness)).T, 'hmin', 1000., 'hmax', 20000., 'gradation', 3., 'Hessiantype', 1, 'err', np.array([[10., 100.]]))
x7 = md2.mesh.x
y7 = md2.mesh.y

#Fields and tolerances to track changes
field_names = ['x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'x5', 'y5', 'x6', 'y6', 'x7', 'y7']
field_tolerances = [2e-10, 7e-10, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7]
