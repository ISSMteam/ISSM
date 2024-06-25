#Test Name: SquareBamgMesh
import numpy as np
import time
from model import *
from bamg import *


#Simple mesh
md = bamg(model(), 'domain', '../Exp/Square.exp', 'hmax', 100000.)
x1 = md.mesh.x
y1 = md.mesh.y

#hVertices
md = bamg(model(), 'domain', '../Exp/Square.exp', 'hmax', 300000., 'hVertices', np.array([10000., 100000., 400000., 100000.]).reshape(-1, 1))
x2 = md.mesh.x
y2 = md.mesh.y

#big mesh
t0 = time.time()
md = bamg(model(), 'domain', '../Exp/Square.exp', 'hmax', 3000.)
nbelements = md.mesh.numberofelements
elapsedtime = time.time() - t0
if nbelements > 267895 - 50 and nbelements < 267895 + 50:
    nbewithinrange = 1.
else:
    nbewithinrange = 0.

#Fields and tolerances to track changes
field_names = ['x1', 'y1', 'x2', 'y2', 'nbelements', 'elapsed time']
field_tolerances = [2e-9, 2e-9, 1e-13, 1e-13, 1e-13, 8.5]
field_values = [x1, y1, x2, y2, nbewithinrange, elapsedtime]
