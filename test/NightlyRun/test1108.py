#Test Name: ISMIPDFS
import numpy as np
from model import *
from socket import gethostname
from bamg import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from squaremesh import *

from PythonFuncs import *

"""
This test is a test from the ISMP - HOM Intercomparison project.
Pattyn and Payne 2006
"""

#L_list = [5000., 10000., 20000., 40000., 80000., 160000.]
L_list = [80000.]
results = []

for L in L_list:
    nx = 30  #numberof nodes in x direction
    ny = 30
    md = model()
    md = squaremesh(md, L, L, nx, ny)
    md = setmask(md, '', '')  #ice sheet test
    md = parameterize(md, '../Par/ISMIPD.py')
    md.extrude(10, 1.)

    md = setflowequation(md, 'HO', 'all')

#We need one grd on dirichlet: the 4 corners are set to zero
    md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))

    maxx = max(md.mesh.x)
    maxy = max(md.mesh.y)
    posA = np.where(md.mesh.vertexonbase)
    posB = np.unique(np.concatenate((np.where(md.mesh.x == 0.), np.where(md.mesh.x == maxx))))
    posC = np.unique(np.concatenate((np.where(md.mesh.y == 0.), np.where(md.mesh.y == maxy))))
    pos = np.intersect1d(np.intersect1d(posA, posB), posC)

    md.stressbalance.spcvx[pos] = 0.
    md.stressbalance.spcvy[pos] = 0.
    md.stressbalance.spcvz[pos] = 0.

#Create MPCs to have periodic boundary conditions
    posx = np.nonzero(md.mesh.x == 0.)[0]
    posx2 = np.nonzero(md.mesh.x == np.max(md.mesh.x))[0]

    posy = np.intersect1d(np.intersect1d(np.where(md.mesh.y == 0.), np.where(md.mesh.x != 0.)), np.where(md.mesh.x != np.max(md.mesh.x)))[0]  #Don't take the same nodes two times
    posy2 = np.intersect1d(np.intersect1d(np.where(md.mesh.y == np.max(md.mesh.y)), np.where(md.mesh.x != 0.)), np.where(md.mesh.x != np.max(md.mesh.x)))[0]

    md.stressbalance.vertex_pairing = np.vstack((np.hstack((posx.reshape(-1, 1) + 1, posx2.reshape(-1, 1) + 1)), np.hstack((posy.reshape(-1, 1) + 1, posy2.reshape(-1, 1) + 1))))

#Compute the stressbalance
    md.cluster = generic('name', gethostname(), 'np', 8)
    md.verbose = verbose('convergence', True)
    md = solve(md, 'Stressbalance')
    md.stressbalance.reltol = np.nan
    md.stressbalance.abstol = np.nan
    md.stressbalance.vertex_pairing = np.empty((0, 2))
#We need one grid on dirichlet: the 4 corners are set to zero
    md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))
    pos = np.nonzero(np.logical_or.reduce((md.mesh.y == 0., md.mesh.x == 0., md.mesh.x == np.max(md.mesh.x), md.mesh.y == np.max(md.mesh.y))))  #Don't take the same nodes two times
    md.stressbalance.spcvx[pos] = np.squeeze(md.results.StressbalanceSolution.Vx[pos])
    md.stressbalance.spcvy[pos] = np.squeeze(md.results.StressbalanceSolution.Vy[pos])
    md = setflowequation(md, 'FS', 'all')
    md = solve(md, 'Stressbalance')

#Plot the results and save them
    vx = md.results.StressbalanceSolution.Vx
    vy = md.results.StressbalanceSolution.Vy
    vz = md.results.StressbalanceSolution.Vz
    results.append(md.results.StressbalanceSolution)

#       plotmodel(md, 'data', vx, 'data', vy, 'data', vz, 'layer  #all', md.mesh.numberoflayers)

#Fields and tolerances to track changes
field_names = ['Vx80km', 'Vy80km', 'Vz80km']
field_tolerances = [1e-08, 1e-07, 1e-07]
field_values = []
for result in results:
    field_values = field_values + [result.Vx, result.Vy, result.Vz]
